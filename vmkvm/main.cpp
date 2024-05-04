#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <cassert>
#include <cmath>

#include <omp.h>

#include "csr.h"

size_t gen_test_topo(size_t Nx, size_t Ny, CSR<void> &ans) {
    // square grid of nodes of size Nx,Ny
    // every third quad is split into two triangles
    ans.clear();
    ans.ri.push_back(0);
    for (size_t x = 0; x < Nx-1; ++x) {
        for (size_t y = 0; y < Ny-1; ++y) {
            size_t i = x * Ny + y;
            if ((x * (Ny-1) + y) % 3 == 0) {
                size_t a[] = {i, i+1, i+Ny, i+1, i+Ny, i+Ny+1};
                ans.ri.push_back(ans.ri.back() + 3);
                ans.ri.push_back(ans.ri.back() + 3);
                for (int j = 0; j < 6; ++j) ans.d.push_back(a[j]);
            }
            else {
                size_t a[] = {i, i+1, i+Ny+1, i+Ny};
                ans.ri.push_back(ans.ri.back() + 4);
                for (int j = 0; j < 4; ++j) ans.d.push_back(a[j]);
            }
        }
    }
    size_t mem_usage = ans.mem_usage();
    ans.ri.shrink_to_fit();
    ans.d.shrink_to_fit();
    return mem_usage;
}

template<typename T>
void prefix_sum(T *begin, T *end) {
    for (T *i = begin + 1; i < end; ++i) {
        i[0] += i[-1];
    }
}

size_t unique_rows(CSR<void> &c) {
    std::vector<size_t> I(c.ri.size());
    for (size_t i = 0; i < c.ri.size() - 1; ++i) {
        std::sort(&c.d[c.ri[i]], &c.d[c.ri[i+1]]);
        size_t p = std::unique(&c.d[c.ri[i]], &c.d[c.ri[i+1]]) - &c.d[c.ri[i]];
        I[1+i] = I[i] + p; // waat if i remove 'I[i] +', it crashes
    }
    //prefix_sum(&I[0], &I[I.size() - 1]);
    std::vector<size_t> d(I.back());
    #pragma omp barrier
    #pragma omp parallel for
    for (size_t i = 0; i < c.ri.size() - 1; ++i) {
        for (size_t j = 0; j < I[i+1] - I[i]; ++j) {
            d[I[i] + j] = c.d[c.ri[i] + j];
        }
    }
    c.ri = I;
    c.d = d;
    return c.mem_usage();
}

size_t nodes_to_adj(size_t nodes, const CSR<void> &topo, CSR<void> &adj) {
    size_t mem_usage = 0;
    adj.clear();
    adj.ri.resize(nodes + 1);
    mem_usage += adj.mem_usage();
    FOR_CSR_BEGIN(topo, i, k, j)
        // each node topo.d[j] has topo.ri[i+1] - topo.ri[i] - 1 neighbours
        adj.ri[1 + j] += topo.ri[i+1] - topo.ri[i] - 1;
    FOR_CSR_END
    for (size_t i = 1; i < nodes + 1; ++i) {
        adj.ri[i] += adj.ri[i-1];
    }
    adj.d.resize(adj.ri.back());
    std::vector<size_t> indices(1+nodes);
    FOR_CSR_BEGIN(topo, i, j, l)
        // for each node l add all neighbours to its list beginning at adj.ri[l]
        for (size_t k = topo.ri[i]; k < topo.ri[i+1]; ++k) {
            if (k == j) continue;
            adj.d[adj.ri[l] + indices[l]++] = topo.d[k];
        }
    FOR_CSR_END
    return mem_usage + unique_rows(adj);
}

size_t transpose_csr(size_t nodes, const CSR<void> &C, CSR<void> &T) {
    T.clear();
    T.ri.resize(1 + nodes);

    FOR_CSR_BEGIN(C, i, k, j)
        ++T.ri[1 + j];
    FOR_CSR_END

    prefix_sum(&*T.ri.begin(), &*T.ri.end());
    std::vector<size_t> I(nodes, 0);
    T.d.resize(T.ri.back());

    FOR_CSR_BEGIN(C, i, k, j)
        T.d[T.ri[j] + I[j]++] = i;
        assert(I[j] <= T.ri[1+j] - T.ri[j]);
    FOR_CSR_END

    return VEC_MEM_USAGE(I) + T.mem_usage();
}

bool has_edge(size_t a, size_t b, size_t el, const CSR<void> &en) {
    if (a > b) {
        std::swap(a, b);
    }
    for (size_t k = en.ri[el]; k < en.ri[el+1]; ++k) {
        size_t j = en.d[k], nj = en.d[k+1 < en.ri[el+1] ? k+1 : en.ri[el]];
        if (j > nj) {
            if (a == nj && b == j) return true;
        }
        else if (a == j && b == nj) return true;
    }
    return false;
}

size_t en_to_eEe(size_t nodes, const CSR<void> &en, CSR<void> &eEe) {
    CSR<void> T;
    size_t MU = transpose_csr(nodes, en, T);
    eEe.clear();
    eEe.ri.resize(en.ri.size());

    FOR_CSR_BEGIN(en, i, k, j)
        size_t nj = en.d[k+1 < en.ri[i+1] ? k+1 : en.ri[i]];
        for (size_t ek = T.ri[j]; ek < T.ri[j+1]; ++ek) {
            size_t e = T.d[ek];
            if (has_edge(j, nj, e, en)) {
                ++eEe.ri[1 + e];
            }
        }
    FOR_CSR_END

    prefix_sum(&*eEe.ri.begin(), &*eEe.ri.end());

    std::vector<size_t> I(en.ri.size() - 1);
    eEe.d.resize(eEe.ri.back());
    FOR_CSR_BEGIN(en, i, k, j)
        size_t nj = en.d[k+1 < en.ri[i+1] ? k+1 : en.ri[i]];
        for (size_t ek = T.ri[j]; ek < T.ri[j+1]; ++ek) {
            size_t e = T.d[ek];
            if (e == i || has_edge(j, nj, e, en)) {
                eEe.d[eEe.ri[e] + I[e]++] = i;
            }
        }
    FOR_CSR_END

    unique_rows(eEe);

    return MU + VEC_MEM_USAGE(I) + eEe.mem_usage();
}

size_t get_sample_matrix(const CSR<void> &adj, CSR<float> &mat) {
    mat.clear();
    mat.ri = adj.ri;
    mat.d = adj.d;
    mat.a.resize(mat.d.size());
    float sum = 1e300;
    size_t kdiag = SIZE_MAX / 2;
    FOR_CSR_BEGIN(mat, i, k, j)
        if (k == mat.ri[i]) sum = 0;
        if (i == j) kdiag = k;
        else sum += std::abs(mat.a[k] = std::cos(i*j + i + j));
        if (k+1 == mat.ri[i+1]) mat.a[kdiag] = sum * 1.9;
    FOR_CSR_END
    return mat.mem_usage();
}
void get_sample_rhs(std::vector<float> &v) {
    #pragma omp parallel for
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = std::sin(i);
    }
}

double get_time() {
    auto c = std::chrono::high_resolution_clock();
    auto now = c.now().time_since_epoch();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now);
    return ns.count() / 1e9;
}

template<typename f>
void axpby(f a, std::vector<f> &x, f b, const std::vector<f> &y) {
    assert(x.size() == y.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = a * x[i] + b * y[i];
    }
}

template<typename f>
f dot(const std::vector<f> &a, const std::vector<f> &b) {
    assert(a.size() == b.size());
    f sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

template<typename f>
void spmv(const CSR<f> &m, const std::vector<f> &a, std::vector<f> &b) {
    #pragma omp parallel for
    FOR_CSR_BEGIN(m, i, k, j)
        b[i] += m.a[k] * a[j];
    FOR_CSR_END
}

template<typename f>
void fill(f c, std::vector<f> &v) {
    #pragma omp parallel for
    for (size_t i = 0; i < v.size(); ++i) v[i] = c;
}

template<typename f>
void CG(const CSR<f> &m, const std::vector<f> &b, std::vector<f> &x, int n) {
    using vecf = std::vector<f>;
    assert(x.size() == b.size());
    vecf r = b, z(b.size()), q(b.size()), p(b.size());
    f rho = 0;
    for (int i = 0; i < n; ++i) {
        f rho_2 = dot(r, r);
        if (i == 0) p = r;
        else {
            f bet = rho_2 / rho;
            axpby(bet, p, (f)1., r);
        }
        fill((f)0, q);
        spmv(m, p, q);
        rho = rho_2;
        f alp = rho / dot(p, q);
        axpby((f)1., x, alp, p);
        axpby((f)1., r, -alp, q);
    }
}

int handle_args(int argc, char **args, bool &print_topo, bool &print, size_t N[2], long &ntr, long &its) {
    if (argc >= 3) {
        int j = 0;
        for (int i = 1; i < argc; ++i) {
            if (std::string(args[i]) == "--topo") {
                print_topo = true;
            }
            else if (std::string(args[i]) == "--help") {
                goto HELP;
            }
            else if (std::string(args[i]) == "--no-print") {
                print = false;
            }
            else if (j < 2) {
                N[j++] = std::stoll(args[i]);
            }
            else if (std::string(args[i]) == "--omp-ntr") {
                assert(i < argc - 1);
                ntr = std::stol(args[i+1]);
                omp_set_num_threads(ntr);
            }
            else if (std::string(args[i]) == "--its") {
                assert(i < argc - 1);
                its = std::stol(args[i+1]);
            }
        }
    }
    else {
HELP:
        std::cerr << "Usage: " << args[0] << " Nx Ny [--topo] [--no-print]\n";
        return 1;
    }
    return 0;
}


int main(int argc, char **args) {
    double T = get_time(), T2;
#define MT(s) T2 = get_time(); std::cerr << s << ": " << T2-T << "s\n"; T = get_time();
    bool print_topo = false, print = true;
    size_t N[2];
    long ntr = 1, its = 10;
    if (handle_args(argc, args, print_topo, print, N, ntr, its)) return 1;
    MT("Input handling");

    size_t nodes = N[0] * N[1];
    CSR<void> topo, nen, eEe;
    CSR<float> mat;

    std::cout << "Generating test topology used " << gen_test_topo(N[0], N[1], topo)
        << " bytes of RAM\n";
    MT("Generating test topology");

    std::cout << "Creating nen adjacency used " << nodes_to_adj(nodes, topo, nen)
        << " bytes of RAM\n";
    MT("Creating node adjacency from topology");

    std::cout << "Creating eEe adjacency used " << en_to_eEe(nodes, topo, eEe)
        << " bytes of RAM\n";
    MT("Creating element adjacency from topology");

    std::cout << "Creating sample matrix used " << get_sample_matrix(eEe, mat)
        << " bytes of RAM\n";
    MT("Creating matrix");

    std::vector<float> rhs(mat.ri.size() - 1),
        x(mat.ri.size() - 1), r(mat.ri.size() - 1, 0);

    get_sample_rhs(rhs);
    std::cout << "Generating sample rhs used " << VEC_MEM_USAGE(rhs) << " bytes\n";
    MT("RHS took");

    spmv(mat, x, r);
    axpby(1.f, r, -1.f, rhs);
    std::cout << "\033[31mResidual\033[0m: " << dot(r, r) << "\n";
    r.clear(); r.resize(x.size());
    MT("res");
    CG(mat, rhs, x, its);
    MT("solving.");
    spmv(mat, x, r);
    axpby(1.f, r, -1.f, rhs);
    MT("Calc residual");
    std::cout << "\033[31mResidual\033[0m: " << dot(r, r) << "\n";
    MT("Dot&print residual");


    if (print) {
        if (print_topo) {
            CSR<void> tt;
            transpose_csr(N[0] * N[1], topo, tt);
            print_csr(topo, "topo");
            print_csr(tt, "transposed topo");
            MT("Printing topology");
        }
        print_csr(nen, "nen");
        print_csr(eEe, "eEe");
        print_csr(mat, "eEe mat");
        MT("Printing adjacency");
    }
}
