#ifndef CSR_HEADER
#define CSR_HEADER

#include <vector>
#include <algorithm>
#include <iostream>

template<typename data_type>
struct CSR {
    std::vector<size_t> ri; // row index
    std::vector<size_t> d;
    std::vector<data_type> a;
    void clear() {
        ri.clear();
        d.clear();
        a.clear();
    }
#define VEC_MEM_USAGE(v) (v.capacity() * sizeof(v[0]))
    inline size_t mem_usage() {
        return VEC_MEM_USAGE(ri) + VEC_MEM_USAGE(d) + VEC_MEM_USAGE(a);
    }
};

template<>
struct CSR<void> {
    std::vector<size_t> ri; // row index
    std::vector<size_t> d;
    void clear() {
        ri.clear();
        d.clear();
    }
    inline size_t mem_usage() {
        return VEC_MEM_USAGE(ri) + VEC_MEM_USAGE(d);
    }
};

#define FOR_CSR_BEGIN(csr, i, k, j) for (size_t i = 0; i < csr.ri.size() - 1; ++i) for (size_t k = csr.ri[i]; k < csr.ri[i+1]; ++k) { size_t j = csr.d[k];
#define FOR_CSR_END }

template<typename datat>
void print_csr(CSR<datat> c, const std::string &name) {
    std::cout << "CSR " << name << ":\n";
    FOR_CSR_BEGIN(c, i, k, j)
        if (k == c.ri[i]) std::cout << "    ";
        std::cout << j << ',' << c.a[k] << " ";
        if (k+1 == c.ri[i+1]) std::cout << "\n";
    FOR_CSR_END
}
template<>
void print_csr(CSR<void> c, const std::string &name) {
    std::cout << "CSR " << name << ":\n";
    FOR_CSR_BEGIN(c, i, k, j)
        if (k == c.ri[i]) std::cout << "    ";
        std::cout << j << " ";
        if (k+1 == c.ri[i+1]) std::cout << "\n";
    FOR_CSR_END
}

template<typename T>
void prefix_sum(T begin, T end) {
    auto prev = *begin;
    ++begin; // begin = a+1
    // here prev = prefix sum of a in range [0, 1)
    //      begin = a + 1
    while (begin < end) {
        // here begin = a+i
        //      prev = prefix sum of a in range [0, i)
        *begin += prev; // we add a[i] to prev and store in a[i]
        prev = *begin;  // make prev = prefix sum of a in range [0, i+1)
        ++begin; // begin = a + i+1
    }
}

size_t transpose_csr(size_t nodes, const CSR<void> &C, CSR<void> &T) {
    T.clear();
    T.ri.resize(1 + nodes);

    FOR_CSR_BEGIN(C, i, k, j)
        ++T.ri[1 + j];
    FOR_CSR_END

    prefix_sum(T.ri.begin(), T.ri.end());
    std::vector<size_t> I(nodes, 0);
    T.d.resize(T.ri.back());

    FOR_CSR_BEGIN(C, i, k, j)
        T.d[T.ri[j] + I[j]++] = i;
        assert(I[j] <= T.ri[1+j] - T.ri[j]);
    FOR_CSR_END

    return VEC_MEM_USAGE(I) + T.mem_usage();
}

size_t unique_rows(CSR<void> &c) {
    std::vector<size_t> I(c.ri.size());
    for (size_t i = 0; i < c.ri.size() - 1; ++i) {
        std::sort(&c.d[c.ri[i]], &c.d[c.ri[i+1]]);
        size_t p = std::unique(&c.d[c.ri[i]], &c.d[c.ri[i+1]]) - &c.d[c.ri[i]];
        I[i+1] = p;
    }
    I[0] = 0;
    prefix_sum(I.begin(), I.end());

    std::vector<size_t> d(I.back());
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

#endif
