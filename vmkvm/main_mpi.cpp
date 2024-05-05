#include <mpi.h>

#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <omp.h>

#define assert(x) do { if(x) continue; std::fprintf(stderr,\
        "\033[41m Failed assert at " __FILE__ ":%d \033[0m\n\t\033[31m" # x \
        "\033[0m\n", (int)__LINE__); MPI_Abort(MPI_COMM_WORLD, 17); } while(0)
#define SUCCESS(x) assert((x) == MPI_SUCCESS)

#include "csr.h"
#include "adj.h"

void mpi_exit(int code) {
    SUCCESS(MPI_Finalize());
    std::exit(code);
}

struct cli_params {
    bool print_topo, print;
    size_t N[4];
    long ntr, its;
    cli_params() : print_topo(false), print(true), ntr(1), its(10) {
        std::memset(N, 0, sizeof(N));
    }
    void handle_args(int rank, int argc, char **args) {
        if (argc >= 3) {
            int j = 0;
            for (int i = 1; i < argc; ++i) {
                std::string s(args[i]);
                if (s == "--help") goto HELP;
                else if (j * sizeof(*N) < sizeof(N)) N[j++] = std::stoll(args[i]);
                else if (s == "--topo") print_topo = true;
                else if (s == "--no-print") print = false;
                else if (s == "--omp-ntr") {
                    assert(i < argc - 1);
                    ntr = std::stol(args[i+1]);
                }
                else if (s == "--its") {
                    assert(i < argc - 1);
                    its = std::stol(args[i+1]);
                }
            }
        }
        else {
HELP:
            if (rank == 0) {
                std::cerr << "Usage: " << args[0]
                    << " Nx Ny Px Py [--topo] [--no-print] [--omp-ntr N] [--its N]\n";
            }
            mpi_exit(1);
        }
    }
};

/*
   Grid Nx * Ny with (Nx+1) * (Ny+1) nodes is split into Px * Py parts
   */

inline size_t ibeg(size_t i, size_t n, size_t p) {
    // split array of size n to p pieces so that all have similar size
    // max difference in size is one
    // return beginning of i-th piece
    size_t rem = n % p;
    return i * (n / p) + (i < rem ? i : rem);
    // first rem pieces get n/p + 1 elements, last n-rem get n/p elements
    // ibeg(p) == n, ibeg(0) = 0
}

size_t gen_test_topo(size_t Nx, size_t Ny, size_t px, size_t py, size_t Px, size_t Py, CSR<void> &ans) {
    // square grid of nodes of size Nx,Ny (grid of elements of size (Nx-1)(Ny-1))
    // every third quad is split into two triangles
    assert(Px < Nx && Py < Ny && px < Px && py < Py);
    size_t xbeg = ibeg(px, Nx-1, Px), xend = ibeg(px+1, Nx-1, Px);
    size_t ybeg = ibeg(py, Ny-1, Py), yend = ibeg(py+1, Ny-1, Py);
#if 1
    if (xbeg > 0) --xbeg;
    if (xend < Nx-1) ++xend;
    if (ybeg > 0) --ybeg;
    if (yend < Ny-1) ++yend; // increments and decrements to capture interface elements
#endif
    ans.clear();
    ans.ri.push_back(0);
    for (size_t x = xbeg; x < xend; ++x) {
        for (size_t y = ybeg; y < yend; ++y) {
            size_t i = x * Ny + y;
            if ((x * (Ny-1) + y) % 3 == 0) {
                size_t A[] = {i, i+1, i+Ny, i+1, i+Ny, i+Ny+1};
                for (size_t *a = A; sizeof(*A) * (a-A) < sizeof(A); a += 3) {
                    int c = 0;
                    for (int j = 0; j < 3; ++j) {
                        ++c;
                        ans.d.push_back(a[j]);
                    }
                    ans.ri.push_back(ans.ri.back() + c);
                }
            }
            else {
                size_t a[] = {i, i+1, i+Ny+1, i+Ny};
                int c = 0;
                for (int j = 0; j < 4; ++j) {
                    ++c;
                    ans.d.push_back(a[j]);
                }
                ans.ri.push_back(ans.ri.back() + c);
            }
        }
    }
    size_t mem_usage = ans.mem_usage();
    ans.ri.shrink_to_fit();
    ans.d.shrink_to_fit();
    return mem_usage;
}

double get_time() {
    auto c = std::chrono::high_resolution_clock();
    auto now = c.now().time_since_epoch();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now);
    return ns.count() / 1e9;
}

int main(int argc, char **args) {
    double T = get_time(), T2;
#define MT(s) T2 = get_time(); std::cout << s << ": " << T2-T << "s\n"; T = get_time();

    int world_size = -1, my_rank = -1;

    SUCCESS(MPI_Init(NULL, NULL));
    SUCCESS(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
    SUCCESS(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));
    MT("MPI Init");

    struct cli_params par;
    par.handle_args(my_rank, argc, args);
    size_t Nx = par.N[0], Ny = par.N[1];
    size_t Px = par.N[2], Py = par.N[3];
    size_t px = my_rank / Py, py = my_rank % Py;
    assert(px < Px && py < Py);
    assert(Px * Py == (size_t)world_size && world_size > 0);
    assert(Px <= Nx+1 && Py <= Ny+1);
    omp_set_num_threads(par.ntr);
    MT("Input handling");

    CSR<void> topo, nen, eEe;
    CSR<float> mat;

    std::cout << "Generating test topology used " << gen_test_topo(Nx, Ny, px, py, Px, Py, topo)
        << " bytes of RAM\n";
    MT("Generating test topology");

    std::cout << "Creating nen adjacency used " << nodes_to_adj(Nx * Ny, topo, nen)
        << " bytes of RAM\n";
    MT("Creating node adjacency from topology");

    std::cout << "Creating eEe adjacency used " << en_to_eEe(Nx * Ny, topo, eEe)
        << " bytes of RAM\n";
    MT("Creating element adjacency from topology");
    /*if (my_rank == 0) {
        for (int i = 0; sizeof(par.N[0]) * i < sizeof(par.N); ++i) {
            printf("\t%d, ", (int)par.N[i]);
        }
    }
    printf("\n- %d %d\n", world_size, my_rank);*/

    for (int i = 0; i < world_size; ++i) {
        SUCCESS(MPI_Barrier(MPI_COMM_WORLD));
        if (i != my_rank) continue;
        std::cout << "\033[45mRank " << my_rank << "\033[0m\n";
        if (par.print) {
            if (par.print_topo) {
                CSR<void> tt;
                transpose_csr(Nx * Ny, topo, tt);
                print_csr(topo, "topo");
                print_csr(tt, "transposed topo");
                MT("Printing topology");
            }
            print_csr(nen, "nen");
            print_csr(eEe, "eEe");
            //print_csr(mat, "eEe mat");
            MT("Printing adjacency");
        }
        std::cout.flush();
        std::cerr.flush();
    }

    mpi_exit(0);
}
