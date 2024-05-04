#include <mpi.h>

#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <omp.h>

#include "csr.h"

#define SUCCESS(x) assert((x) == MPI_SUCCESS)
#define ASSERT(x) do { if(x) continue; std::fprintf(stderr, "\033[41m Failed assert at " __FILE__ ":%d \033[0m\n", (int)__LINE__); MPI_Abort(MPI_COMM_WORLD, 1); } while(0)

void mpi_exit(int code) {
    SUCCESS(MPI_Finalize());
    std::exit(code);
}

struct cli_params {
    bool print_topo, print;
    size_t N[4];
    long ntr, its;
    cli_params() : print_topo(false), print(true), ntr(1), its(10) {
        std::memset(N, 0xff, sizeof(N));
    }
};

void handle_args(int rank, int argc, char **args, struct cli_params &par) {
    if (argc >= 3) {
        int j = 0;
        for (int i = 1; i < argc; ++i) {
            std::string s(args[i]);
            if (s == "--help") goto HELP;
            else if (s == "--topo") par.print_topo = true;
            else if (s == "--no-print") par.print = false;
            else if (j * sizeof(par.N[0]) < sizeof(par.N)) par.N[j++] = std::stoll(args[i]);
            else if (s == "--omp-ntr") {
                ASSERT(i < argc - 1);
                par.ntr = std::stol(args[i+1]);
            }
            else if (s == "--its") {
                ASSERT(i < argc - 1);
                par.its = std::stol(args[i+1]);
            }
        }
    }
    else {
HELP:
        if (rank == 0) std::cerr << "Usage: " << args[0] << " Nx Ny [--topo] [--no-print]\n";
        mpi_exit(1);
    }
}

int main(int argc, char **args) {
    int world_size = -1, my_rank = -1;

    SUCCESS(MPI_Init(NULL, NULL));
    SUCCESS(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
    SUCCESS(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));

    struct cli_params par;
    handle_args(my_rank, argc, args, par);
    omp_set_num_threads(par.ntr);

    /*if (my_rank == 0) {
        for (int i = 0; sizeof(par.N[0]) * i < sizeof(par.N); ++i) {
            printf("\t%d, ", (int)par.N[i]);
        }
    }
    printf("\n- %d %d\n", world_size, my_rank);*/

    mpi_exit(0);
}
