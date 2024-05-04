#include <mpi.h>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <cassert>
#include <cmath>

#include <omp.h>

#include "csr.h"

#define SUCCESS(x) assert((x) == MPI_SUCCESS)

int main(int argc, char **args) {
    int world_size = -1, my_rank = -1;

    SUCCESS(MPI_Init(NULL, NULL));
    SUCCESS(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
    SUCCESS(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));

    printf("- %d %d\n", world_size, my_rank);

    SUCCESS(MPI_Finalize());
}
