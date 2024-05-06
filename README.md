# numerical-practicum
Some numerical tasks

To build & run execute the following in a directory with code:
```. build.sh && ./a.out Nx Ny```

To build & run MPI application execute (P=Px*Py):
```. mpi_build.sh && mpiexec -n P ./a.out Nx Ny Px Py```

# Tests
Tested single-threaded application on my laptop (CPU: Intel(R) Core(TM) i5-8300H CPU @ 2.30GHz, 4 cores, 8 threads; RAM: SODIMM DDR4 Synchronous 2667 MHz, width 64 bits, size 8 GiB); calculated memory bandwidth 2667 MHz * 64 bits should be equal to 21 GB/s, but for some reason actual memory bandwidth measured with benchmark is 25 GB/s (multi-threaded write) or 19 GB/s (single-threaded `memset`). Number of iterations of CG solver is always 10.
| Size of grid (Nx*Ny) | Time of solving stage | CG memory throuhgput | Time of entire application run |
|----------------------|-----------------------|----------------------|--------------------------------|
| 1000 * 1000          | 0.335261 s            | 4.38231 GB/s         | 4.23418 s                      |
| 1414 * 1414          | 0.757708 s            | 3.88945 GB/s         | 8.2245 s                       |
| 2000 * 2000          | 1.5642 s              | 3.76783 GB/s         | 16.6548 s                      |
| 2828 * 2828          | 2.95907 s             | 3.98652 GB/s         | 33.0256 s                      |
