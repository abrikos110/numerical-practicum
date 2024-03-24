#include <iostream>
#include <vector>
#include <string>

#include "csr.h"

size_t gen_test_topo(size_t Nx, size_t Ny, CSR<size_t> &ans) {
    // square grid of nodes of size Nx,Ny
    // every third quad is split into two triangles
    ans.clear();
    ans.ri.push_back(0);
    for (size_t x = 0; x < Nx-1; ++x) {
        for (size_t y = 0; y < Ny-1; ++y) {
            size_t i = x * Ny + y;
            if (i % 3 == 0) {
                size_t a[] = {i, i+1, i+Ny, i+1, i+Ny, i+Ny+1};
                ans.ri.push_back(ans.ri.back() + 3);
                ans.ri.push_back(ans.ri.back() + 3);
                for (int j = 0; j < 6; ++j) ans.d.push_back(a[j]);
            }
            else {
                size_t a[] = {i, i+1, i+Ny, i+Ny+1};
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

size_t nodes_to_adj(size_t num_nodes, const CSR<size_t> &topo, CSR<size_t> &adj) {
    size_t mem_usage = 0;
    adj.clear();
    adj.ri.resize(num_nodes + 1);
    mem_usage += adj.mem_usage();
    for (size_t i = 0; i < topo.ri.size() - 1; ++i) {
        for (size_t j = topo.ri[i]; j < topo.ri[i+1]; ++j) {
            // each node topo.d[j] has topo.ri[i+1] - topo.ri[i] - 1 neighbours
            adj.ri[1 + topo.d[j]] += topo.ri[i+1] - topo.ri[i] - 1;
        }
    }
    for (size_t i = 1; i < num_nodes + 1; ++i) {
        adj.ri[i] += adj.ri[i-1];
    }
    adj.d.resize(adj.ri.back());
    std::vector<size_t> indices(1+num_nodes);
    for (size_t i = 0; i < topo.ri.size() - 1; ++i) {
        for (size_t j = topo.ri[i]; j < topo.ri[i+1]; ++j) {
            // for each node l add all neighbours to its list beginning at adj.ri[l]
            size_t l = topo.d[j];
            for (size_t k = topo.ri[i]; k < topo.ri[i+1]; ++k) {
                if (k == j) continue;
                adj.d[adj.ri[l] + indices[l]] = topo.d[k];
                ++indices[l];
            }
        }
    }
    indices[0] = 0; // new row index
    for (size_t i = 0; i < adj.ri.size() - 1; ++i) {
        std::sort(&adj.d[adj.ri[i]], &adj.d[adj.ri[i+1]]);
        size_t p = std::unique(&adj.d[adj.ri[i]], &adj.d[adj.ri[i+1]])
            - &adj.d[adj.ri[i]];
        indices[1+i] = indices[i] + p;
    }
    std::vector<size_t> new_d(indices[num_nodes]);
    for (size_t i = 0; i < adj.ri.size() - 1; ++i) {
        for (size_t j = 0; j < indices[i+1] - indices[i]; ++j) {
            new_d[indices[i] + j] = adj.d[adj.ri[i] + j];
        }
    }
    adj.ri = indices;
    adj.d = new_d;
    return mem_usage + adj.mem_usage();
}

double get_time() {
    auto c = std::chrono::high_resolution_clock();
    auto now = c.now().time_since_epoch();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now);
    return ns.count() / 1e9;
}


int main(int argc, char **args) {
    double T = get_time();
#define MT(s) std::cerr << s << ": " << get_time() - T << "s\n"; T = get_time();
    bool print_topo = false, print = true;
    size_t N[2];
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
        }
    }
    else {
HELP:
        std::cerr << "Usage: " << args[0] << " Nx Ny [--topo] [--no-print]\n";
        return 1;
    }
    MT("Input handling");

    CSR<size_t> topo, adj;
    std::cout << "Generating test topology used " << gen_test_topo(N[0], N[1], topo)
        << " bytes of RAM\n";
    MT("Generating test topology");
    std::cout << "Creating adjacency used " << nodes_to_adj(N[0] * N[1], topo, adj)
        << " bytes of RAM\n";
    MT("Creating node adjacency from topology");

    if (print) {
        if (print_topo) {
            print_csr(topo, "topo");
            MT("Printing topology");
        }
        print_csr(adj, "adj");
        MT("Printing adjacency");
    }
}
