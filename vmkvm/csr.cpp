#include <iostream>
#include <vector>
#include <string>

template<typename data_type>
struct CSR {
    std::vector<size_t> ri; // row index
    std::vector<data_type> d;
    void clear() {
        ri.clear();
        d.clear();
    }
};

void gen_test_topo(size_t Nx, size_t Ny, CSR<size_t> &ans) {
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
    ans.ri.shrink_to_fit();
    ans.d.shrink_to_fit();
}

void nodes_to_adj(size_t num_nodes, const CSR<size_t> &topo, CSR<size_t> &adj) {
    adj.clear();
    adj.ri.resize(num_nodes + 1);
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
    std::vector<size_t> indices(num_nodes);
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
}

template<typename datat>
void print_csr(CSR<datat> c, const std::string &name) {
    std::cout << "CSR " << name << ":\n";
    for (size_t i = 0; i < c.ri.size() - 1; ++i) {
        std::cout << "    ";
        for (size_t j = c.ri[i]; j < c.ri[i+1]; ++j) {
            std::cout << c.d[j] << " ";
        }
        std::cout << "\n";
    }
}

int main(int argc, char **args) {
    if (argc != 3) {
        std::cerr << "Usage: " << args[0] << " Nx Ny\n";
        return 1;
    }
    size_t Nx = std::stoll(args[1]), Ny = std::stoll(args[2]);
    CSR<size_t> topo, adj;
    gen_test_topo(Nx, Ny, topo);
    nodes_to_adj(Nx * Ny, topo, adj);

    print_csr(topo, "topo");
    print_csr(adj, "adj");
}
