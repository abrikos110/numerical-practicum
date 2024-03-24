#ifndef CSR_HEADER
#define CSR_HEADER

#include <vector>


template<typename data_type>
struct CSR {
    std::vector<size_t> ri; // row index
    std::vector<data_type> d;
    void clear() {
        ri.clear();
        d.clear();
    }
#define VEC_MEM_USAGE(v) (v.capacity() * sizeof(v[0]))
    inline size_t mem_usage() {
        return VEC_MEM_USAGE(ri) + VEC_MEM_USAGE(d);
    }
};

#define FOR_EACH_CSR(csr, i, k) for (size_t i = 0; i < csr.ri.size() - 1; ++i) for (size_t j = c.ri[i]; j < c.ri[i+1]; ++j)

template<typename datat>
void print_csr(CSR<datat> c, const std::string &name) {
    std::cout << "CSR " << name << ":\n";
    FOR_EACH_CSR(c, i, j) {
        if (j == c.ri[i]) std::cout << "    ";
        std::cout << c.d[j] << " ";
        if (j+1 == c.ri[i+1]) std::cout << "\n";
    }
}

#endif
