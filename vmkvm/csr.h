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

#define FOR_CSR_BEGIN(csr, i, k, j) for (size_t i = 0; i < csr.ri.size() - 1; ++i) for (size_t k = csr.ri[i]; k < csr.ri[i+1]; ++k) { size_t j = csr.d[k];
#define FOR_CSR_END }

template<typename datat>
void print_csr(CSR<datat> c, const std::string &name) {
    std::cout << "CSR " << name << ":\n";
    FOR_CSR_BEGIN(c, i, k, j)
        if (k == c.ri[i]) std::cout << "    ";
        std::cout << j << " ";
        if (k+1 == c.ri[i+1]) std::cout << "\n";
    FOR_CSR_END
}

#endif
