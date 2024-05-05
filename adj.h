#ifndef ADJACENCY_HEADER
#define ADJACENCY_HEADER

#include "csr.h"

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

#endif
