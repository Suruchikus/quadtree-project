#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include "io.h"
#include "quadtree_p.h"
#include "membership.h"

double run_membership_test(
    const MXQuadtreeBits& qt,
    const std::vector<Point>& queries,
    uint64_t repeats,
    uint64_t& found_out
) {
    using clock = std::chrono::high_resolution_clock;

    double total_ms = 0.0;
    uint64_t total_found = 0;

    for (uint64_t r = 0; r < repeats; r++) {
        uint64_t found = 0;

        auto start = clock::now();

        for (const auto& q : queries) {
            if (qt.membership(q)) {
                found++;
            }
        }

        auto end = clock::now();

        std::chrono::duration<double, std::milli> elapsed = end - start;

        total_ms += elapsed.count();
        total_found += found;
    }

    found_out = total_found / repeats;

    double avg_ms = total_ms / (double)repeats;

    return avg_ms;
}

static inline int compute_N_pow2(const std::vector<Point>& pts) {
    int maxCoord = 0;
    for (const auto& p : pts) {
        maxCoord = std::max(maxCoord, std::max(p.x, p.y));
    }

    int N = 1;
    while (N <= maxCoord) {
        N <<= 1;
    }
    return N;
}

static inline int compute_D_from_N(int N) {
    int D = 0;
    while (N > 1) {
        N >>= 1;
        D++;
    }
    return D;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./run_experiments <path_to_points_xy>\n";
        return 1;
    }

    std::string path = argv[1];
    std::vector<Point> pts = load_points_xy(path);

    //int N = compute_N_pow2(pts);
    int N = 0;
    if (path.find("gis_sparse") != std::string::npos) {
        N = 1 << 26;
    } else if (path.find("gis_med") != std::string::npos) {
        N = 1 << 22;
    } else if (path.find("gis_dense") != std::string::npos) {
        N = 1 << 19;
    } else if (path.find("rdf_") != std::string::npos) {
        N = 1 << 26;  // ALL RDF datasets use same grid
    } else {
        N = compute_N_pow2(pts);
    }
    int D = compute_D_from_N(N);

    Rect region{0, 0, N, N};

    MXQuadtreeBits qt;
    MXQuadtreeBits::Params params;
    params.D = D;

    qt.build(region, pts, params);
    const auto& st = qt.stats();

    std::cout << "Points = " << st.points << "\n";
    std::cout << "Grid N = " << st.N << "\n";
    std::cout << "Depth D = " << st.D << "\n\n";

    std::cout << "Bits: T=" << st.T_bits
              << " EX=" << st.EX_bits
              << " UL=" << st.UL_bits
              << " ULD=" << st.ULD_bits << "\n";

    std::cout << "Bits per point (bpp)=" << st.bpp << "\n\n";
    std::cout << "Unary to leaf=" << st.unary_to_leaf_nodes << "\n\n";
    std::cout << "leaf=" << st.leaf_nodes << "\n\n";
    std::cout << "fullblock=" << st.fullblock_nodes << "\n\n";
    std::cout << "internal nodes=" << st.internal_nodes << "\n\n";

    double pct_fullblock =
        100.0 * st.fullblock_nodes / (double)st.total_nodes;

    double pct_unary_to_leaf =
        100.0 * st.unary_to_leaf_nodes / (double)st.total_nodes;

    std::cout << "Fullblock % = " << pct_fullblock << "%\n";
    std::cout << "Unary-to-leaf % = " << pct_unary_to_leaf << "%\n\n";

    std::cout << "Build completed.\n";

    const uint64_t Q = 100000;

    MembershipQueryGenerator qgen(region, pts, 1);

    QuerySet filled_qs = qgen.random_filled(Q);
    QuerySet empty_qs = qgen.random_empty(Q);
    QuerySet isolated_qs = qgen.isolated_filled(Q);

   const uint64_t repeats = 100;

    uint64_t found_filled = 0;
    uint64_t found_empty = 0;

    std::cout << "\nRunning membership tests...\n";

    // Filled queries
    double filled_ms = run_membership_test(
        qt,
        filled_qs.queries,
        repeats,
        found_filled
    );

    // Empty queries
    double empty_ms = run_membership_test(
        qt,
        empty_qs.queries,
        repeats,
        found_empty
    );

    uint64_t found_isolated = 0;

    double isolated_ms = run_membership_test(
        qt,
        isolated_qs.queries,
        repeats,
        found_isolated
    );

    double isolated_us_per_query =
        (isolated_ms * 1000.0) / (double)isolated_qs.queries.size();

    // Convert to microseconds per query
    double filled_us_per_query =
        (filled_ms * 1000.0) / (double)filled_qs.queries.size();

    double empty_us_per_query =
        (empty_ms * 1000.0) / (double)empty_qs.queries.size();

    std::cout << "\nMembership Results\n";
    std::cout << "------------------\n";

    std::cout << "Filled queries:\n";
    std::cout << "  avg time = " << filled_us_per_query << " us/query\n";
    std::cout << "  found    = " << found_filled << "\n";

    std::cout << "Empty queries:\n";
    std::cout << "  avg time = " << empty_us_per_query << " us/query\n";
    std::cout << "  found    = " << found_empty << "\n";
    
    std::cout << "Isolated filled queries:\n";
    std::cout << "  avg time = " << isolated_us_per_query << " us/query\n";
    std::cout << "  found    = " << found_isolated << "\n";

    return 0;
}