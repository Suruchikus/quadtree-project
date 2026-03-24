#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "io.h"
#include "quadtree.h"

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

    int N = compute_N_pow2(pts);
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
              << " UM=" << st.UM_bits
              << " UL=" << st.UL_bits
              << " ULD=" << st.ULD_bits
              << " UML=" << st.UML_bits
              << " UMD=" << st.UMD_bits << "\n";

    std::cout << "Bits per point (bpp)=" << st.bpp << "\n\n";
    std::cout << "Normal mixed nodes=" << st.mixed_internal<< "\n\n";


    std::cout << "Build completed.\n";

    qt.print_logical_quadtree_stats();
    return 0;
}