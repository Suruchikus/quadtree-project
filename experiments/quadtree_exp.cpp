#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "io.h"
#include "quadtree.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./singleton_exp <path_to_points_xy>\n";
        return 1;
    }

    std::string path = argv[1];
    std::vector<Point> pts = load_points_xy(path);

    if (pts.empty()) {
        std::cerr << "No points loaded.\n";
        return 1;
    }

    // Compute padded grid size N (power of two)
    int maxCoord = 0;
    for (const auto& p : pts) {
        maxCoord = std::max(maxCoord, std::max(p.x, p.y));
    }
    int N = 1;
    while (N <= maxCoord) N <<= 1;

    Rect region{0, 0, N, N};

    std::cout << "Loaded points: " << pts.size() << "\n";
    std::cout << "Region: [" << region.xmin << "," << region.ymin
              << " - " << region.xmax << "," << region.ymax << ")\n";
    std::cout << "Grid N=" << N << " (D=";
    {
        int D = 0, t = N;
        while (t > 1) { t >>= 1; D++; }
        std::cout << D << ")\n\n";
    }

    // Query rect (same as your experiments)
    Rect Q{0, 0, 3000, 3000};

    SingletonStopQuadTreeV2 tree(region);
    tree.build(pts);
    auto st = tree.stats();

    std::cout << "SingletonV2 nodes=" << st.nodes
              << " internal=" << st.internal_nodes
              << " leaves=" << st.leaves
              << " singletons=" << st.singleton_leaves
              << " depth=" << st.max_depth << "\n";

    std::vector<Point> out;
    tree.range_query(Q, out);
    std::cout << "SingletonV2 query returned " << out.size() << " points\n";

    std::cout << "SingletonV2 bytes(total allocated)=" << tree.bytes_used() << "\n";
    double bpp = (8.0 * (double)tree.bytes_used()) / (double)pts.size();
    std::cout << "SingletonV2 BitsPerPoint=" << bpp << "\n";

    std::cout << "  T bytes=" << tree.bytes_T() << "\n";
    std::cout << "  rankT bytes=" << tree.bytes_rankT() << "\n";
    std::cout << "  S(levels) bytes=" << tree.bytes_S_levels() << "\n";
    std::cout << "  rankS(levels) bytes=" << tree.bytes_rankS_levels() << "\n";
    std::cout << "  P1(levels) bytes=" << tree.bytes_P1_levels() << "\n";

    return 0;
}
