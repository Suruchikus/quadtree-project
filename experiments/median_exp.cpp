#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

#include "io.h"
#include "median_quadtree.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./median_exp <path_to_points_csv_or_xy>\n";
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
    std::cout << "Grid N=" << N << "\n\n";

    // ---- Median-based succinct quadtree ----
    MedianBitQuadTree med(region);
    med.build(pts);
    auto stM = med.stats();

    std::cout << "Median internal=" << stM.internal_nodes
              << " leaves=" << stM.leaves
              << " depth=" << stM.max_depth << "\n";

    // Example query
    Rect Q{0, 0, 3000, 3000};
    std::vector<Point> outM;
    med.range_query(Q, outM);
    std::cout << "Median query returned " << outM.size() << " points\n";

    std::cout << "Median bytes(total allocated)=" << med.bytes_used() << "\n";
    double bppM = (8.0 * (double)med.bytes_used()) / (double)pts.size();
    std::cout << "Median BitsPerPoint=" << bppM << "\n";
    std::cout << "  T bytes=" << med.bytes_T() << "\n";
    std::cout << "  rank bytes=" << med.bytes_rank() << "\n";
    std::cout << "  splits bytes=" << med.bytes_splits() << "\n";

    std::cout.flush();
    return 0;
}
