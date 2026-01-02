#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

#include "io.h"
#include "timer.h"
#include "quadtree.h"


int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./run_experiments <path_to_points_csv_or_xy>\n";
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

    // QuadTree params 
    int leaf_size = 1;


    std::cout << "Loaded points: " << pts.size() << "\n";
    std::cout << "Region: [" << region.xmin << "," << region.ymin << " - " << region.xmax << "," << region.ymax << ")\n";
    std::cout << "Grid N=" << N << "\n";
    std::cout << "QuadTree params: leaf_size=" << leaf_size << "\n\n";

    std::cout.flush();

    int max_level = 0;
    int temp = N;
    while (temp > 1) {
        temp >>= 1;
        max_level++;
    }


    QuadTree ptr(region, leaf_size);
    ptr.set_debug(false);   // turn off per-level printing for benchmark runs
    ptr.set_unit_cell_mode(true);
    ptr.build(pts);
    auto st = ptr.compute_stats();

    std::cout << "Tree nodes=" << st.nodes << "\n";
    std::cout << "Internal nodes=" << st.internal_nodes << "\n";
    std::cout << "Leaves=" << st.leaves << "\n";

    std::cout << "Max depth=" << st.max_depth << "\n";
    std::cout << "Levels=" << (st.max_depth + 1) << "\n";

    std::cout << "Normal splits=" << st.normal_splits << "\n";
    std::cout << "Heavy splits=" << st.heavy_splits << "\n";

    std::cout << "Pointer QuadTree build DONE\n";
    Rect testQ{0, 0, 3000, 3000};  //[0,0)x[3000x3000)
    std::vector<Point> out;
    ptr.range_query(testQ, out);
    std::cout << "Test query [0,0)-(100,100) returned " << out.size() << " points\n";

    std::cout.flush();
}
