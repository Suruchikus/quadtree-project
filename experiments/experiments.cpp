#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

#include "io.h"
#include "timer.h"
#include "array_quadtree.h"
#include "bit_quadtree.h"


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


    // BitQuadTree sq(region);
    // sq.build(pts);
    // auto sta = sq.stats();
    // std::cout << "Succinct internal=" << sta.internal_nodes << " leaves=" << sta.leaves
    //         << " depth=" << sta.max_depth << "\n";
    Rect Q{0,0,3000,3000};
    // std::vector<Point> out1;
    // sq.range_query(Q, out1);
    // std::cout << "Succinct query returned " << out1.size() << " points\n";

    // std::cout << "Succinct bytes(total allocated)=" << sq.bytes_used() << "\n";

    // double bits_per_point = (8.0 * (double)sq.bytes_used()) / (double)pts.size();
    // std::cout << "Succinct BitsPerPoint(impl)=" << bits_per_point << "\n";

    // double bits_per_node = (8.0 * (double)sq.bytes_used()) / (double)sta.internal_nodes; 
    // std::cout << "Succinct BitsPerInternalNode(impl)=" << bits_per_node << "\n";
    // std::cout << "  T bytes=" << sq.bytes_T() << "\n";
    // std::cout << "  meta bytes=" << sq.bytes_meta() << "\n";
    // std::cout << "  rank bytes=" << sq.bytes_rank() << "\n";

    // auto hp = sq.analyze_heavy_paths();

    // std::cout << "\n--- Heavy-Path Decomposition Indicators ---\n";
    // std::cout << "nodes=" << hp.nodes
    //         << " internal=" << hp.internal_nodes
    //         << " leaves=" << hp.leaves << "\n";

    // std::cout << "deg1=" << hp.deg_hist[1]
    //         << " deg2=" << hp.deg_hist[2]
    //         << " deg3=" << hp.deg_hist[3]
    //         << " deg4=" << hp.deg_hist[4] << "\n";

    // std::cout << "light_edges=" << hp.light_edges
    //         << " rho(light/internal)=" << hp.rho << "\n";

    // std::cout << "avg_heavy_share=" << hp.avg_heavy_share << "\n";

    // std::cout << "heavy_path_len(avg)=" << hp.avg_heavy_path_len
    //         << " median=" << hp.median_heavy_path_len
    //         << " p90=" << hp.p90_heavy_path_len
    //         << " max=" << hp.max_heavy_path_len << "\n";

    // ---- Simple geo-split succinct quadtree baseline ----
SimpleBitQuadTree simple(region);
simple.build(pts);
auto staS = simple.stats();

std::cout << "SimpleGeo internal=" << staS.internal_nodes
          << " leaves=" << staS.leaves
          << " depth=" << staS.max_depth << "\n";

std::vector<Point> outS;
simple.range_query(Q, outS);
std::cout << "SimpleGeo query returned " << outS.size() << " points\n";

std::cout << "SimpleGeo bytes(total allocated)=" << simple.bytes_used() << "\n";
double bppS = (8.0 * (double)simple.bytes_used()) / (double)pts.size();
std::cout << "SimpleGeo BitsPerPoint=" << bppS << "\n";
std::cout << "  T bytes=" << simple.bytes_T() << "\n";
std::cout << "  rank bytes=" << simple.bytes_rank() << "\n";

    std::cout.flush();
}
