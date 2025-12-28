// experiments.cpp  (QuadTree-only benchmark, same output format)
//
// Usage:
//   ./run_experiments <path_to_points_csv_or_xy>
//
// Output CSV columns:
//   WindowSize, avg_us_per_query, avg_points_per_query
//
// Notes:
// - Uses centered square range queries.
// - Builds ONE QuadTree once and queries it for each window size.
// - Uses g_sink to prevent over-optimization.

static volatile long long g_sink = 0;

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

#include "io.h"
#include "timer.h"
#include "modified_quadtree.h"

// ---------- helpers ----------
static std::vector<Rect> make_square_range_queries(const Rect& region, int qcount, int side, uint64_t seed) {
    // Simple xorshift RNG (local)
    uint64_t s = seed ? seed : 88172645463325252ull;
    auto next_u64 = [&]() {
        uint64_t x = s;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        s = x;
        return x;
    };
    auto next_int = [&](int lo, int hi) { // [lo, hi)
        return lo + (int)(next_u64() % (uint64_t)(hi - lo));
    };

    std::vector<Rect> qs;
    qs.reserve(qcount);

    // ensure side fits (clamp if needed)
    int wx = region.xmax - region.xmin;
    int wy = region.ymax - region.ymin;
    int max_side = std::max(1, std::min(wx, wy) - 1);
    if (side > max_side) side = max_side;

    int max_x0 = region.xmax - side;
    int max_y0 = region.ymax - side;

    // guard
    if (max_x0 <= region.xmin) max_x0 = region.xmin + 1;
    if (max_y0 <= region.ymin) max_y0 = region.ymin + 1;

    for (int i = 0; i < qcount; i++) {
        int x0 = next_int(region.xmin, max_x0);
        int y0 = next_int(region.ymin, max_y0);
        qs.push_back(Rect{x0, y0, x0 + side, y0 + side});
    }
    return qs;
}

static std::vector<Rect> make_centered_square_range_queries(const Rect& region,
                                                            const std::vector<Point>& pts,
                                                            int qcount,
                                                            int side,
                                                            uint64_t seed) {
    // Simple xorshift RNG (local)
    uint64_t s = seed ? seed : 88172645463325252ull;
    auto next_u64 = [&]() {
        uint64_t x = s;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        s = x;
        return x;
    };
    auto next_int = [&](int lo, int hi) { // [lo, hi)
        return lo + (int)(next_u64() % (uint64_t)(hi - lo));
    };

    std::vector<Rect> qs;
    qs.reserve(qcount);

    // ensure side fits (clamp if needed)
    int wx = region.xmax - region.xmin;
    int wy = region.ymax - region.ymin;
    int max_side = std::max(1, std::min(wx, wy) - 1);
    if (side > max_side) side = max_side;

    for (int i = 0; i < qcount; i++) {
        // pick a random existing point as the center
        const Point& c = pts[next_int(0, (int)pts.size())];

        int half = side / 2;

        // propose top-left
        int x0 = c.x - half;
        int y0 = c.y - half;

        // clamp so [x0, x0+side) fits in region
        if (x0 < region.xmin) x0 = region.xmin;
        if (y0 < region.ymin) y0 = region.ymin;
        if (x0 + side > region.xmax) x0 = region.xmax - side;
        if (y0 + side > region.ymax) y0 = region.ymax - side;

        qs.push_back(Rect{x0, y0, x0 + side, y0 + side});
    }

    return qs;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./run_experiments <path_to_points_csv_or_xy>\n";
        return 1;
    }

    std::string path = argv[1];
    std::vector<Point> pts = load_points_xy(path);

    // Compute padded grid size N (power of two)
    int maxCoord = 0;
    for (const auto& p : pts) {
        maxCoord = std::max(maxCoord, std::max(p.x, p.y));
    }
    int N = 1;
    while (N <= maxCoord) N <<= 1;

    Rect region{0, 0, N, N};

    // QuadTree params (match your previous style)
    int leaf_size = 1;

    // Paper-style window sizes + paper-style query count
    std::vector<int> window_sizes = {4, 16, 64, 256, 1024};
    const int Q = 100000;
    const uint64_t seed_base = 123456789ull;

    std::cout << "Loaded points: " << pts.size() << "\n";
    std::cout << "Region: [" << region.xmin << "," << region.ymin
              << " - " << region.xmax << "," << region.ymax << ")\n";
    std::cout << "Grid N=" << N << "\n";
    std::cout << "QuadTree params: leaf_size=" << leaf_size << "\n\n";

    // Build QuadTree once
    QuadTree qt(region, pts, leaf_size);

    std::cout << "WindowSize, avg_us_per_query, avg_points_per_query\n";

    for (int side : window_sizes) {
        // Centered queries (same as your current run)
        auto qs = make_centered_square_range_queries(
            region, pts, Q, side,
            seed_base + (uint64_t)side * 1315423911ull
        );

        std::vector<Point> out;
        out.reserve(128);

        long long total_reported = 0;

        Timer t;
        for (const auto& q : qs) {
            out.clear();
            qt.range_query(q, out);
            total_reported += (long long)out.size();
        }
        double ms_total = t.ms();

        // Prevent optimization
        g_sink += total_reported;

        double avg_us_per_query = (ms_total * 1000.0) / (double)Q;
        double avg_points_per_query = (double)total_reported / (double)Q;

        std::cout << side << ", " << avg_us_per_query << ", " << avg_points_per_query << "\n";
    }

    std::cout << "sink=" << g_sink << "\n";
    return 0;
}