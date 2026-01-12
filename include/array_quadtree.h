#pragma once

#include <cstdint>
#include <vector>
#include <array>

#include "point.h"
#include "rect.h"


class ArrayQuadTree {
public:

// Node stored in array: children are indices, -1 means absent
    struct Node {
        uint8_t meta = 0;     // bit0: useHeavy, bits1-2: heavy quadrant under split A
        int32_t child[4];     // indices into nodes_ or -1

        Node() {
            child[0] = child[1] = child[2] = child[3] = -1;
        }
    };
    
    struct Stats {
        size_t nodes = 0;
        size_t internal_nodes = 0;
        size_t leaves = 0;
        size_t normal_splits = 0;
        size_t heavy_splits = 0;
        int max_depth = 0;
    };

    size_t bytes_used() const;
    double bits_per_node() const;
    double bits_per_point(size_t n_points) const;


    explicit ArrayQuadTree(const Rect& world);

    void set_debug(bool on) { debug_ = on; }

    // Build from points (unique, in bounds). After build, points are not kept.
    void build(const std::vector<Point>& points);

    // Range query on half-open rectangle q: [xmin,xmax) x [ymin,ymax)
    // Returns all points in q. In Option 1, each leaf outputs exactly one point.
    void range_query(const Rect& q, std::vector<Point>& out) const;

    Stats compute_stats() const;

    size_t node_count() const { return nodes_.size(); }

private:

    // Build task (stack item) used only during build
    struct BuildTask {
        int32_t node_id;
        Rect region;
        uint32_t start;
        uint32_t count;
        int level;
    };

    Rect world_;
    bool debug_ = false;

    // Array of nodes (root is node 0)
    std::vector<Node> nodes_;

    // Build-time point storage (freed after build)
    std::vector<Point> pts_;
    std::vector<Point> scratch_;

    // ---- Geometry helpers ----
    static void geo_split(const Rect& r, int& mx, int& my);
    static int quadrant_of(int mx, int my, const Point& p);
    static Rect quadrant_rect(const Rect& r, int mx, int my, int q);
    static int balance_score(const std::array<int,4>& counts, int n);

    // Choose split A vs heavy split B, encode node meta, and output chosen (mx,my)
    void choose_split_and_encode(int32_t node_id,
                                 const Rect& region,
                                 uint32_t start, uint32_t count,
                                 int& out_mx, int& out_my,
                                 std::array<int,4>& out_countsA,
                                 std::array<int,4>& out_countsB);

    // Reconstruct chosen split (mx,my) from region + meta (used during queries)
    static void chosen_split(const Rect& region, uint8_t meta, int& mx, int& my);

    // Iterative build implementation
    void build_iterative();
};
