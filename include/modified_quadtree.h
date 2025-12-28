#pragma once

#include <cstdint>
#include <vector>
#include "tree_base.h"

class QuadTree : public BaseTree {
public:
    QuadTree(const Rect& region, const std::vector<Point>& points, int leaf_size);

    bool member(const Point& p) const override;
    void range_query(const Rect& r, std::vector<Point>& out) const override;
    Stats stats() const override;

private:
    struct Node {
        Rect region;
        std::vector<Point> points;   // used only for leaves
        Node* child[4];

        // Split encoding:
        // mode bit0 => use/store mx_store; otherwise use geometric mx
        // mode bit1 => use/store my_store; otherwise use geometric my
        uint8_t mode = 0;
        int mx_store = 0;
        int my_store = 0;

        Node(const Rect& r) : region(r) {
            child[0] = child[1] = child[2] = child[3] = nullptr;
        }

        bool is_leaf() const {
            return child[0] == nullptr &&
                   child[1] == nullptr &&
                   child[2] == nullptr &&
                   child[3] == nullptr;
        }
    };

    struct BuildTask {
        Node* node;
        std::vector<Point> points;
    };

    Node* root = nullptr;
    int leaf_size = 0;

private:
    // Uses node's (possibly stored) split lines
    static int quadrant(const Node* node, const Point& p);

    // split line helpers
    static int geo_mx(const Rect& r) { return (r.xmin + r.xmax) / 2; }
    static int geo_my(const Rect& r) { return (r.ymin + r.ymax) / 2; }

    static int node_mx(const Node* node);
    static int node_my(const Node* node);

    static bool region_can_split(const Rect& r);

    // median helpers (do not fully sort)
    static int median_x(std::vector<Point>& pts);
    static int median_y(std::vector<Point>& pts);

    // choose mode + optional stored medians by comparing imbalance
    static void choose_split_mode(Node* node, std::vector<Point>& pts);

    void build_iterative(const Rect& region, const std::vector<Point>& points);

    // optional: clean up (not required by your base, but good to have)
    static void destroy(Node* n);
};