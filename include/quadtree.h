#pragma once
#include <cstdint>  //for bits representation
#include <vector>

#include "point.h"
#include "rect.h"

struct Node{
    // meta bits:
    // bit0: useHeavy (0 = normal split, 1 = heavy split)
    // bit1-bit2: heavyQuadrant (0..3) valid only if useHeavy==1
    uint8_t meta = 0; 

    // Children pointers (only created for non-empty quadrants)
    Node* child[4] = {nullptr, nullptr, nullptr, nullptr};

    // Leaf points are stored as a range in one global points array:
    // points live in pts[start .. start+count)
    // I am testing with 1 point per leaf, but we can optimize and if the number of points are less than height at that node we can conver to leaf and store all points there. Will think a bit more about this.
    uint32_t start = 0;
    uint32_t count = 0;

    Node() = default;
    ~Node();
};

//This is for building the tree iteratively. Basically, instead of using recursion we have our own Task stack so that we don't get stack overflow problems
struct BuildTask {
    Node* node;
    Rect region;
    uint32_t start;
    uint32_t count;
    int level;
};


class QuadTree {
    public:
        QuadTree(const Rect& world, int leaf_size = 1);

        //Bild tree from list of points from input
        void build(const std::vector<Point>& input);

        void set_debug (bool on) {debug_ =on;}

        Node* root() const {return root_;}

        const std::vector<Point>& points() const {return pts_;}

        struct Stats {
            size_t nodes = 0;
            size_t internal_nodes = 0;
            size_t leaves = 0;
            size_t normal_splits = 0; // internal nodes only
            size_t heavy_splits  = 0; // internal nodes only
            int max_depth = 0;        // root depth = 0
        };
        void range_query(const Rect& q, std::vector<Point>& out) const;
        Stats compute_stats() const;

        // We build the tree until we reach 1x1 grid. Helps eliminate the need to store points
        void set_unit_cell_mode(bool on) { unit_cell_mode_ = on; }



    private:
        Rect world_;  //The entire region of input
        int leaf_size_ = 1;  //Can tune this
        bool debug_ = false;

        Node* root_=nullptr;

        bool unit_cell_mode_ = false;

        std::vector<Point> pts_;  //All points will be stored here globally, and during build this will be reordered

        //Helper functions for build, search, range queries etc
        static int balance_score(const std::array<int,4>& counts, int n);

        static void geo_split(const Rect& r, int& mx, int& my);

        static int quadrant_of(int mx, int my, const Point& p);

        static Rect quadrant_rect(const Rect& r, int mx, int my, int q);

        void choose_split_and_encode(Node* node, const Rect& region, uint32_t start, uint32_t count, int& out_mx, int& out_my, std::array<int,4>& out_countsA, std::array<int,4>& out_countsB) const;

        void build_iterative();

        // Reconstruct the chosen split lines (mx,my) for this node from region + meta
        static void chosen_split(const Rect& region, uint8_t meta, int& mx, int& my);

        void debug_print_node(int level, const Rect& region, uint32_t start, uint32_t count, int mxA, int myA, const std::array<int,4>& cA, int scoreA, int mxB, int myB, const std::array<int,4>& cB, int scoreB, int mxChosen, int myChosen, uint8_t meta) const;

};


