#pragma once

#include<cstdint>
#include<vector>
#include<array>

#include "point.h"
#include "rect.h"

class BitQuadTree{
    public:
        struct Stats{
            size_t nodes=0;
            size_t internal_nodes = 0;
            size_t leaves = 0;
            int max_depth = 0;
        };

        struct HPStats {
    // Degree histogram for internal nodes only
    // deg_hist[d] = number of internal nodes with degree d (d in {1,2,3,4})
    uint64_t deg_hist[5] = {0,0,0,0,0};

    uint64_t internal_nodes = 0;
    uint64_t nodes = 0;
    uint64_t leaves = 0;

    uint64_t heavy_edges = 0;   // internal_nodes (each internal has 1 heavy edge)
    uint64_t light_edges = 0;   // sum(deg-1) over internal nodes
    double rho = 0.0;           // light_edges / internal_nodes

    // Heavy child dominance
    double avg_heavy_share = 0.0; // mean over internal nodes of maxChildSize / subtreeSize

    // Heavy path length distribution (in edges or nodes)
    double avg_heavy_path_len = 0.0; // mean length in nodes
    int median_heavy_path_len = 0;
    int p90_heavy_path_len = 0;
    int max_heavy_path_len = 0;
};

HPStats analyze_heavy_paths() const;


        size_t bytes_T() const;
        size_t bytes_meta() const;
        size_t bytes_rank() const;


        explicit BitQuadTree(const Rect& world);

        void set_debug(bool on) {debug_ = on;}

        void build(const std::vector<Point>& points);

        void range_query(const Rect& q, std::vector<Point>& out) const;

        size_t bytes_used() const;

        Stats stats() const {return stats_;}

        size_t debug_T_bits() const { return T_bits_; }


    private:
        Rect world_;
        bool debug_ = false;

        std::vector<Point> pts_;
        std::vector<Point> scratch_;

        // T stores 4 bits per internal node, packed into 64-bit words.
        // Bit positions: child bit of internal node i and quadrant q is at pos = 4*i + q.
        std::vector<uint64_t> T_;
        size_t T_bits_ = 0; // number of valid bits currently in T_

        std::vector<uint64_t> metaBits_;
        size_t meta_bits_ = 0; // 3 bits per node   

        std::vector<uint32_t> rank_super_;

        Stats stats_;

        static void geo_split(const Rect& r, int& mx, int& my);
        static int quadrant_of(int mx, int my, const Point& p);
        static Rect quadrant_rect(const Rect& r, int mx, int my, int q);
        static int balance_score(const std::array<int,4>& counts, int n);

        void choose_split_and_encode(const Rect& region,
                                    uint32_t start, uint32_t count,
                                    uint8_t& out_meta,
                                    int& out_mx, int& out_my);

        static void chosen_split(const Rect& region, uint8_t meta, int& mx, int& my);
        
        void push_bit(bool b);
        bool get_bit(size_t pos) const;
        void push_meta3(uint8_t m);
        uint8_t get_meta3(size_t node_i) const;
        uint8_t get_mask(size_t internal_i) const;
        uint32_t rank1(size_t pos) const;
        void build_rank();
        void build_bfs();
};


class SimpleBitQuadTree {
public:
    struct Stats{
        size_t nodes=0;
        size_t internal_nodes=0;
        size_t leaves=0;
        int max_depth=0;
    };

    explicit SimpleBitQuadTree(const Rect& world) : world_(world) {}

    void build(const std::vector<Point>& points);
    void range_query(const Rect& q, std::vector<Point>& out) const;

    // memory accounting
    size_t bytes_T() const { return T_.capacity() * sizeof(uint64_t); }
    size_t bytes_rank() const { return rank_super_.capacity() * sizeof(uint32_t); }
    size_t bytes_used() const { return bytes_T() + bytes_rank(); }

    Stats stats() const { return stats_; }

private:
    Rect world_;
    std::vector<Point> pts_;
    std::vector<Point> scratch_;

    std::vector<uint64_t> T_;
    size_t T_bits_ = 0;

    std::vector<uint32_t> rank_super_;
    Stats stats_;

    static void geo_split(const Rect& r, int& mx, int& my);
    static int quadrant_of(int mx, int my, const Point& p);
    static Rect quadrant_rect(const Rect& r, int mx, int my, int q);

    void push_bit(bool b);
    bool get_bit(size_t pos) const;
    uint8_t get_mask(size_t node_i) const;

    uint32_t rank1(size_t pos) const;
    void build_rank();
    void build_bfs();
};
