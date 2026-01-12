#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>

#include "point.h"
#include "rect.h"

// Stop at count==1, store point as (x_off,y_off) inside region using 2*(D-depth) bits.
// Uses BFS-level layout: nodes at same depth are contiguous, so we store S + P1 per level.
class SingletonStopQuadTreeV2 {
public:
    struct Stats{
        size_t nodes = 0;
        size_t internal_nodes = 0;
        size_t leaves = 0;             // unit leaves + singleton leaves
        size_t singleton_leaves = 0;   // count==1 and region > 1x1
        int max_depth = 0;
    };

    explicit SingletonStopQuadTreeV2(const Rect& world);

    void build(const std::vector<Point>& points);
    void range_query(const Rect& q, std::vector<Point>& out) const;

    // memory accounting
    size_t bytes_T() const { return T_.capacity() * sizeof(uint64_t); }
    size_t bytes_rankT() const { return rankT_.capacity() * sizeof(uint32_t); }

    size_t bytes_S_levels() const;
    size_t bytes_rankS_levels() const;
    size_t bytes_P1_levels() const;

    size_t bytes_used() const {
        return bytes_T() + bytes_rankT() + bytes_S_levels() + bytes_rankS_levels() + bytes_P1_levels();
    }

    Stats stats() const { return stats_; }
    int D() const { return D_; }

private:
    Rect world_;
    int D_ = 0; // N=2^D

    std::vector<Point> pts_;
    std::vector<Point> scratch_;

    // T: 4 bits per node (BFS over all nodes)
    std::vector<uint64_t> T_;
    size_t T_bits_ = 0;
    std::vector<uint32_t> rankT_;

    // BFS level boundaries: node indices [level_start[d], level_start[d+1]) are depth d
    std::vector<uint32_t> level_start_; // size max_depth+2 (last is nodes)

    // Per-level singleton markers and payloads
    // Slev[d]: 1 bit per node at depth d (1 => singleton leaf)
    std::vector<std::vector<uint64_t>> Slev_;
    std::vector<size_t> Slev_bits_;
    std::vector<std::vector<uint32_t>> rankSlev_;

    // P1lev[d]: packed stream of offsets. Each singleton uses 2*(D-d) bits.
    std::vector<std::vector<uint64_t>> P1lev_;
    std::vector<size_t> P1lev_bits_;

    Stats stats_;

    // geometry
    static void geo_split(const Rect& r, int& mx, int& my);
    static int quadrant_of(int mx, int my, const Point& p);
    static Rect quadrant_rect(const Rect& r, int mx, int my, int q);

    // bit ops
    void push_bit(std::vector<uint64_t>& bv, size_t& bits_used, bool b);
    bool get_bit(const std::vector<uint64_t>& bv, size_t pos) const;

    void push_mask4(uint8_t mask4);
    uint8_t get_mask4(size_t node_i) const;

    // rank over generic bitvector
    uint32_t rank1(const std::vector<uint64_t>& bv,
                   const std::vector<uint32_t>& sup,
                   size_t bits_used,
                   size_t pos) const;

    void build_rank(const std::vector<uint64_t>& bv, size_t bits_used, std::vector<uint32_t>& out_sup);

    uint32_t rank1_T(size_t pos) const { return rank1(T_, rankT_, T_bits_, pos); }

    uint32_t rank1_Slev(int d, size_t pos) const {
        return rank1(Slev_[d], rankSlev_[d], Slev_bits_[d], pos);
    }

    // payload pack/unpack (LSB-first)
    void push_bits(std::vector<uint64_t>& bv, size_t& bits_used, uint64_t value, int nbits);
    uint64_t get_bits(const std::vector<uint64_t>& bv, size_t bitpos, int nbits) const;

    void push_singleton_payload(int depth, const Rect& region, const Point& p);
    Point read_singleton_payload(int depth, const Rect& region, uint32_t singleton_idx) const;

    // build
    void build_bfs();
    void ensure_level_storage(int depth);
};
