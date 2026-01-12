#pragma once

#include <cstdint>
#include <vector>
#include <array>

#include "point.h"
#include "rect.h"

class MedianBitQuadTree {
public:
    struct Stats{
        size_t nodes = 0;
        size_t internal_nodes = 0;
        size_t leaves = 0;
        int max_depth = 0;
    };

    explicit MedianBitQuadTree(const Rect& world) : world_(world) {}

    void build(const std::vector<Point>& points);
    void range_query(const Rect& q, std::vector<Point>& out) const;

    // memory accounting (match your style)
    size_t bytes_T() const { return T_.capacity() * sizeof(uint64_t); }
    size_t bytes_rank() const { return rank_super_.capacity() * sizeof(uint32_t); }
    size_t bytes_splits() const {
        return dx_.capacity() * sizeof(uint32_t) + dy_.capacity() * sizeof(uint32_t);
    }
    size_t bytes_used() const { return bytes_T() + bytes_rank() + bytes_splits(); }

    Stats stats() const { return stats_; }

private:
    Rect world_;

    std::vector<Point> pts_;
    std::vector<Point> scratch_;

    // 4 bits per node (INCL leaves), packed into 64-bit words.
    std::vector<uint64_t> T_;
    size_t T_bits_ = 0;

    // Split offsets for each node in BFS "all nodes stream".
    // For leaves or mask==0, dx/dy = 0.
    std::vector<uint32_t> dx_;
    std::vector<uint32_t> dy_;

    std::vector<uint32_t> rank_super_;
    Stats stats_;

    static int quadrant_of(int mx, int my, const Point& p);
    static Rect quadrant_rect(const Rect& r, int mx, int my, int q);

    // Compute median split; must guarantee progress (mx in (xmin,xmax), my in (ymin,ymax))
    void median_split(const Rect& r, uint32_t start, uint32_t count, int& mx, int& my);

    void push_bit(bool b);
    bool get_bit(size_t pos) const;
    uint8_t get_mask(size_t node_i) const;

    uint32_t rank1(size_t pos) const;
    void build_rank();
    void build_bfs();
};
