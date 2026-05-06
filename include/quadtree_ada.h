#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include <string>

#include "point.h"
#include "rect.h"
#include "rank_support.h"

class MXQuadtreeBits {
public:
    struct Params {
        int D = 0;
    };

    struct Stats {
        uint64_t points = 0;
        uint64_t N = 0;
        int D = 0;

        uint64_t T_bits = 0;
        uint64_t EX_bits = 0;
        uint64_t UL_bits = 0;
        uint64_t ULL_bits = 0;
        uint64_t ULD_bits = 0;

        uint64_t AD_bits = 0;
        uint64_t AC_bits = 0;

        uint64_t total_nodes = 0;
        uint64_t unary_to_leaf_nodes = 0;
        uint64_t internal_nodes = 0;
        uint64_t fullblock_nodes = 0;
        uint64_t leaf_nodes = 0;

        uint64_t adaptive_nodes = 0;

        uint64_t rank_T_bits = 0;
        uint64_t rank_EX_bits = 0;
        uint64_t rank_UL_bits = 0;
        uint64_t rank_AD_bits = 0;
        uint64_t rank_bits = 0;

        double bpp = 0.0;
    };

    MXQuadtreeBits() = default;

    void build(const Rect& region, const std::vector<Point>& pts, const Params& params);
    bool membership(const Point& q) const;

    const Stats& stats() const {
        return stats_;
    }

private:
    static constexpr int AC_WIDTH = 2;   // 4 adaptive options
    static constexpr int ULL_WIDTH = 6;  // supports unary length up to 63

    struct Node {
        Rect r;
        int depth = 0;
        std::vector<int> ids;
    };

    struct NodeAnalysis {
        bool is_unit_leaf = false;
        bool is_fullblock = false;
        bool expandable = false;

        bool adaptive = false;
        uint8_t split_code = 0;

        int xcut = 0;
        int ycut = 0;

        uint8_t childMask = 0;
        uint8_t nonempty_children = 0;

        std::array<std::vector<int>, 4> child_ids;
        std::array<Rect, 4> child_rects;
    };

    struct QueryState {
        Rect r;
        int depth = 0;
        uint64_t t_pos = 0;
    };

    Node make_root() const;
    Node make_child(const NodeAnalysis& a, const Node& parent, int child_idx) const;

    static inline int midpoint(int a, int b) {
        return (a + b) >> 1;
    }

    static inline bool is_unit_rect(const Rect& r) {
        return ((r.xmax - r.xmin) == 1) &&
               ((r.ymax - r.ymin) == 1);
    }

    static inline int quadrant_index(const Point& p, int xm, int ym) {
        const bool east = (p.x >= xm);
        const bool north = (p.y >= ym);

        if (!east && !north) return 0;
        if (east && !north) return 1;
        if (!east && north) return 2;
        return 3;
    }

    NodeAnalysis analyze_node(const Rect& r, int depth, const std::vector<int>& ids) const;

    // Midpoint-only analyzer used only for unary-to-leaf paths.
    NodeAnalysis analyze_node_midpoint(const Rect& r, int depth, const std::vector<int>& ids) const;

    Params params_;
    Stats stats_;
    Rect region_{0, 0, 0, 0};
    std::vector<Point> points_;
    bool root_is_fullblock_ = false;
    bool root_is_unitleaf_ = false;

    std::vector<uint64_t> T_;
    std::vector<uint64_t> EX_;
    std::vector<uint64_t> UL_;
    std::vector<uint64_t> ULL_;
    std::vector<uint64_t> ULD_;

    std::vector<uint64_t> AD_;
    std::vector<uint64_t> AC_;

    uint64_t T_len_ = 0;
    uint64_t EX_len_ = 0;
    uint64_t UL_len_ = 0;
    uint64_t ULL_len_ = 0;
    uint64_t ULD_len_ = 0;
    uint64_t AD_len_ = 0;
    uint64_t AC_len_ = 0;

    static void push_bits(std::vector<uint64_t>& dst, uint64_t& bit_len, uint64_t value, int width);
    static uint64_t get_bit(const std::vector<uint64_t>& src, uint64_t bit_pos);
    static uint64_t read_bits(const std::vector<uint64_t>& src, uint64_t bit_pos, int width);

    Rect child_rect(const Rect& r, int child_idx) const;

    void split_rects_for_code(
        const Rect& r,
        bool adaptive,
        uint8_t split_code,
        std::array<Rect, 4>& out_rects,
        int& xcut,
        int& ycut
    ) const;

    Rect child_rect_with_split(
        const Rect& r,
        int child_idx,
        bool adaptive,
        uint8_t split_code
    ) const;

    uint8_t choose_adaptive_split(
        const Rect& r,
        int depth,
        const std::vector<int>& ids,
        bool& use_adaptive
    ) const;

    struct UnarySkipResult {
        uint8_t L = 0;
        std::vector<uint8_t> dirs;
        Node endpoint;
        NodeAnalysis endpoint_analysis;
    };

    // Follows unary path using midpoint splits only.
    UnarySkipResult follow_unary_chain_midpoint(const Node& start) const;

    void build_bfs_contracted();

    uint64_t rank1_T(uint64_t bit_pos) const;
    uint64_t rank1_EX(uint64_t bit_pos) const;
    uint64_t rank0_EX(uint64_t bit_pos) const;
    uint64_t rank1_UL(uint64_t bit_pos) const;
    uint64_t rank1_AD(uint64_t bit_pos) const;

    RankSupport64 rank_T_;
    RankSupport64 rank_EX_;
    RankSupport64 rank_UL_;
    RankSupport64 rank_AD_;
};