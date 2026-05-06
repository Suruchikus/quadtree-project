#include "quadtree_test.h"

#include <algorithm>

void MXQuadtreeBits::push_bits(
    std::vector<uint64_t>& dst,
    uint64_t& bit_len,
    uint64_t v,
    int width
) {
    if (width <= 0) return;

    while (width > 0) {
        const uint64_t word_idx = bit_len >> 6;
        const uint64_t bit_off = bit_len & 63ULL;

        if (word_idx >= dst.size()) {
            dst.push_back(0ULL);
        }

        const int space = 64 - (int)bit_off;
        const int take = std::min(width, space);

        const uint64_t mask = (take == 64)
            ? ~0ULL
            : ((1ULL << take) - 1ULL);

        const uint64_t chunk = (v & mask) << bit_off;
        dst[word_idx] |= chunk;

        v >>= take;
        width -= take;
        bit_len += (uint64_t)take;
    }
}

uint64_t MXQuadtreeBits::get_bit(const std::vector<uint64_t>& src, uint64_t bit_pos) {
    const uint64_t w = bit_pos >> 6;
    const uint64_t b = bit_pos & 63ULL;

    if (w >= src.size()) return 0ULL;

    return (src[w] >> b) & 1ULL;
}

void MXQuadtreeBits::split_rects_for_code(
    const Rect& r,
    bool adaptive,
    uint8_t split_code,
    std::array<Rect, 4>& out_rects,
    int& xcut,
    int& ycut
) const {
    const int w = r.xmax - r.xmin;
    const int h = r.ymax - r.ymin;

    if (!adaptive) {
        xcut = midpoint(r.xmin, r.xmax);
        ycut = midpoint(r.ymin, r.ymax);
    } else {
        const bool x_high = ((split_code >> 1) & 1u) != 0;
        const bool y_high = (split_code & 1u) != 0;

        xcut = r.xmin + (x_high ? (3 * w) / 4 : w / 4);
        ycut = r.ymin + (y_high ? (3 * h) / 4 : h / 4);

        if (xcut <= r.xmin) xcut = r.xmin + 1;
        if (xcut >= r.xmax) xcut = r.xmax - 1;

        if (ycut <= r.ymin) ycut = r.ymin + 1;
        if (ycut >= r.ymax) ycut = r.ymax - 1;
    }

    out_rects[0] = Rect{r.xmin, r.ymin, xcut, ycut};
    out_rects[1] = Rect{xcut, r.ymin, r.xmax, ycut};
    out_rects[2] = Rect{r.xmin, ycut, xcut, r.ymax};
    out_rects[3] = Rect{xcut, ycut, r.xmax, r.ymax};
}

uint8_t MXQuadtreeBits::choose_adaptive_split(
    const Rect& r,
    int depth,
    const std::vector<int>& ids,
    bool& use_adaptive
) const {
    use_adaptive = false;

    if (depth >= params_.D) return 0;

    const int w = r.xmax - r.xmin;
    const int h = r.ymax - r.ymin;

    if (w < 4 || h < 4) return 0;

    const int mx = midpoint(r.xmin, r.xmax);
    const int my = midpoint(r.ymin, r.ymax);

    uint8_t mid_mask = 0;

    for (int id : ids) {
        const Point& p = points_[id];
        const int q = quadrant_index(p, mx, my);
        mid_mask |= (uint8_t)(1u << q);
    }

    int mid_nonempty = 0;
    for (int i = 0; i < 4; i++) {
        if ((mid_mask >> i) & 1u) {
            mid_nonempty++;
        }
    }

    uint8_t best_code = 0;
    int best_nonempty = mid_nonempty;
    uint64_t best_largest_empty_area = 0;

    for (uint8_t code = 0; code < 4; code++) {
        std::array<Rect, 4> rects;
        int xcut = 0;
        int ycut = 0;

        split_rects_for_code(r, true, code, rects, xcut, ycut);

        std::array<uint64_t, 4> counts{0, 0, 0, 0};

        for (int id : ids) {
            const Point& p = points_[id];
            const int q = quadrant_index(p, xcut, ycut);
            counts[q]++;
        }

        int nonempty = 0;
        uint64_t largest_empty_area = 0;

        for (int i = 0; i < 4; i++) {
            if (counts[i] > 0) {
                nonempty++;
            } else {
                const uint64_t rw = (uint64_t)(rects[i].xmax - rects[i].xmin);
                const uint64_t rh = (uint64_t)(rects[i].ymax - rects[i].ymin);
                const uint64_t area = rw * rh;
                largest_empty_area = std::max(largest_empty_area, area);
            }
        }

        if (nonempty < best_nonempty ||
            (nonempty == best_nonempty &&
             nonempty < mid_nonempty &&
             largest_empty_area > best_largest_empty_area)) {
            best_nonempty = nonempty;
            best_largest_empty_area = largest_empty_area;
            best_code = code;
        }
    }

    if (best_nonempty < mid_nonempty) {
        use_adaptive = true;
        return best_code;
    }

    return 0;
}

MXQuadtreeBits::NodeAnalysis
MXQuadtreeBits::analyze_node(
    const Rect& r,
    int depth,
    const std::vector<int>& ids
) const {
    NodeAnalysis a;

    const uint64_t w = (uint64_t)(r.xmax - r.xmin);
    const uint64_t h = (uint64_t)(r.ymax - r.ymin);
    const uint64_t area = w * h;

    a.is_unit_leaf = depth >= params_.D;
    a.is_fullblock = !a.is_unit_leaf && (uint64_t)ids.size() == area;
    a.expandable = !a.is_unit_leaf && !a.is_fullblock;

    if (!a.expandable) {
        return a;
    }

    bool use_adaptive = false;
    const uint8_t split_code =
        choose_adaptive_split(r, depth, ids, use_adaptive);

    a.adaptive = use_adaptive;
    a.split_code = split_code;

    split_rects_for_code(
        r,
        a.adaptive,
        a.split_code,
        a.child_rects,
        a.xcut,
        a.ycut
    );

    for (int i = 0; i < 4; i++) {
        a.child_ids[i].clear();
    }

    for (int id : ids) {
        const Point& p = points_[id];
        const int q = quadrant_index(p, a.xcut, a.ycut);
        a.child_ids[q].push_back(id);
    }

    uint8_t mask = 0;
    uint8_t count = 0;

    for (int i = 0; i < 4; i++) {
        if (!a.child_ids[i].empty()) {
            mask |= (uint8_t)(1u << i);
            count++;
        }
    }

    a.childMask = mask;
    a.nonempty_children = count;

    return a;
}

MXQuadtreeBits::Node MXQuadtreeBits::make_root() const {
    Node root;
    root.r = region_;
    root.depth = 0;
    root.ids.reserve(points_.size());

    for (int i = 0; i < (int)points_.size(); i++) {
        root.ids.push_back(i);
    }

    return root;
}

MXQuadtreeBits::Node
MXQuadtreeBits::make_child(
    const NodeAnalysis& a,
    const Node& parent,
    int child_idx
) const {
    Node child;
    child.r = a.child_rects[child_idx];
    child.depth = parent.depth + 1;
    child.ids = a.child_ids[child_idx];
    return child;
}

void MXQuadtreeBits::build(
    const Rect& region,
    const std::vector<Point>& pts,
    const Params& params
) {
    params_ = params;
    region_ = region;
    points_ = pts;

    root_is_fullblock_ = false;
    root_is_unitleaf_ = false;

    stats_ = Stats{};
    stats_.points = points_.size();
    stats_.N = (uint64_t)(region_.xmax - region_.xmin);
    stats_.D = params_.D;

    build_bfs();
}

void MXQuadtreeBits::build_bfs() {
    T_.clear();
    EX_.clear();
    AD_.clear();
    AC_.clear();

    T_len_ = 0;
    EX_len_ = 0;
    AD_len_ = 0;
    AC_len_ = 0;

    if (points_.empty()) {
        return;
    }

    Node root = make_root();
    NodeAnalysis root_a = analyze_node(root.r, root.depth, root.ids);

    if (!root_a.expandable) {
        root_is_fullblock_ = root_a.is_fullblock;
        root_is_unitleaf_ = root_a.is_unit_leaf;

        stats_.total_nodes = 1;

        if (root_a.is_fullblock) stats_.fullblock_nodes = 1;
        if (root_a.is_unit_leaf) stats_.leaf_nodes = 1;

        return;
    }

    std::vector<Node> curr;
    curr.push_back(std::move(root));

    while (!curr.empty()) {
        std::vector<Node> next;

        for (const Node& node : curr) {
            NodeAnalysis a = analyze_node(node.r, node.depth, node.ids);

            push_bits(AD_, AD_len_, a.adaptive ? 1ULL : 0ULL, 1);
            if (a.adaptive) {
                push_bits(AC_, AC_len_, (uint64_t)(a.split_code & 3u), AC_WIDTH);
                stats_.adaptive_nodes++;
            }

            push_bits(T_, T_len_, (uint64_t)a.childMask, 4);
            stats_.internal_nodes++;

            for (int i = 0; i < 4; i++) {
                if (((a.childMask >> i) & 1u) == 0u) {
                    continue;
                }

                Node child = make_child(a, node, i);

                if (child.depth == params_.D) {
                    stats_.leaf_nodes++;
                    continue;
                }

                NodeAnalysis ca = analyze_node(child.r, child.depth, child.ids);

                if (ca.is_fullblock) {
                    push_bits(EX_, EX_len_, 1ULL, 1);
                    stats_.fullblock_nodes++;
                } else {
                    push_bits(EX_, EX_len_, 0ULL, 1);
                    next.push_back(std::move(child));
                }
            }
        }

        curr.swap(next);
    }

    rank_T_.build(T_, T_len_, 8);
    rank_EX_.build(EX_, EX_len_, 8);
    rank_AD_.build(AD_, AD_len_, 8);

    stats_.rank_T_bits = rank_T_.space_in_bits();
    stats_.rank_EX_bits = rank_EX_.space_in_bits();
    stats_.rank_AD_bits = rank_AD_.space_in_bits();

    stats_.rank_bits =
        stats_.rank_T_bits +
        stats_.rank_EX_bits +
        stats_.rank_AD_bits;

    uint64_t ones_T = 0;
    for (uint64_t i = 0; i < T_len_; i++) {
        ones_T += get_bit(T_, i);
    }

    stats_.total_nodes = 1 + ones_T;

    stats_.T_bits = T_len_;
    stats_.EX_bits = EX_len_;
    stats_.AD_bits = AD_len_;
    stats_.AC_bits = AC_len_;

    const uint64_t total_bits =
        T_len_ +
        EX_len_ +
        AD_len_ +
        AC_len_ +
        stats_.rank_bits;

    stats_.bpp = stats_.points == 0
        ? 0.0
        : (double)total_bits / (double)stats_.points;
}