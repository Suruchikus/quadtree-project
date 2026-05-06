#include "quadtree_ada.h"
#include <iostream>
#include <limits>

// ------------------------------------------------------------
// Rank helpers
// ------------------------------------------------------------

uint64_t MXQuadtreeBits::rank1_T(uint64_t bit_pos) const {
    return rank_T_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank1_EX(uint64_t bit_pos) const {
    return rank_EX_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank1_UL(uint64_t bit_pos) const {
    return rank_UL_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank0_EX(uint64_t bit_pos) const {
    return rank_EX_.rank0_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank1_AD(uint64_t bit_pos) const {
    return rank_AD_.rank1_before(bit_pos);
}

// ------------------------------------------------------------
// Adaptive split helpers
// ------------------------------------------------------------

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

Rect MXQuadtreeBits::child_rect_with_split(
    const Rect& r,
    int child_idx,
    bool adaptive,
    uint8_t split_code
) const {
    std::array<Rect, 4> rects;
    int xcut = 0;
    int ycut = 0;

    split_rects_for_code(r, adaptive, split_code, rects, xcut, ycut);
    return rects[child_idx];
}

Rect MXQuadtreeBits::child_rect(const Rect& r, int child_idx) const {
    return child_rect_with_split(r, child_idx, false, 0);
}

uint8_t MXQuadtreeBits::choose_adaptive_split(
    const Rect& r,
    int depth,
    const std::vector<int>& ids,
    bool& use_adaptive
) const {
    use_adaptive = false;

    if (is_unit_rect(r)) return 0;

    const int w = r.xmax - r.xmin;
    const int h = r.ymax - r.ymin;

    if (w < 4 || h < 4) return 0;

    const int mx = midpoint(r.xmin, r.xmax);
    const int my = midpoint(r.ymin, r.ymax);

    uint8_t mid_mask = 0;
    for (int id : ids) {
        const Point& p = points_[id];
        const int qi = quadrant_index(p, mx, my);
        mid_mask |= (uint8_t)(1u << qi);
    }

    int mid_nonempty = 0;
    for (int i = 0; i < 4; i++) {
        if ((mid_mask >> i) & 1u) mid_nonempty++;
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
            const int qi = quadrant_index(p, xcut, ycut);
            counts[qi]++;
        }

        int nonempty = 0;
        uint64_t largest_empty_area = 0;

        for (int i = 0; i < 4; i++) {
            if (counts[i] > 0) {
                nonempty++;
            } else {
                const uint64_t rw = (uint64_t)(rects[i].xmax - rects[i].xmin);
                const uint64_t rh = (uint64_t)(rects[i].ymax - rects[i].ymin);
                largest_empty_area = std::max(largest_empty_area, rw * rh);
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

// ------------------------------------------------------------
// Node analyzers
// ------------------------------------------------------------

MXQuadtreeBits::NodeAnalysis
MXQuadtreeBits::analyze_node_midpoint(
    const Rect& r,
    int depth,
    const std::vector<int>& ids
) const {
    NodeAnalysis a;

    const uint64_t w = (uint64_t)(r.xmax - r.xmin);
    const uint64_t h = (uint64_t)(r.ymax - r.ymin);
    const uint64_t area = w * h;

    a.is_unit_leaf = is_unit_rect(r);
    a.is_fullblock = (!a.is_unit_leaf && (uint64_t)ids.size() == area);
    a.expandable = (!a.is_unit_leaf && !a.is_fullblock);

    if (!a.expandable) {
        return a;
    }

    a.adaptive = false;
    a.split_code = 0;

    split_rects_for_code(
        r,
        false,
        0,
        a.child_rects,
        a.xcut,
        a.ycut
    );

    for (int i = 0; i < 4; i++) {
        a.child_ids[i].clear();
    }

    for (int id : ids) {
        const Point& p = points_[id];
        const int qi = quadrant_index(p, a.xcut, a.ycut);
        a.child_ids[qi].push_back(id);
    }

    uint8_t mask = 0;
    uint8_t cnt = 0;

    for (int i = 0; i < 4; i++) {
        if (!a.child_ids[i].empty()) {
            mask |= (uint8_t)(1u << i);
            cnt++;
        }
    }

    a.childMask = mask;
    a.nonempty_children = cnt;

    return a;
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

    a.is_unit_leaf = is_unit_rect(r);
    a.is_fullblock = (!a.is_unit_leaf && (uint64_t)ids.size() == area);
    a.expandable = (!a.is_unit_leaf && !a.is_fullblock);

    if (!a.expandable) {
        return a;
    }

    bool use_adaptive = false;
    uint8_t split_code = choose_adaptive_split(r, depth, ids, use_adaptive);

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
        const int qi = quadrant_index(p, a.xcut, a.ycut);
        a.child_ids[qi].push_back(id);
    }

    uint8_t mask = 0;
    uint8_t cnt = 0;

    for (int i = 0; i < 4; i++) {
        if (!a.child_ids[i].empty()) {
            mask |= (uint8_t)(1u << i);
            cnt++;
        }
    }

    a.childMask = mask;
    a.nonempty_children = cnt;

    return a;
}

// ------------------------------------------------------------
// Basic helpers
// ------------------------------------------------------------

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

    build_bfs_contracted();
}

MXQuadtreeBits::Node
MXQuadtreeBits::make_root() const {
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

// ------------------------------------------------------------
// Bit helpers
// ------------------------------------------------------------

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
        const int take = (width < space) ? width : space;

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

uint64_t MXQuadtreeBits::get_bit(const std::vector<uint64_t>& src, uint64_t i) {
    const uint64_t w = i >> 6;
    const uint64_t b = i & 63ULL;

    if (w >= src.size()) return 0;

    return (src[w] >> b) & 1ULL;
}

uint64_t MXQuadtreeBits::read_bits(
    const std::vector<uint64_t>& src,
    uint64_t bit_pos,
    int width
) {
    if (width <= 0) return 0ULL;

    uint64_t result = 0;
    int written = 0;

    while (width > 0) {
        const uint64_t word_idx = bit_pos >> 6;
        const uint64_t bit_off = bit_pos & 63ULL;

        if (word_idx >= src.size()) break;

        const int available = 64 - (int)bit_off;
        const int take = (width < available) ? width : available;

        const uint64_t mask = (take == 64)
            ? ~0ULL
            : ((1ULL << take) - 1ULL);

        const uint64_t chunk = (src[word_idx] >> bit_off) & mask;

        result |= (chunk << written);

        bit_pos += (uint64_t)take;
        width -= take;
        written += take;
    }

    return result;
}

// ------------------------------------------------------------
// Build
// ------------------------------------------------------------

void MXQuadtreeBits::build_bfs_contracted() {
    T_.clear();
    EX_.clear();
    UL_.clear();
    ULL_.clear();
    ULD_.clear();
    AD_.clear();
    AC_.clear();

    T_len_ = 0;
    EX_len_ = 0;
    UL_len_ = 0;
    ULL_len_ = 0;
    ULD_len_ = 0;
    AD_len_ = 0;
    AC_len_ = 0;

    if (points_.empty()) return;

    Node root = make_root();
    NodeAnalysis ra = analyze_node(root.r, root.depth, root.ids);

    if (!ra.expandable) {
        root_is_fullblock_ = ra.is_fullblock;
        root_is_unitleaf_ = ra.is_unit_leaf;

        stats_.total_nodes = 1;

        if (ra.is_fullblock) stats_.fullblock_nodes = 1;
        if (ra.is_unit_leaf) stats_.leaf_nodes = 1;

        return;
    }

    std::vector<Node> curr;
    curr.push_back(std::move(root));

    while (!curr.empty()) {
        std::vector<Node> next;

        for (const Node& node : curr) {
            NodeAnalysis a = analyze_node(node.r, node.depth, node.ids);

            push_bits(T_, T_len_, (uint64_t)a.childMask, 4);
            stats_.internal_nodes++;

            for (int i = 0; i < 4; i++) {
                if (((a.childMask >> i) & 1u) == 0u) {
                    continue;
                }

                Node child = make_child(a, node, i);
                NodeAnalysis ca = analyze_node(child.r, child.depth, child.ids);

                // Terminal unit leaf: EX=1, UL=0.
                if (ca.is_unit_leaf) {
                    push_bits(EX_, EX_len_, 1ULL, 1);
                    push_bits(UL_, UL_len_, 0ULL, 1);
                    stats_.leaf_nodes++;
                    continue;
                }

                // Terminal fullblock: EX=1, UL=0.
                if (ca.is_fullblock) {
                    push_bits(EX_, EX_len_, 1ULL, 1);
                    push_bits(UL_, UL_len_, 0ULL, 1);
                    stats_.fullblock_nodes++;
                    continue;
                }

                // Try midpoint-only unary-to-leaf.
                NodeAnalysis ca_mid =
                    analyze_node_midpoint(child.r, child.depth, child.ids);

                const bool is_mid_unary =
                    ca_mid.expandable && ca_mid.nonempty_children == 1;

                if (is_mid_unary) {
                    UnarySkipResult sk =
                        follow_unary_chain_midpoint(child);

                    if (sk.endpoint_analysis.is_unit_leaf) {
                        push_bits(EX_, EX_len_, 1ULL, 1);
                        push_bits(UL_, UL_len_, 1ULL, 1);

                        push_bits(ULL_, ULL_len_, (uint64_t)sk.L, ULL_WIDTH);

                        for (uint8_t dir : sk.dirs) {
                            push_bits(ULD_, ULD_len_, (uint64_t)(dir & 3u), 2);
                        }

                        stats_.unary_to_leaf_nodes++;
                        continue;
                    }
                }

                // Otherwise explicit child: EX=0, then AD/AC for this child node.
                push_bits(EX_, EX_len_, 0ULL, 1);

                push_bits(AD_, AD_len_, ca.adaptive ? 1ULL : 0ULL, 1);

                if (ca.adaptive) {
                    push_bits(AC_, AC_len_, (uint64_t)(ca.split_code & 3u), AC_WIDTH);
                    stats_.adaptive_nodes++;
                }

                next.push_back(std::move(child));
            }
        }

        curr.swap(next);
    }

    rank_T_.build(T_, T_len_, 8);
    rank_EX_.build(EX_, EX_len_, 8);
    rank_UL_.build(UL_, UL_len_, 8);
    rank_AD_.build(AD_, AD_len_, 8);

    stats_.rank_T_bits  = rank_T_.space_in_bits();
    stats_.rank_EX_bits = rank_EX_.space_in_bits();
    stats_.rank_UL_bits = rank_UL_.space_in_bits();
    stats_.rank_AD_bits = rank_AD_.space_in_bits();

    stats_.rank_bits =
        stats_.rank_T_bits +
        stats_.rank_EX_bits +
        stats_.rank_UL_bits +
        stats_.rank_AD_bits;

    uint64_t ones_T = 0;
    for (uint64_t i = 0; i < T_len_; i++) {
        ones_T += get_bit(T_, i);
    }

    stats_.total_nodes = 1 + ones_T;

    stats_.T_bits = T_len_;
    stats_.EX_bits = EX_len_;
    stats_.UL_bits = UL_len_;
    stats_.ULL_bits = ULL_len_;
    stats_.ULD_bits = ULD_len_;
    stats_.AD_bits = AD_len_;
    stats_.AC_bits = AC_len_;

    const uint64_t total_bits =
        T_len_ +
        EX_len_ +
        UL_len_ +
        ULL_len_ +
        ULD_len_ +
        AD_len_ +
        AC_len_ +
        stats_.rank_bits;

    stats_.bpp = stats_.points == 0
        ? 0.0
        : (double)total_bits / (double)stats_.points;
}

// ------------------------------------------------------------
// Unary skipping: midpoint only
// ------------------------------------------------------------

MXQuadtreeBits::UnarySkipResult
MXQuadtreeBits::follow_unary_chain_midpoint(const Node& start) const {
    UnarySkipResult res;

    Node cur = start;
    NodeAnalysis ca = analyze_node_midpoint(cur.r, cur.depth, cur.ids);

    if (!(ca.expandable && ca.nonempty_children == 1)) {
        res.L = 0;
        res.endpoint = cur;
        res.endpoint_analysis = ca;
        return res;
    }

    while (true) {
        int only_child = -1;

        for (int i = 0; i < 4; i++) {
            if (((ca.childMask >> i) & 1u) != 0u) {
                only_child = i;
                break;
            }
        }

        if (only_child < 0) {
            res.endpoint = cur;
            res.endpoint_analysis = ca;
            return res;
        }

        res.dirs.push_back((uint8_t)only_child);
        res.L++;

        Node nxt = make_child(ca, cur, only_child);
        NodeAnalysis na =
            analyze_node_midpoint(nxt.r, nxt.depth, nxt.ids);

        if (!(na.expandable && na.nonempty_children == 1)) {
            res.endpoint = std::move(nxt);
            res.endpoint_analysis = na;
            return res;
        }

        cur = std::move(nxt);
        ca = na;

        if (res.L >= ((1u << ULL_WIDTH) - 1u)) {
            res.endpoint = cur;
            res.endpoint_analysis = ca;
            return res;
        }
    }
}

// ------------------------------------------------------------
// Membership
// ------------------------------------------------------------

bool MXQuadtreeBits::membership(const Point& q) const {
    if (q.x < region_.xmin || q.x >= region_.xmax ||
        q.y < region_.ymin || q.y >= region_.ymax) {
        return false;
    }

    if (points_.empty()) return false;
    if (root_is_fullblock_) return true;
    if (root_is_unitleaf_) return true;

    Rect curr_r = region_;
    uint64_t curr_t_pos = 0;

    bool curr_adaptive = false;
    uint8_t curr_split_code = 0;

    while (true) {
        const uint64_t childMask = read_bits(T_, curr_t_pos, 4);

        std::array<Rect, 4> rects;
        int xcut = 0;
        int ycut = 0;

        split_rects_for_code(
            curr_r,
            curr_adaptive,
            curr_split_code,
            rects,
            xcut,
            ycut
        );

        const int child_idx = quadrant_index(q, xcut, ycut);

        if (((childMask >> child_idx) & 1ULL) == 0ULL) {
            return false;
        }

        Rect child_r = rects[child_idx];

        const uint64_t child_t_bit_pos =
            curr_t_pos + (uint64_t)child_idx;

        const uint64_t ex_pos = rank1_T(child_t_bit_pos);
        const uint64_t ex_bit = get_bit(EX_, ex_pos);

        if (ex_bit == 0ULL) {
            const uint64_t explicit_idx = rank0_EX(ex_pos);

            const uint64_t ad_bit = get_bit(AD_, explicit_idx);

            bool next_adaptive = false;
            uint8_t next_split_code = 0;

            if (ad_bit == 1ULL) {
                const uint64_t ac_idx = rank1_AD(explicit_idx);
                const uint64_t ac_pos = ac_idx * (uint64_t)AC_WIDTH;

                next_adaptive = true;
                next_split_code =
                    (uint8_t)read_bits(AC_, ac_pos, AC_WIDTH);
            }

            const uint64_t next_node_idx = 1ULL + explicit_idx;

            curr_t_pos = 4ULL * next_node_idx;
            curr_r = child_r;
            curr_adaptive = next_adaptive;
            curr_split_code = next_split_code;

            continue;
        }

        const uint64_t ul_pos = rank1_EX(ex_pos);
        const uint64_t ul_bit = get_bit(UL_, ul_pos);

        if (ul_bit == 0ULL) {
            return true;
        }

        const uint64_t unary_idx = rank1_UL(ul_pos);
        const uint64_t len_pos = unary_idx * (uint64_t)ULL_WIDTH;
        const uint64_t L = read_bits(ULL_, len_pos, ULL_WIDTH);

        uint64_t dir_start = 0;

        for (uint64_t j = 0; j < unary_idx; j++) {
            const uint64_t prev_L =
                read_bits(ULL_, j * (uint64_t)ULL_WIDTH, ULL_WIDTH);

            dir_start += prev_L * 2ULL;
        }

        Rect r = child_r;

        for (uint64_t step = 0; step < L; step++) {
            const uint64_t got =
                read_bits(ULD_, dir_start + 2ULL * step, 2);

            const int xm = midpoint(r.xmin, r.xmax);
            const int ym = midpoint(r.ymin, r.ymax);

            const uint64_t want =
                (uint64_t)quadrant_index(q, xm, ym);

            if (want != got) {
                return false;
            }

            r = child_rect(r, (int)got);
        }

        return true;
    }

    return false;
}