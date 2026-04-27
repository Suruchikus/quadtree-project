#include "quadtree_v2.h"
#include <iostream>


// Rank/Select Support functions
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



// Helper functions for the main build
MXQuadtreeBits::NodeAnalysis
MXQuadtreeBits::analyze_node(const Rect& r, int depth, const std::vector<int>& ids) const {
    NodeAnalysis a;

    const uint64_t w = (uint64_t)(r.xmax - r.xmin);
    const uint64_t h = (uint64_t)(r.ymax - r.ymin);
    const uint64_t area = w * h;

    a.is_unit_leaf = (depth >= params_.D);
    a.is_fullblock = (!a.is_unit_leaf && (uint64_t)ids.size() == area);
    a.expandable = (!a.is_unit_leaf && !a.is_fullblock);

    const int xm = midpoint(r.xmin, r.xmax);
    const int ym = midpoint(r.ymin, r.ymax);

    a.child_rects[0] = Rect{r.xmin, r.ymin, xm, ym};
    a.child_rects[1] = Rect{xm, r.ymin, r.xmax, ym};
    a.child_rects[2] = Rect{r.xmin, ym, xm, r.ymax};
    a.child_rects[3] = Rect{xm, ym, r.xmax, r.ymax};

    for (int i = 0; i < 4; i++) {
        a.child_ids[i].clear();
    }

    for (int id : ids) {
        const Point& p = points_[id];
        const int qi = quadrant_index(p, xm, ym);
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

void MXQuadtreeBits::build(const Rect& region, const std::vector<Point>& pts, const Params& params) {
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
MXQuadtreeBits::make_child(const NodeAnalysis& a, const Node& parent, int child_idx) const {
    Node child;
    child.r = a.child_rects[child_idx];
    child.depth = parent.depth + 1;
    child.ids = a.child_ids[child_idx];
    return child;
}

Rect MXQuadtreeBits::child_rect(const Rect& r, int child_idx) const {
    const int xm = midpoint(r.xmin, r.xmax);
    const int ym = midpoint(r.ymin, r.ymax);

    if (child_idx == 0) return Rect{r.xmin, r.ymin, xm, ym};
    if (child_idx == 1) return Rect{xm, r.ymin, r.xmax, ym};
    if (child_idx == 2) return Rect{r.xmin, ym, xm, r.ymax};
    return Rect{xm, ym, r.xmax, r.ymax};
}

// Helpers of bit quadtree
void MXQuadtreeBits::push_bits(std::vector<uint64_t>& dst, uint64_t& bit_len, uint64_t v, int width) {
    if (width <= 0) return;

    while (width > 0) {
        const uint64_t word_idx = bit_len >> 6;
        const uint64_t bit_off = bit_len & 63ULL;
        if (word_idx >= dst.size()) dst.push_back(0ULL);

        const int space = 64 - (int)bit_off;
        const int take = (width < space) ? width : space;

        const uint64_t mask = (take == 64) ? ~0ULL : ((1ULL << take) - 1ULL);
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

uint64_t MXQuadtreeBits::read_bits(const std::vector<uint64_t>& src, uint64_t bit_pos, int width) {
    if (width <= 0) return 0ULL;

    uint64_t result = 0;
    int written = 0;

    while (width > 0) {
        const uint64_t word_idx = bit_pos >> 6;
        const uint64_t bit_off = bit_pos & 63ULL;

        if (word_idx >= src.size()) break;

        const int available = 64 - (int)bit_off;
        const int take = (width < available) ? width : available;

        const uint64_t mask = (take == 64) ? ~0ULL : ((1ULL << take) - 1ULL);
        const uint64_t chunk = (src[word_idx] >> bit_off) & mask;

        result |= (chunk << written);

        bit_pos += (uint64_t)take;
        width -= take;
        written += take;
    }

    return result;
}


void MXQuadtreeBits::build_bfs_contracted() {
    T_.clear();
    EX_.clear();
    UL_.clear();
    ULD_.clear();

    T_len_ = 0;
    EX_len_ = 0;
    UL_len_ = 0;
    ULD_len_ = 0;

    uleaf_count_by_depth_.assign(params_.D + 1, 0);

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

            // Store only internal explicit nodes in T.
            // Since depth-D leaves are never pushed, every node here has depth < D.
            push_bits(T_, T_len_, (uint64_t)a.childMask, 4);
            stats_.internal_nodes++;

            for (int i = 0; i < 4; i++) {
                if (((a.childMask >> i) & 1u) == 0u) continue;

                Node child = make_child(a, node, i);
                NodeAnalysis ca = analyze_node(child.r, child.depth, child.ids);

                // Case 1: direct unit leaf.
                // The 1-bit in T is enough.
                // No EX, no UL, no child pushed.
                if (child.depth == params_.D) {
                    stats_.leaf_nodes++;
                    continue;
                }

                const bool is_unary =
                    ca.expandable && ca.nonempty_children == 1;

                if (is_unary) {
                    UnarySkipResult sk = follow_unary_chain(child);

                    if (sk.endpoint_analysis.is_unit_leaf) {
                        // Case 2: unary-to-leaf.
                        // EX=1 says terminal shortcut.
                        // UL=1 says unary-to-leaf.
                        push_bits(EX_, EX_len_, 1ULL, 1);
                        push_bits(UL_, UL_len_, 1ULL, 1);

                        for (uint8_t dir : sk.dirs) {
                            push_bits(ULD_, ULD_len_, (uint64_t)(dir & 3u), 2);
                        }

                        uleaf_count_by_depth_[child.depth]++;

                        stats_.unary_to_leaf_nodes++;
                    } else {
                        // Case 3: unary-to-mixed or unary-to-fullblock.
                        // Do not compress. Store immediate child as explicit.
                        push_bits(EX_, EX_len_, 0ULL, 1);

                        next.push_back(std::move(child));
                        //stats_.internal_nodes++;
                    }
                } else {
                    if (ca.is_fullblock) {
                        // Case 4: fullblock terminal before depth D.
                        // EX=1 terminal, UL=0 fullblock/direct stop.
                        push_bits(EX_, EX_len_, 1ULL, 1);
                        push_bits(UL_, UL_len_, 0ULL, 1);

                        stats_.fullblock_nodes++;
                    } else {
                        // Case 5: normal explicit child.
                        push_bits(EX_, EX_len_, 0ULL, 1);

                        next.push_back(std::move(child));
                        //stats_.internal_nodes++;
                    }
                }
            }
        }

        curr.swap(next);
    }

    rank_T_.build(T_, T_len_, 8);
    rank_EX_.build(EX_, EX_len_, 8);
    rank_UL_.build(UL_, UL_len_, 8);

    stats_.rank_T_bits  = rank_T_.space_in_bits();
    stats_.rank_EX_bits = rank_EX_.space_in_bits();
    stats_.rank_UL_bits = rank_UL_.space_in_bits();

    stats_.rank_bits =
        stats_.rank_T_bits +
        stats_.rank_EX_bits +
        stats_.rank_UL_bits;

    uint64_t ones_T = 0;
    for (uint64_t i = 0; i < T_len_; i++) {
        ones_T += get_bit(T_, i);
    }

    stats_.total_nodes = 1 + ones_T;

    stats_.T_bits = T_len_;
    stats_.EX_bits = EX_len_;
    stats_.UL_bits = UL_len_;
    stats_.ULD_bits = ULD_len_;

    const uint64_t uleaf_count_bits =
    (uint64_t)uleaf_count_by_depth_.size() * 64ULL;

    const uint64_t total_bits =
        T_len_ +
        EX_len_ +
        UL_len_ +
        ULD_len_ +
        uleaf_count_bits +
        stats_.rank_bits;

    stats_.bpp = stats_.points == 0
        ? 0.0
        : (double)total_bits / (double)stats_.points;
}

// Unary Skipping
MXQuadtreeBits::UnarySkipResult
MXQuadtreeBits::follow_unary_chain(const Node& start) const {
    UnarySkipResult res;
    Node cur = start;
    NodeAnalysis ca = analyze_node(cur.r, cur.depth, cur.ids);

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
        NodeAnalysis na = analyze_node(nxt.r, nxt.depth, nxt.ids);

        if (!(na.expandable && na.nonempty_children == 1)) {
            res.endpoint = std::move(nxt);
            res.endpoint_analysis = na;
            return res;
        }

        cur = std::move(nxt);
        ca = na;

        if (res.L >= 31) {
            res.endpoint = cur;
            res.endpoint_analysis = ca;
            return res;
        }
    }
}

// Membership Queries
bool MXQuadtreeBits::membership(const Point& q) const {
    if (q.x < region_.xmin || q.x >= region_.xmax ||
        q.y < region_.ymin || q.y >= region_.ymax) {
        return false;
    }

    if (points_.empty()) {
        return false;
    }

    if (root_is_fullblock_) {
        return true;
    }

    if (root_is_unitleaf_) {
        return true;
    }

    Rect curr_r = region_;
    int curr_depth = 0;
    uint64_t curr_t_pos = 0;

    while (true) {
        const uint64_t childMask = read_bits(T_, curr_t_pos, 4);

        const int xm = midpoint(curr_r.xmin, curr_r.xmax);
        const int ym = midpoint(curr_r.ymin, curr_r.ymax);
        const int child_idx = quadrant_index(q, xm, ym);

        if (((childMask >> child_idx) & 1ULL) == 0ULL) {
            return false;
        }

        Rect child_r = child_rect(curr_r, child_idx);
        const int child_depth = curr_depth + 1;

        // Direct unit leaf: no EX/UL exists for this edge.
        if (child_depth == params_.D) {
            return true;
        }

        const uint64_t child_t_bit_pos =
            curr_t_pos + (uint64_t)child_idx;

        const uint64_t ex_pos = rank1_T(child_t_bit_pos);
        const uint64_t ex_bit = get_bit(EX_, ex_pos);

        if (ex_bit == 0ULL) {
            // Normal explicit child.
            const uint64_t next_node_idx = 1ULL + rank0_EX(ex_pos);

            curr_t_pos = 4ULL * next_node_idx;
            curr_r = child_r;
            curr_depth = child_depth;
            continue;
        }

        // EX=1 means terminal shortcut.
        const uint64_t ul_pos = rank1_EX(ex_pos);
        const uint64_t ul_bit = get_bit(UL_, ul_pos);

        if (ul_bit == 0ULL) {
            // Fullblock terminal.
            return true;
        }

        // UL=1 means unary-to-leaf.
        const uint64_t total_uleaf_before_pos = rank1_UL(ul_pos);

        uint64_t dir_start = 0;
        uint64_t count_before_level = 0;

        for (int d = 0; d < child_depth; d++) {
            const uint64_t rem_len_d =
                (uint64_t)(params_.D - d);

            dir_start +=
                uleaf_count_by_depth_[d] * rem_len_d * 2ULL;

            count_before_level += uleaf_count_by_depth_[d];
        }

        const uint64_t local_idx =
            total_uleaf_before_pos - count_before_level;

        const uint64_t rem_len =
            (uint64_t)(params_.D - child_depth);

        dir_start += local_idx * rem_len * 2ULL;

        Rect r = child_r;

        for (uint64_t step = 0; step < rem_len; step++) {
            const uint64_t got =
                read_bits(ULD_, dir_start + 2ULL * step, 2);

            const int xm2 = midpoint(r.xmin, r.xmax);
            const int ym2 = midpoint(r.ymin, r.ymax);

            const uint64_t want =
                (uint64_t)quadrant_index(q, xm2, ym2);

            if (want != got) {
                return false;
            }

            r = child_rect(r, (int)got);
        }

        return true;
    }

    return false;
}