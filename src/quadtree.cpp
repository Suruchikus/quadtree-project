#include "quadtree.h"
#include <iostream>


// Rank/Select Support functions
uint64_t MXQuadtreeBits::rank1_T(uint64_t bit_pos) const {
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < bit_pos && i < T_len_; i++) {
        cnt += get_bit(T_, i);
    }
    return cnt;
}

uint64_t MXQuadtreeBits::rank1_EX(uint64_t bit_pos) const {
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < bit_pos && i < EX_len_; i++) {
        cnt += get_bit(EX_, i);
    }
    return cnt;
}

uint64_t MXQuadtreeBits::rank1_UM(uint64_t bit_pos) const {
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < bit_pos && i < UM_len_; i++) {
        cnt += get_bit(UM_, i);
    }
    return cnt;
}

uint64_t MXQuadtreeBits::select1_UML(uint64_t k) const {
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < UML_len_; i++) {
        if (get_bit(UML_, i)) {
            if (cnt == k) return i;
            cnt++;
        }
    }
    return UML_len_;
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

    if (!a.expandable) {
        return a;
    }

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

// Level-wise helpers
void MXQuadtreeBits::push_bits_level(std::vector<std::vector<uint64_t>>& dst_levels,
                                     std::vector<uint64_t>& len_levels,
                                     int level,
                                     uint64_t value,
                                     int width) {
    if (level < 0) return;
    if ((size_t)level >= dst_levels.size()) {
        dst_levels.resize((size_t)level + 1);
        len_levels.resize((size_t)level + 1, 0);
    }
    push_bits(dst_levels[(size_t)level], len_levels[(size_t)level], value, width);
}

uint64_t MXQuadtreeBits::get_bit_level(const std::vector<std::vector<uint64_t>>& src_levels,
                                       uint64_t bit_pos,
                                       int level) const {
    if (level < 0 || (size_t)level >= src_levels.size()) return 0ULL;
    return get_bit(src_levels[(size_t)level], bit_pos);
}

uint64_t MXQuadtreeBits::read_bits_level(const std::vector<std::vector<uint64_t>>& src_levels,
                                         uint64_t bit_pos,
                                         int width,
                                         int level) const {
    if (level < 0 || (size_t)level >= src_levels.size()) return 0ULL;
    return read_bits(src_levels[(size_t)level], bit_pos, width);
}

uint64_t MXQuadtreeBits::sum_level_bits(const std::vector<uint64_t>& len_levels) const {
    uint64_t total = 0;
    for (uint64_t x : len_levels) total += x;
    return total;
}

// Build for the contracted tree
void MXQuadtreeBits::build_bfs_contracted() {
    T_.clear();
    EX_.clear();
    UM_.clear();
    UML_.clear();
    UMD_.clear();

    UL_level_.clear();
    ULD_level_.clear();
    UL_level_len_.clear();
    ULD_level_len_.clear();

    T_len_ = EX_len_ = UM_len_ = UML_len_ = UMD_len_ = 0;

    // Stop streams are only needed for parent depths 0..D-2.
    const size_t stop_levels = (params_.D >= 1) ? (size_t)params_.D - 1 : 0;

    UL_level_.resize(stop_levels);
    ULD_level_.resize(stop_levels);
    UL_level_len_.assign(stop_levels, 0);
    ULD_level_len_.assign(stop_levels, 0);

    if (points_.empty()) return;

    Node root = make_root();
    NodeAnalysis ra = analyze_node(root.r, root.depth, root.ids);

    if (!ra.expandable) {
        root_is_fullblock_ = ra.is_fullblock;
        root_is_unitleaf_ = ra.is_unit_leaf;

        stats_.total_nodes = 1;

        if (ra.is_fullblock) {
            stats_.fullblock_nodes = 1;
        }
        if (ra.is_unit_leaf) {
            stats_.leaf_nodes = 1;
        }
        return;
    }

    std::vector<Node> curr;
    curr.push_back(std::move(root));

    while (!curr.empty()) {
        std::vector<Node> next;

        for (const Node& node : curr) {
            NodeAnalysis a = analyze_node(node.r, node.depth, node.ids);

            // write topology bits for this explicit node
            push_bits(T_, T_len_, (uint64_t)a.childMask, 4);
            stats_.internal_nodes++;

            // process each non-empty child
            for (int i = 0; i < 4; i++) {
                if (((a.childMask >> i) & 1u) == 0u) continue;

                Node child = make_child(a, node, i);
                NodeAnalysis ca = analyze_node(child.r, child.depth, child.ids);

                const bool is_unary = (ca.expandable && ca.nonempty_children == 1);

                if (is_unary) {
                    UnarySkipResult sk = follow_unary_chain(child);
                    const bool terminal = (sk.endpoint.depth == params_.D);

                    if (terminal) {
                        // EX=1 => no next explicit node in T
                        // unary-to-leaf stop:
                        //   UL_level_[d] = 0
                        //   ULD_level_[d] stores directions
                        // for parent depth D-1, no UL/ULD is needed because the stop is inferable
                        push_bits(EX_, EX_len_, 1ULL, 1);

                        if (node.depth <= params_.D - 2) {
                            push_bits_level(UL_level_, UL_level_len_, node.depth, 0ULL, 1);

                            for (size_t j = 0; j < sk.dirs.size(); j++) {
                                push_bits_level(ULD_level_, ULD_level_len_, node.depth,
                                                (uint64_t)(sk.dirs[j] & 3u), 2);
                            }
                        }

                        stats_.unary_to_leaf_nodes++;
                    } else {
                        // EX=0 => there is a next explicit node in T
                        // UM=0 => unary-to-mixed
                        push_bits(EX_, EX_len_, 0ULL, 1);
                        push_bits(UM_, UM_len_, 0ULL, 1);

                        for (size_t j = 0; j < sk.dirs.size(); j++) {
                            push_bits(UMD_, UMD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                            push_bits(UML_, UML_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                        }

                        stats_.unary_to_mixed_nodes++;
                        next.push_back(std::move(sk.endpoint));
                    }
                } else {
                    if (ca.is_fullblock) {
                        // EX=1 => stop at full block
                        // UL_level_[d] = 1, except at parent depth D-1 where no UL is needed
                        push_bits(EX_, EX_len_, 1ULL, 1);

                        if (node.depth <= params_.D - 2) {
                            push_bits_level(UL_level_, UL_level_len_, node.depth, 1ULL, 1);
                        }

                        stats_.fullblock_nodes++;
                    } else if (child.depth == params_.D) {
                        // EX=1 => direct leaf at depth D
                        // no UL/ULD needed at parent depth D-1
                        push_bits(EX_, EX_len_, 1ULL, 1);

                        if (node.depth <= params_.D - 2) {
                            // This branch should normally not happen for a correct MX quadtree build,
                            // because child.depth == D implies node.depth == D-1.
                            // Kept guarded for safety/minimal disruption.
                            push_bits_level(UL_level_, UL_level_len_, node.depth, 1ULL, 1);
                        }

                        stats_.leaf_nodes++;
                    } else {
                        // EX=0 => there is a next explicit node in T
                        // UM=1 => normal mixed continue
                        push_bits(EX_, EX_len_, 0ULL, 1);
                        push_bits(UM_, UM_len_, 1ULL, 1);

                        next.push_back(std::move(child));
                        stats_.mixed_internal++;
                    }
                }
            }
        }

        curr.swap(next);
    }

    uint64_t ones_T = 0;
    for (uint64_t i = 0; i < T_len_; i++) {
        ones_T += get_bit(T_, i);
    }

    stats_.total_nodes = 1 + ones_T;

    stats_.T_bits = T_len_;
    stats_.EX_bits = EX_len_;
    stats_.UM_bits = UM_len_;
    stats_.UL_bits = sum_level_bits(UL_level_len_);
    stats_.ULD_bits = sum_level_bits(ULD_level_len_);
    stats_.UML_bits = UML_len_;
    stats_.UMD_bits = UMD_len_;

    const uint64_t total_bits =
        T_len_ + EX_len_ + UM_len_ +
        stats_.UL_bits + stats_.ULD_bits +
        UML_len_ + UMD_len_;

    stats_.bpp = (stats_.points == 0) ? 0.0 : (double)total_bits / (double)stats_.points;
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
bool MXQuadtreeBits::contains(const Point& q) const {
    (void)q;
    // Not implemented yet for the new level-wise UL/ULD storage.
    // Later, queries can use:
    // - parent depth d to choose UL_level_[d] and ULD_level_[d]
    // - rank among EX=1 stops at that same depth
    // - fixed unary-to-leaf length = D - (d+1)
    return false;
}