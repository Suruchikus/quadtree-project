#include "quadtree.h"
#include <iostream>


// Rank/Select Support functions
// uint64_t MXQuadtreeBits::rank1_T(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < bit_pos && i < T_len_; i++) {
//         cnt += get_bit(T_, i);
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::rank1_EX(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < bit_pos && i < EX_len_; i++) {
//         cnt += get_bit(EX_, i);
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::rank1_UM(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < bit_pos && i < UM_len_; i++) {
//         cnt += get_bit(UM_, i);
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::rank1_UL(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < bit_pos && i < UL_len_; i++) {
//         cnt += get_bit(UL_, i);
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::select1_UML(uint64_t k) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < UML_len_; i++) {
//         if (get_bit(UML_, i)) {
//             if (cnt == k) return i;
//             cnt++;
//         }
//     }
//     return UML_len_;
// }

// uint64_t MXQuadtreeBits::rank0_UM(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < bit_pos && i < UM_len_; i++) {
//         cnt += (get_bit(UM_, i) == 0ULL);
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::rank0_EX(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < bit_pos && i < EX_len_; i++) {
//         cnt += (get_bit(EX_, i) == 0ULL);
//     }
//     return cnt;
// }

uint64_t MXQuadtreeBits::rank1_T(uint64_t bit_pos) const {
    return rank_T_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank1_EX(uint64_t bit_pos) const {
    return rank_EX_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank1_UM(uint64_t bit_pos) const {
    return rank_UM_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank1_UL(uint64_t bit_pos) const {
    return rank_UL_.rank1_before(bit_pos);
}

uint64_t MXQuadtreeBits::rank0_UM(uint64_t bit_pos) const {
    return rank_UM_.rank0_before(bit_pos);
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


// Build for the contracted tree
void MXQuadtreeBits::build_bfs_contracted() {
    T_.clear();
    EX_.clear();
    UL_.clear();
    UM_.clear();
    UML_.clear();
    UMD_.clear();
    ULD_.clear();
    ULL_.clear();

    uint64_t x=0;

    T_len_ = EX_len_ = UL_len_ = ULD_len_ = ULL_len_ = UM_len_ = UML_len_ = UMD_len_ = 0;


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

    int bfs_level = 0;

    while (!curr.empty()) {
        std::vector<Node> next;

        for (const Node& node : curr) {
            NodeAnalysis a = analyze_node(node.r, node.depth, node.ids);

            uint64_t mask_to_write = 0;
            if (a.is_fullblock) {
                mask_to_write = 0b1111ULL;
            } else {
                mask_to_write = (uint64_t)a.childMask;
            }

            push_bits(T_, T_len_, mask_to_write, 4);
            stats_.internal_nodes++;

            for (int i = 0; i < 4; i++) {
                if (((a.childMask >> i) & 1u) == 0u) continue;

                Node child = make_child(a, node, i);
                NodeAnalysis ca = analyze_node(child.r, child.depth, child.ids);

                const bool is_unary = (ca.expandable && ca.nonempty_children == 1);

                // if (is_unary) {
                //     UnarySkipResult sk = follow_unary_chain(child);

                //     if (sk.endpoint_analysis.is_fullblock) {
                //         //Unary to Fullblock but we process it as unary to mixed leaf, EX=0, UM=0
                //         push_bits(EX_, EX_len_, 0ULL, 1);
                //         push_bits(UM_, UM_len_, 0ULL, 1);

                //         for (size_t j = 0; j < sk.dirs.size(); j++) {
                //             push_bits(UMD_, UMD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                //             push_bits(UML_, UML_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                //         }
                //         stats_.unary_to_mixed_nodes++;
                //         next.push_back(std::move(sk.endpoint));

                //     } else if (sk.endpoint_analysis.is_unit_leaf) {
                //         // unary-to-leaf which means EX=1 and UL=1
                //         push_bits(EX_, EX_len_, 1ULL, 1);
                //         push_bits(UL_, UL_len_, 1ULL, 1);

                //         const uint64_t real_len = (uint64_t)sk.dirs.size();
                //         const uint64_t uniform_len = (uint64_t)(params_.D - (bfs_level + 1));

                //         if(real_len<uniform_len){
                //             x++;
                //         }

                //         for (size_t j = 0; j < sk.dirs.size(); j++) {
                //             push_bits(ULD_, ULD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                //             push_bits(ULL_, ULL_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                //         }

                //         stats_.unary_to_leaf_nodes++;

                //     } else {
                //         // Unary to Mixed, EX=0 and UM=0
                //         push_bits(EX_, EX_len_, 0ULL, 1);
                //         push_bits(UM_, UM_len_, 0ULL, 1);

                //         for (size_t j = 0; j < sk.dirs.size(); j++) {
                //             push_bits(UMD_, UMD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                //             push_bits(UML_, UML_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                //         }

                //         stats_.unary_to_mixed_nodes++;
                //         next.push_back(std::move(sk.endpoint));
                //     }
                // }
                if (is_unary) {
                    UnarySkipResult sk = follow_unary_chain(child);

                    if (sk.endpoint_analysis.is_unit_leaf) {
                        // unary-to-leaf: keep compressed
                        push_bits(EX_, EX_len_, 1ULL, 1);
                        push_bits(UL_, UL_len_, 1ULL, 1);

                        for (size_t j = 0; j < sk.dirs.size(); j++) {
                            push_bits(ULD_, ULD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                            push_bits(ULL_, ULL_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                        }

                        stats_.unary_to_leaf_nodes++;

                    } else {
                        // unary-to-mixed or unary-to-fullblock:
                        // DO NOT skip. Store the immediate child as normal explicit node.
                        push_bits(EX_, EX_len_, 0ULL, 1);

                        next.push_back(std::move(child));
                        stats_.mixed_internal++;
                    }
                } else {
                    //Non Unary cases
                    if (ca.is_fullblock) {
                        // EX=1 => full block
                        push_bits(EX_, EX_len_, 1ULL, 1);
                        push_bits(UL_, UL_len_, 0ULL, 1);

                        stats_.fullblock_nodes++;

                    } else if (child.depth == params_.D) {
                        // EX=1 => direct leaf at depth D
                        push_bits(EX_, EX_len_, 1ULL, 1);
                        push_bits(UL_, UL_len_, 0ULL, 1);

                        stats_.leaf_nodes++;

                    } else {
                        // EX=0 => normal mixed continue
                        push_bits(EX_, EX_len_, 0ULL, 1);
                        //push_bits(UM_, UM_len_, 1ULL, 1);

                        next.push_back(std::move(child));
                        stats_.mixed_internal++;
                    }
                }
            }
        }

        curr.swap(next);
        bfs_level++;
    }

    rank_T_.build(T_, T_len_, 8);    // 8 words = 512 bits
    rank_EX_.build(EX_, EX_len_, 8);
    //rank_UM_.build(UM_, UM_len_, 8);
    rank_UL_.build(UL_, UL_len_, 8);

    stats_.rank_T_bits  = rank_T_.space_in_bits();
    stats_.rank_EX_bits = rank_EX_.space_in_bits();
    //stats_.rank_UM_bits = rank_UM_.space_in_bits();
    stats_.rank_UL_bits = rank_UL_.space_in_bits();

    stats_.rank_bits =
        stats_.rank_T_bits +
        stats_.rank_EX_bits +
        //stats_.rank_UM_bits +
        stats_.rank_UL_bits;

    uint64_t ones_T = 0;
    for (uint64_t i = 0; i < T_len_; i++) {
        ones_T += get_bit(T_, i);
    }

    stats_.total_nodes = 1 + ones_T;

    stats_.T_bits = T_len_;
    stats_.EX_bits = EX_len_;
    stats_.UM_bits = UM_len_;
    stats_.UL_bits = UL_len_;
    stats_.UML_bits = UML_len_;
    stats_.UMD_bits = UMD_len_;
    stats_.ULD_bits = ULD_len_;
    stats_.ULL_bits = ULL_len_;

    // const uint64_t total_bits =
    //     T_len_ + EX_len_ + UM_len_ +
    //     UL_len_ + ULL_len_ + stats_.ULD_bits +
    //     UML_len_ + UMD_len_ + stats_.rank_bits;

        const uint64_t total_bits =
        T_len_ + EX_len_ +
        UL_len_ + ULD_len_ +
         stats_.rank_bits;

    stats_.bpp = (stats_.points == 0)
        ? 0.0
        : (double)total_bits / (double)stats_.points;
    std::cout << "leaf of mix=" << x << "\n\n";

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
    // 1. Outside Grid
    if (q.x < region_.xmin || q.x >= region_.xmax ||
        q.y < region_.ymin || q.y >= region_.ymax) {
        return false;
    }

    // 2. Empty Grid
    if (points_.empty()) {
        return false;
    }

    // 3. Root special cases
    if (root_is_fullblock_) {
        return true;
    }

    if (root_is_unitleaf_) {
        return true;
    }

    // 4. Start navigation at the root explicit node
    Rect curr_r = region_;
    int curr_depth = 0;

    // Root is the first explicit node, so its 4-bit block starts at bit 0 in T.
    uint64_t curr_t_pos = 0;

    while (true) {
        // Read current node's 4-bit child mask from T
        const uint64_t childMask = read_bits(T_, curr_t_pos, 4);

        // Find which child quadrant contains q
        const int xm = midpoint(curr_r.xmin, curr_r.xmax);
        const int ym = midpoint(curr_r.ymin, curr_r.ymax);
        const int child_idx = quadrant_index(q, xm, ym);

        // If that child bit is 0, query point is absent
        if (((childMask >> child_idx) & 1ULL) == 0ULL) {
            return false;
        }

        // Map this chosen 1-child to its EX position.
        // EX is aligned with the 1-bits of T.
        const uint64_t child_t_bit_pos = curr_t_pos + (uint64_t)child_idx;
        const uint64_t ex_pos = rank1_T(child_t_bit_pos);

        // Read EX for this child
        const uint64_t ex_bit = get_bit(EX_, ex_pos);

        if (ex_bit == 1ULL) {
            // This child is a stopping case.
            // UL is aligned with EX=1 entries.
            const uint64_t ul_pos = rank1_EX(ex_pos);
            const uint64_t ul_bit = get_bit(UL_, ul_pos);

            Rect child_r = child_rect(curr_r, child_idx);
            const int child_depth = curr_depth + 1;

            if (ul_bit == 0ULL) {
                // direct fullblock OR direct leaf
                return true;
            } else {
                // unary-to-leaf

                const uint64_t prev_uleaf = rank1_UL(ul_pos);

                // find start/end direction index in ULD/ULL
                uint64_t dir_start = 0;
                uint64_t seen_paths = 0;

                while (seen_paths < prev_uleaf && dir_start < ULL_len_) {
                    if (get_bit(ULL_, dir_start) == 1ULL) {
                        seen_paths++;
                    }
                    dir_start++;
                }

                uint64_t dir_end = dir_start;
                while (dir_end < ULL_len_ && get_bit(ULL_, dir_end) == 0ULL) {
                    dir_end++;
                }
                if (dir_end < ULL_len_) dir_end++;

                Rect r = child_rect(curr_r, child_idx);

                for (uint64_t j = dir_start; j < dir_end; j++) {
                    const uint64_t got = read_bits(ULD_, 2 * j, 2);

                    const int xm2 = midpoint(r.xmin, r.xmax);
                    const int ym2 = midpoint(r.ymin, r.ymax);
                    const uint64_t want = (uint64_t)quadrant_index(q, xm2, ym2);

                    if (want != got) {
                        return false;
                    }

                    r = child_rect(r, (int)got);
                }

                return true;
            }
        } else {
            // continuation case: use UM aligned with EX=0
            const uint64_t um_pos = rank0_EX(ex_pos);
            const uint64_t um_bit = get_bit(UM_, um_pos);

            Rect child_r = child_rect(curr_r, child_idx);
            const int child_depth = curr_depth + 1;

            if (um_bit == 1ULL) {
                const uint64_t next_node_idx = 1 + ex_pos - rank1_EX(ex_pos);
                curr_t_pos = 4 * next_node_idx;
                curr_r = child_r;
                curr_depth = child_depth;
                continue;
            } else {
                // unary continuation
                const uint64_t prev_unary = rank0_UM(um_pos);

                // find start/end direction index in UMD/UML
                uint64_t dir_start = 0;
                uint64_t seen_paths = 0;

                while (seen_paths < prev_unary && dir_start < UML_len_) {
                    if (get_bit(UML_, dir_start) == 1ULL) {
                        seen_paths++;
                    }
                    dir_start++;
                }

                uint64_t dir_end = dir_start;
                while (dir_end < UML_len_ && get_bit(UML_, dir_end) == 0ULL) {
                    dir_end++;
                }
                // include the final 1 bit
                if (dir_end < UML_len_) dir_end++;

                Rect r = child_r;
                for (uint64_t j = dir_start; j < dir_end; j++) {
                    const uint64_t got = read_bits(UMD_, 2 * j, 2);    

                    const int xm2 = midpoint(r.xmin, r.xmax);
                    const int ym2 = midpoint(r.ymin, r.ymax);
                    const uint64_t want = (uint64_t)quadrant_index(q, xm2, ym2);

                    if (want != got) {
                        return false;
                    }

                    r = child_rect(r, (int)got);
                }

                // after unary path, land at endpoint explicit node
                const uint64_t next_node_idx = 1 + ex_pos - rank1_EX(ex_pos);
                curr_t_pos = 4 * next_node_idx;
                curr_r = r;
                curr_depth = child_depth + (int)(dir_end - dir_start);
                continue;
            }
        }   
    }

    return false;
}