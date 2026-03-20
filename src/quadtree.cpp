#include "quadtree.h"
#include <iostream>

//Rank/Select Support functions
// uint64_t MXQuadtreeBits::rank1_T(uint64_t bit_pos) const {
//     uint64_t cnt = 0;
//     for(uint64_t i=0;i<bit_pos && i<T_len_;i++){
//         cnt+=get_bit(T_,i);
//     }
//     return cnt;
// }


// uint64_t MXQuadtreeBits::rank_unary_S(uint64_t child_index) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < child_index; i++) {
//         uint64_t s = state_at_S(i);
//         if (s == 2ULL || s == 3ULL) cnt++;
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::rank_continue_S(uint64_t child_index) const {
//     uint64_t cnt = 0;
//     for (uint64_t i = 0; i < child_index; i++) {
//         uint64_t s = state_at_S(i);
//         if (s == 0ULL || s == 3ULL) cnt++;
//     }
//     return cnt;
// }

// uint64_t MXQuadtreeBits::next_node_id_for_child(uint64_t child_bit_pos) const {
//     const uint64_t child_index = rank1_T(child_bit_pos);
//     return 1ULL + rank_continue_S(child_index);
// }


//Helper functions for the main build
MXQuadtreeBits::NodeAnalysis
MXQuadtreeBits::analyze_node(const Rect& r, int depth, const std::vector<int>& ids) const {
    NodeAnalysis a; //Creates the result object that we will fill

    //Compute the width, height and area. Area will tell us if a region is full block or not because if number of points = area then its fullblock
    const uint64_t w=(uint64_t)(r.xmax-r.xmin);
    const uint64_t h=(uint64_t)(r.ymax-r.ymin);
    const uint64_t area = w*h;

    //Classifying the nodes
    a.is_unit_leaf = (depth >= params_.D);
    a.is_fullblock = (!a.is_unit_leaf && (uint64_t)ids.size() == area);
    a.expandable = (!a.is_unit_leaf && !a.is_fullblock);

    if (!a.expandable) {
        return a;
    }

    //We compute the 4 child rectangles if a is expandable and will contribute to some kind of child
    const int xm = midpoint(r.xmin, r.xmax);
    const int ym = midpoint(r.ymin, r.ymax);

    a.child_rects[0] = Rect{r.xmin, r.ymin, xm, ym};     // SW
    a.child_rects[1] = Rect{xm, r.ymin, r.xmax, ym};     // SE
    a.child_rects[2] = Rect{r.xmin, ym, xm, r.ymax};     // NW
    a.child_rects[3] = Rect{xm, ym, r.xmax, r.ymax};     // NE


    //Distribute the points into 4 children
    for(int i=0;i<4;i++){
        a.child_ids[i].clear();
    }

    for(int id:ids){
        const Point& p = points_[id];
        const int qi = quadrant_index(p,xm,ym);
        a.child_ids[qi].push_back(id);
    }

    //Compute child mask and non-empty children
    uint8_t mask=0;
    uint8_t cnt=0;

    for(int i=0;i<4;i++){
        if(!a.child_ids[i].empty()){
            mask |= (uint8_t)(1u << i);
            cnt++;
        }
    }

    a.childMask = mask;
    a.nonempty_children = cnt;

    return a;

}

void MXQuadtreeBits::build(const Rect& region, const std::vector<Point>& pts, const Params& params){
    params_ = params;
    region_ = region;
    points_ = pts;
    root_is_fullblock_=false;
    root_is_unitleaf_=false;

    stats_ = Stats{};
    stats_.points=points_.size();
    stats_.N=(uint64_t)(region_.xmax-region_.xmin);
    stats_.D=params_.D;

    build_bfs_contracted();
}

//Make the root of the tree
MXQuadtreeBits::Node 
MXQuadtreeBits::make_root() const{
    Node root;
    root.r = region_;
    root.depth = 0;
    root.ids.reserve(points_.size());

    for(int i=0; i<(int)points_.size();i++){
        root.ids.push_back(i);
    }

    return root;
}

//Make the child of the tree of each node
MXQuadtreeBits::Node 
MXQuadtreeBits::make_child(const NodeAnalysis& a, const Node& parent, int child_idx) const {
    Node child;
    child.r=a.child_rects[child_idx];
    child.depth = parent.depth +1;
    child.ids=a.child_ids[child_idx];
    return child;
}

//Calculate the child rectangle during query time
Rect MXQuadtreeBits::child_rect(const Rect& r, int child_idx) const {
    const int xm = midpoint(r.xmin, r.xmax);
    const int ym = midpoint(r.ymin, r.ymax);

    if (child_idx == 0) return Rect{r.xmin, r.ymin, xm, ym};     // SW
    if (child_idx == 1) return Rect{xm, r.ymin, r.xmax, ym};     // SE
    if (child_idx == 2) return Rect{r.xmin, ym, xm, r.ymax};     // NW
    return Rect{xm, ym, r.xmax, r.ymax};                         // NE
}

//Helpers of bit quadtree
void MXQuadtreeBits::push_bits(std::vector<uint64_t>& dst, uint64_t& bit_len, uint64_t v, int width) {
    if (width <= 0) return;
    while (width > 0) {
        const uint64_t word_idx = bit_len >> 6;
        const uint64_t bit_off  = bit_len & 63ULL;
        if (word_idx >= dst.size()) dst.push_back(0ULL);

        const int space = 64 - (int)bit_off;
        const int take  = (width < space) ? width : space;

        const uint64_t mask  = (take == 64) ? ~0ULL : ((1ULL << take) - 1ULL);
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

uint64_t MXQuadtreeBits::read_bits(const std::vector<uint64_t>& src, uint64_t bit_pos, int width){
    if(width<=0) return 0ULL;

    uint64_t result = 0;
    int written =0;

    while(width>0){
        const uint64_t word_idx = bit_pos>>6;
        const uint64_t bit_off = bit_pos & 63ULL;

        if(word_idx>=src.size()) break;

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

//Build for the contracted tree
void MXQuadtreeBits::build_bfs_contracted() {
    T_.clear();
    UM_.clear();
    ULL_.clear();
    ULD_.clear();
    UML_.clear();
    UMD_.clear();

    T_len_=UM_len_=ULL_len_=ULD_len_=UML_len_=UMD_len_=0;

    if(points_.empty()) return;

    Node root = make_root();  //Create the root of the quadtree
    NodeAnalysis ra = analyze_node(root.r,root.depth,root.ids);   //Analyse the root

    //If root is not expandable nothing to do
    if(!ra.expandable){
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

    while(!curr.empty()){
        std::vector<Node> next;

        for(const Node& node: curr){
            NodeAnalysis a = analyze_node(node.r,node.depth,node.ids);

            //write in T (4 bits child mask)
            push_bits(T_,T_len_,(uint64_t)a.childMask,4);
            stats_.internal_nodes++;

            //process each of the child
            for(int i=0;i<4;i++){
                if(((a.childMask>>i) & 1u) == 0u) continue; //if the child is 0 we do nothing just exit

                Node child = make_child(a,node,i);
                NodeAnalysis ca = analyze_node(child.r,child.depth,child.ids);

                bool is_unary = (ca.expandable && ca.nonempty_children == 1);

            if (is_unary) {
                UnarySkipResult sk = follow_unary_chain(child);

                const uint8_t L = sk.L;
                const bool terminal = (sk.endpoint.depth == params_.D);

                // 11 = unary leaf
                if(terminal){
                    push_bits(UM_, UM_len_, 3ULL, 2);
                    for (size_t j = 0; j < sk.dirs.size(); j++) {
                        push_bits(ULD_, ULD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                        push_bits(ULL_, ULL_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                    }
                }
                //10 = unary to mixed
                else{
                    push_bits(UM_, UM_len_, 2ULL, 2);
                    for (size_t j = 0; j < sk.dirs.size(); j++) {
                        push_bits(UMD_, UMD_len_, (uint64_t)(sk.dirs[j] & 3u), 2);
                        push_bits(UML_, UML_len_, (j + 1 == sk.dirs.size()) ? 1ULL : 0ULL, 1);
                    }
                    stats_.unary_to_mixed_nodes++;
                    //
                    next.push_back(std::move(sk.endpoint));
                }

            } else {
                // non-unary child
                if (ca.is_fullblock) {
                    // 01 = non-unary full block stop
                    push_bits(UM_, UM_len_, 1ULL, 2);
                    stats_.fullblock_nodes++;
                } else if (child.depth == params_.D) {
                    // 01 = non-unary leaf stop
                    push_bits(UM_, UM_len_, 1ULL, 2);
                    stats_.leaf_nodes++;
                } else {
                    // 00 = non-unary continue
                    push_bits(UM_, UM_len_, 0ULL, 2);
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
    stats_.UM_bits = UM_len_;
    stats_.ULL_bits = ULL_len_;
    stats_.ULD_bits = ULD_len_;
    stats_.UML_bits = UML_len_;
    stats_.UMD_bits = UMD_len_;

    uint64_t total_bits = T_len_ + UM_len_ + ULL_len_ + ULD_len_ + UML_len_ + UMD_len_;
    stats_.bpp = (stats_.points == 0) ? 0.0 : (double)total_bits / (double)stats_.points;
}

//Unary Skipping
MXQuadtreeBits::UnarySkipResult
MXQuadtreeBits::follow_unary_chain(const Node& start) const{
    UnarySkipResult res;
    Node cur = start;
    NodeAnalysis ca=analyze_node(cur.r,cur.depth,cur.ids);

    //check if by chance start is not unary expandable return len 0
    if(!(ca.expandable && ca.nonempty_children==1)) {
        res.L=0;
        res.endpoint=cur;
        res.endpoint_analysis=ca;
        return res;
    }

    while(true){
        int only_child = -1;
        for(int i=0;i<4;i++){
            if(((ca.childMask>>i) & 1u) != 0u){
                only_child =i;
                break;
            }
        }

        if(only_child <0){
            res.endpoint = cur;
            res.endpoint_analysis = ca;
            return res;
        }

        res.dirs.push_back((uint8_t)only_child);
        res.L++;

        Node nxt = make_child(ca,cur,only_child);
        NodeAnalysis na = analyze_node(nxt.r,nxt.depth,nxt.ids);

        //stop at first non unary node
        //Case 1: more than 1 child (2,3,4)
        if(!(na.expandable && na.nonempty_children == 1)){
            res.endpoint = std::move(nxt);
            res.endpoint_analysis = na;
            return res;
        }

        cur = std::move(nxt);
        ca=na;

        //5 bit cap because length cant be more than 20, but have to make this variable
        if(res.L>=31){
            res.endpoint=cur;
            res.endpoint_analysis = ca;
            return res;
        }
    }
}

//Membership Queries
// bool MXQuadtreeBits::contains(const Point& q) const{
//     if(points_.empty()) return false;
//     if(!region_.contains(q)) return false;

//     if(T_len_==0){
//         if(root_is_fullblock_){
//             return true;
//         }
//         if(root_is_unitleaf_){
//             for(const Point& p:points_){
//                 if(p.x==q.x && p.y==q.y){
//                     return true;
//                 }
//             }
//             return false;
//         }
//         return false;
//     }

//     QueryState cur;
//     cur.r=region_;
//     cur.depth=0;
//     cur.t_pos=0;

//     while (true) {
//         const uint64_t mask = read_bits(T_, cur.t_pos, 4);

//         const int xm = midpoint(cur.r.xmin, cur.r.xmax);
//         const int ym = midpoint(cur.r.ymin, cur.r.ymax);
//         const int qi = quadrant_index(q, xm, ym);

//         if (((mask >> qi) & 1ULL) == 0ULL) {
//             return false;
//         }

//         const Rect child_r = child_rect(cur.r, qi);
//         const int child_depth = cur.depth + 1;

//         const uint64_t child_bit_pos = cur.t_pos + (uint64_t)qi;
//         const uint64_t child_index = rank1_T(child_bit_pos);
//         const uint64_t s = state_at_S(child_index);

//         if (s == 1ULL) {
//             // 01 = non unary stop
//             return true;
//         }

//         if (s == 2ULL || s == 3ULL) {
//             // unary terminal or unary continue
//             const uint64_t unary_index = rank_unary_S(child_index);
//             const uint64_t L = read_bits(U_, unary_index * 5ULL, 5);
//             const uint64_t p_off = unary_path_offset(unary_index);

//             Rect walk_r = child_r;
//             int walk_depth = child_depth;

//             for (uint64_t step = 0; step < L; step++) {
//                 const int wxm = midpoint(walk_r.xmin, walk_r.xmax);
//                 const int wym = midpoint(walk_r.ymin, walk_r.ymax);

//                 const int qdir = quadrant_index(q, wxm, wym);
//                 const int sdir = (int)read_bits(P_, p_off + 2ULL * step, 2);

//                 if (qdir != sdir) {
//                     return false;
//                 }

//                 walk_r = child_rect(walk_r, sdir);
//                 walk_depth++;
//             }

//             if (s == 2ULL) {
//                 // 10 = unary terminal
//                 return true;
//             }

//             // 11 = unary continue
//             const uint64_t next_node_id = next_node_id_for_child(child_bit_pos);
//             cur.r = walk_r;
//             cur.depth = walk_depth;
//             cur.t_pos = 4ULL * next_node_id;
//             continue;
//         }

//         // s == 0: non-unary continue
//         if (child_depth == params_.D) {
//             // membership stop, but explicit in topology
//             return true;
//         }

//         const uint64_t next_node_id = next_node_id_for_child(child_bit_pos);
//         cur.r = child_r;
//         cur.depth = child_depth;
//         cur.t_pos = 4ULL * next_node_id;
//     }
    
// }

