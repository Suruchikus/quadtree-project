#include "quadtree.h"
#include <cassert>
#include <array>
#include <cstdlib>
#include <iostream>



//Destructor, when tree is deleted this function will delete all childrent and free the space
Node::~Node() {
    for(int i=0;i<4;i++) {
        delete child[i];
        child[i]=nullptr;
    }
}

//Constructor, basically initalizes
QuadTree:: QuadTree (const Rect& world, int leaf_size):
    world_(world), leaf_size_(leaf_size) {
        root_ = nullptr;
    }

void QuadTree::debug_print_node(int level,
                                const Rect& region,
                                uint32_t start, uint32_t count,
                                int mxA, int myA, const std::array<int,4>& cA, int scoreA,
                                int mxB, int myB, const std::array<int,4>& cB, int scoreB,
                                int mxChosen, int myChosen,
                                uint8_t meta) const
{
    // Decode meta
    const bool useHeavy = (meta & 0x1) != 0;
    const int heavyQ = (meta >> 1) & 0x3;

    // Indent by level (visual tree shape)
    for (int i = 0; i < level; i++) std::cout << "  ";

    std::cout << "L" << level
              << " region=[" << region.xmin << "," << region.ymin
              << " -> " << region.xmax << "," << region.ymax << ")"
              << " slice=[" << start << "," << (start + count) << ")"
              << " n=" << count
              << "\n";

    // Print A
    for (int i = 0; i < level; i++) std::cout << "  ";
    std::cout << "  A: (mx,my)=(" << mxA << "," << myA << ") "
              << "counts=[NW " << cA[0] << ", NE " << cA[1]
              << ", SW " << cA[2] << ", SE " << cA[3] << "] "
              << "score=" << scoreA << "\n";

    // Print B
    for (int i = 0; i < level; i++) std::cout << "  ";
    std::cout << "  B: (mx,my)=(" << mxB << "," << myB << ") "
              << "counts=[NW " << cB[0] << ", NE " << cB[1]
              << ", SW " << cB[2] << ", SE " << cB[3] << "] "
              << "score=" << scoreB
              << "  (heavyQ from A = " << heavyQ << ")"
              << "\n";

    // Print chosen
    for (int i = 0; i < level; i++) std::cout << "  ";
    std::cout << "  CHOSEN: (mx,my)=(" << mxChosen << "," << myChosen << ") "
              << (useHeavy ? "[HEAVY]" : "[NORMAL]")
              << "\n";
}


void QuadTree::geo_split(const Rect& r, int& mx, int & my){
    mx = r.xmin + (r.xmax - r.xmin) / 2;
    my = r.ymin + (r.ymax - r.ymin) / 2;
}

int QuadTree::quadrant_of(int mx, int my, const Point& p){
    // east if x >= mx, north if y >= my (half-open boundaries)
    const bool east  = (p.x >= mx);
    const bool north = (p.y >= my);

    // Quadrant indices:
    // 0: NW (west, north)
    // 1: NE (east, north)
    // 2: SW (west, south)
    // 3: SE (east, south)
    if (!east && north) return 0; //NW
    if (east && north) return 1; //NE
    if(!east && !north) return 2; //SW
    return 3;
}

Rect QuadTree::quadrant_rect(const Rect& r, int mx, int my, int q){
    Rect out = r;

    switch (q) {
        case 0: // NW: [xmin, mx) x [my, ymax)
            out.xmax = mx;
            out.ymin = my;
            break;
        case 1: // NE: [mx, xmax) x [my, ymax)
            out.xmin = mx;
            out.ymin = my;
            break;
        case 2: // SW: [xmin, mx) x [ymin, my)
            out.xmax = mx;
            out.ymax = my;
            break;
        case 3: // SE: [mx, xmax) x [ymin, my)
            out.xmin = mx;
            out.ymax = my;
            break;
        default:
            // Should never happen; return full region as fallback
            break;
    }

    return out;

}

int QuadTree::balance_score(const std::array<int,4>& counts, int n){
    //This is the ideal score and we wanna calculate how far of is the number of points in each quadrant from the ideal
    const int target = n/4;

    int score = 0;
    for(int i=0;i<4;i++){
        score += std::abs(counts[i] - target);
    }

    return score;
}

void QuadTree::choose_split_and_encode(Node* node, const Rect& region, uint32_t start, uint32_t count, int& out_mx, int& out_my, std::array<int,4>& out_countsA, std::array<int,4>& out_countsB) const {
    //First we try a normal quadtree geometric split
    int mxA, myA;
    geo_split(region, mxA, myA);

    //Count how many points will go in each quadrant with above split
    out_countsA = {0,0,0,0};
    for (uint32_t i=0; i< count;i++){
        const Point& p = pts_[start + i];
        int q = quadrant_of(mxA, myA, p);
        out_countsA[q]++;
    }

    //Look for the quadrant with the highest number of points
    int heavyQ = 0;
    for (int q=1; q<4; q++){
        if (out_countsA[q]>out_countsA[heavyQ]) heavyQ=q;
    }

    //Second we will compute the split based on heavy quadrant
    Rect heavyRegion = quadrant_rect(region, mxA, myA, heavyQ);

    int mxB,myB;
    geo_split(heavyRegion,mxB,myB);

    //Count distribution for second split strategy
    out_countsB = {0,0,0,0};
    for (uint32_t i = 0;i<count;i++){
        const Point& p = pts_[start + i];
        int q = quadrant_of(mxB,myB,p);
        out_countsB[q]++;
    }

    //Find the balance score
    const int n = static_cast<int> (count);
    int scoreA = balance_score(out_countsA,n);
    int scoreB = balance_score(out_countsB,n);

    //Is the second split better?
    bool useHeavy = (scoreB<scoreA);

    //Encode the decision
    uint8_t meta = 0;
    if(useHeavy) {
        meta |= 0x1; //bit0=1
        meta |= (uint8_t(heavyQ)<< 1); //bits1-2=heavyQ
        out_mx = mxB;
        out_my = myB;
    } else {
        //if useHeavy = 0 it means we use normal split
        out_mx = mxA;
        out_my = myA;
    }

    node->meta = meta;
}

void QuadTree::chosen_split(const Rect& region, uint8_t meta, int& mx, int& my) {
    // Split A: midpoint of region
    int mxA, myA;
    geo_split(region, mxA, myA);

    const bool useHeavy = (meta & 0x1) != 0;
    if (!useHeavy) {
        mx = mxA;
        my = myA;
        return;
    }

    // Heavy quadrant index stored in bits1-2
    const int heavyQ = (meta >> 1) & 0x3;

    // Heavy region is defined wrt split A
    Rect heavyRegion = quadrant_rect(region, mxA, myA, heavyQ);

    // Split B: midpoint of heavy region
    geo_split(heavyRegion, mx, my);
}




void QuadTree::build (const std::vector<Point>& input) {
    //Clear any old tree 
    delete root_;
    root_ = nullptr;

    pts_ = input; //Copy the input points to our internal vector where we will reorder too

    //Create root node
    root_ = new Node();
    root_->start = 0;
    root_->count = static_cast<uint32_t>(pts_.size());

    build_iterative(); //Calling the function that will build the tree

}

void QuadTree::build_iterative() {
    //No points
    if (!root_ || root_->count == 0) return;

    //temp vector to hold the reordering of points
    std::vector<Point> scratch;

    //We are using stack to avoid recursion and stackoverflow issues
    std::vector<BuildTask> st;
    st.push_back(BuildTask{root_,world_,root_->start,root_->count,0});

    while(!st.empty()){
        BuildTask t = st.back();
        st.pop_back();

        Node* node = t.node;
        const Rect region = t.region;
        const uint32_t start = t.start;
        const uint32_t count = t.count;
        const int level = t.level;

        //Store the start and end of the global array that this node holds
        node->start = start;
        node->count = count;

        //Check teh cndition where the leaf cannot split
        const int w = (region.xmax - region.xmin);
        const int h = (region.ymax - region.ymin);
        const bool cannot_split = (w < 2) || (h < 2); // i.e., 1x? or ?x1
        const bool is_unit = (w == 1 && h == 1);

        if (unit_cell_mode_) {
            // ONLY stop when we cannot split further (unit cell)
            if (is_unit) {
                // With unique grid points, count should be exactly 1 here.
                continue;
            }
        } else {
            // Normal mode (Option 0): stop by leaf_size OR cannot split. This was my first experiment keeping here for now
            if (count <= static_cast<uint32_t>(leaf_size_) || cannot_split) {
                continue;
            }
        }


        //Choose the split and encode the meta data
        int mxChosen, myChosen;
        std::array<int,4> countsA, countsB;

        choose_split_and_encode(node, region, start, count, mxChosen,myChosen,countsA,countsB);

        // For Debugging purposes
        const int n = static_cast<int>(count);
        const int scoreA = balance_score(countsA, n);
        const int scoreB = balance_score(countsB, n);

        // For debug printing, also compute mxA/myA and mxB/myB deterministically
        int mxA, myA;
        geo_split(region, mxA, myA);

        // heavy quadrant under A (same tie-breaking as in choose_split_and_encode)
        int heavyQ = 0;
        for (int q = 1; q < 4; q++) {
            if (countsA[q] > countsA[heavyQ]) heavyQ = q;
        }
        Rect heavyRegion = quadrant_rect(region, mxA, myA, heavyQ);
        int mxB, myB;
        geo_split(heavyRegion, mxB, myB);

        if (debug_) {
            debug_print_node(level, region, start, count,mxA, myA, countsA, scoreA,mxB, myB, countsB, scoreB,mxChosen, myChosen,node->meta);
        }

        //Rearrange the points of global array according to the split
        std::array<int,4> countsChosen{0,0,0,0};

        //Count of how many points go in each quadrant
        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            int q = quadrant_of(mxChosen, myChosen, p);
            countsChosen[q]++;
        }

        //The starting indices of each quadrant 
        std::array<uint32_t,4> off;
        off[0] = 0;
        off[1] = off[0] + static_cast<uint32_t>(countsChosen[0]);
        off[2] = off[1] + static_cast<uint32_t>(countsChosen[1]);
        off[3] = off[2] + static_cast<uint32_t>(countsChosen[2]);

        // Temp array to write into scratch[0..count) and then copy back
        if (scratch.size() < count) scratch.resize(count);

        // Current write positions (relative to 0..count)
        std::array<uint32_t,4> pos = off;

        // Second pass- scatter into scratch grouped by quadrant
        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            int q = quadrant_of(mxChosen, myChosen, p);
            scratch[pos[q]++] = p;
        }

        // Copy back into the global array slice
        for (uint32_t i = 0; i < count; i++) {
            pts_[start + i] = scratch[i];
        }

        //We create only non-empty children and push tasks for next iteration
        for (int q = 0; q < 4; q++) {
            const int cq = countsChosen[q];
            if (cq == 0) {
                node->child[q] = nullptr;
                continue;
            }

            const uint32_t chStart = start + off[q];
            const uint32_t chCount = static_cast<uint32_t>(cq);

            Rect chRegion = quadrant_rect(region, mxChosen, myChosen, q);

            // --- Guard 1: skip empty regions (shouldn't happen for valid children)
            const int cw = chRegion.xmax - chRegion.xmin;
            const int ch = chRegion.ymax - chRegion.ymin;
            if (cw <= 0 || ch <= 0) {
                // This would mean our geometry produced a zero-area region.
                // With cq>0, that indicates inconsistency; better to skip to avoid infinite loops.
                node->child[q] = nullptr;
                continue;
            }

            // --- Guard 2: ensure progress (child region must be strictly smaller than parent region)
            if (cw == (region.xmax - region.xmin) && ch == (region.ymax - region.ymin)) {
                // No geometric progress; avoid infinite loop
                node->child[q] = nullptr;
                continue;
            }

            Node* child = new Node();
            node->child[q] = child;

            st.push_back(BuildTask{child, chRegion, chStart, chCount, level + 1});
        }

    }
}

void QuadTree::range_query(const Rect& q, std::vector<Point>& out) const {
    out.clear();
    if (!root_) return;

    // Stack holds: node pointer + its region
    struct QTask {
        Node* node;
        Rect region;
    };

    std::vector<QTask> st;
    st.push_back(QTask{root_, world_});

    while (!st.empty()) {
        QTask t = st.back();
        st.pop_back();

        Node* node = t.node;
        const Rect region = t.region;

        // If query doesn't intersect this node region, skip
        if (!region.intersects(q)) continue;

        // Check if leaf (no children)
        bool hasChild = false;
        for (int i = 0; i < 4; i++) {
            if (node->child[i]) { hasChild = true; break; }
        }

        if (!hasChild) {
            // leaf is a unit cell, point is implicit
            // The leaf region is [x,x+1) x [y,y+1)
            // It intersects q, so output the point (x,y)
            out.push_back(Point{region.xmin, region.ymin});
            continue;
        }

        // Internal node: reconstruct chosen split, then push children
        int mx, my;
        chosen_split(region, node->meta, mx, my);

        // Push children that exist (and their regions)
        for (int qi = 0; qi < 4; qi++) {
            Node* ch = node->child[qi];
            if (!ch) continue;

            Rect cr = quadrant_rect(region, mx, my, qi);
            // safety: skip empty regions
            if (cr.xmin >= cr.xmax || cr.ymin >= cr.ymax) continue;

            st.push_back(QTask{ch, cr});
        }
    }
}


QuadTree::Stats QuadTree::compute_stats() const {
    Stats s;
    if (!root_) return s;

    std::vector<std::pair<Node*, int>> st;
    st.push_back({root_, 0});

    while (!st.empty()) {
        auto [cur, depth] = st.back();
        st.pop_back();

        s.nodes++;
        s.max_depth = std::max(s.max_depth, depth);

        bool hasChild = false;
        for (int q = 0; q < 4; q++) {
            if (cur->child[q]) {
                hasChild = true;
                st.push_back({cur->child[q], depth + 1});
            }
        }

        if (!hasChild) {
            s.leaves++;
        } else {
            s.internal_nodes++;
            bool useHeavy = (cur->meta & 0x1) != 0;
            if (useHeavy) s.heavy_splits++;
            else s.normal_splits++;
        }
    }

    return s;
}
