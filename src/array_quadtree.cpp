#include "array_quadtree.h"

#include <algorithm> // std::max
#include <cstdlib>   // std::abs

size_t ArrayQuadTree::bytes_used() const {
    // actual allocated bytes for nodes_ (capacity, not size)
    return nodes_.capacity() * sizeof(Node);
}

double ArrayQuadTree::bits_per_node() const {
    if (nodes_.empty()) return 0.0;
    return (8.0 * (double)bytes_used()) / (double)nodes_.size();
}

double ArrayQuadTree::bits_per_point(size_t n_points) const {
    if (n_points == 0) return 0.0;
    return (8.0 * (double)bytes_used()) / (double)n_points;
}


//Constructor
ArrayQuadTree::ArrayQuadTree(const Rect& world)
    : world_(world)
{
    // nodes_ will be allocated in build()
}

//Helpers
void ArrayQuadTree::geo_split(const Rect& r, int& mx, int& my) {
    mx = r.xmin + (r.xmax - r.xmin) / 2;
    my = r.ymin + (r.ymax - r.ymin) / 2;
}

int ArrayQuadTree::quadrant_of(int mx, int my, const Point& p) {
    // east if x >= mx, north if y >= my (half-open boundaries)
    const bool east  = (p.x >= mx);
    const bool north = (p.y >= my);

    // 0: NW, 1: NE, 2: SW, 3: SE
    if (!east && north) return 0; // NW
    if ( east && north) return 1; // NE
    if (!east && !north) return 2; // SW
    return 3; // SE
}

Rect ArrayQuadTree::quadrant_rect(const Rect& r, int mx, int my, int q) {
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
            // invalid quadrant id; return full region as fallback
            break;
    }

    return out;
}

int ArrayQuadTree::balance_score(const std::array<int,4>& counts, int n) {
    const int target = n / 4;
    int score = 0;
    for (int i = 0; i < 4; i++) {
        score += std::abs(counts[i] - target);
    }
    return score;
}

void ArrayQuadTree::choose_split_and_encode(int32_t node_id,
                                            const Rect& region,
                                            uint32_t start, uint32_t count,
                                            int& out_mx, int& out_my,
                                            std::array<int,4>& out_countsA,
                                            std::array<int,4>& out_countsB)
{
    // --- Split A: normal midpoint split of the current region ---
    int mxA, myA;
    geo_split(region, mxA, myA);

    out_countsA = {0,0,0,0};
    for (uint32_t i = 0; i < count; i++) {
        const Point& p = pts_[start + i];
        const int q = quadrant_of(mxA, myA, p);
        out_countsA[q]++;
    }

    // find the quadrant with max points (tie-breaking: smallest index)
    int heavyQ = 0;
    for (int q = 1; q < 4; q++) {
        if (out_countsA[q] > out_countsA[heavyQ]) heavyQ = q;
    }

    // --- Split B: split inside the heavy quadrant ---
    Rect heavyRegion = quadrant_rect(region, mxA, myA, heavyQ);

    int mxB, myB;
    geo_split(heavyRegion, mxB, myB);

    out_countsB = {0,0,0,0};
    for (uint32_t i = 0; i < count; i++) {
        const Point& p = pts_[start + i];
        const int q = quadrant_of(mxB, myB, p);
        out_countsB[q]++;
    }

    // --- Choose the more balanced split ---
    const int n = static_cast<int>(count);
    const int scoreA = balance_score(out_countsA, n);
    const int scoreB = balance_score(out_countsB, n);

    const bool useHeavy = (scoreB < scoreA);

    // --- Encode decision in meta (3 bits logically) ---
    uint8_t meta = 0;
    if (useHeavy) {
        meta |= 0x1;                 // bit0 = 1 => heavy split chosen
        meta |= (uint8_t(heavyQ) << 1); // bits1-2 = heavy quadrant id
        out_mx = mxB;
        out_my = myB;
    } else {
        // bit0 = 0 => normal split chosen (heavyQ bits don't matter)
        out_mx = mxA;
        out_my = myA;
    }

    nodes_[node_id].meta = meta;
}

void ArrayQuadTree::chosen_split(const Rect& region, uint8_t meta, int& mx, int& my) {
    // Split A: midpoint of region
    int mxA, myA;
    geo_split(region, mxA, myA);

    const bool useHeavy = (meta & 0x1) != 0;
    if (!useHeavy) {
        mx = mxA;
        my = myA;
        return;
    }

    const int heavyQ = (meta >> 1) & 0x3;
    Rect heavyRegion = quadrant_rect(region, mxA, myA, heavyQ);

    // Split B: midpoint of heavy region
    geo_split(heavyRegion, mx, my);
}

void ArrayQuadTree::build(const std::vector<Point>& points) {
    // reset old structure
    nodes_.clear();

    // copy points into build-time array (we will reorder it)
    pts_ = points;
    scratch_.clear();

    // empty input => empty tree
    if (pts_.empty()) return;

    // create root node at index 0
    nodes_.push_back(Node{});

    // build iteratively from root slice [0, pts_.size())
    build_iterative();

    // Option 1: points are implicit after build, so we can free payload storage
    pts_.clear();
    pts_.shrink_to_fit();
    scratch_.clear();
    scratch_.shrink_to_fit();
}

void ArrayQuadTree::build_iterative() {
    // Root covers entire world and owns all points initially
    std::vector<BuildTask> st;
    st.push_back(BuildTask{0, world_, 0u, static_cast<uint32_t>(pts_.size()), 0});

    while (!st.empty()) {
        BuildTask t = st.back();
        st.pop_back();

        const int32_t node_id = t.node_id;
        const Rect region = t.region;
        const uint32_t start = t.start;
        const uint32_t count = t.count;
        const int level = t.level;

        // Option 1: stop ONLY at unit cell
        const int w = region.xmax - region.xmin;
        const int h = region.ymax - region.ymin;
        const bool is_unit = (w == 1 && h == 1);

        if (is_unit) {
            // leaf: no children to create
            continue;
        }

        // If there are no points in this region, it should not exist as a node in our build,
        // because we only create children for non-empty quadrants.
        // (Root may still get count==0 if input empty, handled earlier.)
        if (count == 0) {
            continue;
        }

        // Choose split and store meta
        int mxChosen, myChosen;
        std::array<int,4> countsA, countsB;
        choose_split_and_encode(node_id, region, start, count,
                                mxChosen, myChosen, countsA, countsB);

        // Compute counts per quadrant under the CHOSEN split
        std::array<int,4> countsChosen{0,0,0,0};
        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            const int q = quadrant_of(mxChosen, myChosen, p);
            countsChosen[q]++;
        }

        // Compute offsets (prefix sums) for stable partition into 4 groups
        std::array<uint32_t,4> off;
        off[0] = 0;
        off[1] = off[0] + static_cast<uint32_t>(countsChosen[0]);
        off[2] = off[1] + static_cast<uint32_t>(countsChosen[1]);
        off[3] = off[2] + static_cast<uint32_t>(countsChosen[2]);

        // Ensure scratch can hold this slice
        if (scratch_.size() < count) scratch_.resize(count);

        // Scatter points into scratch grouped by quadrant
        std::array<uint32_t,4> pos = off;
        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            const int q = quadrant_of(mxChosen, myChosen, p);
            scratch_[pos[q]++] = p;
        }

        // Copy back into pts_ slice
        for (uint32_t i = 0; i < count; i++) {
            pts_[start + i] = scratch_[i];
        }

        // Create only non-empty children (indices, not pointers)
        // First clear existing children to -1
        nodes_[node_id].child[0] = nodes_[node_id].child[1] =
        nodes_[node_id].child[2] = nodes_[node_id].child[3] = -1;

        for (int q = 0; q < 4; q++) {
            const int cq = countsChosen[q];
            if (cq == 0) {
                continue; // no child
            }

            const uint32_t chStart = start + off[q];
            const uint32_t chCount = static_cast<uint32_t>(cq);

            Rect chRegion = quadrant_rect(region, mxChosen, myChosen, q);

            // Guard: skip empty region (can happen in degenerate splits)
            const int cw = chRegion.xmax - chRegion.xmin;
            const int ch = chRegion.ymax - chRegion.ymin;
            if (cw <= 0 || ch <= 0) {
                continue;
            }

            // Guard: ensure geometric progress (avoid infinite loop)
            if (cw == w && ch == h) {
                continue;
            }

            // Allocate child node
            const int32_t child_id = static_cast<int32_t>(nodes_.size());
            nodes_.push_back(Node{});
            nodes_[node_id].child[q] = child_id;

            // Push task
            st.push_back(BuildTask{child_id, chRegion, chStart, chCount, level + 1});
        }
    }
}

ArrayQuadTree::Stats ArrayQuadTree::compute_stats() const {
    Stats s;
    if (nodes_.empty()) return s;

    struct STTask {
        int32_t node_id;
        int depth;
    };

    std::vector<STTask> st;
    st.push_back(STTask{0, 0});

    while (!st.empty()) {
        STTask t = st.back();
        st.pop_back();

        const int32_t id = t.node_id;
        const int depth = t.depth;

        s.nodes++;
        if (depth > s.max_depth) s.max_depth = depth;

        bool hasChild = false;
        for (int q = 0; q < 4; q++) {
            const int32_t ch = nodes_[id].child[q];
            if (ch != -1) {
                hasChild = true;
                st.push_back(STTask{ch, depth + 1});
            }
        }

        if (!hasChild) {
            s.leaves++;
        } else {
            s.internal_nodes++;
            const bool useHeavy = (nodes_[id].meta & 0x1) != 0;
            if (useHeavy) s.heavy_splits++;
            else s.normal_splits++;
        }
    }

    return s;
}

void ArrayQuadTree::range_query(const Rect& q, std::vector<Point>& out) const {
    out.clear();
    if (nodes_.empty()) return;

    struct QTask {
        int32_t node_id;
        Rect region;
    };

    std::vector<QTask> st;
    st.push_back(QTask{0, world_});

    while (!st.empty()) {
        QTask t = st.back();
        st.pop_back();

        const int32_t id = t.node_id;
        const Rect region = t.region;

        // Prune if this region doesn't intersect the query
        if (!region.intersects(q)) continue;

        // Leaf check: no children
        bool hasChild = false;
        for (int qi = 0; qi < 4; qi++) {
            if (nodes_[id].child[qi] != -1) {
                hasChild = true;
                break;
            }
        }

        if (!hasChild) {
            // Option 1: leaf should be 1x1 cell => implicit point
            // It intersects q, so report (xmin, ymin)
            out.push_back(Point{region.xmin, region.ymin});
            continue;
        }

        // Internal node: reconstruct chosen split from meta
        int mx, my;
        chosen_split(region, nodes_[id].meta, mx, my);

        // Push existing children with their regions
        for (int qi = 0; qi < 4; qi++) {
            const int32_t ch = nodes_[id].child[qi];
            if (ch == -1) continue;

            Rect cr = quadrant_rect(region, mx, my, qi);
            if (cr.xmin >= cr.xmax || cr.ymin >= cr.ymax) continue; // safety

            st.push_back(QTask{ch, cr});
        }
    }
}





