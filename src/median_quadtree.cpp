#include "median_quadtree.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>

static constexpr size_t SUPER = 512; // bits per superblock (must be multiple of 64)

static inline uint32_t popcount_u64(uint64_t x) {
    return (uint32_t)__builtin_popcountll(x);
}

// Same quadrant convention as your existing code:
// 0: NW, 1: NE, 2: SW, 3: SE
int MedianBitQuadTree::quadrant_of(int mx, int my, const Point& p) {
    const bool east  = (p.x >= mx);
    const bool north = (p.y >= my);

    if (!east && north) return 0; // NW
    if ( east && north) return 1; // NE
    if (!east && !north) return 2; // SW
    return 3; // SE
}

Rect MedianBitQuadTree::quadrant_rect(const Rect& r, int mx, int my, int q) {
    Rect out = r;
    switch (q) {
        case 0: // NW
            out.xmax = mx;
            out.ymin = my;
            break;
        case 1: // NE
            out.xmin = mx;
            out.ymin = my;
            break;
        case 2: // SW
            out.xmax = mx;
            out.ymax = my;
            break;
        case 3: // SE
            out.xmin = mx;
            out.ymax = my;
            break;
        default:
            break;
    }
    return out;
}

// --- Bitvector ops (copied style from your SimpleBitQuadTree) ---
void MedianBitQuadTree::push_bit(bool b) {
    const size_t word = T_bits_ >> 6;
    const size_t bit  = T_bits_ & 63;
    if (word >= T_.size()) T_.push_back(0ULL);
    if (b) T_[word] |= (1ULL << bit);
    T_bits_++;
}

bool MedianBitQuadTree::get_bit(size_t pos) const {
    const size_t word = pos >> 6;
    const size_t bit  = pos & 63;
    return (T_[word] >> bit) & 1ULL;
}

uint8_t MedianBitQuadTree::get_mask(size_t node_i) const {
    const size_t base = 4 * node_i;
    uint8_t m = 0;
    m |= (uint8_t)get_bit(base + 0) << 0;
    m |= (uint8_t)get_bit(base + 1) << 1;
    m |= (uint8_t)get_bit(base + 2) << 2;
    m |= (uint8_t)get_bit(base + 3) << 3;
    return m; // 0..15
}

uint32_t MedianBitQuadTree::rank1(size_t pos) const {
    if (pos == 0) return 0;
    if (pos > T_bits_) pos = T_bits_;

    const size_t sb = pos / SUPER;
    uint32_t res = rank_super_[sb];

    const size_t sb_start = sb * SUPER;
    size_t cur = sb_start;

    while (cur + 64 <= pos) {
        const size_t w = cur >> 6;
        res += popcount_u64(T_[w]);
        cur += 64;
    }

    if (cur < pos) {
        const size_t w = cur >> 6;
        const size_t rem = pos - cur; // 1..63
        const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
        res += popcount_u64(T_[w] & mask);
    }

    return res;
}

void MedianBitQuadTree::build_rank() {
    rank_super_.clear();
    const size_t num_super = (T_bits_ + SUPER - 1) / SUPER + 1;
    rank_super_.resize(num_super, 0);

    uint32_t running = 0;
    for (size_t sb = 0; sb + 1 < num_super; sb++) {
        rank_super_[sb] = running;

        const size_t sb_start = sb * SUPER;
        const size_t sb_end = std::min(sb_start + SUPER, T_bits_);

        size_t pos = sb_start;
        while (pos + 64 <= sb_end) {
            const size_t w = pos >> 6;
            running += popcount_u64(T_[w]);
            pos += 64;
        }

        if (pos < sb_end) {
            const size_t w = pos >> 6;
            const size_t rem = sb_end - pos;
            const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
            running += popcount_u64(T_[w] & mask);
        }
    }

    rank_super_[num_super - 1] = running;
}

// --- Median split (safe, guarantees region shrinks) ---
void MedianBitQuadTree::median_split(const Rect& r, uint32_t start, uint32_t count,
                                     int& mx, int& my) {
    const int w = r.xmax - r.xmin;
    const int h = r.ymax - r.ymin;

    // fallback to geometric if degenerate or tiny
    auto geo_split = [&](int& gx, int& gy) {
        gx = r.xmin + w / 2;
        gy = r.ymin + h / 2;
    };

    if (count == 0 || w <= 1 || h <= 1) {
        geo_split(mx, my);
        mx = std::min(std::max(mx, r.xmin + 1), r.xmax - 1);
        my = std::min(std::max(my, r.ymin + 1), r.ymax - 1);
        return;
    }

    // collect coords
    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(count);
    ys.reserve(count);

    for (uint32_t i = 0; i < count; i++) {
        const Point& p = pts_[start + i];
        xs.push_back(p.x);
        ys.push_back(p.y);
    }

    const size_t mid = xs.size() / 2;

    std::nth_element(xs.begin(), xs.begin() + mid, xs.end());
    std::nth_element(ys.begin(), ys.begin() + mid, ys.end());

    mx = xs[mid];
    my = ys[mid];

    // Ensure splits are strictly inside region to guarantee progress.
    // If median hits boundary or would not shrink, fallback to geo split.
    if (mx <= r.xmin || mx >= r.xmax) {
        int gx, gy;
        geo_split(gx, gy);
        mx = gx;
    }
    if (my <= r.ymin || my >= r.ymax) {
        int gx, gy;
        geo_split(gx, gy);
        my = gy;
    }

    // Clamp (still must be strictly inside)
    mx = std::min(std::max(mx, r.xmin + 1), r.xmax - 1);
    my = std::min(std::max(my, r.ymin + 1), r.ymax - 1);
}

// --- Public build ---
void MedianBitQuadTree::build(const std::vector<Point>& points) {
    pts_ = points;
    scratch_.clear();

    T_.clear();
    T_bits_ = 0;
    dx_.clear();
    dy_.clear();
    rank_super_.clear();
    stats_ = Stats{};

    if (pts_.empty()) return;

    build_bfs();
    build_rank();

    // points are implicit after build
    pts_.clear();
    pts_.shrink_to_fit();
    scratch_.clear();
    scratch_.shrink_to_fit();
}

void MedianBitQuadTree::build_bfs() {
    struct Task {
        Rect region;
        uint32_t start;
        uint32_t count;
        int depth;
    };

    std::vector<Task> q;
    q.reserve(1024);
    q.push_back(Task{world_, 0u, (uint32_t)pts_.size(), 0});

    size_t head = 0;
    while (head < q.size()) {
        Task t = q[head++];
        const Rect region = t.region;
        const uint32_t start = t.start;
        const uint32_t count = t.count;
        const int depth = t.depth;

        stats_.max_depth = std::max(stats_.max_depth, depth);

        const int w = region.xmax - region.xmin;
        const int h = region.ymax - region.ymin;
        const bool is_unit = (w == 1 && h == 1);

        if (is_unit) {
            // leaf node exists in BFS order: mask=0000, dx/dy=0
            stats_.nodes++;
            stats_.leaves++;

            dx_.push_back(0);
            dy_.push_back(0);

            push_bit(false); push_bit(false); push_bit(false); push_bit(false);
            continue;
        }

        if (count == 0) {
            // should not happen (we only enqueue non-empty children)
            continue;
        }

        int mx, my;
        median_split(region, start, count, mx, my);

        // child counts
        std::array<int,4> counts{0,0,0,0};
        for (uint32_t i = 0; i < count; i++) {
            const int qi = quadrant_of(mx, my, pts_[start + i]);
            counts[qi]++;
        }

        // offsets into slice
        std::array<uint32_t,4> off;
        off[0] = 0;
        off[1] = off[0] + (uint32_t)counts[0];
        off[2] = off[1] + (uint32_t)counts[1];
        off[3] = off[2] + (uint32_t)counts[2];

        // stable partition into scratch_
        if (scratch_.size() < count) scratch_.resize(count);
        std::array<uint32_t,4> pos = off;

        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            const int qi = quadrant_of(mx, my, p);
            scratch_[pos[qi]++] = p;
        }
        for (uint32_t i = 0; i < count; i++) pts_[start + i] = scratch_[i];

        // internal node
        stats_.internal_nodes++;
        stats_.nodes++;

        // store split offsets aligned with node index
        dx_.push_back((uint32_t)(mx - region.xmin));
        dy_.push_back((uint32_t)(my - region.ymin));

        // topology bits
        for (int qi = 0; qi < 4; qi++) push_bit(counts[qi] > 0);

        // enqueue non-empty children in quadrant order (NW,NE,SW,SE)
        for (int qi = 0; qi < 4; qi++) {
            if (counts[qi] == 0) continue;

            Rect cr = quadrant_rect(region, mx, my, qi);
            const int cw = cr.xmax - cr.xmin;
            const int ch = cr.ymax - cr.ymin;

            // guard: skip invalid/no-progress regions
            if (cw <= 0 || ch <= 0) continue;
            if (cw == w && ch == h) continue;

            const uint32_t chStart = start + off[qi];
            const uint32_t chCount = (uint32_t)counts[qi];

            q.push_back(Task{cr, chStart, chCount, depth + 1});
        }
    }

    assert(T_bits_ == 4 * stats_.nodes);
    assert(dx_.size() == stats_.nodes);
    assert(dy_.size() == stats_.nodes);

    T_.shrink_to_fit();
    dx_.shrink_to_fit();
    dy_.shrink_to_fit();
    rank_super_.shrink_to_fit();
}

void MedianBitQuadTree::range_query(const Rect& q, std::vector<Point>& out) const {
    out.clear();
    if (stats_.nodes == 0) return;

    struct Task {
        size_t node_i; // BFS index in all-nodes stream
        Rect region;
        int depth;
    };

    std::vector<Task> st;
    st.push_back(Task{0, world_, 0});

    while (!st.empty()) {
        Task t = st.back();
        st.pop_back();

        const size_t i = t.node_i;
        const Rect region = t.region;

        if (!region.intersects(q)) continue;

        const int w = region.xmax - region.xmin;
        const int h = region.ymax - region.ymin;

        const uint8_t mask = get_mask(i);

        if ((w == 1 && h == 1) || mask == 0) {
            out.push_back(Point{region.xmin, region.ymin});
            continue;
        }

        // reconstruct split
        const int mx = region.xmin + (int)dx_[i];
        const int my = region.ymin + (int)dy_[i];

        const size_t baseBit = 4 * i;

        for (int qi = 0; qi < 4; qi++) {
            if ((mask & (1u << qi)) == 0) continue;

            const size_t child = 1 + (size_t)rank1(baseBit + (size_t)qi);

            Rect cr = quadrant_rect(region, mx, my, qi);
            if (cr.xmin >= cr.xmax || cr.ymin >= cr.ymax) continue;

            st.push_back(Task{child, cr, t.depth + 1});
        }
    }
}
