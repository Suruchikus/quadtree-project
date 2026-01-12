#include "quadtree.h"

#include <cassert>

static constexpr size_t SUPER = 512;

static inline uint32_t popcount_u64(uint64_t x) {
    return (uint32_t)__builtin_popcountll(x);
}

SingletonStopQuadTreeV2::SingletonStopQuadTreeV2(const Rect& world)
    : world_(world)
{
    const int N = world_.xmax - world_.xmin;
    int d = 0;
    int t = N;
    while (t > 1) { t >>= 1; d++; }
    D_ = d;
}

size_t SingletonStopQuadTreeV2::bytes_S_levels() const {
    size_t s = 0;
    for (const auto& v : Slev_) s += v.capacity() * sizeof(uint64_t);
    return s;
}
size_t SingletonStopQuadTreeV2::bytes_rankS_levels() const {
    size_t s = 0;
    for (const auto& v : rankSlev_) s += v.capacity() * sizeof(uint32_t);
    return s;
}
size_t SingletonStopQuadTreeV2::bytes_P1_levels() const {
    size_t s = 0;
    for (const auto& v : P1lev_) s += v.capacity() * sizeof(uint64_t);
    return s;
}

// ---------------- geometry ----------------
void SingletonStopQuadTreeV2::geo_split(const Rect& r, int& mx, int& my) {
    mx = r.xmin + (r.xmax - r.xmin) / 2;
    my = r.ymin + (r.ymax - r.ymin) / 2;
}
int SingletonStopQuadTreeV2::quadrant_of(int mx, int my, const Point& p) {
    const bool east  = (p.x >= mx);
    const bool north = (p.y >= my);
    if (!east && north) return 0; // NW
    if ( east && north) return 1; // NE
    if (!east && !north) return 2; // SW
    return 3; // SE
}
Rect SingletonStopQuadTreeV2::quadrant_rect(const Rect& r, int mx, int my, int q) {
    Rect out = r;
    switch (q) {
        case 0: out.xmax = mx; out.ymin = my; break;
        case 1: out.xmin = mx; out.ymin = my; break;
        case 2: out.xmax = mx; out.ymax = my; break;
        case 3: out.xmin = mx; out.ymax = my; break;
        default: break;
    }
    return out;
}

// ---------------- bits ----------------
void SingletonStopQuadTreeV2::push_bit(std::vector<uint64_t>& bv, size_t& bits_used, bool b) {
    const size_t word = bits_used >> 6;
    const size_t bit  = bits_used & 63;
    if (word >= bv.size()) bv.push_back(0ULL);
    if (b) bv[word] |= (1ULL << bit);
    bits_used++;
}
bool SingletonStopQuadTreeV2::get_bit(const std::vector<uint64_t>& bv, size_t pos) const {
    const size_t word = pos >> 6;
    const size_t bit  = pos & 63;
    return (bv[word] >> bit) & 1ULL;
}
void SingletonStopQuadTreeV2::push_mask4(uint8_t mask4) {
    for (int b = 0; b < 4; b++) push_bit(T_, T_bits_, ((mask4 >> b) & 1u) != 0);
}
uint8_t SingletonStopQuadTreeV2::get_mask4(size_t node_i) const {
    const size_t base = 4 * node_i;
    uint8_t m = 0;
    m |= (uint8_t)get_bit(T_, base + 0) << 0;
    m |= (uint8_t)get_bit(T_, base + 1) << 1;
    m |= (uint8_t)get_bit(T_, base + 2) << 2;
    m |= (uint8_t)get_bit(T_, base + 3) << 3;
    return m;
}

// ---------------- rank ----------------
uint32_t SingletonStopQuadTreeV2::rank1(const std::vector<uint64_t>& bv,
                                       const std::vector<uint32_t>& sup,
                                       size_t bits_used,
                                       size_t pos) const
{
    if (pos == 0) return 0;
    if (pos > bits_used) pos = bits_used;

    const size_t sb = pos / SUPER;
    uint32_t res = sup[sb];

    const size_t sb_start = sb * SUPER;
    size_t cur = sb_start;

    while (cur + 64 <= pos) {
        const size_t w = cur >> 6;
        res += popcount_u64(bv[w]);
        cur += 64;
    }
    if (cur < pos) {
        const size_t w = cur >> 6;
        const size_t rem = pos - cur;
        const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
        res += popcount_u64(bv[w] & mask);
    }
    return res;
}

void SingletonStopQuadTreeV2::build_rank(const std::vector<uint64_t>& bv, size_t bits_used, std::vector<uint32_t>& out_sup) {
    out_sup.clear();
    const size_t num_super = (bits_used + SUPER - 1) / SUPER + 1;
    out_sup.resize(num_super, 0);

    uint32_t running = 0;
    for (size_t sb = 0; sb + 1 < num_super; sb++) {
        out_sup[sb] = running;

        const size_t sb_start = sb * SUPER;
        const size_t sb_end = std::min(sb_start + SUPER, bits_used);

        size_t p = sb_start;
        while (p + 64 <= sb_end) {
            const size_t w = p >> 6;
            running += popcount_u64(bv[w]);
            p += 64;
        }
        if (p < sb_end) {
            const size_t w = p >> 6;
            const size_t rem = sb_end - p;
            const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
            running += popcount_u64(bv[w] & mask);
        }
    }
    out_sup[num_super - 1] = running;
}

// ---------------- payload ----------------
void SingletonStopQuadTreeV2::push_bits(std::vector<uint64_t>& bv, size_t& bits_used, uint64_t value, int nbits) {
    for (int i = 0; i < nbits; i++) push_bit(bv, bits_used, ((value >> i) & 1ULL) != 0);
}
uint64_t SingletonStopQuadTreeV2::get_bits(const std::vector<uint64_t>& bv, size_t bitpos, int nbits) const {
    uint64_t v = 0;
    for (int i = 0; i < nbits; i++) {
        if (get_bit(bv, bitpos + (size_t)i)) v |= (1ULL << i);
    }
    return v;
}

void SingletonStopQuadTreeV2::ensure_level_storage(int depth) {
    if ((int)Slev_.size() > depth) return;
    const int need = depth + 1;
    Slev_.resize(need);
    Slev_bits_.resize(need, 0);
    rankSlev_.resize(need);

    P1lev_.resize(need);
    P1lev_bits_.resize(need, 0);
}

void SingletonStopQuadTreeV2::push_singleton_payload(int depth, const Rect& region, const Point& p) {
    const int k = D_ - depth; // remaining bits per axis
    // store offsets within region
    const uint32_t xoff = (uint32_t)(p.x - region.xmin);
    const uint32_t yoff = (uint32_t)(p.y - region.ymin);
    push_bits(P1lev_[depth], P1lev_bits_[depth], xoff, k);
    push_bits(P1lev_[depth], P1lev_bits_[depth], yoff, k);
}

Point SingletonStopQuadTreeV2::read_singleton_payload(int depth, const Rect& region, uint32_t singleton_idx) const {
    const int k = D_ - depth;
    const size_t stride = (size_t)(2 * k);
    const size_t base = (size_t)singleton_idx * stride;
    const uint32_t xoff = (uint32_t)get_bits(P1lev_[depth], base, k);
    const uint32_t yoff = (uint32_t)get_bits(P1lev_[depth], base + (size_t)k, k);
    return Point{ region.xmin + (int)xoff, region.ymin + (int)yoff };
}

// ---------------- build ----------------
void SingletonStopQuadTreeV2::build(const std::vector<Point>& points) {
    pts_ = points;
    scratch_.clear();

    T_.clear(); T_bits_ = 0;
    rankT_.clear();

    level_start_.clear();

    Slev_.clear(); Slev_bits_.clear(); rankSlev_.clear();
    P1lev_.clear(); P1lev_bits_.clear();

    stats_ = Stats{};

    if (pts_.empty()) return;

    build_bfs();
    build_rank(T_, T_bits_, rankT_);

    // build rank per level for Slev
    for (int d = 0; d < (int)Slev_.size(); d++) {
        build_rank(Slev_[d], Slev_bits_[d], rankSlev_[d]);
    }

    // make points implicit
    pts_.clear(); pts_.shrink_to_fit();
    scratch_.clear(); scratch_.shrink_to_fit();
}

void SingletonStopQuadTreeV2::build_bfs() {
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

    // record level starts in BFS node index space
    int current_depth = -1;

    while (head < q.size()) {
        Task t = q[head++];

        if (t.depth != current_depth) {
            current_depth = t.depth;
            // node index about to be emitted = stats_.nodes
            if ((int)level_start_.size() <= current_depth) level_start_.resize(current_depth + 1);
            level_start_[current_depth] = (uint32_t)stats_.nodes;
            ensure_level_storage(current_depth);
        }

        const Rect region = t.region;
        const uint32_t start = t.start;
        const uint32_t count = t.count;
        const int depth = t.depth;

        stats_.max_depth = std::max(stats_.max_depth, depth);

        const int w = region.xmax - region.xmin;
        const int h = region.ymax - region.ymin;
        const bool is_unit = (w == 1 && h == 1);

        // Emit node index (BFS)
        const size_t node_i = stats_.nodes;
        (void)node_i;

        if (is_unit) {
            stats_.nodes++;
            stats_.leaves++;

            // Slev bit for this node in this depth
            push_bit(Slev_[depth], Slev_bits_[depth], false);
            push_mask4(0u);
            continue;
        }

        // stop at singleton
        if (count == 1) {
            stats_.nodes++;
            stats_.leaves++;
            stats_.singleton_leaves++;

            push_bit(Slev_[depth], Slev_bits_[depth], true);
            push_singleton_payload(depth, region, pts_[start]);

            push_mask4(0u);
            continue;
        }

        // internal
        if (count == 0) continue;

        int mx, my;
        geo_split(region, mx, my);

        std::array<int,4> counts{0,0,0,0};
        for (uint32_t i = 0; i < count; i++) {
            counts[quadrant_of(mx, my, pts_[start + i])] ++;
        }

        // partition offsets
        std::array<uint32_t,4> off;
        off[0] = 0;
        off[1] = off[0] + (uint32_t)counts[0];
        off[2] = off[1] + (uint32_t)counts[1];
        off[3] = off[2] + (uint32_t)counts[2];

        if (scratch_.size() < count) scratch_.resize(count);
        std::array<uint32_t,4> pos = off;

        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            const int qi = quadrant_of(mx, my, p);
            scratch_[pos[qi]++] = p;
        }
        for (uint32_t i = 0; i < count; i++) pts_[start + i] = scratch_[i];

        stats_.nodes++;
        stats_.internal_nodes++;

        push_bit(Slev_[depth], Slev_bits_[depth], false);

        uint8_t mask = 0;
        for (int qi = 0; qi < 4; qi++) if (counts[qi] > 0) mask |= (uint8_t)(1u << qi);
        push_mask4(mask);

        // enqueue children
        for (int qi = 0; qi < 4; qi++) {
            if (counts[qi] == 0) continue;

            Rect cr = quadrant_rect(region, mx, my, qi);
            const int cw = cr.xmax - cr.xmin;
            const int ch = cr.ymax - cr.ymin;
            if (cw <= 0 || ch <= 0) continue;
            if (cw == w && ch == h) continue;

            q.push_back(Task{cr, start + off[qi], (uint32_t)counts[qi], depth + 1});
        }
    }

    // append final end marker for levels (optional but handy)
    level_start_.push_back((uint32_t)stats_.nodes);

    assert(T_bits_ == 4 * stats_.nodes);
}
 
// ---------------- query ----------------
void SingletonStopQuadTreeV2::range_query(const Rect& qrect, std::vector<Point>& out) const {
    out.clear();
    if (stats_.nodes == 0) return;

    struct QTask {
        size_t node_i;
        Rect region;
        int depth;
    };

    std::vector<QTask> st;
    st.push_back(QTask{0, world_, 0});

    while (!st.empty()) {
        QTask t = st.back();
        st.pop_back();

        if (!t.region.intersects(qrect)) continue;

        const uint8_t mask = get_mask4(t.node_i);
        const int w = t.region.xmax - t.region.xmin;
        const int h = t.region.ymax - t.region.ymin;

        if (mask == 0) {
            // leaf: either unit-cell leaf or singleton leaf
            if (w == 1 && h == 1) {
                Point p{t.region.xmin, t.region.ymin};
                if (qrect.contains(p)) out.push_back(p);
                continue;
            }

            const int d = t.depth;
            if (d >= (int)level_start_.size()-1) continue;

            const uint32_t level_start = level_start_[d];
            const uint32_t pos_in_level = (uint32_t)t.node_i - level_start;

            // is singleton?
            if (!get_bit(Slev_[d], pos_in_level)) {
                // should not happen (non-unit leaf must be singleton), but be safe
                continue;
            }

            const uint32_t singleton_idx = rank1_Slev(d, (size_t)pos_in_level + 1) - 1;
            Point p = read_singleton_payload(d, t.region, singleton_idx);
            if (qrect.contains(p)) out.push_back(p);
            continue;
        }

        // internal
        int mx, my;
        geo_split(t.region, mx, my);

        const size_t baseBit = 4 * t.node_i;

        for (int qi = 0; qi < 4; qi++) {
            if ((mask & (1u << qi)) == 0) continue;

            const size_t child = 1 + (size_t)rank1_T(baseBit + (size_t)qi);
            Rect cr = quadrant_rect(t.region, mx, my, qi);
            if (cr.xmin >= cr.xmax || cr.ymin >= cr.ymax) continue;

            st.push_back(QTask{child, cr, t.depth + 1});
        }
    }
}
