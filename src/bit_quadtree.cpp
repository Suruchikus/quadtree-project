#include "bit_quadtree.h"

#include<algorithm>
#include<cstdlib>
#include<cassert>

static constexpr size_t SUPER = 512; // bits per superblock (must be multiple of 64)

static inline uint32_t popcount_u64(uint64_t x) {
    return (uint32_t)__builtin_popcountll(x);
}

//Constructor
BitQuadTree::BitQuadTree(const Rect& world)
    : world_(world)
{
}

size_t BitQuadTree::bytes_T() const { return T_.size() * sizeof(uint64_t); }
size_t BitQuadTree::bytes_meta() const { return ((meta_bits_ + 63) / 64) * sizeof(uint64_t); }
size_t BitQuadTree::bytes_rank() const { return rank_super_.size() * sizeof(uint32_t); }

size_t BitQuadTree::bytes_used() const {
    return bytes_T() + bytes_meta() + bytes_rank();
}

uint32_t BitQuadTree::rank1(size_t pos) const {
    // rank1 over [0, pos)
    if (pos == 0) return 0;
    if (pos > T_bits_) pos = T_bits_;

    const size_t sb = pos / SUPER;
    uint32_t res = rank_super_[sb];

    const size_t sb_start = sb * SUPER;
    size_t cur = sb_start;

    // count full words from sb_start to pos
    while (cur + 64 <= pos) {
        const size_t w = cur >> 6;
        res += popcount_u64(T_[w]);
        cur += 64;
    }

    // leftover bits
    if (cur < pos) {
        const size_t w = cur >> 6;
        const size_t rem = pos - cur; // 1..63
        const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
        res += popcount_u64(T_[w] & mask);
    }

    return res;
}

BitQuadTree::HPStats BitQuadTree::analyze_heavy_paths() const {
    HPStats out;

    // total nodes in your current "all-nodes stored" representation:
    // meta_bits_ stores 3 bits per node (even leaves have meta=0)
    const uint64_t n = (uint64_t)(meta_bits_ / 3);
    out.nodes = n;

    if (n == 0) return out;

    // subtree leaf counts (each leaf contributes 1)
    // fits in uint32_t since leaves ~ 6.7M
    std::vector<uint32_t> sub(n, 0);

    // heavy path length from node (in nodes)
    // depth is <= ~30 so uint8_t is enough
    std::vector<uint8_t> hlen(n, 1);

    // Because children indices are always > parent index in BFS order,
    // we can compute bottom-up by iterating i = n-1..0.
    for (int64_t i = (int64_t)n - 1; i >= 0; --i) {
        const uint8_t mask = get_mask((size_t)i);

        // degree = number of existing children
        int deg = 0;
        for (int q = 0; q < 4; q++) deg += ((mask >> q) & 1u);

        if (deg == 0) {
            // leaf => one point
            sub[(size_t)i] = 1;
            hlen[(size_t)i] = 1;
            continue;
        }

        // internal node
        out.internal_nodes++;

        if (deg >= 1 && deg <= 4) out.deg_hist[deg]++;

        out.heavy_edges += 1;
        out.light_edges += (uint64_t)(deg - 1);

        const size_t baseBit = 4 * (size_t)i;

        uint32_t sum = 0;
        uint32_t best = 0;
        size_t bestChild = (size_t)-1;

        // iterate children in fixed quadrant order
        for (int q = 0; q < 4; q++) {
            if (((mask >> q) & 1u) == 0) continue;

            // child node index in the all-nodes BFS stream:
            // each 1-bit corresponds to one node entry
            const size_t child = 1 + (size_t)rank1(baseBit + (size_t)q);

            const uint32_t cs = sub[child];
            sum += cs;

            if (cs > best) {
                best = cs;
                bestChild = child;
            }
        }

        sub[(size_t)i] = sum;

        // heavy share = best / sum
        if (sum > 0) {
            out.avg_heavy_share += (double)best / (double)sum;
        }

        // heavy path length from i = 1 + heavy path length from heavy child
        // (since bestChild > i, hlen[bestChild] already computed)
        if (bestChild != (size_t)-1) {
            hlen[(size_t)i] = (uint8_t)(1 + hlen[bestChild]);
        } else {
            hlen[(size_t)i] = 1;
        }

        if ((int)hlen[(size_t)i] > out.max_heavy_path_len)
            out.max_heavy_path_len = (int)hlen[(size_t)i];
    }

    // leaves count = sub[root] (since each leaf is one point)
    out.leaves = sub[0];

    if (out.internal_nodes > 0) {
        out.rho = (double)out.light_edges / (double)out.internal_nodes;
        out.avg_heavy_share /= (double)out.internal_nodes;
    }

    // Heavy path length stats (histogram lengths, depth <= ~64 is safe)
    // We compute distribution over INTERNAL nodes only (more meaningful).
    std::vector<uint64_t> lenHist(128, 0);
    uint64_t internalCount = 0;
    uint64_t sumLen = 0;

    for (uint64_t i = 0; i < n; i++) {
        const uint8_t mask = get_mask((size_t)i);
        if (mask == 0) continue; // leaf
        const int L = (int)hlen[(size_t)i];
        if (L >= (int)lenHist.size()) continue;
        lenHist[(size_t)L]++;
        sumLen += (uint64_t)L;
        internalCount++;
    }

    if (internalCount > 0) {
        out.avg_heavy_path_len = (double)sumLen / (double)internalCount;

        auto percentile_from_hist = [&](double frac) -> int {
            const uint64_t target = (uint64_t)std::ceil(frac * (double)internalCount);
            uint64_t acc = 0;
            for (size_t L = 0; L < lenHist.size(); L++) {
                acc += lenHist[L];
                if (acc >= target) return (int)L;
            }
            return (int)lenHist.size() - 1;
        };

        out.median_heavy_path_len = percentile_from_hist(0.50);
        out.p90_heavy_path_len = percentile_from_hist(0.90);
    }

    return out;
}



void BitQuadTree::build_rank() {
    rank_super_.clear();

    // number of superblocks needed to cover [0, T_bits_]
    const size_t num_super = (T_bits_ + SUPER - 1) / SUPER + 1;
    rank_super_.resize(num_super, 0);

    uint32_t running = 0;
    for (size_t sb = 0; sb + 1 < num_super; sb++) {
        rank_super_[sb] = running;

        const size_t sb_start = sb * SUPER;
        const size_t sb_end = std::min(sb_start + SUPER, T_bits_);

        // count bits in [sb_start, sb_end)
        size_t pos = sb_start;
        while (pos + 64 <= sb_end) {
            const size_t w = pos >> 6;
            running += popcount_u64(T_[w]);
            pos += 64;
        }

        // handle leftover bits (not multiple of 64)
        if (pos < sb_end) {
            const size_t w = pos >> 6;
            const size_t rem = sb_end - pos; // 1..63
            const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
            running += popcount_u64(T_[w] & mask);
        }
    }

    // last entry
    rank_super_[num_super - 1] = running;
}


void BitQuadTree::push_bit(bool b) {
    const size_t word = T_bits_ >> 6;        // / 64
    const size_t bit  = T_bits_ & 63;        // % 64

    if (word >= T_.size()) {
        T_.push_back(0ULL);
    }
    if (b) {
        T_[word] |= (1ULL << bit);
    }
    T_bits_++;
}

bool BitQuadTree::get_bit(size_t pos) const{
    const size_t word = pos>>6;
    const size_t bit = pos & 63;
    return (T_[word]>>bit) & 1ULL;
}

uint8_t BitQuadTree::get_mask(size_t internal_i) const {
    const size_t base = 4 * internal_i;
    uint8_t m = 0;
    m |= (uint8_t)get_bit(base + 0) << 0;
    m |= (uint8_t)get_bit(base + 1) << 1;
    m |= (uint8_t)get_bit(base + 2) << 2;
    m |= (uint8_t)get_bit(base + 3) << 3;
    return m; // in [0,15]
}

void BitQuadTree::push_meta3(uint8_t m) {
    // we only use the lowest 3 bits
    m &= 0x7;

    // append 3 bits, LSB first (same style as push_bit)
    for (int b = 0; b < 3; b++) {
        const bool bitval = ((m >> b) & 1u) != 0;

        const size_t word = meta_bits_ >> 6; // /64
        const size_t bit  = meta_bits_ & 63; // %64

        if (word >= metaBits_.size()) {
            metaBits_.push_back(0ULL);
        }
        if (bitval) {
            metaBits_[word] |= (1ULL << bit);
        }
        meta_bits_++;
    }
}

uint8_t BitQuadTree::get_meta3(size_t node_i) const {
    const size_t base = 3 * node_i;
    uint8_t m = 0;

    for (int b = 0; b < 3; b++) {
        const size_t pos = base + (size_t)b;
        const size_t word = pos >> 6;
        const size_t bit  = pos & 63;
        const uint8_t bitval = (uint8_t)((metaBits_[word] >> bit) & 1ULL);
        m |= (bitval << b);
    }

    return m; // 0..7
}





void BitQuadTree::geo_split(const Rect& r, int& mx, int& my) {
    mx = r.xmin + (r.xmax - r.xmin) / 2;
    my = r.ymin + (r.ymax - r.ymin) / 2;
}

int BitQuadTree::quadrant_of(int mx, int my, const Point& p) {
    const bool east  = (p.x >= mx);
    const bool north = (p.y >= my);

    // 0: NW, 1: NE, 2: SW, 3: SE
    if (!east && north) return 0;
    if ( east && north) return 1;
    if (!east && !north) return 2;
    return 3;
}

Rect BitQuadTree::quadrant_rect(const Rect& r, int mx, int my, int q) {
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

int BitQuadTree::balance_score(const std::array<int,4>& counts, int n) {
    const int target = n / 4;
    int score = 0;
    for (int i = 0; i < 4; i++) {
        score += std::abs(counts[i] - target);
    }
    return score;
}

void BitQuadTree::chosen_split(const Rect& region, uint8_t meta, int& mx, int& my) {
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
    geo_split(heavyRegion, mx, my);
}

void BitQuadTree::choose_split_and_encode(const Rect& region,uint32_t start, uint32_t count,uint8_t& out_meta,int& out_mx, int& out_my)
{
    // --- Split A ---
    int mxA, myA;
    geo_split(region, mxA, myA);

    std::array<int,4> cA{0,0,0,0};
    for (uint32_t i = 0; i < count; i++) {
        const Point& p = pts_[start + i];
        const int q = quadrant_of(mxA, myA, p);
        cA[q]++;
    }

    // heavy quadrant under A
    int heavyQ = 0;
    for (int q = 1; q < 4; q++) {
        if (cA[q] > cA[heavyQ]) heavyQ = q;
    }

    // --- Split B: split heavy region of A ---
    Rect heavyRegion = quadrant_rect(region, mxA, myA, heavyQ);

    int mxB, myB;
    geo_split(heavyRegion, mxB, myB);

    std::array<int,4> cB{0,0,0,0};
    for (uint32_t i = 0; i < count; i++) {
        const Point& p = pts_[start + i];
        const int q = quadrant_of(mxB, myB, p);
        cB[q]++;
    }

    const int n = static_cast<int>(count);
    const int scoreA = balance_score(cA, n);
    const int scoreB = balance_score(cB, n);

    const bool useHeavy = (scoreB < scoreA);

    uint8_t meta = 0;
    if (useHeavy) {
        meta |= 0x1;
        meta |= (uint8_t(heavyQ) << 1);
        out_mx = mxB;
        out_my = myB;
    } else {
        out_mx = mxA;
        out_my = myA;
    }

    out_meta = meta;
}

void BitQuadTree::build(const std::vector<Point>& points) {
    // reset
    pts_ = points;
    scratch_.clear();

    T_.clear();
    T_bits_ = 0;
    metaBits_.clear();
    meta_bits_ = 0;
    rank_super_.clear();
    stats_ = Stats{};

    if (pts_.empty()) return;

    build_bfs();

    // build rank support (we'll implement next)
    build_rank();

    // points are implicit after build
    pts_.clear();
    pts_.shrink_to_fit();
    scratch_.clear();
    scratch_.shrink_to_fit();
}

void BitQuadTree::build_bfs() {
    struct Task {
        Rect region;
        uint32_t start;
        uint32_t count;
        int depth;
    };

    // BFS queue
    std::vector<Task> q;
    q.reserve(1024);
    q.push_back(Task{world_, 0u, static_cast<uint32_t>(pts_.size()), 0});

    // We'll process tasks in FIFO style using an index
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
            // Leaf node exists in BFS order: store mask=0000 and meta=0
            stats_.nodes++;
            stats_.leaves++;

            push_meta3(0);

            // push 4 zero bits for mask
            push_bit(false);
            push_bit(false);
            push_bit(false);
            push_bit(false);

            continue;
        }


        // This should not happen since we only enqueue non-empty children
        if (count == 0) {
            continue;
        }

        // Compute split choice + meta + chosen (mx,my)
        uint8_t meta = 0;
        int mx, my;
        choose_split_and_encode(region, start, count, meta, mx, my);

        // Determine child counts under chosen split
        std::array<int,4> counts{0,0,0,0};
        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            const int qi = quadrant_of(mx, my, p);
            counts[qi]++;
        }

        // Compute offsets into the slice [start, start+count)
        std::array<uint32_t,4> off;
        off[0] = 0;
        off[1] = off[0] + static_cast<uint32_t>(counts[0]);
        off[2] = off[1] + static_cast<uint32_t>(counts[1]);
        off[3] = off[2] + static_cast<uint32_t>(counts[2]);

        // Stable partition into scratch_
        if (scratch_.size() < count) scratch_.resize(count);
        std::array<uint32_t,4> pos = off;

        for (uint32_t i = 0; i < count; i++) {
            const Point& p = pts_[start + i];
            const int qi = quadrant_of(mx, my, p);
            scratch_[pos[qi]++] = p;
        }

        for (uint32_t i = 0; i < count; i++) {
            pts_[start + i] = scratch_[i];
        }

        // This node is internal: append its 4 child-existence bits to T_
        // and store its meta in meta_
        stats_.internal_nodes++;
        stats_.nodes++;
        push_meta3(meta);

        for (int qi = 0; qi < 4; qi++) {
            const bool exists = (counts[qi] > 0);
            push_bit(exists);
        }

        // Enqueue children in quadrant order (NW,NE,SW,SE) so BFS order is consistent
        for (int qi = 0; qi < 4; qi++) {
            if (counts[qi] == 0) continue;

            Rect cr = quadrant_rect(region, mx, my, qi);
            const int cw = cr.xmax - cr.xmin;
            const int ch = cr.ymax - cr.ymin;

            // guard: skip invalid / no-progress regions
            if (cw <= 0 || ch <= 0) continue;
            if (cw == w && ch == h) continue;

            const uint32_t chStart = start + off[qi];
            const uint32_t chCount = static_cast<uint32_t>(counts[qi]);

            q.push_back(Task{cr, chStart, chCount, depth + 1});
        }
    }

    assert(meta_bits_ == 3 * stats_.nodes);
    assert(T_bits_ == 4 * stats_.nodes);
    T_.shrink_to_fit();
    metaBits_.shrink_to_fit();
    rank_super_.shrink_to_fit();


}

void BitQuadTree::range_query(const Rect& q, std::vector<Point>& out) const {
    out.clear();
    if (stats_.leaves == 0 && stats_.internal_nodes == 0) return;

    struct Task {
        // internal node index (BFS order)
        size_t internal_i;
        Rect region;
        int depth;
    };

    std::vector<Task> st;
    st.push_back(Task{0, world_, 0});

    while (!st.empty()) {
        Task t = st.back();
        st.pop_back();

        const size_t i = t.internal_i;
        const Rect region = t.region;

        if (!region.intersects(q)) continue;

        const int w = region.xmax - region.xmin;
        const int h = region.ymax - region.ymin;

        const uint8_t mask = get_mask(i);

        if (w == 1 && h == 1 || mask == 0) {
            out.push_back(Point{region.xmin, region.ymin});
            continue;
        }

        // Reconstruct chosen split from meta
        int mx, my;
        const uint8_t meta = get_meta3(i);
        chosen_split(region, meta, mx, my);


        const size_t baseBit = 4 * i;

        // Traverse children that exist
        for (int qi = 0; qi < 4; qi++) {
            if ((mask & (1u << qi)) == 0) continue;

            // internal index of this child (BFS order among internal nodes)
            const size_t child_internal = 1+(size_t)rank1(baseBit + (size_t)qi);

            Rect cr = quadrant_rect(region, mx, my, qi);
            if (cr.xmin >= cr.xmax || cr.ymin >= cr.ymax) continue;

            st.push_back(Task{child_internal, cr, t.depth + 1});
        }
    }
}



// ===============================
// SimpleBitQuadTree implementation
// ===============================

static inline uint32_t popcount_u64_simple(uint64_t x) {
    return (uint32_t)__builtin_popcountll(x);
}

void SimpleBitQuadTree::geo_split(const Rect& r, int& mx, int& my) {
    mx = r.xmin + (r.xmax - r.xmin) / 2;
    my = r.ymin + (r.ymax - r.ymin) / 2;
}

int SimpleBitQuadTree::quadrant_of(int mx, int my, const Point& p) {
    const bool east  = (p.x >= mx);
    const bool north = (p.y >= my);

    // 0: NW, 1: NE, 2: SW, 3: SE
    if (!east && north) return 0;
    if ( east && north) return 1;
    if (!east && !north) return 2;
    return 3;
}

Rect SimpleBitQuadTree::quadrant_rect(const Rect& r, int mx, int my, int q) {
    Rect out = r;
    switch (q) {
        case 0: out.xmax = mx; out.ymin = my; break; // NW
        case 1: out.xmin = mx; out.ymin = my; break; // NE
        case 2: out.xmax = mx; out.ymax = my; break; // SW
        case 3: out.xmin = mx; out.ymax = my; break; // SE
        default: break;
    }
    return out;
}

void SimpleBitQuadTree::push_bit(bool b) {
    const size_t word = T_bits_ >> 6;
    const size_t bit  = T_bits_ & 63;
    if (word >= T_.size()) T_.push_back(0ULL);
    if (b) T_[word] |= (1ULL << bit);
    T_bits_++;
}

bool SimpleBitQuadTree::get_bit(size_t pos) const {
    const size_t word = pos >> 6;
    const size_t bit  = pos & 63;
    return (T_[word] >> bit) & 1ULL;
}

uint8_t SimpleBitQuadTree::get_mask(size_t node_i) const {
    const size_t base = 4 * node_i;
    uint8_t m = 0;
    m |= (uint8_t)get_bit(base + 0) << 0;
    m |= (uint8_t)get_bit(base + 1) << 1;
    m |= (uint8_t)get_bit(base + 2) << 2;
    m |= (uint8_t)get_bit(base + 3) << 3;
    return m;
}

uint32_t SimpleBitQuadTree::rank1(size_t pos) const {
    if (pos == 0) return 0;
    if (pos > T_bits_) pos = T_bits_;

    const size_t sb = pos / SUPER;
    uint32_t res = rank_super_[sb];

    const size_t sb_start = sb * SUPER;
    size_t cur = sb_start;

    while (cur + 64 <= pos) {
        const size_t w = cur >> 6;
        res += popcount_u64_simple(T_[w]);
        cur += 64;
    }

    if (cur < pos) {
        const size_t w = cur >> 6;
        const size_t rem = pos - cur; // 1..63
        const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
        res += popcount_u64_simple(T_[w] & mask);
    }

    return res;
}

void SimpleBitQuadTree::build_rank() {
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
            running += popcount_u64_simple(T_[w]);
            pos += 64;
        }

        if (pos < sb_end) {
            const size_t w = pos >> 6;
            const size_t rem = sb_end - pos;
            const uint64_t mask = (rem == 64) ? ~0ULL : ((1ULL << rem) - 1ULL);
            running += popcount_u64_simple(T_[w] & mask);
        }
    }

    rank_super_[num_super - 1] = running;
}

void SimpleBitQuadTree::build(const std::vector<Point>& points) {
    pts_ = points;
    scratch_.clear();

    T_.clear();
    T_bits_ = 0;
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

void SimpleBitQuadTree::build_bfs() {
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
            stats_.nodes++;
            stats_.leaves++;

            // leaf => push 0000 mask
            push_bit(false); push_bit(false); push_bit(false); push_bit(false);
            continue;
        }

        if (count == 0) continue;

        int mx, my;
        geo_split(region, mx, my);

        std::array<int,4> counts{0,0,0,0};
        for (uint32_t i = 0; i < count; i++) {
            const int qi = quadrant_of(mx, my, pts_[start + i]);
            counts[qi]++;
        }

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

        // internal node
        stats_.internal_nodes++;
        stats_.nodes++;

        for (int qi = 0; qi < 4; qi++) push_bit(counts[qi] > 0);

        // enqueue non-empty children in quadrant order
        for (int qi = 0; qi < 4; qi++) {
            if (counts[qi] == 0) continue;

            Rect cr = quadrant_rect(region, mx, my, qi);
            const int cw = cr.xmax - cr.xmin;
            const int ch = cr.ymax - cr.ymin;

            if (cw <= 0 || ch <= 0) continue;
            if (cw == w && ch == h) continue;

            const uint32_t chStart = start + off[qi];
            const uint32_t chCount = (uint32_t)counts[qi];
            q.push_back(Task{cr, chStart, chCount, depth + 1});
        }
    }

    // For simple tree: 4 bits per node
    assert(T_bits_ == 4 * stats_.nodes);
    T_.shrink_to_fit();
    rank_super_.shrink_to_fit();
}

void SimpleBitQuadTree::range_query(const Rect& q, std::vector<Point>& out) const {
    out.clear();
    if (stats_.nodes == 0) return;

    struct Task {
        size_t node_i; // BFS index in "all nodes stream"
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

        int mx, my;
        geo_split(region, mx, my);

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










