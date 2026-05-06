#include "membership.h"

#include <algorithm>
#include <random>
#include <iostream>
#include <numeric>
#include <limits>
#include <cmath>
#include <functional>

uint64_t MembershipQueryGenerator::encode_point(const Point& p) {
    return ((uint64_t)(uint32_t)p.x << 32) | (uint32_t)p.y;
}

MembershipQueryGenerator::MembershipQueryGenerator(
    const Rect& region,
    const std::vector<Point>& points,
    uint64_t seed
)
    : region_(region), points_(points), seed_(seed) {

    encoded_points_.reserve(points_.size());

    for (const Point& p : points_) {
        encoded_points_.push_back(encode_point(p));
    }

    std::sort(encoded_points_.begin(), encoded_points_.end());
}

bool MembershipQueryGenerator::contains_point(const Point& p) const {
    const uint64_t key = encode_point(p);

    return std::binary_search(
        encoded_points_.begin(),
        encoded_points_.end(),
        key
    );
}

QuerySet MembershipQueryGenerator::random_filled(uint64_t count) const {
    QuerySet qs;
    qs.name = "filled_random";
    qs.queries.reserve(count);

    if (points_.empty()) {
        return qs;
    }

    std::mt19937_64 rng(seed_);
    std::uniform_int_distribution<uint64_t> dist(0, points_.size() - 1);

    for (uint64_t i = 0; i < count; i++) {
        qs.queries.push_back(points_[dist(rng)]);
    }

    return qs;
}

QuerySet MembershipQueryGenerator::random_empty(uint64_t count) const {
    QuerySet qs;
    qs.name = "empty_random";
    qs.queries.reserve(count);

    std::mt19937_64 rng(seed_ + 100);

    std::uniform_int_distribution<int> xdist(region_.xmin, region_.xmax - 1);
    std::uniform_int_distribution<int> ydist(region_.ymin, region_.ymax - 1);

    while (qs.queries.size() < count) {
        Point p{xdist(rng), ydist(rng)};

        if (!contains_point(p)) {
            qs.queries.push_back(p);
        }
    }

    return qs;
}

QuerySet MembershipQueryGenerator::isolated_filled(uint64_t count) const {
    QuerySet qs;
    qs.name = "isolated_filled";

    const uint64_t n = points_.size();

    if (n == 0) {
        return qs;
    }

    qs.queries.reserve(std::min<uint64_t>(count, n));

    std::cout << "Building KD-tree for isolated filled queries...\n";

    std::vector<uint32_t> idx(n);
    for (uint64_t i = 0; i < n; i++) {
        idx[i] = (uint32_t)i;
    }

    auto coord = [&](uint32_t id, int axis) -> int {
        return axis == 0 ? points_[id].x : points_[id].y;
    };

    std::function<void(int, int, int)> build_kd =
        [&](int l, int r, int depth) {
            if (l >= r) return;

            int axis = depth & 1;
            int m = l + (r - l) / 2;

            std::nth_element(
                idx.begin() + l,
                idx.begin() + m,
                idx.begin() + r,
                [&](uint32_t a, uint32_t b) {
                    if (coord(a, axis) != coord(b, axis)) {
                        return coord(a, axis) < coord(b, axis);
                    }
                    return coord(a, 1 - axis) < coord(b, 1 - axis);
                }
            );

            build_kd(l, m, depth + 1);
            build_kd(m + 1, r, depth + 1);
        };

    build_kd(0, (int)n, 0);

    auto dist2 = [&](uint32_t a, uint32_t b) -> uint64_t {
        int64_t dx = (int64_t)points_[a].x - (int64_t)points_[b].x;
        int64_t dy = (int64_t)points_[a].y - (int64_t)points_[b].y;
        return (uint64_t)(dx * dx + dy * dy);
    };

    std::function<void(int, int, int, uint32_t, uint64_t&)> nearest_rec =
        [&](int l, int r, int depth, uint32_t target, uint64_t& best) {
            if (l >= r) return;

            int axis = depth & 1;
            int m = l + (r - l) / 2;

            uint32_t cur = idx[m];

            if (cur != target) {
                uint64_t d = dist2(target, cur);
                if (d < best) best = d;
            }

            int target_coord = coord(target, axis);
            int cur_coord = coord(cur, axis);

            int near_l, near_r, far_l, far_r;

            if (target_coord < cur_coord) {
                near_l = l;
                near_r = m;
                far_l = m + 1;
                far_r = r;
            } else {
                near_l = m + 1;
                near_r = r;
                far_l = l;
                far_r = m;
            }

            nearest_rec(near_l, near_r, depth + 1, target, best);

            int64_t diff = (int64_t)target_coord - (int64_t)cur_coord;
            uint64_t diff2 = (uint64_t)(diff * diff);

            if (diff2 < best) {
                nearest_rec(far_l, far_r, depth + 1, target, best);
            }
        };

    struct ScoredPoint {
        uint64_t nearest_dist2;
        uint32_t point_idx;
    };

    std::vector<ScoredPoint> scored;
    scored.reserve(n);

    std::cout << "Computing nearest-neighbor distance for each point...\n";

    for (uint64_t i = 0; i < n; i++) {
        uint64_t best = std::numeric_limits<uint64_t>::max();

        nearest_rec(0, (int)n, 0, (uint32_t)i, best);

        scored.push_back(ScoredPoint{best, (uint32_t)i});

        if ((i + 1) % 1000000 == 0) {
            std::cout << "  processed " << (i + 1) << " / " << n << "\n";
        }
    }

    uint64_t take = std::min<uint64_t>(count, scored.size());

    std::nth_element(
        scored.begin(),
        scored.begin() + take,
        scored.end(),
        [](const ScoredPoint& a, const ScoredPoint& b) {
            return a.nearest_dist2 > b.nearest_dist2;
        }
    );

    std::sort(
        scored.begin(),
        scored.begin() + take,
        [](const ScoredPoint& a, const ScoredPoint& b) {
            return a.nearest_dist2 > b.nearest_dist2;
        }
    );

    for (uint64_t i = 0; i < take; i++) {
        qs.queries.push_back(points_[scored[i].point_idx]);
    }

    std::cout << "Isolated filled queries generated = "
              << qs.queries.size() << "\n";

    if (!scored.empty()) {
        std::cout << "Largest nearest-neighbor distance squared = "
                  << scored[0].nearest_dist2 << "\n";
    }

    return qs;
}

