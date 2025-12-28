#include "modified_quadtree.h"

#include <algorithm>
#include <array>
#include <limits>

QuadTree::QuadTree(const Rect& region, const std::vector<Point>& points, int leaf_size)
    : root(nullptr), leaf_size(leaf_size) {
    build_iterative(region, points);
}

int QuadTree::node_mx(const Node* node) {
    return (node->mode & 0x1) ? node->mx_store : geo_mx(node->region);
}

int QuadTree::node_my(const Node* node) {
    return (node->mode & 0x2) ? node->my_store : geo_my(node->region);
}

int QuadTree::quadrant(const Node* node, const Point& p) {
    const int mx = node_mx(node);
    const int my = node_my(node);

    // y < my => bottom; else => top
    // x < mx => left;  else => right
    if (p.y < my) {
        return (p.x < mx) ? 2 : 3;
    } else {
        return (p.x < mx) ? 0 : 1;
    }
}

bool QuadTree::region_can_split(const Rect& r) {
    // This prevents infinite splitting when midpoints stop changing.
    return (r.xmax - r.xmin) > 1 || (r.ymax - r.ymin) > 1;
}

int QuadTree::median_x(std::vector<Point>& pts) {
    const size_t n = pts.size();
    const size_t k = n / 2;
    std::nth_element(pts.begin(), pts.begin() + k, pts.end(),[](const Point& a, const Point& b) { return a.x < b.x; });
    return pts[k].x;
}

int QuadTree::median_y(std::vector<Point>& pts) {
    const size_t n = pts.size();
    const size_t k = n / 2;
    std::nth_element(pts.begin(), pts.begin() + k, pts.end(), [](const Point& a, const Point& b) { return a.y < b.y; });
    return pts[k].y;
}

static inline int imbalance_score(const std::array<int,4>& c, int n) {
    const int t = n / 4;
    auto dev = [&](int v){ return std::abs(v - t); };
    return std::max({dev(c[0]), dev(c[1]), dev(c[2]), dev(c[3])});
}

static inline std::array<int,4> counts_for_split(const std::vector<Point>& pts, int mx, int my){
    std::array<int,4> c{0,0,0,0};
    for (const auto& p : pts) {
        // 0 TL, 1 TR, 2 BL, 3 BR
        if (p.y < my) {
            if (p.x < mx) c[2]++; else c[3]++;
        } else {
            if (p.x < mx) c[0]++; else c[1]++;
        }
    }
    return c;
}

static inline bool split_is_degenerate(const std::array<int,4>& c, int n) {
    // If all points go to one child (or effectively no progress), treat as degenerate.
    return (c[0] == n) || (c[1] == n) || (c[2] == n) || (c[3] == n);
}

void QuadTree::choose_split_mode(Node* node, std::vector<Point>& pts) {
    const Rect& r = node->region;
    const int n = (int)pts.size();

    // Default geometric split
    const int mx_g = geo_mx(r);
    const int my_g = geo_my(r);

    // If region boundaries make geometric split unusable, we still handle it below.
    // Candidate medians (computed on copies to avoid messing with pts ordering too much)
    int mx_m = mx_g;
    int my_m = my_g;

    {
        std::vector<Point> tmp = pts;
        mx_m = median_x(tmp);
    }
    {
        std::vector<Point> tmp = pts;
        my_m = median_y(tmp);
    }

    // Clamp medians to be inside the region (otherwise they create empty rectangles)
    // If clamping ruins the split, we’ll fall back to geometric automatically.
    auto clamp_inside_x = [&](int x){
        // Want xmin < mx < xmax ideally.
        if (x <= r.xmin) return r.xmin + 1;
        if (x >= r.xmax) return r.xmax - 1;
        return x;
    };
    auto clamp_inside_y = [&](int y){
        if (y <= r.ymin) return r.ymin + 1;
        if (y >= r.ymax) return r.ymax - 1;
        return y;
    };

    // Only clamp when there is room; otherwise keep as-is
    if ((r.xmax - r.xmin) > 1) mx_m = clamp_inside_x(mx_m);
    if ((r.ymax - r.ymin) > 1) my_m = clamp_inside_y(my_m);

    struct Opt { uint8_t mode; int mx; int my; };
    // mode bits: 1 => store/use mx, 2 => store/use my
    const Opt options[4] = {
        {0,  mx_g, my_g},          // 00: geo, geo
        {1,  mx_m, my_g},          // 01: medx, geo
        {2,  mx_g, my_m},          // 10: geo, medy
        {3,  mx_m, my_m},          // 11: medx, medy
    };

    int bestScore = std::numeric_limits<int>::max();
    Opt best = options[0];

    for (const auto& opt : options) {
        // If split line equals boundary (no room), treat as poor candidate.
        if ((r.xmax - r.xmin) > 1) {
            if (opt.mx <= r.xmin || opt.mx >= r.xmax) continue;
        }
        if ((r.ymax - r.ymin) > 1) {
            if (opt.my <= r.ymin || opt.my >= r.ymax) continue;
        }

        auto c = counts_for_split(pts, opt.mx, opt.my);
        if (split_is_degenerate(c, n)) continue;

        const int score = imbalance_score(c, n);
        if (score < bestScore) {
            bestScore = score;
            best = opt;
        }
    }

    node->mode = best.mode;
    if (node->mode & 0x1) node->mx_store = best.mx;
    if (node->mode & 0x2) node->my_store = best.my;
}

void QuadTree::build_iterative(const Rect& region, const std::vector<Point>& points) {
    root = new Node(region);

    std::vector<BuildTask> stack;
    stack.push_back({root, points});

    while (!stack.empty()) {
        BuildTask task = std::move(stack.back());
        stack.pop_back();

        Node* node = task.node;
        std::vector<Point>& pts = task.points;
        const Rect& r = node->region;

        // Leaf / stop conditions
        if ((int)pts.size() <= leaf_size || !region_can_split(r)) {
            node->points = std::move(pts);
            continue;
        }

        // Choose split mode for this node (geometric vs median-x/y/both)
        choose_split_mode(node, pts);

        const int mx = node_mx(node);
        const int my = node_my(node);

        // Build 4 child rectangles using chosen cross (mx,my)
        // (0 TL, 1 TR, 2 BL, 3 BR)
        Rect sub[4] = {
            {r.xmin, my,  mx,    r.ymax}, // 0
            {mx,     my,  r.xmax, r.ymax}, // 1
            {r.xmin, r.ymin, mx,  my},     // 2
            {mx,     r.ymin, r.xmax, my}   // 3
        };

        std::vector<Point> buckets[4];
        buckets[0].reserve(pts.size()/4 + 1);
        buckets[1].reserve(pts.size()/4 + 1);
        buckets[2].reserve(pts.size()/4 + 1);
        buckets[3].reserve(pts.size()/4 + 1);

        for (const auto& p : pts) {
            int q;
            // Use chosen mx,my directly for speed
            if (p.y < my) q = (p.x < mx) ? 2 : 3;
            else          q = (p.x < mx) ? 0 : 1;
            buckets[q].push_back(p);
        }

        // If no progress (all points stayed in one bucket), stop splitting to avoid infinite depth
        int nonEmpty = 0;
        int maxBucket = 0;
        for (int i = 0; i < 4; i++) {
            if (!buckets[i].empty()) nonEmpty++;
            maxBucket = std::max(maxBucket, (int)buckets[i].size());
        }
        if (nonEmpty <= 1 || maxBucket == (int)pts.size()) {
            node->points = std::move(pts);
            continue;
        }

        for (int i = 0; i < 4; i++) {
            if (!buckets[i].empty()) {
                node->child[i] = new Node(sub[i]);
                stack.push_back({node->child[i], std::move(buckets[i])});
            }
        }
    }
}

bool QuadTree::member(const Point& p) const {
    Node* cur = root;

    while (cur != nullptr && !cur->is_leaf()) {
        int q = quadrant(cur, p);
        cur = cur->child[q];
    }

    if (cur == nullptr) return false;

    for (const auto& pt : cur->points) {
        if (pt.x == p.x && pt.y == p.y) return true;
    }
    return false;
}

void QuadTree::range_query(const Rect& q, std::vector<Point>& out) const {
    if (root == nullptr) return;

    std::vector<Node*> stack;
    stack.push_back(root);

    while (!stack.empty()) {
        Node* node = stack.back();
        stack.pop_back();

        if (!node->region.intersects(q)) continue;

        if (node->is_leaf()) {
            for (const auto& p : node->points) {
                if (q.contains(p)) out.push_back(p);
            }
        } else {
            for (int i = 0; i < 4; i++) {
                if (node->child[i] != nullptr) stack.push_back(node->child[i]);
            }
        }
    }
}

Stats QuadTree::stats() const {
    return Stats{};
}

void QuadTree::destroy(Node* n) {
    if (!n) return;
    for (int i = 0; i < 4; i++) destroy(n->child[i]);
    delete n;
}







