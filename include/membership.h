#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include "point.h"
#include "rect.h"
#include "quadtree_v2.h"

struct QuerySet {
    std::string name;
    std::vector<Point> queries;
};

struct MembershipResult {
    std::string dataset_name;
    std::string query_type;

    uint64_t num_queries = 0;
    uint64_t repeats = 0;

    double avg_total_ms = 0.0;
    double avg_ns_per_query = 0.0;

    uint64_t found_count = 0;
};

class MembershipExperiment {
public:
    MembershipExperiment(
        const Rect& region,
        const std::vector<Point>& points,
        uint64_t seed = 1
    );

    QuerySet make_random_filled_queries(uint64_t count);
    QuerySet make_random_empty_queries(uint64_t count);

    // Simple reusable version first.
    // Later we can optimize this for exact farthest-neighbor computation.
    QuerySet make_isolated_filled_queries(uint64_t count);

    MembershipResult run(
        const std::string& dataset_name,
        const QuerySet& query_set,
        const MXQuadtreeBits& tree,
        uint64_t repeats
    ) const;

    static void write_csv(
        const std::string& out_path,
        const std::vector<MembershipResult>& results
    );

private:
    Rect region_;
    std::vector<Point> points_;
    uint64_t seed_;

    std::vector<uint64_t> encoded_points_;

    static uint64_t encode_point(const Point& p);
    bool contains_point(const Point& p) const;
};