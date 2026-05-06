#pragma once

#include <vector>
#include <string>
#include <cstdint>

#include "point.h"
#include "rect.h"

struct QuerySet {
    std::string name;
    std::vector<Point> queries;
};

class MembershipQueryGenerator {
public:
    MembershipQueryGenerator(
        const Rect& region,
        const std::vector<Point>& points,
        uint64_t seed = 1
    );

    QuerySet random_filled(uint64_t count) const;
    QuerySet random_empty(uint64_t count) const;
    QuerySet isolated_filled(uint64_t count) const;

private:
    Rect region_;
    std::vector<Point> points_;
    std::vector<uint64_t> encoded_points_;
    uint64_t seed_;

    static uint64_t encode_point(const Point& p);

    bool contains_point(const Point& p) const;
};