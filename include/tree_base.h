#pragma once

#include <vector>
#include "point.h"
#include "rect.h"
#include "stats.h"

class BaseTree {
    public:
        virtual ~BaseTree() = default;

        virtual bool member(const Point& p) const = 0;

        virtual void range_query(const Rect& r, std::vector<Point>& out) const = 0;

        virtual Stats stats() const = 0;
};