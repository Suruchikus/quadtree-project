//We are using Half-open intervals to ensure that every point belongs to exactly one region, every split is clean, and no boundary cases need special handling.

#pragma once

#include "point.h"

struct Rect {
    int xmin;
    int ymin;
    int xmax;
    int ymax;

    inline bool contains(const Point& p) const {
        return p.x >= xmin && p.x < xmax && p.y >= ymin && p.y<ymax;
    } 

    inline bool intersects(const Rect& r) const {
        return !(r.xmax <= xmin || r.xmin >= xmax || r.ymax <= ymin || r.ymin >= ymax);
    }
};