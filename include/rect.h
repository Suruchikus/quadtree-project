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