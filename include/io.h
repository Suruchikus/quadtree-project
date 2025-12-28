#pragma once
#include <string>
#include <vector>
#include "point.h"
#include "rect.h"

std::vector<Point> load_points_xy(const std::string& path);
Rect bounding_box_region(const std::vector<Point>& pts, int padding = 1);
