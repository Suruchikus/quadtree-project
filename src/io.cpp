#include "io.h"
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <sstream>

Rect bounding_box_region(const std::vector<Point>& pts, int padding) {
    int xmin = pts[0].x, xmax = pts[0].x;
    int ymin = pts[0].y, ymax = pts[0].y;

    for (const auto& p : pts) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }

    Rect r;
    r.xmin = xmin - padding;
    r.ymin = ymin - padding;
    r.xmax = xmax + 1 + padding;
    r.ymax = ymax + 1 + padding;
    return r;
}


static bool parse_point_line(const std::string& line, Point& out) {
    // Skip empty lines
    if (line.empty()) return false;

    // Try comma-separated: x,y
    {
        size_t comma = line.find(',');
        if (comma != std::string::npos) {
            std::string xs = line.substr(0, comma);
            std::string ys = line.substr(comma + 1);

            // Allow trailing spaces
            try {
                out.x = std::stoi(xs);
                out.y = std::stoi(ys);
                return true;
            } catch (...) {
                return false;
            }
        }
    }

    // Try whitespace-separated: x y
    {
        int x, y;
        // simple manual parse using stringstream
        std::istringstream iss(line);
        if (iss >> x >> y) {
            out.x = x;
            out.y = y;
            return true;
        }
    }

    return false;
}

std::vector<Point> load_points_xy(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open file: " + path);

    std::vector<Point> pts;
    pts.reserve(1 << 20);

    std::string line;
    while (std::getline(in, line)) {
        Point p;
        if (parse_point_line(line, p)) {
            pts.push_back(p);
        }
    }

    if (pts.empty()) throw std::runtime_error("No points loaded from: " + path);
    return pts;
}