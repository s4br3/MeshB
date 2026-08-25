#include "geom_2d.hpp"
ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin) {
    ProjectionFrame frame;
    frame.origin = origin;
    Vec3 absNormal = {std::abs(normal[0]), std::abs(normal[1]), std::abs(normal[2])};
    Vec3 helper = {0, 0, 0};
    if (absNormal[0] <= absNormal[1] && absNormal[0] <= absNormal[2]) {
        helper[0] = 1.0; 
    } else if (absNormal[1] <= absNormal[0] && absNormal[1] <= absNormal[2]) {
        helper[1] = 1.0;
    } else {
        helper[2] = 1.0;
    }
    frame.u = cross(helper, normal).normalize();
    frame.v = cross(normal, frame.u).normalize();
    return frame;
}
void intersect2DAllPoints(
    const Edge& e1, const Edge& e2, SpatialGrid2D& grid,
    std::vector<size_t>& outs,
    double eps)
{
    outs.clear();
    Vec2 v12 = e1.second - e1.first;
    Vec2 v34 = e2.second - e2.first;
    Vec2 v13 = e2.first - e1.first;
    double denom = cross2D(v12, v34);
    double cross1334 = cross2D(v13, v34);
    double cross1312 = cross2D(v13, v12);
    auto addIdUnique = [&](size_t id) {
        for (auto existing : outs) {
            if (existing == id) return;
        }
        outs.push_back(id);
    };
    if (std::abs(denom) >= eps) {
        double t = cross1334 / denom;
        double u = cross1312 / denom;
        if (t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps) {
            if (std::abs(t) <= eps) {
                addIdUnique(grid.getOrAdd(e1.first));
            } else if (std::abs(t - 1.0) <= eps) {
                addIdUnique(grid.getOrAdd(e1.second));
            }
            else if (std::abs(u) <= eps) {
                addIdUnique(grid.getOrAdd(e2.first));
            } else if (std::abs(u - 1.0) <= eps) {
                addIdUnique(grid.getOrAdd(e2.second));
            }
            else {
                addIdUnique(grid.getOrAdd(e1.first + v12 * t));
            }
        }
        return;
    }
    if (std::abs(cross1312) > eps) return;
    if (pointOnSegment(e2.first, e1, eps)) addIdUnique(grid.getOrAdd(e2.first));
    if (pointOnSegment(e2.second, e1, eps)) addIdUnique(grid.getOrAdd(e2.second));
    if (pointOnSegment(e1.first, e2, eps)) addIdUnique(grid.getOrAdd(e1.first));
    if (pointOnSegment(e1.second, e2, eps)) addIdUnique(grid.getOrAdd(e1.second));
}
bool pointOnSegment(const Vec2& p, const Edge& e, double eps)
{
    Vec2 v = e.second - e.first;
    Vec2 w = p - e.first;
    double cross = cross2D(v, w);
    if (std::abs(cross) > eps) return false;
    auto inRange = [eps](double v, double a, double b) {
        return v >= std::min(a,b) - eps && v <= std::max(a,b) + eps;
    };
    return inRange(p.x, e.first.x, e.second.x) && inRange(p.y, e.first.y, e.second.y);
}
bool isPointInsidePolygon(const Vec2& pt, const std::vector<Edge>& edges) {
    bool inside = false;
    for (const auto& edge : edges) {
        const Vec2& v1 = edge.first;
        const Vec2& v2 = edge.second;
        if (((v1.y > pt.y) != (v2.y > pt.y)) &&
            (pt.x < (v2.x - v1.x) * (pt.y - v1.y) / (v2.y - v1.y) + v1.x)) 
        {
            inside = !inside;
        }
    }
    return inside;
}