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
    const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D,
    SpatialGrid2D& grid,
    std::vector<size_t>& outs,
    double eps)
{
    outs.clear();
    double x12 = B.x - A.x, y12 = B.y - A.y;
    double x34 = D.x - C.x, y34 = D.y - C.y;
    double x13 = C.x - A.x, y13 = C.y - A.y;
    double denom = x12 * y34 - y12 * x34;
    auto addIdUnique = [&](size_t id) {
        for (auto existing : outs) {
            if (existing == id) return;
        }
        outs.push_back(id);
    };
    if (std::abs(denom) >= eps) {
        double t = (x13 * y34 - y13 * x34) / denom;
        double u = (x13 * y12 - y13 * x12) / denom;
        if (t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps) {
            if (std::abs(t) <= eps) {
                addIdUnique(grid.getOrAdd(A));
            } else if (std::abs(t - 1.0) <= eps) {
                addIdUnique(grid.getOrAdd(B));
            }
            else if (std::abs(u) <= eps) {
                addIdUnique(grid.getOrAdd(C));
            } else if (std::abs(u - 1.0) <= eps) {
                addIdUnique(grid.getOrAdd(D));
            }
            else {
                addIdUnique(grid.getOrAdd( {A.x + t * x12, A.y + t * y12} ));
            }
        }
        return;
    }
    if (std::abs(x13 * y12 - y13 * x12) > eps) return;
    if (pointOnSegment(C, A, B, eps)) addIdUnique(grid.getOrAdd(C));
    if (pointOnSegment(D, A, B, eps)) addIdUnique(grid.getOrAdd(D));
    if (pointOnSegment(A, C, D, eps)) addIdUnique(grid.getOrAdd(A));
    if (pointOnSegment(B, C, D, eps)) addIdUnique(grid.getOrAdd(B));
}
bool pointOnSegment(
    const Vec2& P,
    const Vec2& A,
    const Vec2& B,
    double eps)
{
    double vx = B.x - A.x, vy = B.y - A.y;
    double wx = P.x - A.x, wy = P.y - A.y;
    double cross = vx * wy - vy * wx;
    if (std::abs(cross) > eps) return false;
    auto inRange = [eps](double v, double a, double b) {
        return v >= std::min(a,b) - eps && v <= std::max(a,b) + eps;
    };
    return inRange(P.x, A.x, B.x) && inRange(P.y, A.y, B.y);
}
bool isPointInsidePolygon(const Vec2& pt, const std::vector<std::pair<Vec2, Vec2>>& edges) {
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