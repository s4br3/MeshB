#include "geom_3d.hpp"
#include "geom_2d.hpp"
bool pointInTriangle(const Vec3& p, const TriVerts& vertices, double eps) {
    Vec3 v0 = vertices[2] - vertices[0];
    Vec3 v1 = vertices[1] - vertices[0];
    Vec3 v2 = p - vertices[0];
    double dot00 = dot(v0, v0);
    double dot01 = dot(v0, v1);
    double dot02 = dot(v0, v2);
    double dot11 = dot(v1, v1);
    double dot12 = dot(v1, v2);
    double denom = (dot00 * dot11 - dot01 * dot01);
    if (std::fabs(denom) <= eps) return false;
    double u = (dot11 * dot02 - dot01 * dot12)/denom;
    double v = (dot00 * dot12 - dot01 * dot02)/denom;
    return (u >= -eps) && (v >= -eps) && (u + v <= 1 + eps);
}
std::array<double, 3> distancesToPlane(const TriVerts& vertices, const Vec3& n2, const Vec3& centre2) {
    const double p2Dist = dot(n2, centre2);
    return {dot(n2, vertices[0]) - p2Dist,
            dot(n2, vertices[1]) - p2Dist,
            dot(n2, vertices[2]) - p2Dist};
}
bool boundingBoxOverlap(const BBox& box1, const BBox& box2, double eps) {
    auto intervalsOverlap = [eps](double a0, double a1, double b0, double b1) {
        return !(a1 < b0 - eps || b1 < a0 - eps);
    };
    return (intervalsOverlap(box1.min[0], box1.max[0], box2.min[0], box2.max[0]) &&
            intervalsOverlap(box1.min[1], box1.max[1], box2.min[1], box2.max[1]) &&
            intervalsOverlap(box1.min[2], box1.max[2], box2.min[2], box2.max[2]));
}
bool coplanar(const Vec3& n1, const Vec3& centre1,
              const Vec3& n2, const Vec3& centre2, double eps) {
    constexpr double ANGLE_EPS = 1e-6; 
    return std::abs(dot(n1, n2)) >= 1.0 - ANGLE_EPS && 
           std::abs(dot(n1, centre1) - dot(n1, centre2)) < eps;
}
std::optional<std::pair<Vec3, Vec3>> segmentAOnB(
    const TriVerts& vertices,
    const std::array<double, 3>& distances,
    const double eps)
{
    std::vector<Vec3> pts;
    pts.reserve(3);
    for (int i = 0; i < 3; ++i) {
        int j = (i + 1) % 3;
        double di = distances[i];
        double dj = distances[j];
        if (std::abs(di) <= eps) {
            pts.push_back(vertices[i]);
        } 
        else if ((di > eps && dj < -eps) || (di < -eps && dj > eps)) {
            double t = di / (di - dj);
            pts.push_back(vertices[i] + (vertices[j] - vertices[i]) * t);
        }
    }
    if (pts.size() >= 2) {
        return std::make_pair(pts[0], pts[1]);
    } else if (pts.size() == 1) {
        return std::make_pair(pts[0], pts[0]);
    }
    return std::nullopt;
}
void findCycles(const PolyLine& edges, SpatialGrid3D& nodeGrid, std::vector<PolyLine>& outLoops) {
    if (edges.empty()) return;
    std::unordered_map<size_t, std::vector<size_t>> adj;
    for (const auto& seg : edges) {
        size_t u = nodeGrid.getOrAdd(seg.first);
        size_t v = nodeGrid.getOrAdd(seg.second);
        if (u != v) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<size_t> toRemove;
        for (const auto& [node, neighbors] : adj) {
            if (neighbors.size() <= 1) toRemove.push_back(node);
        }
        for (size_t node : toRemove) {
            auto it = adj.find(node);
            if (it == adj.end()) continue;
            for (size_t nb : it->second) {
                auto nbIt = adj.find(nb);
                if (nbIt != adj.end()) {
                    auto& vec = nbIt->second;
                    vec.erase(std::remove(vec.begin(), vec.end(), node), vec.end());
                }
            }
            adj.erase(it);
            changed = true;
        }
    }
    const auto& pts = nodeGrid.getUniquePoints();
    while (!adj.empty()) {
        size_t start = adj.begin()->first;
        std::vector<size_t> path;
        size_t curr = start;
        while (true) {
            path.push_back(curr);
            if (adj.find(curr) == adj.end() || adj[curr].empty()) break;
            size_t next = adj[curr].back();
            adj[curr].pop_back();
            auto& nextAdj = adj[next];
            auto it = std::find(nextAdj.begin(), nextAdj.end(), curr);
            if (it != nextAdj.end()) nextAdj.erase(it);
            if (adj[curr].empty()) adj.erase(curr);
            if (nextAdj.empty()) adj.erase(next);
            auto pathIt = std::find(path.begin(), path.end(), next);
            if (pathIt != path.end()) {
                std::vector<size_t> cycle(pathIt, path.end());
                if (cycle.size() >= 3) {
                    PolyLine loop;
                    for (size_t i = 0; i < cycle.size(); ++i) {
                        loop.push_back({pts[cycle[i]], pts[cycle[(i + 1) % cycle.size()]]});
                    }
                    outLoops.push_back(loop);
                }
                path.erase(pathIt, path.end());
                if (path.empty()) break;
                curr = path.back(); 
                path.pop_back(); 
            } else {
                curr = next;
            }
        }
    }
}
bool isCentroidInHole(
    const Vec3& centre,
    const std::vector<PolyLine>& allHoles, 
    const ProjectionFrame& frame) 
{
    auto p = frame.to2D(centre);
    Vec2 pt = {p.x, p.y};
    
    std::vector<std::pair<Vec2, Vec2>> edges2D;
    for (const PolyLine& hole : allHoles) {
        for (const auto& edge : hole) {
            auto p1 = frame.to2D(edge.first);
            auto p2 = frame.to2D(edge.second);
            edges2D.push_back({ {p1.x, p1.y}, {p2.x, p2.y} });
        }
    }
    return isPointInsidePolygon(pt, edges2D);
}