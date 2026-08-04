#include "bvh_collisions.hpp"
#include <cstddef>
bool pointInTriangle(const Vec3& p, const Triangle& t, const std::vector<Vec3>& nodes, double eps) {
    Vec3 v0 = nodes[t.v[2]] - nodes[t.v[0]];
    Vec3 v1 = nodes[t.v[1]] - nodes[t.v[0]];
    Vec3 v2 = p - nodes[t.v[0]];
    double dot00 = BVH::dot(v0, v0);
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
std::array<double, 3> distancesToPlane(const Triangle& t1, const std::vector<Vec3>& nodes1, const Vec3& n2, const Vec3& centre2) {
    const double p2Dist = BVH::dot(n2, centre2);
    return {BVH::dot(n2, nodes1[t1.v[0]]) - p2Dist,
            BVH::dot(n2, nodes1[t1.v[1]]) - p2Dist,
            BVH::dot(n2, nodes1[t1.v[2]]) - p2Dist};
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
    return std::abs(BVH::dot(n1, n2)) >= 1.0 - ANGLE_EPS && 
           std::abs(BVH::dot(n1, centre1) - BVH::dot(n1, centre2)) < eps;
}
std::pair<Vec3, Vec3> segmentAOnB(const Triangle& t1, const std::vector<Vec3>& nodes1, const std::array<double, 3>& distances, const double eps)
{
    int unique = uniqueSignIndex(distances, eps);
    if (unique < 0) {
        return { nodes1[t1.v[0]], nodes1[t1.v[1]] };
    }
    int a = (unique + 1) % 3;
    int b = (unique + 2) % 3;
    auto intersectOnEdge = [&](int iOther) -> Vec3 {
        double dOther = distances[iOther];
        double dUnique = distances[unique];
        double denom = (dOther - dUnique);
        if (std::fabs(denom) <= eps) {
            return nodes1[t1.v[unique]];
        }
        double t = dOther / denom;
        return nodes1[t1.v[iOther]] + (nodes1[t1.v[unique]] - nodes1[t1.v[iOther]]) * t;
    };
    Vec3 A = intersectOnEdge(a);
    Vec3 B = intersectOnEdge(b);
    return {A, B};
}
std::optional<std::pair<Vec3, Vec3>> findIntersectionPointsNC(
    const Triangle& t1, const std::vector<Vec3>& nodes1, const Vec3& n1, const Vec3& c1,
    const Triangle& t2, const std::vector<Vec3>& nodes2, const Vec3& n2, const Vec3& c2,
    const double eps)
{
    std::array<double, 3> distAB = distancesToPlane(t1, nodes1, n2, c2);
    std::array<double, 3> distBA = distancesToPlane(t2, nodes2, n1, c1);
    Vec3 direction = BVH::cross(n1, n2);
    if (BVH::dot(direction, direction) <= eps * eps) {
        return std::nullopt;
    }
    size_t largestAxis = direction.get_largest_axis();
    auto [A, B] = segmentAOnB(t1, nodes1, distAB, eps);
    auto [C, D] = segmentAOnB(t2, nodes2, distBA, eps);
    double line1A = A[largestAxis], line1B = B[largestAxis];
    double line2A = C[largestAxis], line2B = D[largestAxis];
    if (line1A > line1B) { std::swap(line1A, line1B); std::swap(A, B); }
    if (line2A > line2B) { std::swap(line2A, line2B); std::swap(C, D); }
    double startVal = std::max(line1A, line2A);
    double endVal   = std::min(line1B, line2B);
    if (startVal > endVal + eps) {
        return std::nullopt;
    }
    Vec3 start = (line1A >= line2A) ? A : C;
    Vec3 end   = (line1B <= line2B) ? B : D;
    return std::make_pair(start, end);
}
std::vector<Vec3> findIntersectionPointsC(
    const Triangle& t1, const std::vector<Vec3>& nodes1,
    const Triangle& t2, const std::vector<Vec3>& nodes2,
    const Vec3& n2, const Vec3& c2,
    const double eps)
{
    std::vector<Vec3> poly, out;
    poly.reserve(6);
    out.reserve(6);
    poly = { nodes1[t1.v[0]], nodes1[t1.v[1]], nodes1[t1.v[2]] };
    for (int i = 0; i < 3; ++i) {
        if (poly.empty()) break;
        const Vec3& e0 = nodes2[t2.v[i]];
        const Vec3& e1 = nodes2[t2.v[(i + 1) % 3]];
        Vec3 edgeDir = e1 - e0;
        Vec3 inwardNormal = BVH::cross(n2, edgeDir);
        double distC = BVH::dot(c2 - e0, inwardNormal);
        int s = sign(distC, eps);
        out.clear();
        out.reserve(6);
        Vec3 S = poly.back();
        for (size_t j = 0; j < poly.size(); ++j) {
            const Vec3& E = poly[j];
            double distE = BVH::dot(E - e0, inwardNormal);
            double distS = BVH::dot(S - e0, inwardNormal);
            bool EInside = (distE * s) >= -eps;
            bool SInside = (distS * s) >= -eps;
            if (EInside) {
                if (!SInside) {
                    double denom = (distS - distE);
                    if (std::fabs(denom) > eps) {
                        double t = distS / denom;
                        out.push_back(S + (E - S) * t);
                    }
                }
                out.push_back(E);
            } else if (SInside) {
                double denom = (distS - distE);
                if (std::fabs(denom) > eps) {
                    double t = distS / denom;
                    out.push_back(S + (E - S) * t);
                }
            }
            S = E;
        }
        poly.swap(out);
    }
    return poly;
}
Bvh buildMeshBVH(const MeshData& mesh) {
    if (mesh.triangles.empty()) return Bvh();
    std::vector<BBox> bboxes(mesh.triangles.size());
    for (size_t i = 0; i < mesh.triangles.size(); i++) {
        bboxes[i] = mesh.triangles[i].bounds(mesh.nodes);
    }
    typename BVH::DefaultBuilder<Node>::Config config;
    config.quality = BVH::DefaultBuilder<Node>::Quality::High;
    return BVH::DefaultBuilder<Node>::build(bboxes, mesh.centres, config);
}
void findAllCollisions(
    const Bvh& bvh1, const Bvh& bvh2,
    const MeshData& mesh1, 
    const MeshData& mesh2,
    std::unordered_map<size_t, PolyLine>& NCAcoords,
    std::unordered_map<size_t, std::vector<PolyLine>>& CAcoords,
    std::unordered_map<size_t, std::vector<size_t>>& Atris,
    std::unordered_map<size_t, PolyLine>& NCBcoords,
    std::unordered_map<size_t, std::vector<PolyLine>>& CBcoords,
    std::unordered_map<size_t, std::vector<size_t>>& Btris,
    const double eps)
{
    std::vector<std::pair<size_t, size_t>> stack;
    stack.reserve(128);
    stack.emplace_back(0, 0);
    while (!stack.empty()) {
        auto [node1Idx, node2Idx] = stack.back();
        stack.pop_back();
        const auto& node1 = bvh1.nodes[node1Idx];
        const auto& node2 = bvh2.nodes[node2Idx];
        
        if (!boundingBoxOverlap(node1.get_bbox(), node2.get_bbox(), eps)) {
            continue;
        }
        
        if (node1.is_leaf() && node2.is_leaf()) {
            size_t prim1Id = bvh1.prim_ids[node1.index.first_id()];
            size_t prim2Id = bvh2.prim_ids[node2.index.first_id()];
            const Triangle& t1 = mesh1.triangles[prim1Id];
            const Triangle& t2 = mesh2.triangles[prim2Id];
            const Vec3& n1 = t1.normal(mesh1.nodes);
            const Vec3& n2 = t2.normal(mesh2.nodes);
            const Vec3& c1 = mesh1.centres[prim1Id];
            const Vec3& c2 = mesh2.centres[prim2Id];            
            
            if (coplanar(t1, n1, c1, t2, n2, c2, eps)) {
                std::vector<Vec3> res = findIntersectionPointsC(t1, mesh1.nodes, t2, mesh2.nodes, n2, c2, eps);
                if (res.size() > 1) {
                    Atris.try_emplace(prim1Id, std::vector<size_t>{}).first->second.push_back(prim2Id);
                    Btris.try_emplace(prim2Id, std::vector<size_t>{}).first->second.push_back(prim1Id);
                    PolyLine temp;
                    for (size_t i = 0; i < res.size(); i++) {
                        temp.push_back({res[i], res[(i+1) % res.size()]});                        
                    }
                    CAcoords.try_emplace(prim1Id, std::vector<PolyLine>{}).first->second.push_back(temp);
                    CBcoords.try_emplace(prim2Id, std::vector<PolyLine>{}).first->second.push_back(temp);
                }
            } else {
                if (auto optResult = findIntersectionPointsNC(t1, mesh1.nodes, n1, c1, t2, mesh2.nodes, n2, c2, eps)){
                    auto [startPt, endPt] = *optResult;
                    Atris[prim1Id].push_back(prim2Id);
                    Btris[prim2Id].push_back(prim1Id);
                    Vec3 diff = startPt - endPt;
                    if (BVH::dot(diff, diff) > eps * eps) {
                    NCAcoords[prim1Id].push_back({startPt, endPt});
                    NCBcoords[prim2Id].push_back({startPt, endPt});
                }
                    else {
                        NCAcoords[prim1Id].push_back({startPt, startPt});
                        NCBcoords[prim2Id].push_back({startPt, startPt});
                    }
                }
            }
            continue;
        }
        
        if (node1.is_leaf()) {
            size_t l2 = node2.index.first_id();
            stack.emplace_back(node1Idx, l2);
            stack.emplace_back(node1Idx, l2 + 1);
        } else if (node2.is_leaf()) {
            size_t l1 = node1.index.first_id();
            stack.emplace_back(l1, node2Idx);
            stack.emplace_back(l1 + 1, node2Idx);
        } else {
            size_t l1 = node1.index.first_id(); size_t r1 = l1 + 1;
            size_t l2 = node2.index.first_id(); size_t r2 = l2 + 1;
            stack.emplace_back(l1, l2);
            stack.emplace_back(l1, r2);
            stack.emplace_back(r1, l2);
            stack.emplace_back(r1, r2);
        }
    }
}
CollisionContext detectCollisions(MeshData meshA, MeshData meshB) {
    CollisionContext ctx = {meshA, meshB};
    Bvh bvhA = buildMeshBVH(ctx.meshDataA);
    Bvh bvhB = buildMeshBVH(ctx.meshDataB);
    ctx.eps = 1e-7; 
    if (!bvhA.nodes.empty() && !bvhB.nodes.empty()) {
        ctx.eps = computeMeshEpsilon(bvhA.nodes[0].get_bbox(), bvhB.nodes[0].get_bbox());
        findAllCollisions(
            bvhA, bvhB, ctx.meshDataA, ctx.meshDataB,
            ctx.NCAcoords, ctx.CAcoords, ctx.Atris,
            ctx.NCBcoords, ctx.CBcoords, ctx.Btris,
            ctx.eps
        );
    }
    return ctx;
}
