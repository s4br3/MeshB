#include "bvh_collisions.hpp"
#include <cstddef>
bool pointInTriangle(const Vec3& p, const TriVerts& vertices, double eps) {
    Vec3 v0 = vertices[2] - vertices[0];
    Vec3 v1 = vertices[1] - vertices[0];
    Vec3 v2 = p - vertices[0];
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
std::array<double, 3> distancesToPlane(const TriVerts& vertices, const Vec3& n2, const Vec3& centre2) {
    const double p2Dist = BVH::dot(n2, centre2);
    return {BVH::dot(n2, vertices[0]) - p2Dist,
            BVH::dot(n2, vertices[1]) - p2Dist,
            BVH::dot(n2, vertices[2]) - p2Dist};
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
std::optional<std::pair<Vec3, Vec3>> segmentAOnB(
    const TriVerts& vertices,
    const std::array<double, 3>& distances,
    const double eps)
{
    int unique = uniqueSignIndex(distances, eps);
    if (unique < 0) {
        return std::nullopt;
    }
    int a = (unique + 1) % 3;
    int b = (unique + 2) % 3;
    auto intersectOnEdge = [&](int iOther) -> Vec3 {
        double dOther = distances[iOther];
        double dUnique = distances[unique];
        double denom = (dOther - dUnique);
        if (std::fabs(denom) <= eps) {
            return vertices[unique];
        }
        double t = dOther / denom;
        return vertices[iOther] + (vertices[unique] - vertices[iOther]) * t;
    };
    Vec3 A = intersectOnEdge(a);
    Vec3 B = intersectOnEdge(b);
    return {{A, B}};
}
std::optional<std::pair<Vec3, Vec3>> findIntersectionPointsNC(
    const TriVerts& vertices1, const Vec3& n1, const Vec3& c1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps)
{
    std::array<double, 3> distAB = distancesToPlane(vertices1, n2, c2);
    std::array<double, 3> distBA = distancesToPlane(vertices2, n1, c1);
    Vec3 direction = BVH::cross(n1, n2);
    if (BVH::dot(direction, direction) <= eps * eps) {
        return std::nullopt;
    }
    size_t largestAxis = direction.get_largest_axis();
    std::optional<std::pair<Vec3, Vec3>> temp1 = segmentAOnB(vertices1, distAB, eps);
    std::optional<std::pair<Vec3, Vec3>> temp2 = segmentAOnB(vertices2, distBA, eps);
    if (!temp1 || !temp2){
        return std::nullopt;
    }
    auto [A, B] = *temp1;
    auto [C, D] = *temp2;
    double line1A = A[largestAxis], line1B = B[largestAxis];
    double line2A = C[largestAxis], line2B = D[largestAxis];
    if (line1A > line1B) { std::swap(line1A, line1B); std::swap(A, B); }
    if (line2A > line2B) { std::swap(line2A, line2B); std::swap(C, D); }
    Vec3 start = (line1A > line2A + eps) ? A : ((line2A > line1A + eps) ? C : ((A + C) * 0.5));
    Vec3 end   = (line1B < line2B - eps) ? B : ((line2B < line1B - eps) ? D : ((B + D) * 0.5));
    return std::make_pair(start, end);
}
std::vector<Vec3> findIntersectionPointsC(
    const TriVerts& vertices1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps)
{
    std::vector<Vec3> poly, out;
    poly.reserve(6);
    out.reserve(6);
    poly = {vertices1[0], vertices1[1], vertices1[2]};
    for (int i = 0; i < 3; ++i) {
        if (poly.empty()) break;
        const Vec3& e0 = vertices2[i];
        const Vec3& e1 = vertices2[(i + 1) % 3];
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
    CollisionContext& ctx)
{
    std::vector<std::pair<size_t, size_t>> stack;
    stack.reserve(128);
    stack.emplace_back(0, 0);
    while (!stack.empty()) {
        auto [node1Idx, node2Idx] = stack.back();
        stack.pop_back();
        const Node& node1 = bvh1.nodes[node1Idx];
        const Node& node2 = bvh2.nodes[node2Idx];
        
        if (!boundingBoxOverlap(node1.get_bbox(), node2.get_bbox(), ctx.eps)) {
            continue;
        }
        
        if (node1.is_leaf() && node2.is_leaf()) {
            size_t prim1Id = bvh1.prim_ids[node1.index.first_id()];
            size_t prim2Id = bvh2.prim_ids[node2.index.first_id()];
            const Triangle& t1 = ctx.meshDataA.triangles[prim1Id];
            const Triangle& t2 = ctx.meshDataB.triangles[prim2Id];
            const TriVerts& vs1 = t1.getVertices(ctx.meshDataA.nodes);
            const TriVerts& vs2 = t2.getVertices(ctx.meshDataB.nodes);
            const Vec3& n1 = ctx.meshDataA.normals[prim1Id];
            const Vec3& n2 = ctx.meshDataB.normals[prim2Id];
            const Vec3& c1 = ctx.meshDataA.centres[prim1Id];
            const Vec3& c2 = ctx.meshDataB.centres[prim2Id];            
            
            if (coplanar(n1, c1, n2, c2, ctx.eps)) {
                std::vector<Vec3> res = findIntersectionPointsC(vs1, vs2, n2, c2, ctx.eps);
                if (res.size() > 1) {
                    ctx.Atris[prim1Id].push_back(prim2Id);
                    ctx.Btris[prim2Id].push_back(prim1Id);
                    PolyLine temp;
                    for (size_t i = 0; i < res.size(); i++) {
                        temp.push_back({res[i], res[(i+1) % res.size()]});                        
                    }
                    ctx.CAcoords[prim1Id].push_back(temp);
                    ctx.CBcoords[prim2Id].push_back(temp);
                }
            } else {
                if (std::optional<std::pair<Vec3, Vec3>> optResult = findIntersectionPointsNC(vs1, n1, c1, vs2, n2, c2, ctx.eps)){
                    auto [startPt, endPt] = *optResult;
                    ctx.Atris[prim1Id].push_back(prim2Id);
                    ctx.Btris[prim2Id].push_back(prim1Id);
                    Vec3 diff = startPt - endPt;
                    if (BVH::dot(diff, diff) > ctx.eps * ctx.eps) {
                    ctx.NCAcoords[prim1Id].push_back({startPt, endPt});
                    ctx.NCBcoords[prim2Id].push_back({startPt, endPt});
                }
                    else {
                        ctx.NCAcoords[prim1Id].push_back({startPt, startPt});
                        ctx.NCBcoords[prim2Id].push_back({startPt, startPt});
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
CollisionContext detectCollisions(const MeshData& meshA, const MeshData& meshB) {
    CollisionContext ctx = {meshA, meshB};
    Bvh bvhA = buildMeshBVH(ctx.meshDataA);
    Bvh bvhB = buildMeshBVH(ctx.meshDataB);
    ctx.eps = 1e-7; 
    if (!bvhA.nodes.empty() && !bvhB.nodes.empty()) {
        ctx.eps = computeMeshEpsilon(bvhA.nodes[0].get_bbox(), bvhB.nodes[0].get_bbox());
        findAllCollisions(bvhA, bvhB, ctx);
    }
    return ctx;
}
