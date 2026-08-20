#include "bvh_collisions.hpp"
#include "math_utils.hpp"
#include "geom_3d.hpp"
#include <cstddef>
thread_local static std::vector<Vec3> s_polyBuffer;
thread_local static std::vector<Vec3> s_outBuffer;
thread_local static std::vector<Vec3> s_uniquePolyBuffer;
std::optional<std::pair<Vec3, Vec3>> findIntersectionPointsNC(
    const TriVerts& vertices1, const Vec3& n1, const Vec3& c1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps)
{
    std::array<double, 3> distAB = distancesToPlane(vertices1, n2, c2);
    std::array<double, 3> distBA = distancesToPlane(vertices2, n1, c1);
    Vec3 direction = cross(n1, n2);
    int largestAxis = getLargestAxis(direction);
    std::optional<std::pair<Vec3, Vec3>> temp1 = segmentAOnB(vertices1, distAB, eps);
    std::optional<std::pair<Vec3, Vec3>> temp2 = segmentAOnB(vertices2, distBA, eps);
    if (!temp1 || !temp2) {
        return std::nullopt;
    }
    auto [A, B] = *temp1;
    auto [C, D] = *temp2;
    double line1A = A[largestAxis], line1B = B[largestAxis];
    double line2A = C[largestAxis], line2B = D[largestAxis];
    if (line1A > line1B) { std::swap(line1A, line1B); std::swap(A, B); }
    if (line2A > line2B) { std::swap(line2A, line2B); std::swap(C, D); }
    if (line1B < line2A - eps || line2B < line1A - eps) {
        return std::nullopt;
    }
    Vec3 start = (line1A > line2A + eps) ? A : ((line2A > line1A + eps) ? C : ((A + C) * 0.5));
    Vec3 end   = (line1B < line2B - eps) ? B : ((line2B < line1B - eps) ? D : ((B + D) * 0.5));
    return std::make_pair(start, end);
}
const std::vector<Vec3>& findIntersectionPointsC(
    const TriVerts& vertices1,
    const TriVerts& vertices2, const Vec3& n2, const Vec3& c2,
    const double eps)
{
    s_polyBuffer.clear();
    s_outBuffer.clear();
    s_uniquePolyBuffer.clear();
    s_polyBuffer.assign(vertices1.begin(), vertices1.end());
    for (int i = 0; i < 3; ++i) {
        if (s_polyBuffer.empty()) break;
        const Vec3& e0 = vertices2[i];
        const Vec3& e1 = vertices2[(i + 1) % 3];
        Vec3 edgeDir = e1 - e0;
        Vec3 inwardNormal = cross(n2, edgeDir);
        double distC = dot(c2 - e0, inwardNormal);
        int s = sign(distC, eps);
        s_outBuffer.clear();
        Vec3 S = s_polyBuffer.back();
        double distS = dot(S - e0, inwardNormal);
        bool SInside = (distS * s) >= -eps;
        for (const Vec3& E : s_polyBuffer) {
            double distE = dot(E - e0, inwardNormal);
            bool EInside = (distE * s) >= -eps;
            if (EInside != SInside) {
                double denom = (distS - distE);
                if (std::fabs(denom) > eps) {
                    double t = distS / denom;
                    s_outBuffer.push_back(S + (E - S) * t);
                }
            }
            if (EInside) s_outBuffer.push_back(E);
            S = E;
            distS = distE;
            SInside = EInside;
        }
        s_polyBuffer.swap(s_outBuffer);
    }
    const double sqrEps = eps * eps;
    for (const auto& pt : s_polyBuffer) {
        bool isDuplicate = false;
        for (const auto& existing : s_uniquePolyBuffer) {
            Vec3 diff = pt - existing;
            if (dot(diff, diff) <= sqrEps) {
                isDuplicate = true;
                break;
            }
        }
        if (!isDuplicate) {
            s_uniquePolyBuffer.push_back(pt);
        }
    }
    return s_uniquePolyBuffer;
}
void findAllCollisions(
    const BVH& bvh1, const BVH& bvh2,
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
        if (!boundingBoxOverlap(node1.getBbox(), node2.getBbox(), ctx.eps)) {
            continue;
        }
        if (node1.isLeaf() && node2.isLeaf()) {
            for (size_t i = node1.firstId; i < node1.firstId + node1.primCount; ++i) {
                size_t prim1Id = bvh1.primIds[i];
                for (size_t j = node2.firstId; j < node2.firstId + node2.primCount; ++j) {
                    size_t prim2Id = bvh2.primIds[j];
                    const TriVerts& t1 = ctx.meshDataA.triangles[prim1Id].getVertices(ctx.meshDataA.nodes);
                    const TriVerts& t2 = ctx.meshDataB.triangles[prim2Id].getVertices(ctx.meshDataB.nodes);
                    const Vec3& n1 = ctx.meshDataA.normals[prim1Id];
                    const Vec3& n2 = ctx.meshDataB.normals[prim2Id];
                    const Vec3& c1 = ctx.meshDataA.centres[prim1Id];
                    const Vec3& c2 = ctx.meshDataB.centres[prim2Id];
                    if (coplanar(n1, c1, n2, c2, ctx.eps)) {
                        std::vector<Vec3> res = findIntersectionPointsC(t1, t2, n2, c2, ctx.eps);
                        if (res.size() > 1) {
                            ctx.Atris[prim1Id].push_back(prim2Id);
                            ctx.Btris[prim2Id].push_back(prim1Id);
                            PolyLine temp;
                            for (size_t k = 0; k < res.size(); k++) {
                                temp.push_back({res[k], res[(k + 1) % res.size()]});                        
                            }
                            ctx.CAcoords[prim1Id].push_back(temp);
                            ctx.CBcoords[prim2Id].push_back(temp);
                        }
                    } else {
                        if (std::optional<std::pair<Vec3, Vec3>> optResult = findIntersectionPointsNC(t1, n1, c1, t2, n2, c2, ctx.eps)) {
                            auto [startPt, endPt] = *optResult;
                            ctx.Atris[prim1Id].push_back(prim2Id);
                            ctx.Btris[prim2Id].push_back(prim1Id);
                            Vec3 diff = startPt - endPt;
                            if (dot(diff, diff) > ctx.eps * ctx.eps) {
                                ctx.NCAcoords[prim1Id].push_back({startPt, endPt});
                                ctx.NCBcoords[prim2Id].push_back({startPt, endPt});
                            }
                            else {
                                ctx.NCAcoords[prim1Id].push_back({startPt, startPt});
                                ctx.NCBcoords[prim2Id].push_back({startPt, startPt});
                            }
                        }
                    }
                }
            }
            continue;
        }
        if (node1.isLeaf()) {
            size_t l2 = node2.firstId;
            stack.emplace_back(node1Idx, l2);
            stack.emplace_back(node1Idx, l2 + 1);
        } else if (node2.isLeaf()) {
            size_t l1 = node1.firstId;
            stack.emplace_back(l1, node2Idx);
            stack.emplace_back(l1 + 1, node2Idx);
        } else {
            size_t l1 = node1.firstId; size_t r1 = l1 + 1;
            size_t l2 = node2.firstId; size_t r2 = l2 + 1;
            stack.emplace_back(l1, l2);
            stack.emplace_back(l1, r2);
            stack.emplace_back(r1, l2);
            stack.emplace_back(r1, r2);
        }
    }
}
BVH buildMeshBVH(const MeshData &mesh, double eps){
    std::vector<BBox> boxes;
    for (const Triangle& tri : mesh.triangles){
        boxes.push_back(tri.bounds(mesh.nodes));
    }
    return BVH(boxes, mesh.centres, eps);
}
CollisionContext detectCollisions(const MeshData& meshA, const MeshData& meshB) {
    CollisionContext ctx = {meshA, meshB, 0};
    BVH bvhA = buildMeshBVH(ctx.meshDataA, ctx.eps);
    BVH bvhB = buildMeshBVH(ctx.meshDataB, ctx.eps);
    if (!bvhA.nodes.empty() && !bvhB.nodes.empty()) {
        ctx.eps = computeMeshesEpsilon(bvhA.nodes[0].getBbox(), bvhB.nodes[0].getBbox());
        BVH bvhA = buildMeshBVH(ctx.meshDataA, ctx.eps);
        BVH bvhB = buildMeshBVH(ctx.meshDataB, ctx.eps);
        findAllCollisions(bvhA, bvhB, ctx);
    }
    return ctx;
}
