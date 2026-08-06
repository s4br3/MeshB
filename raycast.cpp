#include "raycast.hpp"
#include "bvh_collisions.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <array>
static bool intersectRayTri(const Vec3& orig, const Vec3& dir, const Triangle& t, const std::vector<Vec3>& nodes, double eps) {
    Vec3 e1 = nodes[t.v[1]] - nodes[t.v[0]];
    Vec3 e2 = nodes[t.v[2]] - nodes[t.v[0]];
    Vec3 pvec = BVH::cross(dir, e2);
    double det = BVH::dot(e1, pvec);
    if (std::abs(det) < eps) return false;
    double inv_det = 1.0 / det;
    Vec3 tvec = orig - nodes[t.v[0]];
    double u = BVH::dot(tvec, pvec) * inv_det;
    if (u < eps || u > 1.0 + eps) return false;
    Vec3 qvec = BVH::cross(tvec, e1);
    double v = BVH::dot(dir, qvec) * inv_det;
    if (v < eps || u + v > 1.0 + eps) return false;
    double t_hit = BVH::dot(e2, qvec) * inv_det;
    return t_hit > -eps; 
}
static bool rayBoxIntersect(const Vec3& orig, const Vec3& invDir, const BBox& box) {
    double tmin = -std::numeric_limits<double>::infinity();
    double tmax = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 3; ++i) {
        double t1 = (box.min[i] - orig[i]) * invDir[i];
        double t2 = (box.max[i] - orig[i]) * invDir[i];
        tmin = std::max(tmin, std::min(t1, t2));
        tmax = std::min(tmax, std::max(t1, t2));
    }
    return tmax >= std::max(0.0, tmin);
}
int countRayIntersections(
    const Bvh& bvh, const MeshData& mesh, 
    const Vec3& orig, const Vec3& direction, 
    double eps) 
{
    if (bvh.nodes.empty()) return 0;
    Vec3 dir = BVH::normalize(direction);
    Vec3 invDir = {1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]};
    int hitCount = 0;
    std::vector<size_t> stack;
    stack.reserve(32);
    stack.push_back(0);
    while (!stack.empty()) {
        size_t nodeIdx = stack.back();
        stack.pop_back();
        const Node& node = bvh.nodes[nodeIdx];
        if (!rayBoxIntersect(orig, invDir, node.get_bbox())) {
            continue;
        }
        if (node.is_leaf()) {
            size_t prim_id = bvh.prim_ids[node.index.first_id()]; 
            const Triangle& tri = mesh.triangles[prim_id];
            if (intersectRayTri(orig, dir, tri, mesh.nodes, eps)) {
                hitCount++;
            }
        } else {
            stack.push_back(node.index.first_id());
            stack.push_back(node.index.first_id() + 1);
        }
    }
    return hitCount;
}
bool isInsideMesh(const Bvh& bvh, const MeshData& mesh, const Vec3& point, double eps) {
    static const std::array<Vec3, 3> perturbedRayDirs = {{
        {1.0,      0.001337, 0.002468},
        {0.003579, 1.0,      0.004680},
        {0.005791, 0.006802, 1.0     }
    }};
    int insideVotes = 0;
    for (const Vec3& dir : perturbedRayDirs) {
        int count = countRayIntersections(bvh, mesh, point, dir, eps);
        if ((count % 2) != 0) {
            insideVotes++;
        }
    }
    return insideVotes >= 2;
}
FaceClass classifyFace(
    const Bvh& targetBvh,
    const MeshData& targetMesh,
    const Vec3& triCentre,
    const Vec3& triNormal,
    double eps)
{
    if (!targetBvh.nodes.empty()) {
        std::vector<size_t> stack;
        stack.reserve(32);
        stack.push_back(0);
        BBox ptBox = BBox::make_empty();
        ptBox.extend(triCentre - Vec3{eps, eps, eps});
        ptBox.extend(triCentre + Vec3{eps, eps, eps});
        while (!stack.empty()) {
            size_t nodeIdx = stack.back();
            stack.pop_back();
            const Node& node = targetBvh.nodes[nodeIdx];
            if (!boundingBoxOverlap(ptBox, node.get_bbox(), eps)) continue;
            if (node.is_leaf()) {
                size_t primId = targetBvh.prim_ids[node.index.first_id()];
                const Triangle& targetTri = targetMesh.triangles[primId];
                Vec3 targetNorm = targetTri.normal(targetMesh.nodes);
                double distToPlane = std::abs(BVH::dot(targetNorm, triCentre - targetTri.centre(targetMesh.nodes)));
                if (distToPlane <= eps * 10.0) {
                    Vec3 projPt = triCentre - targetNorm * BVH::dot(targetNorm, triCentre - targetMesh.nodes[targetTri.v[0]]);
                    if (pointInTriangle(projPt, targetTri.getVertices(targetMesh.nodes), eps * 10.0)) {
                        double dotN = BVH::dot(triNormal, targetNorm);
                        return (dotN > 0.0) ? FaceClass::CoplanarSame : FaceClass::CoplanarOpp;
                    }
                }
            } else {
                stack.push_back(node.index.first_id());
                stack.push_back(node.index.first_id() + 1);
            }
        }
    }
    bool inside = isInsideMesh(targetBvh, targetMesh, triCentre, eps);
    return inside ? FaceClass::Inside : FaceClass::Outside;
}