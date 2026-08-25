#include "classify.hpp"
#include "geom_3d.hpp"
#include <vector>
#include <algorithm>
bool isMeshClosed(const MeshData& mesh) {
    std::vector<std::pair<size_t, size_t>> edges;
    edges.reserve(mesh.triangles.size() * 3);
    for (const auto& tri : mesh.triangles) {
        for (int i = 0; i < 3; ++i) {
            size_t v1 = tri.v[i];
            size_t v2 = tri.v[(i + 1) % 3];
            edges.push_back({std::min(v1, v2), std::max(v1, v2)});
        }
    }
    std::sort(edges.begin(), edges.end());
    size_t i = 0;
    while (i < edges.size()) {
        size_t count = 1;
        while (i + 1 < edges.size() && edges[i] == edges[i+1]) {
            count++;
            i++;
        }
        if (count != 2) return false; 
        i++;
    }
    return true;
}
double solidAngle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 A = a - p;
    Vec3 B = b - p;
    Vec3 C = c - p;
    double det = dot(A, cross(B, C));
    double al = std::sqrt(dot(A, A));
    double bl = std::sqrt(dot(B, B));
    double cl = std::sqrt(dot(C, C));
    double div = al * bl * cl + dot(A, B) * cl + dot(B, C) * al + dot(C, A) * bl;
    return 2.0 * std::atan2(det, div);
}
std::vector<GWNNodeData> buildGWNData(const MeshData& mesh, const BVH& bvh) {
    std::vector<GWNNodeData> data(bvh.nodes.size());
    for (int i = (int)bvh.nodes.size() - 1; i >= 0; --i) {
        const Node& node = bvh.nodes[i];
        if (node.isLeaf()) {
            Vec3 vArea{0, 0, 0};
            Vec3 center{0, 0, 0};
            double totalArea = 0.0;
            
            for (size_t j = 0; j < node.primCount; ++j) {
                size_t primId = bvh.primIds[node.firstId + j];
                const Triangle& tri = mesh.triangles[primId];
                Vec3 a = mesh.nodes[tri.v[0]];
                Vec3 b = mesh.nodes[tri.v[1]];
                Vec3 c = mesh.nodes[tri.v[2]];
                Vec3 n = cross(b - a, c - a);
                double area = std::sqrt(dot(n, n)) * 0.5;
                vArea = vArea + n * 0.5;
                center = center + (a + b + c) / 3.0 * area;
                totalArea += area;
            }
            data[i].vectorArea = vArea;
            if (totalArea > 1e-12) {
                data[i].center = center / totalArea;
            } else {
                data[i].center = center;
            }
        } else {
            size_t left = node.firstId;
            size_t right = node.firstId + 1;
            
            data[i].vectorArea = data[left].vectorArea + data[right].vectorArea;
            double areaL = std::sqrt(dot(data[left].vectorArea, data[left].vectorArea));
            double areaR = std::sqrt(dot(data[right].vectorArea, data[right].vectorArea));
            double sumArea = areaL + areaR;
            if (sumArea > 1e-12) {
                data[i].center = (data[left].center * areaL + data[right].center * areaR) / sumArea;
            } else {
                data[i].center = (data[left].center + data[right].center) * 0.5;
            }
        }
    }
    return data;
}
double evaluateGWN(const Vec3& p, const MeshData& mesh, const BVH& bvh, const std::vector<GWNNodeData>& gwnData) {
    if (bvh.nodes.empty()) return 0.0;
    double w = 0.0;
    std::vector<size_t> stack;
    stack.reserve(64);
    stack.push_back(0);
    while (!stack.empty()) {
        size_t nodeIdx = stack.back();
        stack.pop_back();
        const Node& node = bvh.nodes[nodeIdx];
        const GWNNodeData& data = gwnData[nodeIdx];
        Vec3 diag = node.boundingBox.max - node.boundingBox.min;
        double size = std::sqrt(dot(diag, diag));
        Vec3 r = data.center - p;
        double dist = r.length();
        if (size / dist < 0.5 && !node.isLeaf()) {
            double r3 = dist * dist * dist;
            if (r3 > 1e-12) {
                w += dot(data.vectorArea, r) / r3;
            }
        } else if (node.isLeaf()) {
            for (size_t j = 0; j < node.primCount; ++j) {
                size_t primId = bvh.primIds[node.firstId + j];
                const Triangle& tri = mesh.triangles[primId];
                Vec3 a = mesh.nodes[tri.v[0]];
                Vec3 b = mesh.nodes[tri.v[1]];
                Vec3 c = mesh.nodes[tri.v[2]];
                w += solidAngle(p, a, b, c);
            }
        } else {
            stack.push_back(node.firstId);
            stack.push_back(node.firstId + 1);
        }
    }
    return w / (4.0 * M_PI);
}
FaceClass classifyFace(
    const BVH& targetBVH,
    const MeshData& targetMesh,
    const Vec3& triCentre,
    const Vec3& triNormal,
    double eps)
{
    if (!targetBVH.nodes.empty()) {
        std::vector<size_t> stack;
        stack.reserve(32);
        stack.push_back(0);
        BBox ptBox;
        ptBox.extend(triCentre - Vec3{eps, eps, eps});
        ptBox.extend(triCentre + Vec3{eps, eps, eps});
        while (!stack.empty()) {
            size_t nodeIdx = stack.back();
            stack.pop_back();
            const Node& node = targetBVH.nodes[nodeIdx];
            
            if (!boundingBoxOverlap(ptBox, node.getBbox(), eps)) continue;
            
            if (node.isLeaf()) {
                for (size_t i = 0; i < node.primCount; ++i) {
                    size_t primId = targetBVH.primIds[node.firstId + i];
                    const Triangle& targetTri = targetMesh.triangles[primId];
                    Vec3 targetNorm = targetTri.normal(targetMesh.nodes);
                    double distToPlane = std::abs(::dot(targetNorm, triCentre - targetTri.centre(targetMesh.nodes)));
                    
                    if (distToPlane <= eps * 10.0) {
                        Vec3 projPt = triCentre - targetNorm * ::dot(targetNorm, triCentre - targetMesh.nodes[targetTri.v[0]]);
                        if (pointInTriangle(projPt, targetTri.getVertices(targetMesh.nodes), eps * 10.0)) {
                            double dotN = ::dot(triNormal, targetNorm);
                            return (dotN > 0.0) ? FaceClass::CoplanarSame : FaceClass::CoplanarOpp;
                        }
                    }
                }
            } else {
                stack.push_back(node.firstId);
                stack.push_back(node.firstId + 1);
            }
        }
    }
    return FaceClass::Outside; 
}