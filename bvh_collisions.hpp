#pragma once
#include "mesh_types.hpp"
#include <bvh/v2/bvh.h>
#include <bvh/v2/node.h>
#include <optional>
#include <bvh/v2/default_builder.h>
using Node = BVH::Node<Scalar, 3>;
using Bvh  = BVH::Bvh<Node>;
bool pointInTriangle(
    const Vec3& p,
    const Triangle& t, const std::vector<Vec3>& nodes,
    double eps);
std::array<double, 3> distancesToPlane(
    const Triangle& t1, const std::vector<Vec3>& nodes,
    const Vec3& n2, const Vec3& centre2);
bool boundingBoxOverlap(const BBox& box1, const BBox& box2, double eps);
bool coplanar(
    const Triangle& t1, const Vec3& n1, const Vec3& centre1,
    const Triangle& t2, const Vec3& n2, const Vec3& centre2,
    double eps);
std::pair<Vec3, Vec3> segmentAOnB(const Triangle& t1, const std::array<double, 3>& distances, double eps);
std::optional<std::pair<Vec3, Vec3>> findIntersectionPointsNC(
    const Triangle& t1, const Vec3& n1, const Vec3& c1,
    const Triangle& t2, const Vec3& n2, const Vec3& c2,
    const double eps);
std::vector<Vec3> findIntersectionPointsC(
    const Triangle& t1,
    const Triangle& t2, const Vec3& n2, const Vec3& c2,
    const double eps);
void findAllCollisions(const Bvh& bvh1, const Bvh& bvh2,
    const MeshData& mesh1, const MeshData& mesh2,
    std::unordered_map<size_t, PolyLine>& NCAcoords, std::unordered_map<size_t, std::vector<PolyLine>>& CAcoords,
    std::unordered_map<size_t, std::vector<size_t>>& Atris,
    std::unordered_map<size_t, PolyLine>& NCBcoords, std::unordered_map<size_t, std::vector<PolyLine>>& CBcoords,
    std::unordered_map<size_t, std::vector<size_t>>& Btris,
    const double eps);
Bvh buildMeshBVH(const MeshData& mesh);
CollisionContext detectCollisions(const MeshData& meshA, const MeshData& meshB);