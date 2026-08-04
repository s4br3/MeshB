#pragma once
#include "mesh_types.hpp"
#include "bvh_collisions.hpp"
static bool intersectRayTri(
    const Vec3& orig, const Vec3& dir,
    const Triangle& t, const std::vector<Vec3>& nodes,
    double eps);
static bool rayBoxIntersect(const Vec3& orig, const Vec3& invDir, const BBox& box);
int countRayIntersections(
    const Bvh& bvh, const MeshData& mesh,
    const Vec3& startPoint, const Vec3& direction,
    double eps);
bool isInsideMesh(const Bvh& bvh, const MeshData& mesh,
    const Vec3& point,
    double eps);
enum class FaceClass {Outside, Inside, CoplanarSame, CoplanarOpp};
FaceClass classifyFace(
    const Bvh& targetBvh, const MeshData& targetMesh,
    const Vec3& triCenter, const Vec3& triNormal,
    double eps);