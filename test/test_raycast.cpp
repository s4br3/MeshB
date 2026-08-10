#include "raycast.hpp"
#include "mesh_types.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
constexpr double EPS = 1e-6;
MeshData createUnitCubeMesh() {
    MeshData mesh;
    mesh.nodes = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
    };
    mesh.triangles = {
        {0, 2, 1}, {0, 3, 2},
        {4, 5, 6}, {4, 6, 7},
        {0, 1, 5}, {0, 5, 4},
        {2, 3, 7}, {2, 7, 6},
        {0, 4, 7}, {0, 7, 3},
        {1, 2, 6}, {1, 6, 5}
    };
    recomputeMeshData(mesh);
    return mesh;
}
void test_raycast_operations() {
    MeshData cube = createUnitCubeMesh();
    std::vector<BBox> boxes;
    for (const auto& tri : cube.triangles) {
        boxes.push_back(tri.bounds(cube.nodes));
    }
    BVH bvh(boxes, cube.centres, EPS);
    int hitsInsideOut = countRayIntersections(bvh, cube, Vec3{0.5, 0.5, 0.5}, Vec3{0, 0, 1}, EPS);
    assert(hitsInsideOut == 1);
    int hitsOutsideThrough = countRayIntersections(bvh, cube, Vec3{0.5, 0.5, -2.0}, Vec3{0, 0, 1}, EPS);
    assert(hitsOutsideThrough == 2);
    int hitsMiss = countRayIntersections(bvh, cube, Vec3{5.0, 5.0, 5.0}, Vec3{0, 0, 1}, EPS);
    assert(hitsMiss == 0);
    assert(isInsideMesh(bvh, cube, Vec3{0.5, 0.5, 0.5}, EPS));
    assert(!isInsideMesh(bvh, cube, Vec3{2.0, 0.5, 0.5}, EPS));
    assert(!isInsideMesh(bvh, cube, Vec3{-1.0, -1.0, -1.0}, EPS));
    FaceClass fcOutside = classifyFace(bvh, cube, Vec3{2.0, 2.0, 2.0}, Vec3{0, 0, 1}, EPS);
    assert(fcOutside == FaceClass::Outside);
    FaceClass fcInside = classifyFace(bvh, cube, Vec3{0.5, 0.5, 0.5}, Vec3{0, 0, 1}, EPS);
    assert(fcInside == FaceClass::Inside);
    FaceClass fcCoplanarSame = classifyFace(bvh, cube, Vec3{0.5, 0.5, 1.0}, Vec3{0, 0, 1}, EPS);
    assert(fcCoplanarSame == FaceClass::CoplanarSame);
    FaceClass fcCoplanarOpp = classifyFace(bvh, cube, Vec3{0.5, 0.5, 1.0}, Vec3{0, 0, -1}, EPS);
    assert(fcCoplanarOpp == FaceClass::CoplanarOpp);
}
int main() {
    test_raycast_operations();
    std::cout << "test_raycast passed successfully.\n";
    return 0;
}