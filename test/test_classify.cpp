#define _USE_MATH_DEFINES
#include "classify.hpp"
#include "bvh_collisions.hpp"
#include "mesh_types.hpp"
#include "vector.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
MeshData createTetrahedron() {
    MeshData mesh;
    mesh.nodes = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    mesh.triangles = {
        {0, 2, 1},
        {0, 1, 3},
        {0, 3, 2},
        {1, 2, 3}
    };
    mesh.tags = {0, 0, 0, 0};
    recomputeMeshData(mesh); 
    return mesh;
}

void testSolidAngle() {
    Vec3 p{0.0, 0.0, 0.0};
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    Vec3 c{0.0, 0.0, 1.0};
    double omega = solidAngle(p, a, b, c);
    assert(std::abs(omega - (M_PI / 2.0)) < 1e-6 && "Solid angle calculation failed!");
    std::cout << "[PASS] testSolidAngle\n";
}

void testGWN() {
    MeshData mesh = createTetrahedron();
    double eps = 1e-6;
    BVH bvh = buildMeshBVH(mesh, eps);
    std::vector<GWNNodeData> gwnData = buildGWNData(mesh, bvh);
    Vec3 pInside{0.2, 0.2, 0.2};
    double wnInside = evaluateGWN(pInside, mesh, bvh, gwnData);
    assert(std::abs(wnInside - 1.0) < 1e-4 && "GWN inside should be ~1.0");
    Vec3 pOutside{2.0, 2.0, 2.0};
    double wnOutside = evaluateGWN(pOutside, mesh, bvh, gwnData);
    assert(std::abs(wnOutside - 0.0) < 1e-4 && "GWN outside should be ~0.0");

    std::cout << "[PASS] testGWN\n";
}

void testClassifyFace() {
    MeshData mesh = createTetrahedron();
    double eps = 1e-5;
    BVH bvh = buildMeshBVH(mesh, eps);
    Vec3 triCentreSame{0.2, 0.2, 0.0};
    Vec3 triNormalSame{0.0, 0.0, -1.0}; 
    FaceClass clsSame = classifyFace(bvh, mesh, triCentreSame, triNormalSame, eps);
    assert(clsSame == FaceClass::CoplanarSame && "Expected CoplanarSame");
    Vec3 triNormalOpp{0.0, 0.0, 1.0};
    FaceClass clsOpp = classifyFace(bvh, mesh, triCentreSame, triNormalOpp, eps);
    assert(clsOpp == FaceClass::CoplanarOpp && "Expected CoplanarOpp");
    Vec3 triCentreOut{5.0, 5.0, 5.0};
    Vec3 triNormalOut{0.0, 1.0, 0.0};
    FaceClass clsOut = classifyFace(bvh, mesh, triCentreOut, triNormalOut, eps);
    assert(clsOut == FaceClass::Outside && "Expected Outside");
    std::cout << "[PASS] testClassifyFace\n";
}

int main() {
    std::cout << "Running tests for classify.cpp...\n";
    testSolidAngle();
    testGWN();
    testClassifyFace();
    std::cout << "All classify tests passed successfully!\n";
    return 0;
}