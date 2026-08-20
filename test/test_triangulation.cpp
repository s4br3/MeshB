#include "triangulation.hpp"
#include "mesh_types.hpp"
#include <cassert>
#include <iostream>
#include <cmath>
constexpr double EPS = 1e-6;
void test_buildSubdividedEdges() {
    std::vector<Vec2> initialPts = {
        {0.0, 0.0},
        {2.0, 2.0},
        {0.0, 2.0},
        {2.0, 0.0}
    };
    std::vector<std::pair<size_t, size_t>> segs = {
        {0, 1},
        {2, 3}
    };
    std::vector<Vec2> outUniquePts;
    std::vector<EdgeKey> outCDTEdges;
    buildSubdividedEdges(initialPts, segs, outUniquePts, outCDTEdges, EPS);
    assert(outUniquePts.size() == 5);
    assert(outCDTEdges.size() == 4);
}
void test_triangulate() {
    SpatialGrid3D nodeGrid(EPS);
    Vec3 normal{0.0, 0.0, 1.0};
    Vec3 origin{0.0, 0.0, 0.0};
    ProjectionFrame frame = computeSharedFrame(normal, origin);
    PolyLine boundary = {
        {Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}},
        {Vec3{1.0, 0.0, 0.0}, Vec3{1.0, 1.0, 0.0}},
        {Vec3{1.0, 1.0, 0.0}, Vec3{0.0, 1.0, 0.0}},
        {Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 0.0}}
    };
    PolyLine cuts = {
        {Vec3{0.5, 0.0, 0.0}, Vec3{0.5, 1.0, 0.0}}
    };
    std::vector<Triangle> result;
    triangulate(boundary, cuts, frame, nodeGrid, EPS, result);
    assert(result.size() >= 4);
    bool foundBottomCut = false;
    bool foundTopCut = false;
    const auto& nodes = nodeGrid.getUniquePoints();
    for (const auto& p : nodes) {
        if (std::abs(p.x - 0.5) < EPS && std::abs(p.y) < EPS && std::abs(p.z) < EPS) {
            foundBottomCut = true;
        }
        if (std::abs(p.x - 0.5) < EPS && std::abs(p.y - 1.0) < EPS && std::abs(p.z) < EPS) {
            foundTopCut = true;
        }
    }
    assert(foundBottomCut);
    assert(foundTopCut);
}
void test_cutTriangles() {
    SpatialGrid3D nodeGrid(EPS);
    Vec3 normal{0.0, 0.0, 1.0};
    Vec3 origin{0.0, 0.0, 0.0};
    ProjectionFrame frame = computeSharedFrame(normal, origin);
    size_t n0 = nodeGrid.getOrAdd({0.0, 0.0, 0.0});
    size_t n1 = nodeGrid.getOrAdd({2.0, 0.0, 0.0});
    size_t n2 = nodeGrid.getOrAdd({0.0, 2.0, 0.0});
    std::vector<Triangle> tris = { Triangle(n0, n1, n2) };
    PolyLine cuts = {
        {Vec3{1.0, -1.0, 0.0}, Vec3{1.0, 3.0, 0.0}}
    };
    std::vector<Triangle> result;
    cutTriangles(tris, nodeGrid.getUniquePoints(), cuts, frame, nodeGrid, EPS, result);
    assert(result.size() >= 3);
    const auto& outNodes = nodeGrid.getUniquePoints();
    size_t nodeCount = outNodes.size();
    for (const Triangle& tri : result) {
        assert(tri.v[0] < nodeCount);
        assert(tri.v[1] < nodeCount);
        assert(tri.v[2] < nodeCount);
    }
}
int main() {
    test_buildSubdividedEdges();
    test_triangulate();
    test_cutTriangles();
    std::cout << "test_triangulation (up to cutTriangles) passed successfully.\n";
    return 0;
}