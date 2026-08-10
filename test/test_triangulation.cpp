#include "triangulation.hpp"
#include "mesh_types.hpp"
#include "math_utils.hpp"
#include <cassert>
#include <iostream>
#include <cmath>
constexpr double EPS = 1e-6;
void test_computeSharedFrame() {
    Vec3 normal{0.0, 0.0, 1.0};
    Vec3 origin{0.0, 0.0, 0.0};
    ProjectionFrame frame = computeSharedFrame(normal, origin);
    assert(std::abs(dot(frame.u, frame.v)) < EPS);
    assert(std::abs(dot(frame.u, normal)) < EPS);
    assert(std::abs(dot(frame.v, normal)) < EPS);
    assert(std::abs(frame.u.length() - 1.0) < EPS);
    assert(std::abs(frame.v.length() - 1.0) < EPS);
}
void test_intersect2DAllPoints() {
    SpatialGrid2D grid(EPS);
    std::vector<CDT::VertInd> outs;
    CDT::V2d<double> A{0.0, 0.0}, B{2.0, 2.0};
    CDT::V2d<double> C{0.0, 2.0}, D{2.0, 0.0};
    intersect2DAllPoints(A, B, C, D, grid, outs, EPS);
    assert(outs.size() == 1);
    auto pts = grid.getUniquePoints();
    assert(std::abs(pts[outs[0]].x - 1.0) < EPS);
    assert(std::abs(pts[outs[0]].y - 1.0) < EPS);
    grid.clear();
    outs.clear();
    CDT::V2d<double> E{0.0, 0.0}, F{2.0, 0.0};
    CDT::V2d<double> G{0.0, 1.0}, H{2.0, 1.0};
    intersect2DAllPoints(E, F, G, H, grid, outs, EPS);
    assert(outs.empty());
}
void test_buildSubdividedEdges() {
    std::vector<CDT::V2d<double>> initialPts = {
        {0.0, 0.0},
        {2.0, 2.0},
        {0.0, 2.0},
        {2.0, 0.0}
    };
    std::vector<std::pair<size_t, size_t>> segs = {
        {0, 1},
        {2, 3}
    };
    std::vector<CDT::V2d<double>> outUniquePts;
    std::vector<CDT::Edge> outCDTEdges;
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
    std::vector<Triangle> result = triangulate(boundary, cuts, frame, nodeGrid, true, EPS);
    assert(result.size() >= 4);
    auto uniquePoints = nodeGrid.getUniquePoints();
    bool foundBottomCut = false;
    bool foundTopCut = false;
    for (const auto& p : uniquePoints) {
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
    std::vector<Vec3> nodes = {
        {0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {0.0, 2.0, 0.0}
    };
    std::vector<Triangle> tris = { Triangle(0, 1, 2) };
    PolyLine cuts = {
        {Vec3{1.0, -1.0, 0.0}, Vec3{1.0, 3.0, 0.0}}
    };
    std::vector<Triangle> result = cutTriangles(tris, nodes, cuts, frame, nodeGrid, true, EPS);
    assert(result.size() >= 3);
    size_t gridNodeCount = nodeGrid.getUniquePoints().size();
    for (const Triangle& tri : result) {
        assert(tri.v[0] < gridNodeCount);
        assert(tri.v[1] < gridNodeCount);
        assert(tri.v[2] < gridNodeCount);
    }
}
int main() {
    test_computeSharedFrame();
    test_triangulate();
    test_cutTriangles();
    std::cout << "test_triangulation (up to cutTriangles) passed successfully.\n";
    return 0;
}