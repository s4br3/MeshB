#include "triangulation.hpp"
#include "mesh_types.hpp"
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
    double EPS = 1e-6;
    SpatialGrid2D grid(EPS);
    std::vector<CDT::VertInd> outs;
    CDT::V2d<double> A{0.0, 0.0}, B{2.0, 2.0};
    CDT::V2d<double> C{0.0, 2.0}, D{2.0, 0.0};
    intersect2DAllPoints(A, B, C, D, grid, outs, EPS);
    auto uniquePts = grid.getUniquePoints();
    assert(uniquePts.size() == 1);
    assert(std::abs(uniquePts[0].x - 1.0) < EPS);
    assert(std::abs(uniquePts[0].y - 1.0) < EPS);
    grid.clear();
    outs.clear();
    CDT::V2d<double> E{0.0, 0.0}, F{2.0, 0.0};
    CDT::V2d<double> G{0.0, 1.0}, H{2.0, 1.0};
    intersect2DAllPoints(E, F, G, H, grid, outs, EPS);
    auto uniquePts2 = grid.getUniquePoints();
    assert(uniquePts2.empty());
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
    
    // Populate nodeGrid and create initial triangle referencing node indices
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
void test_findCycles() {
    SpatialGrid3D nodeGrid(EPS);
    Vec3 p1{0.0, 0.0, 0.0};
    Vec3 p2{1.0, 0.0, 0.0};
    Vec3 p3{1.0, 1.0, 0.0};
    Vec3 p4{0.0, 1.0, 0.0};
    PolyLine edges = {
        {p1, p2},
        {p2, p3},
        {p3, p4},
        {p4, p1}
    };
    std::vector<PolyLine> cycles;
    findCycles(edges, nodeGrid, cycles);
    assert(cycles.size() == 1);
    assert(cycles[0].size() == 4);
    
    Vec3 t1{2.0, 2.0, 0.0};
    Vec3 t2{3.0, 2.0, 0.0};
    Vec3 t3{2.0, 3.0, 0.0};
    edges.push_back({t1, t2});
    edges.push_back({t2, t3});
    edges.push_back({t3, t1});
    
    cycles.clear();
    findCycles(edges, nodeGrid, cycles);
    assert(cycles.size() == 2);
    
    Vec3 d1{-1.0, -1.0, 0.0};
    Vec3 d2{-2.0, -2.0, 0.0};
    edges.push_back({p1, d1}); 
    edges.push_back({d1, d2});
    
    cycles.clear();
    findCycles(edges, nodeGrid, cycles);
    assert(cycles.size() == 2); 
}
void test_isCentroidInHole() {
    Vec3 normal{0.0, 0.0, 1.0};
    Vec3 origin{0.0, 0.0, 0.0};
    ProjectionFrame frame = computeSharedFrame(normal, origin);
    Vec3 p1{-1.0, -1.0, 0.0};
    Vec3 p2{1.0, -1.0, 0.0};
    Vec3 p3{1.0, 1.0, 0.0};
    Vec3 p4{-1.0, 1.0, 0.0};
    PolyLine hole = {
        {p1, p2}, 
        {p2, p3}, 
        {p3, p4}, 
        {p4, p1}
    };
    std::vector<PolyLine> allHoles = {hole};
    Vec3 inside_pt{0.0, 0.0, 0.0};
    assert(isCentroidInHole(inside_pt, allHoles, frame) == true);
    Vec3 outside_pt{2.0, 0.0, 0.0};
    assert(isCentroidInHole(outside_pt, allHoles, frame) == false);
    Vec3 offset_inside_pt{0.5, 0.4, 0.0};
    assert(isCentroidInHole(offset_inside_pt, allHoles, frame) == true);
    Vec3 h2_p1{3.0, 3.0, 0.0};
    Vec3 h2_p2{5.0, 3.0, 0.0};
    Vec3 h2_p3{5.0, 5.0, 0.0};
    Vec3 h2_p4{3.0, 5.0, 0.0};
    PolyLine hole2 = { {h2_p1, h2_p2}, {h2_p2, h2_p3}, {h2_p3, h2_p4}, {h2_p4, h2_p1} };
    allHoles.push_back(hole2);
    Vec3 inside_hole2{4.0, 4.0, 0.0};
    assert(isCentroidInHole(inside_hole2, allHoles, frame) == true);
}
int main() {
    test_computeSharedFrame();
    test_intersect2DAllPoints();
    test_buildSubdividedEdges();
    test_triangulate();
    test_cutTriangles();
    test_findCycles();
    test_isCentroidInHole();
    std::cout << "test_triangulation (up to cutTriangles) passed successfully.\n";
    return 0;
}