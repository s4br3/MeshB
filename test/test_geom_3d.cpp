#include "mesh_types.hpp"
#include "geom_3d.hpp"
#include "geom_2d.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string_view>
constexpr double EPS = 1e-6;
constexpr bool is_near(double a, double b, double eps = EPS) {
    return std::abs(a - b) < eps;
}
void run_test(std::string_view name, void(*test_fn)()) {
    test_fn();
    std::cout << "[PASS] " << name << "\n";
}
void test_bbox() {
    BBox b1;
    b1.extend(Vec3{0.0, 0.0, 0.0});
    b1.extend(Vec3{2.0, 4.0, 4.0});
    assert(is_near(b1.surfaceArea(), 64.0));
    BBox b2({Vec3{1.0, 1.0, 1.0}, Vec3{3.0, 3.0, 3.0}});
    assert(boundingBoxOverlap(b1, b2, EPS));
    BBox b3({Vec3{10.0, 10.0, 10.0}, Vec3{12.0, 12.0, 12.0}});
    assert(!boundingBoxOverlap(b1, b3, EPS));
}
void test_bvh_construction() {
    std::vector<BBox> boxes;
    std::vector<Vec3> centres;
    boxes.reserve(10);
    centres.reserve(10);
    
    for (int i = 0; i < 10; ++i) {
        BBox b;
        Vec3 c{i * 2.0, 0.0, 0.0};
        b.extend(c - Vec3{0.5, 0.5, 0.5});
        b.extend(c + Vec3{0.5, 0.5, 0.5});
        boxes.push_back(b);
        centres.push_back(c);
    }
    BVH bvh(boxes, centres, EPS);
    assert(!bvh.nodes.empty());
    assert(bvh.primIds.size() == 10);
    assert(bvh.nodes[0].getBbox().min.x <= -0.5);
    assert(bvh.nodes[0].getBbox().max.x >= 18.5);
    BVH emptyBVH({}, {}, EPS);
    assert(emptyBVH.nodes.empty());
}
void test_point_in_triangle() {
    TriVerts tri = {Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}};
    assert(pointInTriangle(Vec3{0.5, 0.5, 0.0}, tri, EPS));
    assert(pointInTriangle(Vec3{0.0, 0.0, 0.0}, tri, EPS)); 
    assert(!pointInTriangle(Vec3{2.0, 2.0, 0.0}, tri, EPS));
}
void test_distances_to_plane() {
    TriVerts tri = {Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}};
    Vec3 normal{0.0, 0.0, 1.0};
    Vec3 planeCentre{0.0, 0.0, 1.0};
    auto dists = distancesToPlane(tri, normal, planeCentre);
    for (double d : dists) {
        assert(is_near(d, -1.0));
    }
}
void test_segment_straddle() {
    TriVerts straddleTri = {Vec3{0.0, 0.0, -1.0}, Vec3{2.0, 0.0, 1.0}, Vec3{0.0, 2.0, 1.0}};
    std::array<double, 3> straddleDists = {-1.0, 1.0, 1.0};
    auto segRes = segmentAOnB(straddleTri, straddleDists, EPS);
    assert(segRes.has_value());
    assert(is_near(segRes->first.z, 0.0));
    assert(is_near(segRes->second.z, 0.0));
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

int main(){
    run_test("Bounding Box", test_bbox);
    run_test("BVH Construction", test_bvh_construction);
    run_test("Point in Triangle", test_point_in_triangle);
    run_test("Distances to Plane", test_distances_to_plane);
    run_test("Segment Straddle", test_segment_straddle);
    test_findCycles();
    test_isCentroidInHole();
    return 0;
}