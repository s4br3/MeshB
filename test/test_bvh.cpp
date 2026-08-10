#include "bvh.hpp"
#include "bvh_collisions.hpp"
#include "mesh_types.hpp"
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
void test_non_coplanar_intersection() {
    const std::vector<Vec3> vertexPool = {
        Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 2.0},
        Vec3{0.0, -1.0, 0.5}, Vec3{0.0, 1.0, 0.5}, Vec3{0.0, 0.0, 1.5}
    };
    Triangle t1{0, 1, 2};
    Triangle t2{3, 4, 5};
    auto ncRes = findIntersectionPointsNC(
        t1.getVertices(vertexPool), t1.normal(vertexPool), t1.centre(vertexPool),
        t2.getVertices(vertexPool), t2.normal(vertexPool), t2.centre(vertexPool),
        EPS
    );
    assert(ncRes.has_value());
    const bool validOrder1 = ncRes->first.equals(Vec3{0, 0, 0.5}, EPS) && ncRes->second.equals(Vec3{0, 0, 1.5}, EPS);
    const bool validOrder2 = ncRes->first.equals(Vec3{0, 0, 1.5}, EPS) && ncRes->second.equals(Vec3{0, 0, 0.5}, EPS);
    assert(validOrder1 || validOrder2);
}
void test_coplanar_intersection() {
    TriVerts cp1 = {Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}};
    TriVerts cp2 = {Vec3{0.5, 0.5, 0.0}, Vec3{2.5, 0.5, 0.0}, Vec3{0.5, 2.5, 0.0}};
    auto poly = findIntersectionPointsC(cp1, cp2, Vec3{0.0, 0.0, 1.0}, Vec3{1.0, 1.0, 0.0}, EPS);
    assert(poly.size() >= 3);
}
int main() {
    run_test("Bounding Box", test_bbox);
    run_test("BVH Construction", test_bvh_construction);
    run_test("Point in Triangle", test_point_in_triangle);
    run_test("Distances to Plane", test_distances_to_plane);
    run_test("Segment Straddle", test_segment_straddle);
    run_test("Non-Coplanar Intersection", test_non_coplanar_intersection);
    run_test("Coplanar Intersection", test_coplanar_intersection);
    std::cout << "\nAll BVH tests passed successfully.\n";
    return 0;
}