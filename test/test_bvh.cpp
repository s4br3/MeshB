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
    run_test("Non-Coplanar Intersection", test_non_coplanar_intersection);
    run_test("Coplanar Intersection", test_coplanar_intersection);
    std::cout << "\nAll BVH tests passed successfully.\n";
    return 0;
}