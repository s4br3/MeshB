#include "math_utils.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
constexpr double EPS = 1e-6;
void test_basic_math_helpers() {
    size_t seed = 42;
    size_t h1 = hashCombine(seed, 100);
    size_t h2 = hashCombine(seed, 100);
    size_t h3 = hashCombine(seed, 101);
    assert(h1 == h2);
    assert(h1 != h3);
    BBox b1({Vec3{0,0,0}, Vec3{10,0,0}});
    BBox b2({Vec3{0,0,0}, Vec3{0,20,0}});
    double meshEps = computeMeshEpsilon(b1, b2);
    assert(std::abs(meshEps - 20.0 * 1e-7) < 1e-12);
    std::vector<std::vector<int>> nested = {{1, 2}, {3, 4, 5}, {6}};
    auto flat = flattenVector(nested);
    assert(flat.size() == 6);
    assert(flat[3] == 4);
    assert(sign(0.001, EPS) == 1);
    assert(sign(-0.001, EPS) == -1);
    assert(sign(1e-7, EPS) == 0);
    assert(uniqueSignIndex(std::array<double, 3>{-1.0, 1.0, 1.0}, EPS) == 0);
    assert(uniqueSignIndex(std::array<double, 3>{1.0, -1.0, 1.0}, EPS) == 1);
    assert(uniqueSignIndex(std::array<double, 3>{-1.0, -1.0, 1.0}, EPS) == 2);
    assert(uniqueSignIndex(std::array<double, 3>{1.0, 1.0, 1.0}, EPS) == -1);
}
void test_spatial_grid_3d() {
    SpatialGrid3D grid(0.1);
    size_t idx0 = grid.getOrAdd(Vec3{0.0, 0.0, 0.0});
    size_t idx1 = grid.getOrAdd(Vec3{0.01, 0.01, 0.01});
    size_t idx2 = grid.getOrAdd(Vec3{5.0, 5.0, 5.0});
    assert(idx0 == idx1);
    assert(idx0 != idx2);
    assert(grid.getUniquePoints().size() == 2);
    grid.clear();
    assert(grid.getUniquePoints().empty());
}
void test_spatial_grid_2d() {
    SpatialGrid2D grid(0.1);
    auto idx0 = grid.getOrAdd(CDT::V2d<double>{0.0, 0.0});
    auto idx1 = grid.getOrAdd(CDT::V2d<double>{0.02, 0.02});
    auto idx2 = grid.getOrAdd(CDT::V2d<double>{2.0, 2.0});
    assert(idx0 == idx1);
    assert(idx0 != idx2);
    assert(grid.getUniquePoints().size() == 2);
}
int main() {
    test_basic_math_helpers();
    test_spatial_grid_3d();
    test_spatial_grid_2d();
    std::cout << "test_math_utils passed successfully.\n";
    return 0;
}