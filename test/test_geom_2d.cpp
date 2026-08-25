#include "mesh_types.hpp"
#include "geom_2d.cpp"
#include <cassert>
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
    std::vector<size_t> outs;
    Edge e1{{0.0, 0.0}, {2.0, 2.0}};
    Edge e2{{0.0, 2.0}, {2.0, 0.0}};
    intersect2DAllPoints(e1, e2, grid, outs, EPS);
    auto uniquePts = grid.getUniquePoints();
    assert(uniquePts.size() == 1);
    assert(std::abs(uniquePts[0].x - 1.0) < EPS);
    assert(std::abs(uniquePts[0].y - 1.0) < EPS);
    grid.clear();
    outs.clear();
    Edge e3{{0.0, 0.0}, {2.0, 0.0}};
    Edge e4{{0.0, 1.0}, {2.0, 1.0}};
    intersect2DAllPoints(e3, e4, grid, outs, EPS);
    auto uniquePts2 = grid.getUniquePoints();
    assert(uniquePts2.empty());
}
int main(){
    test_computeSharedFrame();
    test_intersect2DAllPoints();
}