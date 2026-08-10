#include "vector.hpp"
#include <cassert>
#include <cmath>
constexpr double EPS = 1e-6;
int main() {
    Vec3 v1{1.0, 2.0, 3.0};
    Vec3 v2{4.0, 5.0, 6.0};
    assert((v1 + v2).x == 5.0);
    assert((v1 - v2).x == -3.0);
    assert((v2 / 2).equals(Vec3{2.0, 2.5, 3.0},EPS));
    assert(v1.length() == std::sqrt(14.0));
    assert((v1.normalize()).equals(Vec3{1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0), 3.0 / std::sqrt(14.0)}, EPS));
    assert(std::abs(dot(v1, v2) - 32.0) < EPS);
    Vec3 c = cross(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    assert(c.equals(Vec3{0.0, 0.0, 1.0}, EPS));
        assert(std::abs(dot(v1, v2) - 32.0) < EPS);
    assert(v1[0] == 1.0 && v1[1] == 2.0 && v1[2] == 3.0);
    v1[0] = 4.0;
    assert(v1[0] == 4.0);
    assert(getLargestAxis(Vec3{1.0, 5.0, 2.0}) == 1);
    Vec2 v2d1{1.0, 2.0};
    Vec2 v2d2{4.0, 5.0};
    assert((v2d1 + v2d2).equals(Vec2{5.0, 7.0}, EPS));
    assert((v2d1 - v2d2).equals(Vec2{-3.0, -3.0}, EPS));
    assert((v2d2 / 2).equals(Vec2{2.0, 2.5}, EPS));
    assert(v2d1.length() == std::sqrt(5.0));
    assert(v2d1.normalize().equals(Vec2{1.0 / std::sqrt(5.0), 2.0 / std::sqrt(5.0)}, EPS));
    double c2d = cross2D(Vec2{1.0, 0.0}, Vec2{0.0, 1.0});
    assert(std::abs(c2d - 1.0) < EPS);
    return 0;
}