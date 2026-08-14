#include "vector.hpp"
std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}
std::ostream& operator<<(std::ostream& os, const Vec2& v) {
    return os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
}
double dot(const Vec3& lhs, const Vec3& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
double dot(const Vec2& lhs, const Vec2& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y;
}
Vec3 cross(const Vec3& lhs, const Vec3& rhs) {
    return Vec3{
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x
    };
}
double cross2D(const Vec2& lhs, const Vec2& rhs) {
    return lhs.x * rhs.y - lhs.y * rhs.x;
}
int getLargestAxis(const Vec3& v) {
    double ax = std::abs(v.x);
    double ay = std::abs(v.y);
    double az = std::abs(v.z);
    if (ax >= ay && ax >= az) return 0;
    if (ay >= az) return 1;
    return 2;
}
int getLargestAxis(const Vec2& v){
    return std::abs(v.y) > std::abs(v.x);
}
Vec3 abs(const Vec3& v){
    Vec3 r = v;
    r.x = std::abs(r.x);
    r.y = std::abs(r.y);
    r.z = std::abs(r.z);
    return r;
}
Vec2 abs(const Vec2& v){
    Vec2 r = v;
    r.x = std::abs(r.x);
    r.y = std::abs(r.y);
    return r;
}