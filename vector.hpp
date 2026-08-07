#pragma once
#include <iostream>
#include <stdexcept>
#include <array>
#include <cmath>
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    Vec3 normalize() const {
        return *this/length();
    }
    Vec3() = default;
    Vec3(double px, double py, double pz){
        x = px;
        y = py;
        z = pz;
    }
    Vec3(const std::array<double, 3>& p){
        x = p[0];
        y = p[1];
        z = p[2];
    }
    double& operator[](size_t i) {
        switch (i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Vector index out of range");
        }
    }
    const double& operator[](size_t i) const {
        switch (i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Vector index out of range");
        }
    }
    Vec3 operator+(const Vec3& rhs) const {
        Vec3 r;
        r.x = x + rhs.x;
        r.y = y + rhs.y;
        r.z = z + rhs.z;
        return r;
    }
    Vec3 operator-(const Vec3& rhs) const {
        Vec3 r;
        r.x = x - rhs.x;
        r.y = y - rhs.y;
        r.z = z - rhs.z;
        return r;
    }
    Vec3 operator*(const double rhs) const{
        Vec3 r;
        r.x = x * rhs;
        r.y = y * rhs;
        r.z = z * rhs;
        return r;
    }
    Vec3 operator/(const double rhs) const{
        Vec3 r;
        r.x = x / rhs;
        r.y = y / rhs;
        r.z = z / rhs;
        return r;
    }
};
struct Vec2 {
    double x = 0.0, y = 0.0;
    double length() const {
        return std::sqrt(x * x + y * y);
    }
    Vec2 normalize() const {
        return *this/length();
    }
    Vec2() = default;
    Vec2(double px, double py) {
        x = px;
        y = py;
    }
    Vec2(std::array<double, 2> p) {
        x = p[0];
        y = p[1];
    }
    double operator[](size_t i) const {
        switch (i) {
            case 0: return x;
            case 1: return y;
            default: throw std::out_of_range("Vector index out of range");
        }
    }
    Vec2 operator+(const Vec2& rhs) const {
        Vec2 r;
        r.x = x + rhs.x;
        r.y = y + rhs.y;
        return r;
    }
    Vec2 operator-(const Vec2& rhs) const {
        Vec2 r;
        r.x = x - rhs.x;
        r.y = y - rhs.y;
        return r;
    }
    Vec2 operator*(const double rhs) const{
        Vec2 r;
        r.x = x * rhs;
        r.y = y * rhs;
        return r;
    }
    Vec2 operator/(const double rhs) const{
        Vec2 r;
        r.x = x / rhs;
        r.y = y / rhs;
        return r;
    }

};
std::ostream& operator<<(std::ostream& os, const Vec3& v);
std::ostream& operator<<(std::ostream& os, const Vec2& v);
double dot(const Vec3& lhs, const Vec3& rhs);
Vec3 cross(const Vec3& lhs, const Vec3& rhs);
double cross2D(const Vec3& lhs, const Vec3& rhs);
int getLargestAxis(const Vec3& v);
int getLargestAxis(const Vec2& v);