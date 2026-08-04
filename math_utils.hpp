#pragma once
#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>
#include <unordered_map>
#include <bvh/v2/vec.h>
#include <bvh/v2/bbox.h>
#include <CDT.h>
namespace BVH = bvh::v2;
using Scalar = double;
using Vec3 = BVH::Vec<Scalar, 3>;
using BBox = BVH::BBox<Scalar, 3>;
std::ostream& operator<<(std::ostream& os, const Vec3& v);
size_t hashCombine(size_t seed, int64_t v);
double computeMeshEpsilon(const BBox& boxA, const BBox& boxB);
template <class T>
std::vector<T> flattenVector(const std::vector<std::vector<T>>& v) {
    std::vector<T> out;
    size_t total = 0;
    for (const auto& row : v) total += row.size();
    out.reserve(total);
    for (const auto& row : v)
        for (const auto& x : row)
            out.push_back(x);
    return out;
}
template <class T>
int sign(const T& num, double eps) {
    if (num < -eps) return -1;
    if (num >  eps) return  1;
    return 0;
}
template <class T>
int uniqueSignIndex(const std::array<T,3>& dist, double eps) {
    int s0 = sign(dist[0], eps);
    int s1 = sign(dist[1], eps);
    int s2 = sign(dist[2], eps);
    if (s0 != s1 && s0 == s2) return 0;
    if (s1 != s0 && s1 == s2) return 1;
    if (s2 != s0 && s2 == s1) return 2;
    return -1;
}
class SpatialGrid3D {
private:
    struct Key {
        int64_t x, y, z;
        bool operator==(const Key& o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<int64_t>{}(k.x);
            h = hashCombine(h, k.y);
            h = hashCombine(h, k.z);
            return h;
        }
    };
    double cellSize;
    double sqrEps;
    std::unordered_map<Key, std::vector<size_t>, KeyHash> grid;
    std::vector<Vec3> points;
public:
    explicit SpatialGrid3D(double eps) : cellSize(eps), sqrEps(eps * eps) {}
    size_t getOrAdd(const Vec3& p) {
        int64_t cx = static_cast<int64_t>(std::floor(p[0] / cellSize));
        int64_t cy = static_cast<int64_t>(std::floor(p[1] / cellSize));
        int64_t cz = static_cast<int64_t>(std::floor(p[2] / cellSize));
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    Key neighborKey{cx + dx, cy + dy, cz + dz};
                    auto it = grid.find(neighborKey);
                    if (it != grid.end()) {
                        for (size_t idx : it->second) {
                            Vec3 diff = points[idx] - p;
                            double sqrDist = BVH::dot(diff, diff);
                            if (sqrDist < sqrEps) {
                                return idx;
                            }
                        }
                    }
                }
            }
        }
        size_t newIdx = points.size();
        points.push_back(p);
        grid[{cx, cy, cz}].push_back(newIdx);
        return newIdx;
    }
    const std::vector<Vec3>& getUniquePoints() const { return points; }
    void clear() { grid.clear(); points.clear(); }
};
class SpatialGrid2D {
private:
    struct Key {
        int64_t x, y;
        bool operator==(const Key& o) const {
            return x == o.x && y == o.y;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<int64_t>{}(k.x);
            return hashCombine(h, k.y);
        }
    };
    double cellSize;
    double sqrEps;
    std::unordered_map<Key, std::vector<CDT::VertInd>, KeyHash> grid;
    std::vector<CDT::V2d<double>> points;
public:
    explicit SpatialGrid2D(double eps) : cellSize(eps), sqrEps(eps * eps) {}
    CDT::VertInd getOrAdd(const CDT::V2d<double>& p) {
        int64_t cx = static_cast<int64_t>(std::floor(p.x / cellSize));
        int64_t cy = static_cast<int64_t>(std::floor(p.y / cellSize));
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                Key neighborKey{cx + dx, cy + dy};
                auto it = grid.find(neighborKey);
                if (it != grid.end()) {
                    for (CDT::VertInd idx : it->second) {
                        double dx_p = points[idx].x - p.x;
                        double dy_p = points[idx].y - p.y;
                        if ((dx_p * dx_p + dy_p * dy_p) < sqrEps) {
                            return idx;
                        }
                    }
                }
            }
        }
        CDT::VertInd newIdx = static_cast<CDT::VertInd>(points.size());
        points.push_back(p);
        grid[{cx, cy}].push_back(newIdx);
        return newIdx;
    }
    const std::vector<CDT::V2d<double>>& getUniquePoints() const { return points; }
    void clear() { grid.clear(); points.clear(); }
};