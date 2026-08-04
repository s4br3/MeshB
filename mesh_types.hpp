#pragma once
#include "math_utils.hpp"
#include <vector>
#include <array>
#include <unordered_map>
using PolyLine = std::vector<std::pair<Vec3, Vec3>>;
struct Triangle {
    std::array<size_t, 3> v;
    Vec3 normal(const std::vector<Vec3>& nodes) const { 
        return BVH::normalize(BVH::cross(nodes[v[1]] - nodes[v[0]], nodes[v[2]] - nodes[v[0]])); 
    }
    Vec3 centre(const std::vector<Vec3>& nodes) const { 
        return (nodes[v[0]] + nodes[v[1]] + nodes[v[2]]) * (1.0 / 3.0); 
    }
    BBox bounds(const std::vector<Vec3>& nodes) const {
        BBox box = BBox::make_empty();
        box.extend(nodes[v[0]]); 
        box.extend(nodes[v[1]]); 
        box.extend(nodes[v[2]]);
        return box;
    }
};
struct MeshData {
    std::vector<Vec3> nodes;
    std::vector<Triangle> triangles;
    std::vector<size_t> tags;
    std::vector<Vec3> centres;
    std::vector<Vec3> normals;
    MeshData() = default;
    MeshData(const MeshData&) = default;
    MeshData& operator=(const MeshData&) = default;
    MeshData(MeshData&&) noexcept = default;
    MeshData& operator=(MeshData&&) noexcept = default;
};
struct CollisionContext {
    MeshData meshDataA;
    MeshData meshDataB;
    double eps;
    std::unordered_map<size_t, PolyLine> NCAcoords;
    std::unordered_map<size_t, std::vector<PolyLine>> CAcoords;
    std::unordered_map<size_t, std::vector<size_t>> Atris;
    std::unordered_map<size_t, PolyLine> NCBcoords;
    std::unordered_map<size_t, std::vector<PolyLine>> CBcoords;
    std::unordered_map<size_t, std::vector<size_t>> Btris;
};
struct Connection {
    std::unordered_map<size_t, std::vector<Vec3>> meshAIntersections;
    std::unordered_map<size_t, std::vector<Vec3>> meshBIntersections;
    std::unordered_map<size_t, std::vector<size_t>> aToBConnections;
    std::unordered_map<size_t, std::vector<size_t>> bToAConnections;
};
struct ProjectionFrame {
    Vec3 origin;
    Vec3 u, v;
    CDT::V2d<double> to2D(const Vec3& p) const {
        Vec3 rel = p - origin;
        return {BVH::dot(rel, u), BVH::dot(rel, v)};
    }
    Vec3 to3D(const CDT::V2d<double>& p) const {
        return origin + u * p.x + v * p.y;
    }
};