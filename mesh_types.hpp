#pragma once
#include "math_utils.hpp"
#include <vector>
#include <array>
#include <unordered_map>

/**
* @brief Polyline defined as an ordered collection of 3D line segment endpoints.
*/
using PolyLine = std::vector<std::pair<Vec3, Vec3>>;

/**
* @struct Triangle
* @brief Represents a triangular element defined by indices referring to mesh nodes.
*/
struct Triangle {
    std::array<size_t, 3> v; /**< Array of vertex indices in the global node array */

    /**
    * @brief Default constructor for an uninitialized triangle.
    */
    Triangle() {}

    /**
    * @brief Constructs a `Triangle` using an array of node indices.
    * @param[in] vertices - Array of 3 node indices.
    */
    Triangle(const std::array<size_t, 3>& vertices) : v(vertices) {}

    /**
    * @brief Constructs a `Triangle` using 3 node index values.
    * @param[in] a - Index of first vertex.
    * @param[in] b - Index of second vertex.
    * @param[in] c - Index of third vertex.
    */
    Triangle(size_t a, size_t b, size_t c) : v{a, b, c} {}

    /**
    * @brief Computes the surface normal vector of the triangle.
    * @param[in] nodes - Global node list containing vertex positions.
    * @return Normalized surface normal vector (counter-clockwise orientation).
    */
    const Vec3 normal(const std::vector<Vec3>& nodes) const { 
        return BVH::normalize(BVH::cross(nodes[v[1]] - nodes[v[0]], nodes[v[2]] - nodes[v[0]])); 
    }

    /**
    * @brief Computes the centroid coordinate of the triangle.
    * @param[in] nodes - Global node list containing vertex positions.
    * @return Centroid coordinate in 3D space.
    */
    const Vec3 centre(const std::vector<Vec3>& nodes) const { 
        return (nodes[v[0]] + nodes[v[1]] + nodes[v[2]]) * (1.0 / 3.0); 
    }

    /**
    * @brief Calculates the axis-aligned bounding box (AABB) of the triangle.
    * @param[in] nodes - Global node list containing vertex positions.
    * @return Bounding box encompassing all three vertices.
    */
    BBox bounds(const std::vector<Vec3>& nodes) const {
        BBox box = BBox::make_empty();
        box.extend(nodes[v[0]]); 
        box.extend(nodes[v[1]]); 
        box.extend(nodes[v[2]]);
        return box;
    }

    /**
    * @brief Formatted stream output operator for a triangle given
    * @param[in] nodes - Global node list containing vertex positions
    */
};

/**
* @struct MeshData
* @brief Primary data container for 3D triangle surface meshes.
*/
struct MeshData {
    std::vector<Vec3> nodes;           /**< Array of 3D node coordinates */
    std::vector<Triangle> triangles;   /**< Element connectivity table */
    std::vector<size_t> tags;          /**< Element group or boundary tags */
    std::vector<Vec3> centres;         /**< Precomputed element centroid locations */
    std::vector<Vec3> normals;         /**< Precomputed element surface normals */

    MeshData() = default;
    MeshData(const MeshData&) = default;
    MeshData& operator=(const MeshData&) = default;
    MeshData(MeshData&&) noexcept = default;
    MeshData& operator=(MeshData&&) noexcept = default;
};

/**
* @struct CollisionContext
* @brief Stores surface intersection polyline results and element maps between colliding meshes.
*/
struct CollisionContext {
    MeshData meshDataA;                                            /**< First mesh data */
    MeshData meshDataB;                                            /**< Second mesh data */
    double eps;                                                    /**< Distance tolerance threshold */
    std::unordered_map<size_t, PolyLine> NCAcoords;                /**< Mapping of Non-coplanar intersection segments for mesh A */
    std::unordered_map<size_t, std::vector<PolyLine>> CAcoords;    /**< Coplanar intersection polylines for mesh A */
    std::unordered_map<size_t, std::vector<size_t>> Atris;         /**< Mapping of intersected elements in mesh A */
    std::unordered_map<size_t, PolyLine> NCBcoords;                /**< Non-coplanar intersection segments for mesh B */
    std::unordered_map<size_t, std::vector<PolyLine>> CBcoords;    /**< Coplanar intersection polylines for mesh B */
    std::unordered_map<size_t, std::vector<size_t>> Btris;         /**< Mapping of intersected elements in mesh B */
};

/**
* @struct Connection
* @brief Stores vertex mappings and connectivity across non-conformal mesh interfaces.
*/
struct Connection {
    std::unordered_map<size_t, std::vector<Vec3>> meshAIntersections; /**< Intersection points on mesh A elements */
    std::unordered_map<size_t, std::vector<Vec3>> meshBIntersections; /**< Intersection points on mesh B elements */
    std::unordered_map<size_t, std::vector<size_t>> aToBConnections;  /**< Node topology mapping from A to B */
    std::unordered_map<size_t, std::vector<size_t>> bToAConnections;  /**< Node topology mapping from B to A */
};

/**
* @struct ProjectionFrame
* @brief Defines a local 2D coordinate system on a 3D planar surface.
*/
struct ProjectionFrame {
    Vec3 origin; /**< Origin point in 3D space */
    Vec3 u;      /**< Local U-axis orthonormal basis vector */
    Vec3 v;      /**< Local V-axis orthonormal basis vector */

    /**
    * @brief Projects a 3D coordinate point onto the local 2D frame.
    * @param[in] p - 3D Point coordinate.
    * @return Projected 2D point.
    */
    CDT::V2d<double> to2D(const Vec3& p) const {
        Vec3 rel = p - origin;
        return {BVH::dot(rel, u), BVH::dot(rel, v)};
    }

    /**
    * @brief Unprojects a local 2D point back into global 3D space.
    * @param[in] p - Local 2D point.
    * @return Unprojected 3D coordinate.
    */
    Vec3 to3D(const CDT::V2d<double>& p) const {
        return origin + u * p.x + v * p.y;
    }
};