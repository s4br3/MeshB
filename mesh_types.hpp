#pragma once
#include "vector.hpp"
#include "bvh.hpp"
#include <unordered_map>
#include <CDT.h>

using TriVerts = std::array<Vec3, 3>;

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
    Vec3 normal(const std::vector<Vec3>& nodes) const {
        return cross((nodes[v[1]] - nodes[v[0]]), (nodes[v[2]] - nodes[v[0]]));
    }

    /**
    * @brief Computes the centroid coordinate of the triangle.
    * @param[in] nodes - Global node list containing vertex positions.
    * @return Centroid coordinate in 3D space.
    */
    Vec3 centre(const std::vector<Vec3>& nodes) const { 
        return (nodes[v[0]] + nodes[v[1]] + nodes[v[2]]) * (1.0 / 3.0); 
    }

    /**
    * @brief Calculates the axis-aligned bounding box (AABB) of the triangle.
    * @param[in] nodes - Global node list containing vertex positions.
    * @return Bounding box encompassing all three vertices.
    */
    BBox bounds(const std::vector<Vec3>& nodes) const {
        BBox b;
        for (const auto& p : nodes) b.extend(p);
        return b;
    }
    /**
    * @brief Get the coordinates of vertices defining this triangle
    * @param[in] nodes - Global node list containing vertex positions
    * @return Vertex array containing the Vec3s defining this triangle
    */
    std::array<Vec3, 3> getVertices(const std::vector<Vec3>& nodes) const{
        return {nodes[v[0]], nodes[v[1]], nodes[v[2]]};
    }
    /**
    * @brief Compute the area of the triangle
    * @param[in] nodes - Global node list containing vertex positions
    * @return Triangle area
    */
    double getArea(const std::vector<Vec3>& nodes) const{
        return normal(nodes).length()/2;
    }
};
/**
* @brief Formatted stream output operator for a triangle given
* @param[in] os - Output stream
* @param[in] vertices - Vertex array containing the Vec3s defining this triangle
* @return Reference to the output stream
*/
std::ostream& operator<<(std::ostream& os, const std::array<Vec3, 3>& vertices);

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
std::ostream& operator<<(std::ostream& os, const MeshData& mesh);
void recomputeMeshData(MeshData& mesh);
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
        return {dot(rel, u), dot(rel, v)};
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
namespace std {

    /**
    * @brief Hashing for pairs of values (to be used for unordered_map)
    * @param[in] p - Pair of values (of any hashable types)
    * @return The pair's hash
    */
    template <typename T1, typename T2>
    struct hash<pair<T1, T2>> {
        size_t operator()(const pair<T1, T2>& p) const {
            size_t h1 = hash<T1>{}(p.first);
            size_t h2 = hash<T2>{}(p.second);
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}