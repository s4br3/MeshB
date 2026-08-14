#pragma once
#include "vector.hpp"
#include "bvh.hpp"
#include "mesh_types.hpp"
#include <cstddef>
#include <cmath>
#include <functional>
/**
* @enum FaceClass
* @brief Topological classification of a mesh element relative to a target surface.
*/
enum class FaceClass { Outside, Inside, CoplanarSame, CoplanarOpp };

/**
* @struct TupleHash
* @brief Hash functor utility for computing hash values of three-element tuples.
*/
struct TupleHash {
    template <class T1, class T2, class T3>
    std::size_t operator()(const std::tuple<T1, T2, T3>& t) const {
        auto h1 = std::hash<T1>{}(std::get<0>(t));
        auto h2 = std::hash<T2>{}(std::get<1>(t));
        auto h3 = std::hash<T3>{}(std::get<2>(t));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

/**
* @struct GWNNodeData
* @brief Stores precomputed geometric data for a node to accelerate Generalized Winding Number evaluation.
*/
struct GWNNodeData {
    Vec3 vectorArea;    /**< Total surface area of the triangles descendant of this node */
    Vec3 center;        /**< Centroid of the node */
};

/**
* @brief Checks if a given mesh is closed (watertight).
* @param[in] mesh - Input mesh dataset to check.
* @return True if the mesh is closed, false otherwise.
*/
bool isMeshClosed(const MeshData& mesh);

/**
* @brief Computes the solid angle subtended by a triangle (a, b, c) from a point p.
* @param[in] p - Query point coordinate.
* @param[in] a - First triangle vertex.
* @param[in] b - Second triangle vertex.
* @param[in] c - Third triangle vertex.
* @return Computed solid angle value.
*/
double solidAngle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c);

/**
* @brief Builds Generalized Winding Number (GWN) node data for efficient evaluation.
* @param[in] mesh - Input mesh dataset.
* @param[in] bvh - BVH tree constructed from the mesh.
* @return Vector of GWNNodeData containing vector areas and centers.
*/
std::vector<GWNNodeData> buildGWNData(const MeshData& mesh, const BVH& bvh);

/**
* @brief Evaluates the Generalized Winding Number at a given query point.
* @param[in] p - Query point coordinate.
* @param[in] mesh - Input mesh dataset.
* @param[in] bvh - BVH tree of the mesh.
* @param[in] gwnData - Precomputed GWN node data.
* @return Evaluated winding number scalar value.
*/
double evaluateGWN(const Vec3& p, const MeshData& mesh, const BVH& bvh, const std::vector<GWNNodeData>& gwnData);

/**
* @brief Classifies a triangle face as inside, outside, or coplanar relative to a target surface mesh.
* @param[in] targetBVH - BVH tree of target closed domain mesh.
* @param[in] targetMesh - Target mesh dataset.
* @param[in] triCenter - Centroid of query triangle element.
* @param[in] triNormal - Normal vector of query triangle element.
* @param[in] eps - Distance tolerance.
* @return FaceClass classification result.
*/
FaceClass classifyFace(const BVH& targetBVH, const MeshData& targetMesh,
    const Vec3& triCentre, const Vec3& triNormal,
    double eps);