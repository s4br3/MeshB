#pragma once
#include "mesh_types.hpp"
#include "geom_2d.hpp"
#include "math_utils.hpp"


/**
* @brief Subdivides overlapping constraints into a clean, planar constrained segment graph.
* @param[in] initial_pts - Initial set of 2D vertex positions.
* @param[in] segs - Initial input constraint edges.
* @param[out] out_unique_pts - Output deduplicated vertex array.
* @param[out] out_cdt_edges - Output CDT constraint edges.
* @param[in] eps - Distance tolerance.
*/
void buildSubdividedEdges(
    const std::vector<Vec2>& initial_pts, const std::vector<std::pair<size_t, size_t>>& segs,
    std::vector<Vec2>& out_unique_pts, std::vector<EdgeKey>& out_cdt_edges,
    double eps);

/**
* @brief Performs Constrained Delaunay Triangulation on a set of 3D polylines and cut lines projected into 2D space.
* @param[in] polygonSegments - Outer domain boundary polylines.
* @param[in] cuts - Internal cutting/intersection polylines.
* @param[in] frame - Projection frame defining planar mapping.
* @param[in] eps - Spatial tolerance threshold.
* @param[in,out] nodeGrid - Spatial 3D lookup grid for vertex merging.
* @param[out] outTriangles - Output buffer for newly retriangulated 3D Triangle elements.
* @param[in] constrained - Context boolean for constraining to given edges and avoiding insertion of Steiner points.
*/
void triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    double eps, std::vector<Triangle>& outTriangles,
    const bool constrained = true);

/**
* @brief Subdivides and remeshes a list of triangles along intersection polylines.
* @param[in] tris - Vector of triangles to cut.
* @param[in] nodes - Vector of coordinates of vertices.
* @param[in] cuts - Polyline cut paths crossing elements.
* @param[in] frame - Planar projection frame.
* @param[in,out] nodeGrid - Spatial node matching structure.
* @param[in] eps - Tolerance threshold.
* @param[out] outTriangles - Output buffer for subdivided set of triangles.
*/
void cutTriangles(
    const std::vector<Triangle>& tris, const std::vector<Vec3>& nodes, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    const double eps, std::vector<Triangle>& outTriangles);
    
/**
* @brief Retriangulates and cuts an entire surface mesh along non-coplanar and coplanar intersection constraints.
* @param[in] meshData - Source mesh to cut.
* @param[in] NCCuts - Non-coplanar intersection paths indexed by triangle index.
* @param[in] CCuts - Coplanar intersection polylines indexed by triangle index.
* @param[in] eps - Geometric tolerance.
* @param[in] removeTouchingSurfaces - If true, discards generated overlapping patch elements.
* @return Subdivided cut MeshData object.
*/
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCCuts, const std::unordered_map<size_t, std::vector<PolyLine>>& CCuts,
    double eps, bool removeTouchingSurfaces = false);