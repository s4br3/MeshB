#pragma once
#include "mesh_types.hpp"

/**
* @brief Alias for Constrained Delaunay Triangulation structure.
*/
using CDT_Triangulation = CDT::Triangulation<Scalar>;

/**
* @brief Constructs an orthonormal 2D projection frame for a planar surface.
* @param[in] normal - Plane normal vector.
* @param[in] origin - Reference origin point on plane.
* @return Orthonormal ProjectionFrame object.
*/
ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin);

/**
* @brief Computes 2D line-segment intersection point.
* @param[in] A - Start point of segment 1.
* @param[in] B - End point of segment 1.
* @param[in] C - Start point of segment 2.
* @param[in] D - End point of segment 2.
* @param[out] out - Calculated 2D intersection point coordinate.
* @param[in] eps - Distance tolerance.
* @return True if segments intersect within valid bounds.
*/
bool intersect2D(
    const CDT::V2d<double>& A, const CDT::V2d<double>& B,
    const CDT::V2d<double>& C, const CDT::V2d<double>& D,
    CDT::V2d<double>& out, double eps);

/**
* @brief Deduplicates a 2D point against an existing point array.
* @param[in,out] unique_pts - Global list of unique 2D points.
* @param[in] p - Query point coordinate.
* @param[in] eps - Deduplication tolerance distance.
* @return CDT vertex index corresponding to query point.
*/
CDT::VertInd getOrAddUniquePoint(std::vector<CDT::V2d<double>>& unique_pts, const CDT::V2d<double>& p, double eps);

/**
* @brief Subdivides overlapping constraints into a clean, planar constrained segment graph.
* @param[in] initial_pts - Initial set of 2D vertex positions.
* @param[in] segs - Initial input constraint edges.
* @param[out] out_unique_pts - Output deduplicated vertex array.
* @param[out] out_cdt_edges - Output CDT constraint edges.
* @param[in] eps - Distance tolerance.
*/
void buildSubdividedEdges(
    const std::vector<CDT::V2d<double>>& initial_pts, const std::vector<std::pair<size_t, size_t>>& segs,
    std::vector<CDT::V2d<double>>& out_unique_pts, std::vector<CDT::Edge>& out_cdt_edges,
    double eps);

/**
* @brief Performs Constrained Delaunay Triangulation on a set of 3D polylines and cut lines projected into 2D space.
* @param[in] polygonSegments - Outer domain boundary polylines.
* @param[in] cuts - Internal cutting/intersection polylines.
* @param[in] frame - Projection frame defining planar mapping.
* @param[in] eps - Spatial tolerance threshold.
* @param[in,out] nodeGrid - Spatial 3D lookup grid for vertex merging.
* @param[in] coplanar - Flag indicating if triangulation applies to coplanar surface patches.
* @return Array of newly retriangulated 3D Triangle elements.
*/
std::vector<Triangle> triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, double eps, SpatialGrid3D& nodeGrid, bool coplanar);

/**
* @brief Subdivides and remeshes a list of triangles along intersection polylines.
* @param[in] tris - Vector of triangles to cut.
* @param[in] cuts - Polyline cut paths crossing elements.
* @param[in] frame - Planar projection frame.
* @param[in] eps - Tolerance threshold.
* @param[in,out] nodeGrid - Spatial node matching structure.
* @param[in] coplanar - Coplanar surface flag.
* @return Subdivided set of output triangles.
*/
std::vector<Triangle> cutTriangles(
    const std::vector<Triangle>& tris, const PolyLine& cuts,
    const ProjectionFrame& frame, double eps, SpatialGrid3D& nodeGrid, bool coplanar);

/**
* @brief Retriangulates and cuts an entire surface mesh along non-coplanar and coplanar intersection constraints.
* @param[in] meshData - Source mesh to cut.
* @param[in] NCcoords - Non-coplanar intersection paths indexed by triangle index.
* @param[in] Ccoords - Coplanar intersection polylines indexed by triangle index.
* @param[in] eps - Geometric tolerance.
* @param[in] removeInternalSurfaces - If true, discards generated internal/overlapping patch elements.
* @return Subdivided cut MeshData object.
*/
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCcoords, const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
    double eps, bool removeInternalSurfaces = false);