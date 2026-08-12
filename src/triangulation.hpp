#pragma once
#include "mesh_types.hpp"
#include "math_utils.hpp"


/**
* @brief Alias for Constrained Delaunay Triangulation structure.
*/
using CDT_Triangulation = CDT::Triangulation<double>;

/**
* @brief Constructs an orthonormal 2D projection frame for a planar surface.
* @param[in] normal - Plane normal vector.
* @param[in] origin - Reference origin point on plane.
* @return Orthonormal ProjectionFrame object.
*/
ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin);
/**
* @brief Finds whether a point P is on the line segment AB
* @param[in] P - Point to test
* @param[in] A - Start point of segment.
* @param[in] B - End point of segment.
* @param[in] eps - Tolerance value.
* return True if point is on line segment, false otherwise
*/
static bool pointOnSegment(
    const CDT::V2d<double>& P,
    const CDT::V2d<double>& A, const CDT::V2d<double>& B,
    double eps);
/**
* @brief Computes 2D line-segment intersection point.
* @param[in] A - Start point of segment 1.
* @param[in] B - End point of segment 1.
* @param[in] C - Start point of segment 2.
* @param[in] D - End point of segment 2.
* @param [in, out] grid - 2D grid used for vertex deduplication and storing.
* @param[in, out] out - Vector of indices of points in the 2D spatial grid.
* @param[in] eps - Distance tolerance.
* @return True if segments intersect within valid bounds.
*/
void intersect2DAllPoints(
    const CDT::V2d<double>& A, const CDT::V2d<double>& B,
    const CDT::V2d<double>& C, const CDT::V2d<double>& D,
    SpatialGrid2D& grid,
    std::vector<CDT::VertInd>& outs,
    double eps);

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
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    bool coplanar, double eps);

/**
* @brief Subdivides and remeshes a list of triangles along intersection polylines.
* @param[in] tris - Vector of triangles to cut.
* @param[in] nodes - Vector of coordinates of vertices.
* @param[in] cuts - Polyline cut paths crossing elements.
* @param[in] frame - Planar projection frame.
* @param[in,out] nodeGrid - Spatial node matching structure.
* @param[in] coplanar - Coplanar surface flag.
* @param[in] eps - Tolerance threshold.
* @return Subdivided set of output triangles.
*/
std::vector<Triangle> cutTriangles(
    const std::vector<Triangle>& tris, const std::vector<Vec3>& nodes, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    bool coplanar, const double eps);

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
    double eps, bool removeInternalSurfaces = false, bool meshA = true);