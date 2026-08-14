#pragma once
#include "mesh_types.hpp"
#include "math_utils.hpp"
#include <set>
using EdgeKey = std::pair<size_t, size_t>;
using CDT_Triangulation = CDT::Triangulation<double>;

/**
* @brief Standard key of any edge, with smaller index first and larger index second.
* @param [in] u - Index of one endpoint of the edge.
* @param [in] v - Index of the other endpoint of the edge.
* @return The pair of indices in ascending order.
*/
inline EdgeKey makeEdgeKey(size_t u, size_t v) {
    return {std::min(u, v), std::max(u, v)};
}

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
* @param[out] outTriangles - Output buffer for newly retriangulated 3D Triangle elements.
*/
void triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    double eps, std::vector<Triangle>& outTriangles);

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
* @brief Given the adjacencies of a graph, find every cycle.
* @param[in] curr - Current index in the graph.
* @param[in] adjacencies - Adjacency List of the graph. Since the graph is undirected, we want
* an edge from A <-> resulting in A's associated vector containing B and vice versa.
* @param[in] pathStack - Stack containing all the indices traversed (excluding curr) so far.
* @param[in, out] allCycles - Output buffer containing every cycle in the graph.
*/
void dfsEdges(
    size_t curr,
    const std::unordered_map<size_t, std::vector<size_t>>& adjacencies,
    std::set<EdgeKey>& visitedEdges,
    std::vector<size_t>& pathStack,
    std::vector<std::vector<size_t>>& allCycles);

/**
* @brief Find all cycles in a PolyLine soup ("holes").
* @param[in] edges - Group of pairs of vectors (line segments) in arbitrary order.
* @param[in] nodeGrid - Spatial node matching structure.
* @param[out] outLoops - Output buffer containing all valid polyline loops found.
*/
void findCycles(const PolyLine& edges, SpatialGrid3D& nodeGrid, std::vector<PolyLine>& outLoops);

/**
* @brief Determines if a point projected into 2D space falls inside any cycle hole region.
* @param[in] centre - 3D Query coordinate.
* @param[in] allHoles - Vector of closed boundary polyline cycles defining hole regions.
* @param[in] frame - Planar projection frame.
* @return True if point lies inside any hole boundary.
*/
bool isCentroidInHole(
    const Vec3 centre,
    const std::vector<PolyLine>& allHoles, 
    const ProjectionFrame& frame);
    
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