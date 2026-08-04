#pragma once
#include "mesh_types.hpp"
using CDT_Triangulation = CDT::Triangulation<Scalar>;
ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin);
bool intersect2D(
    const CDT::V2d<double>& A, const CDT::V2d<double>& B,
    const CDT::V2d<double>& C, const CDT::V2d<double>& D,
    CDT::V2d<double>& out, double eps);
CDT::VertInd getOrAddUniquePoint(std::vector<CDT::V2d<double>>& unique_pts, const CDT::V2d<double>& p, double eps);
void buildSubdividedEdges(
    const std::vector<CDT::V2d<double>>& initial_pts, const std::vector<std::pair<size_t, size_t>>& segs,
    std::vector<CDT::V2d<double>>& out_unique_pts, std::vector<CDT::Edge>& out_cdt_edges,
    double eps);
std::vector<Triangle> triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, double eps, SpatialGrid3D& nodeGrid, bool coplanar);
std::vector<Triangle> cutTriangles(
    const std::vector<Triangle>& tris, const PolyLine& cuts,
    const ProjectionFrame& frame, double eps, SpatialGrid3D& nodeGrid, bool coplanar);
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCcoords, const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
    double eps, bool removeInternalSurfaces = false);