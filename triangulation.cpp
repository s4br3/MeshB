#include "triangulation.hpp"
#include "CDTUtils.h"
#include "math_utils.hpp"
#include "mesh_types.hpp"
#include <set>

ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin) {
    ProjectionFrame frame;
    frame.origin = origin;
    Vec3 absNormal = {std::abs(normal[0]), std::abs(normal[1]), std::abs(normal[2])};
    Vec3 helper = {0, 0, 0};
    if (absNormal[0] <= absNormal[1] && absNormal[0] <= absNormal[2]) {
        helper[0] = 1.0; 
    } else if (absNormal[1] <= absNormal[0] && absNormal[1] <= absNormal[2]) {
        helper[1] = 1.0;
    } else {
        helper[2] = 1.0;
    }
    frame.u = cross(helper, normal).normalize();
    frame.v = cross(normal, frame.u).normalize();
    return frame;
}
static bool pointOnSegment(
    const CDT::V2d<double>& P,
    const CDT::V2d<double>& A,
    const CDT::V2d<double>& B,
    double eps)
{
    double vx = B.x - A.x, vy = B.y - A.y;
    double wx = P.x - A.x, wy = P.y - A.y;
    double cross = vx * wy - vy * wx;
    if (std::abs(cross) > eps) return false;
    auto inRange = [eps](double v, double a, double b) {
        return v >= std::min(a,b) - eps && v <= std::max(a,b) + eps;
    };
    return inRange(P.x, A.x, B.x) && inRange(P.y, A.y, B.y);
}

static void intersect2DAllPoints(
    const CDT::V2d<double>& A, const CDT::V2d<double>& B,
    const CDT::V2d<double>& C, const CDT::V2d<double>& D,
    SpatialGrid2D& grid,
    std::vector<CDT::VertInd>& outs,
    double eps)
{
    outs.clear();
    double x12 = B.x - A.x, y12 = B.y - A.y;
    double x34 = D.x - C.x, y34 = D.y - C.y;
    double x13 = C.x - A.x, y13 = C.y - A.y;
    double denom = x12 * y34 - y12 * x34;
    auto addIdUnique = [&](CDT::VertInd id) {
        for (auto existing : outs) {
            if (existing == id) return;
        }
        outs.push_back(id);
    };
    if (std::abs(denom) >= eps) {
        double t = (x13 * y34 - y13 * x34) / denom;
        double u = (x13 * y12 - y13 * x12) / denom;
        if (t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps) {
            addIdUnique(grid.getOrAdd({A.x + t * x12, A.y + t * y12}));
        }
        return;
    }
    if (std::abs(x13 * y12 - y13 * x12) > eps) return;
    std::vector<CDT::V2d<double>> candidates;
    if (pointOnSegment(C, A, B, eps)) addIdUnique(grid.getOrAdd(C));
    if (pointOnSegment(D, A, B, eps)) addIdUnique(grid.getOrAdd(D));
    if (pointOnSegment(A, C, D, eps)) addIdUnique(grid.getOrAdd(A));
    if (pointOnSegment(B, C, D, eps)) addIdUnique(grid.getOrAdd(B));
}
void buildSubdividedEdges(
    const std::vector<CDT::V2d<double>>& initialPts,
    const std::vector<std::pair<size_t, size_t>>& segs,
    std::vector<CDT::V2d<double>>& outUniquePts,
    std::vector<CDT::Edge>& outCDTEdges,
    double eps)
{
    outUniquePts.clear();
    outCDTEdges.clear();
    if (segs.empty()) return;
    SpatialGrid2D grid(eps);
    std::vector<std::vector<CDT::VertInd>> pointsOnSegments(segs.size());
    // 1. Populate initial endpoints for ALL segments
    for (size_t i = 0; i < segs.size(); ++i) {
        CDT::V2d<double> A = initialPts[segs[i].first];
        CDT::V2d<double> B = initialPts[segs[i].second];
        pointsOnSegments[i].push_back(grid.getOrAdd(A));
        pointsOnSegments[i].push_back(grid.getOrAdd(B));
    }
    // 2. Compute 2D intersections between all segment pairs
    for (size_t i = 0; i + 1 < segs.size(); ++i) {
        CDT::V2d<double> A = initialPts[segs[i].first];
        CDT::V2d<double> B = initialPts[segs[i].second];
        for (size_t j = i + 1; j < segs.size(); ++j) {
            CDT::V2d<double> C = initialPts[segs[j].first];
            CDT::V2d<double> D = initialPts[segs[j].second];
            std::vector<CDT::VertInd> outs;
            intersect2DAllPoints(A, B, C, D, grid, outs, eps);
            pointsOnSegments[i].insert(pointsOnSegments[i].end(), outs.begin(), outs.end());
            pointsOnSegments[j].insert(pointsOnSegments[j].end(), outs.begin(), outs.end());
        }
    }
    outUniquePts = grid.getUniquePoints();
    for (size_t i = 0; i < segs.size(); ++i) {
        CDT::V2d<double> A = initialPts[segs[i].first];
        CDT::V2d<double> B = initialPts[segs[i].second];
        const double vx = B.x - A.x;
        const double vy = B.y - A.y;
        const double lenSq = vx * vx + vy * vy;
        if (lenSq <= eps * eps) continue;
        std::vector<CDT::VertInd> ids = pointsOnSegments[i];
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        std::vector<std::pair<double, CDT::VertInd>> parametrized;
        parametrized.reserve(ids.size());
        for (CDT::VertInd id : ids) {
            const CDT::V2d<double>& P = outUniquePts[id];
            double t = ((P.x - A.x) * vx + (P.y - A.y) * vy) / lenSq;
            parametrized.push_back({t, id});
        }
        std::sort(parametrized.begin(), parametrized.end());
        std::vector<CDT::VertInd> cleanIds;
        cleanIds.reserve(parametrized.size());
        for (const auto& [t, id] : parametrized) {
            if (cleanIds.empty() || cleanIds.back() != id) {
                cleanIds.push_back(id);
            }
        }
        for (size_t k = 0; k + 1 < cleanIds.size(); ++k) {
            if (cleanIds[k] != cleanIds[k + 1]) {
                outCDTEdges.emplace_back(cleanIds[k], cleanIds[k + 1]);
            }
        }
    }
}
std::vector<Triangle> triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    bool coplanar, const double eps)
{
    if (polygonSegments.empty() && cuts.empty()) return {};
    std::vector<CDT::V2d<double>> initialPts;
    std::vector<std::pair<size_t, size_t>> segs;
    for (const std::pair<Vec3, Vec3>& seg : polygonSegments) {
        size_t idx1 = nodeGrid.getOrAdd(seg.first);
        size_t idx2 = nodeGrid.getOrAdd(seg.second);
        if (idx1 != idx2) segs.push_back({idx1, idx2});
    }
    for (const std::pair<Vec3, Vec3>& cut : cuts) {
        size_t idx1 = nodeGrid.getOrAdd(cut.first);
        size_t idx2 = nodeGrid.getOrAdd(cut.second);
        if (idx1 != idx2) segs.push_back({idx1, idx2});
    }
    for (const Vec3& vec : nodeGrid.getUniquePoints()) initialPts.push_back(frame.to2D(vec));
    std::vector<CDT::V2d<double>> uniquePts;
    std::vector<CDT::Edge> CDTEdges;
    buildSubdividedEdges(initialPts, segs, uniquePts, CDTEdges, eps);
    CDT::RemoveDuplicatesAndRemapEdges(uniquePts, CDTEdges);
    std::set<std::pair<CDT::VertInd, CDT::VertInd>> cleanEdges;
    std::vector<CDT::Edge> deduplicatedCDTEdges;
    for (const CDT::Edge& edge : CDTEdges) {
        CDT::VertInd v1 = edge.v1();
        CDT::VertInd v2 = edge.v2();
        if (v1 == v2) continue;
        std::pair<CDT::VertInd, CDT::VertInd> key = std::make_pair(std::min(v1, v2), std::max(v1, v2));
        if (cleanEdges.insert(key).second) {
            deduplicatedCDTEdges.push_back(CDT::Edge{v1, v2});
        }
    }
    CDTEdges = std::move(deduplicatedCDTEdges);
    CDT_Triangulation cdt;
    cdt.insertVertices(uniquePts);
    cdt.conformToEdges(CDTEdges);
    if (coplanar) {
        cdt.eraseOuterTrianglesAndHoles();
    } else {
        cdt.eraseOuterTriangles();
    }
    std::vector<size_t> localToGlobal(cdt.vertices.size());
    for (size_t i = 0; i < cdt.vertices.size(); ++i) {
        localToGlobal[i] = nodeGrid.getOrAdd(frame.to3D(cdt.vertices[i]));
    }
    std::vector<Triangle> result;
    result.reserve(cdt.triangles.size());
    for (const CDT::Triangle& tri : cdt.triangles) {
        Triangle newTri;
        newTri.v = {
            localToGlobal[tri.vertices[0]],
            localToGlobal[tri.vertices[1]],
            localToGlobal[tri.vertices[2]]
        };
        result.push_back(newTri);
    }
    return result;
}
std::vector<Triangle> cutTriangles(
    const std::vector<Triangle>& tris, const std::vector<Vec3>& nodes, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    bool coplanar, const double eps)
{
    std::unordered_map<std::pair<size_t, size_t>, size_t> edgeCounts;
    std::unordered_map<std::pair<size_t, size_t>, std::pair<Vec3, Vec3>> edgeGeom;
    for (const Triangle& tri : tris) {
        for (int i = 0; i < 3; ++i) {
            size_t u = tri.v[i];
            size_t v = tri.v[(i + 1) % 3];
            std::pair<size_t, size_t> key = {std::min(u, v), std::max(u, v)};
            edgeCounts[key]++;
            if (edgeCounts[key] == 1) {
                edgeGeom[key] = {nodes[u], nodes[v]};
            }
        }
    }
    PolyLine segs;
    for (const auto& [key, count] : edgeCounts) {
        if (count == 1) {
            segs.push_back(edgeGeom[key]);
        }
    }

    return triangulate(segs, cuts, frame, nodeGrid, coplanar, eps);
}
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCcoords,
    const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
    const double eps,
    bool removeInternalSurfaces)
{
    SpatialGrid3D nodeGrid(eps);
    for (const Vec3& node : meshData.nodes) {
        nodeGrid.getOrAdd(node);
    }
    std::unordered_set<size_t> modifiedIndices;
    for (const auto& [idx, _] : NCcoords) modifiedIndices.insert(idx);
    for (const auto& [idx, _] : Ccoords) modifiedIndices.insert(idx);
    size_t tagOffset = 0;
    for (size_t tag : meshData.tags) {
        tagOffset = std::max(tagOffset, tag + 1);
    }
    MeshData outData;
    std::vector<Triangle> extraTriangles;
    for (size_t i = 0; i < meshData.triangles.size(); ++i) {
        if (modifiedIndices.count(i) == 0) {
            outData.triangles.push_back(meshData.triangles[i]);
            outData.tags.push_back(meshData.tags[i]);
            continue;
        }
        ProjectionFrame frame = computeSharedFrame(meshData.normals[i], meshData.centres[i]);
        PolyLine CCuts;
        auto itC = Ccoords.find(i);
        if (itC != Ccoords.end()) {
            PolyLine temp = flattenVector(itC->second);
            CCuts.insert(CCuts.end(), temp.begin(), temp.end());
        }
        std::vector<Triangle> subTris = cutTriangles({meshData.triangles[i]}, meshData.nodes, CCuts, frame, nodeGrid, true, eps);
        auto itNC = NCcoords.find(i);
        PolyLine NCCuts;
        if (itNC != NCcoords.end()) {
            NCCuts.insert(NCCuts.end(), itNC->second.begin(), itNC->second.end());
        }
        subTris = cutTriangles(subTris, nodeGrid.getUniquePoints(), NCCuts, frame, nodeGrid, removeInternalSurfaces, eps);
        if (!subTris.empty()) {
            outData.triangles.push_back(subTris[0]);
            outData.tags.push_back(meshData.tags[i]);
            extraTriangles.insert(extraTriangles.end(), subTris.begin() + 1, subTris.end());
        }
        if (itC != Ccoords.end() && !removeInternalSurfaces) {
            for (const PolyLine& hole : itC->second) {
                std::vector<Triangle> holeTris = triangulate(hole, {}, frame, nodeGrid, false, eps);
                extraTriangles.insert(extraTriangles.end(), std::make_move_iterator(holeTris.begin()), std::make_move_iterator(holeTris.end()));
            }
        }
    }
    for (size_t i = 0; i < extraTriangles.size(); ++i) {
        outData.triangles.push_back(std::move(extraTriangles[i]));
        outData.tags.push_back(tagOffset + i);
    }
    outData.nodes = nodeGrid.getUniquePoints();
    recomputeMeshData(outData);
    return outData;
}