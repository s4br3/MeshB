#include "triangulation.hpp"
ProjectionFrame computeSharedFrame(const Vec3& normal, const Vec3& origin) {
    ProjectionFrame frame;
    frame.origin = origin;
    Vec3 absNormal = {std::abs(normal[0]), std::abs(normal[1]), std::abs(normal[2])};
    if (absNormal[0] >= absNormal[1] && absNormal[0] >= absNormal[2]) {
        frame.u = {0, 1, 0}; frame.v = {0, 0, 1};
    } else if (absNormal[1] >= absNormal[0] && absNormal[1] >= absNormal[2]) {
        frame.u = {1, 0, 0}; frame.v = {0, 0, 1};
    } else {
        frame.u = {1, 0, 0}; frame.v = {0, 1, 0};
    }
    return frame;
}
bool intersect2D(
    const CDT::V2d<double>& A, const CDT::V2d<double>& B,
    const CDT::V2d<double>& C, const CDT::V2d<double>& D,
    CDT::V2d<double>& out, const double eps)
{
    double a1 = B.y - A.y, b1 = A.x - B.x, c1 = a1 * A.x + b1 * A.y;
    double a2 = D.y - C.y, b2 = C.x - D.x, c2 = a2 * C.x + b2 * C.y;
    double det = a1 * b2 - a2 * b1;
    if (std::abs(det) < eps) return false;
    out.x = (b2 * c1 - b1 * c2) / det;
    out.y = (a1 * c2 - a2 * c1) / det;
    auto inRange = [eps](double val, double min_val, double max_val) {
        return val >= std::min(min_val, max_val) - eps && val <= std::max(min_val, max_val) + eps;
    };
    return inRange(out.x, A.x, B.x) && inRange(out.y, A.y, B.y) &&
           inRange(out.x, C.x, D.x) && inRange(out.y, C.y, D.y);
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
    SpatialGrid2D grid(eps);
    std::vector<CDT::V2d<double>> candidate_pts = initialPts;
    for (size_t i = 0; i < segs.size(); ++i) {
        for (size_t j = i + 1; j < segs.size(); ++j) {
            CDT::V2d<double> intersection;
            if (intersect2D(initialPts[segs[i].first], initialPts[segs[i].second],
                            initialPts[segs[j].first], initialPts[segs[j].second], intersection, eps))
            {
                candidate_pts.push_back(intersection);
            }
        }
    }
    for (const auto& p : candidate_pts) {
        grid.getOrAdd(p);
    }
    for (const auto& seg : segs) {
        CDT::V2d<double> start = initialPts[seg.first];
        CDT::V2d<double> end   = initialPts[seg.second];
        double seg_len = std::hypot(end.x - start.x, end.y - start.y);
        if (seg_len < eps) continue;
        std::vector<std::pair<double, CDT::VertInd>> pointsOnSeg;
        for (const auto& p : candidate_pts) {
            double d1 = std::hypot(p.x - start.x, p.y - start.y);
            double d2 = std::hypot(end.x - p.x, end.y - p.y);
            if (std::abs((d1 + d2) - seg_len) < eps) {
                CDT::VertInd idx = grid.getOrAdd(p);
                pointsOnSeg.push_back({d1, idx});
            }
        }
        std::sort(pointsOnSeg.begin(), pointsOnSeg.end());
        for (size_t i = 0; i < pointsOnSeg.size() - 1; ++i) {
            CDT::VertInd v1 = pointsOnSeg[i].second;
            CDT::VertInd v2 = pointsOnSeg[i + 1].second;
            if (v1 != v2) {
                outCDTEdges.push_back(CDT::Edge{v1, v2});
            }
        }
    }
    outUniquePts = grid.getUniquePoints();
}
std::vector<Triangle> triangulate(
    const PolyLine& polygonSegments,
    const PolyLine& cuts,
    const ProjectionFrame& frame,
    const double eps,
    SpatialGrid3D& nodeGrid,
    bool coplanar)
{
    if (polygonSegments.empty() && cuts.empty()) return {};
    std::vector<CDT::V2d<double>> initialPts;
    std::vector<std::pair<size_t, size_t>> segs;
    for (const auto& seg : polygonSegments) {
        size_t idx1 = initialPts.size();
        initialPts.push_back(frame.to2D(seg.first));
        size_t idx2 = initialPts.size();
        initialPts.push_back(frame.to2D(seg.second));
        segs.push_back({idx1, idx2});
    }
    for (const auto& cut : cuts) {
        size_t idx1 = initialPts.size();
        initialPts.push_back(frame.to2D(cut.first));
        size_t idx2 = initialPts.size();
        initialPts.push_back(frame.to2D(cut.second));
        segs.push_back({idx1, idx2});
    }
    std::vector<CDT::V2d<double>> uniquePts;
    std::vector<CDT::Edge> CDTEdges;
    buildSubdividedEdges(initialPts, segs, uniquePts, CDTEdges, eps);
    CDT::RemoveDuplicatesAndRemapEdges(uniquePts, CDTEdges);
    CDT_Triangulation cdt;
    cdt.insertVertices(uniquePts);
    cdt.conformToEdges(CDTEdges);
    if (coplanar){
        cdt.eraseOuterTrianglesAndHoles();
    } else {
        cdt.eraseOuterTriangles();
    }
    std::vector<size_t> localToGlobal(cdt.vertices.size());
    for (size_t i = 0; i < cdt.vertices.size(); ++i) {
        Vec3 p3d = frame.to3D(cdt.vertices[i]);
        localToGlobal[i] = nodeGrid.getOrAdd(p3d);
    }
    std::vector<Triangle> result;
    result.reserve(cdt.triangles.size());
    for (const auto& tri : cdt.triangles) {
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
    const std::vector<Triangle>& tris,
    const std::vector<Vec3>& nodes, 
    const PolyLine& cuts,
    const ProjectionFrame& frame,
    const double eps,
    SpatialGrid3D& nodeGrid,
    bool coplanar)
{
    PolyLine segs;
    for (const Triangle& tri : tris){
        segs.push_back({nodes[tri.v[0]], nodes[tri.v[1]]});
        segs.push_back({nodes[tri.v[1]], nodes[tri.v[2]]});
        segs.push_back({nodes[tri.v[2]], nodes[tri.v[0]]});
    }
    return triangulate(segs, cuts, frame, eps, nodeGrid, coplanar);
}
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCcoords,
    const std::unordered_map<size_t, std::vector<PolyLine>>& Ccoords,
    const double eps,
    bool removeInternalSurfaces)
{
    SpatialGrid3D nodeGrid(eps);
    for (const auto& node : meshData.nodes) {
        nodeGrid.getOrAdd(node);
    }
    std::unordered_set<size_t> modifiedIndices;
    for (const auto& [idx, _] : NCcoords) modifiedIndices.insert(idx);
    for (const auto& [idx, _] : Ccoords) modifiedIndices.insert(idx);
    size_t maxTag = 0;
    for (size_t tag : meshData.tags) {
        maxTag = std::max(maxTag, tag);
    }
    MeshData outData;
    std::vector<Triangle> extraTriangles;
    for (size_t i = 0; i < meshData.triangles.size(); ++i) {
        if (modifiedIndices.count(i) == 0) {
            outData.triangles.push_back(meshData.triangles[i]);
            outData.tags.push_back(meshData.tags[i]);
            continue;
        }
        const auto& tri = meshData.triangles[i];
        const auto& tag = meshData.tags[i];
        ProjectionFrame frame = computeSharedFrame(tri.normal(meshData.nodes), meshData.nodes[tri.v[0]]);
        PolyLine CCuts;
        auto itC = Ccoords.find(i);
        if (itC != Ccoords.end()) {
            auto temp = flattenVector(itC->second);
            CCuts.insert(CCuts.end(), temp.begin(), temp.end());
        }
        std::vector<Triangle> subTris = cutTriangles({tri}, meshData.nodes, CCuts, frame, eps, nodeGrid, true);
        auto itNC = NCcoords.find(i);
        PolyLine NCCuts;
        if (itNC != NCcoords.end()) {
            NCCuts.insert(NCCuts.end(), itNC->second.begin(), itNC->second.end());
        }
        subTris = cutTriangles(subTris, nodeGrid.getUniquePoints(), NCCuts, frame, eps, nodeGrid, false);
        if (!subTris.empty()) {
            outData.triangles.push_back(subTris[0]);
            outData.tags.push_back(tag);
            for (size_t k = 1; k < subTris.size(); ++k) {
                extraTriangles.push_back(subTris[k]);
            }
        }
        if (itC != Ccoords.end() && !removeInternalSurfaces) {
            for (const PolyLine& hole : itC->second) {
                std::vector<Triangle> holeTris;
                if (hole.size() > 3){
                    holeTris = triangulate(hole, hole, frame, eps, nodeGrid, true);
                } else if (hole.size() >= 2) {
                    Vec3 v0 = hole[0].first;
                    Vec3 v1 = hole[0].second;
                    Vec3 v2 = (hole[1].first[0] == v0[0] && hole[1].first[1] == v0[1] && hole[1].first[2] == v0[2]) ? hole[1].second : hole[1].first;
                    size_t idx0 = nodeGrid.getOrAdd(v0);
                    size_t idx1 = nodeGrid.getOrAdd(v1);
                    size_t idx2 = nodeGrid.getOrAdd(v2);
                    Triangle hTri;
                    hTri.v = {idx0, idx1, idx2};
                    holeTris.push_back(hTri);
                }
                extraTriangles.insert(extraTriangles.end(), std::make_move_iterator(holeTris.begin()), std::make_move_iterator(holeTris.end()));
            }
        }
    }
    for (size_t i = 0; i < extraTriangles.size(); ++i) {
        outData.triangles.push_back(std::move(extraTriangles[i]));
        outData.tags.push_back(maxTag + 1 + i);
    }
    outData.nodes = nodeGrid.getUniquePoints();
    const size_t numTris = outData.triangles.size();
    outData.centres.resize(numTris);
    outData.normals.resize(numTris);
    for (size_t i = 0; i < numTris; ++i) {
        outData.centres[i] = outData.triangles[i].centre(outData.nodes);
        outData.normals[i] = outData.triangles[i].normal(outData.nodes);
    }
    return outData;
}