#include "triangulation.hpp"
#include "bvh.hpp"
#include "cdt.hpp"
#include "geom_3d.hpp"
#include "mesh_clean.hpp"
#include <cmath>
#include <unordered_map>
#include <unordered_set>
void buildSubdividedEdges(
    const std::vector<Vec2>& initialPts,
    const std::vector<EdgeKey>& segs,
    std::vector<Vec2>& outUniquePts,
    std::vector<EdgeKey>& outCDTEdges,
    double eps)
{
    outUniquePts.clear();
    outCDTEdges.clear();
    if (segs.empty()) return;
    SpatialGrid2D grid(eps);
    std::vector<BBox> boxes;
    std::vector<Vec3> centres;
    std::vector<std::vector<size_t>> pointsOnSegments(segs.size());
    std::vector<EdgeKey> stack;
    for (size_t i = 0; i < segs.size(); ++i) {
        Vec2 A = initialPts[segs[i].first];
        Vec2 B = initialPts[segs[i].second];
        pointsOnSegments[i].push_back(grid.getOrAdd(A));
        pointsOnSegments[i].push_back(grid.getOrAdd(B));
        BBox box{Vec3{A.x, A.y, 0}, Vec3{B.x, B.y, 0}};
        boxes.push_back(box);
        centres.push_back(Vec3{A.x + B.x, A.y + B.y, 0}/2);
    }
    BVH bvh2D(boxes, centres, eps);
    stack.emplace_back(0, 0);
    while (!stack.empty()){
        auto [node1Idx, node2Idx] = stack.back();
        stack.pop_back();
        const Node& node1 = bvh2D.nodes[node1Idx];
        const Node& node2 = bvh2D.nodes[node2Idx];
        if (!boundingBoxOverlap(node1.getBbox(), node2.getBbox(), eps))continue;
        if (node1.isLeaf() && node2.isLeaf()){
            for(size_t i = node1.firstId; i < node1.firstId + node1.primCount; i++){
                size_t prim1Id = bvh2D.primIds[i];
                for(size_t j = node2.firstId; j < node2.firstId + node2.primCount; j++){
                    size_t prim2Id = bvh2D.primIds[j];
                    if (prim1Id >= prim2Id) continue;
                    Edge e1 = {initialPts[segs[prim1Id].first], initialPts[segs[prim1Id].second]};
                    Edge e2 = {initialPts[segs[prim2Id].first], initialPts[segs[prim2Id].second]};
                    std::vector<size_t> outs;
                    intersect2DAllPoints(e1, e2, grid, outs, eps);
                    pointsOnSegments[prim1Id].insert(pointsOnSegments[prim1Id].end(), outs.begin(), outs.end());
                    pointsOnSegments[prim2Id].insert(pointsOnSegments[prim2Id].end(), outs.begin(), outs.end());
                }
            }
            continue;
        }
        if (node1.isLeaf()){
            size_t l2 = node2.firstId; size_t r2 = l2 + 1;
            stack.emplace_back(node1Idx, l2);
            stack.emplace_back(node1Idx, r2);
        }
        else if (node2.isLeaf()){
            size_t l1 = node1.firstId; size_t r1 = l1 + 1;
            stack.emplace_back(l1, node2Idx);
            stack.emplace_back(r1, node2Idx);
        }
        else{
            size_t l1 = node1.firstId; size_t r1 = l1 + 1;
            size_t l2 = node2.firstId; size_t r2 = l2 + 1;
            stack.emplace_back(l1, l2);
            stack.emplace_back(l1, r2);
            stack.emplace_back(r1, l2);
            stack.emplace_back(r1, r2);
        }
    }
    outUniquePts = grid.getUniquePoints();
    for (size_t i = 0; i < segs.size(); ++i) {
        Vec2 A = initialPts[segs[i].first];
        Vec2 B = initialPts[segs[i].second];
        const double vx = B.x - A.x;
        const double vy = B.y - A.y;
        const double lenSq = vx * vx + vy * vy;
        if (lenSq <= eps * eps) continue;
        std::vector<size_t> ids = pointsOnSegments[i];
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        std::vector<std::pair<double, size_t>> parametrized;
        parametrized.reserve(ids.size());
        for (size_t id : ids) {
            const Vec2& P = outUniquePts[id];
            double t = ((P.x - A.x) * vx + (P.y - A.y) * vy) / lenSq;
            t = std::max(0.0, std::min(1.0, t));
            parametrized.push_back( {t, id} );
        }
        std::sort(parametrized.begin(), parametrized.end());
        std::vector<size_t> cleanIds;
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
void triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    const double eps, std::vector<Triangle>& outTriangles,
    const bool constrained)
{
    std::vector<Vec2> initialPts;
    std::vector<EdgeKey> segs;
    std::vector<Edge> boundary2D;
    SpatialGrid2D grid2D(eps);
    for (const std::pair<Vec3, Vec3>& seg : polygonSegments) {
        Vec2 p1 = frame.to2D(seg.first);
        Vec2 p2 = frame.to2D(seg.second);
        size_t idx1 = grid2D.getOrAdd(p1);
        size_t idx2 = grid2D.getOrAdd(p2);
        if (idx1 != idx2) segs.push_back({idx1, idx2});
        boundary2D.push_back({p1, p2});
    }
    for (const std::pair<Vec3, Vec3>& cut : cuts) {
        Vec2 p1 = frame.to2D(cut.first);
        Vec2 p2 = frame.to2D(cut.second);
        size_t idx1 = grid2D.getOrAdd(p1);
        size_t idx2 = grid2D.getOrAdd(p2);
        segs.push_back({idx1, idx2});
    }
    initialPts = grid2D.getUniquePoints();
    std::vector<Vec2> uniquePts;
    std::vector<EdgeKey> CDTEdges;
    buildSubdividedEdges(initialPts, segs, uniquePts, CDTEdges, eps);
    std::unordered_set<EdgeKey> cleanEdges;
    std::vector<EdgeKey> deduplicatedCDTEdges;
    for (const EdgeKey& edge : CDTEdges) {
        size_t v1 = edge.first;
        size_t v2 = edge.second;
        if (v1 == v2) continue;
        if (cleanEdges.insert(makeEdgeKey(v1, v2)).second) {
            deduplicatedCDTEdges.push_back(EdgeKey {v1, v2} );
        }
    }
    CDTEdges = std::move(deduplicatedCDTEdges);
    if (uniquePts.size() < 3) return;
    std::vector<TriangleCDT> cdtTriangles;
    if (constrained){
        cdtTriangles = calculateCDT(uniquePts, CDTEdges, eps);
    } 
    else{
        cdtTriangles = calculateDelaunay(uniquePts, eps);
    }
    for (const auto& tri : cdtTriangles) {
        Vec2 centroid = {
            (tri.p1.x + tri.p2.x + tri.p3.x) / 3.0,
            (tri.p1.y + tri.p2.y + tri.p3.y) / 3.0
        };
        if (isPointInsidePolygon(centroid, boundary2D)) {
            Triangle newTri;
            newTri.v = {
                nodeGrid.getOrAdd(frame.to3D(tri.p1)),
                nodeGrid.getOrAdd(frame.to3D(tri.p2)),
                nodeGrid.getOrAdd(frame.to3D(tri.p3))
            };
            outTriangles.push_back(newTri);
        }
    }
}
void cutTriangles(
    const std::vector<Triangle>& tris, const std::vector<Vec3>& nodes, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    const double eps, std::vector<Triangle>& outTriangles)
{
    std::unordered_map<EdgeKey, size_t> edgeCounts;
    std::unordered_map<EdgeKey, std::pair<Vec3, Vec3>> edgeGeom;
    for (const Triangle& tri : tris) {
        for (int i = 0; i < 3; ++i) {
            size_t u = tri.v[i];
            size_t v = tri.v[(i + 1) % 3];
            EdgeKey key = makeEdgeKey(u, v);
            edgeCounts[key]++;
            if (edgeCounts[key] == 1) {
                edgeGeom[key] = {nodes[u], nodes[v]};
            }
        }
    }
    PolyLine segs;
    PolyLine allCuts = cuts;
    if (!tris.empty()) {
        const Triangle& originalFace = tris[0];
        for (int i = 0; i < 3; ++i) {
            segs.push_back({nodes[originalFace.v[i]], nodes[originalFace.v[(i + 1) % 3]]});
        }
    }
    for (const auto& [key, count] : edgeCounts) {
        allCuts.push_back(edgeGeom[key]);
    }
    triangulate(segs, allCuts, frame, nodeGrid, eps, outTriangles);
}
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCCuts,
    const std::unordered_map<size_t, std::vector<PolyLine>>& CCuts,
    const double eps, bool removeTouchingSurfaces)
{
    SpatialGrid3D nodeGrid(eps);
    std::vector<size_t> oldToNew(meshData.nodes.size());
    for (size_t i = 0; i <  meshData.nodes.size(); i++) {
        oldToNew[i] = nodeGrid.getOrAdd(meshData.nodes[i]);
    }
    std::unordered_set<size_t> modifiedIndices;
    for (const auto& [idx, _] : NCCuts) modifiedIndices.insert(idx);
    for (const auto& [idx, _] : CCuts) modifiedIndices.insert(idx);
    size_t tagOffset = 0;
    for (size_t tag : meshData.tags) {
        tagOffset = std::max(tagOffset, tag + 1);
    }
    MeshData outData;
    std::vector<Triangle> extraTriangles;
    for (size_t i = 0; i < meshData.triangles.size(); i++) {
        Triangle t = meshData.triangles[i];
            t.v[0] = oldToNew[t.v[0]];
            t.v[1] = oldToNew[t.v[1]];
            t.v[2] = oldToNew[t.v[2]];
        if (modifiedIndices.count(i) == 0) {
            outData.triangles.push_back(t);
            outData.tags.push_back(meshData.tags[i]);
            continue;
        }
        std::vector<Triangle> subTris = {t};
        ProjectionFrame frame = computeSharedFrame(meshData.normals[i], meshData.centres[i]);
        auto itC = CCuts.find(i);
        auto itNC = NCCuts.find(i);
        std::vector<PolyLine> coplanarCuts;
        PolyLine noncoplanarCuts;
        if (itC != CCuts.end()){
            coplanarCuts = itC->second;
        }
        if (itNC != NCCuts.end()){
            noncoplanarCuts = itNC->second;
        }
        std::vector<Triangle> holeTriangles;
        std::vector<Triangle> modifiedTris;
        std::vector<PolyLine> holesNC;
        findCycles(noncoplanarCuts, nodeGrid, holesNC);
        for (const PolyLine& hole : coplanarCuts) {
            if (hole.size() == 3) {
                Triangle holeTri {nodeGrid.getOrAdd(hole[0].first), nodeGrid.getOrAdd(hole[1].first), nodeGrid.getOrAdd(hole[2].first)};
                holeTriangles.push_back(holeTri);
                modifiedTris.push_back(holeTri);
            }
            else {
            std::vector<Triangle> currHoleTris;
            triangulate(hole, {}, frame, nodeGrid, eps, currHoleTris);
            holeTriangles.insert(holeTriangles.end(), currHoleTris.begin(), currHoleTris.end());
            modifiedTris.insert(modifiedTris.end(), currHoleTris.begin(), currHoleTris.end());
            }
        }
        subTris.insert(subTris.end(), holeTriangles.begin(), holeTriangles.end());
        std::vector<Triangle> cutTrisOut;
        cutTriangles(subTris, nodeGrid.getUniquePoints(), noncoplanarCuts, frame, nodeGrid, eps, cutTrisOut);
        if (removeTouchingSurfaces) {
            std::vector<PolyLine> holes = coplanarCuts;
            holes.insert(holes.end(), holesNC.begin(), holesNC.end());
            std::vector<Triangle> goodTris;
            for (const Triangle& tri : cutTrisOut) {
                if (!isCentroidInHole(tri.centre(nodeGrid.getUniquePoints()), holes, frame)) {
                    goodTris.push_back(tri);
                }
            }
            subTris = std::move(goodTris);
        } else {
            subTris = std::move(cutTrisOut);
        }
        
        if (!subTris.empty()) {
            outData.triangles.push_back(subTris[0]);
            outData.tags.push_back(meshData.tags[i]);
            extraTriangles.insert(extraTriangles.end(), subTris.begin() + 1, subTris.end());
        }
    }
    for (size_t i = 0; i < extraTriangles.size(); ++i) {
        outData.triangles.push_back(std::move(extraTriangles[i]));
        outData.tags.push_back(tagOffset + i);
    }
    outData.nodes = std::move(nodeGrid.getUniquePoints());
    cleanMesh(outData, eps);
    return outData;
}