#include "triangulation.hpp"
#include "bvh.hpp"
#include "geom_2d.hpp"
#include "geom_3d.hpp"
#include "mesh_clean.hpp"
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <gmsh.h>
void buildSubdividedEdges(
    const std::vector<Vec2>& initialPts,
    const std::vector<std::pair<size_t, size_t>>& segs,
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
    std::vector<std::pair<size_t, size_t>> stack;
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
                    Vec2 A = initialPts[segs[prim1Id].first], B = initialPts[segs[prim1Id].second];
                    Vec2 C = initialPts[segs[prim2Id].first], D = initialPts[segs[prim2Id].second];
                    std::vector<size_t> outs;
                    intersect2DAllPoints(A, B, C, D, grid, outs, eps);
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
    if (!gmsh::isInitialized()) {
        gmsh::initialize(0, nullptr, false);
        gmsh::option::setNumber("General.Terminal", 0); 
        gmsh::option::setNumber("General.NumThreads", 1);
    }
    gmsh::clear();
    gmsh::model::add("cdt_patch");
    double minX = INFINITY, minY = INFINITY, maxX = -INFINITY, maxY = -INFINITY;
    for (const auto& pt : uniquePts) {
        if (pt.x < minX) minX = pt.x;
        if (pt.y < minY) minY = pt.y;
        if (pt.x > maxX) maxX = pt.x;
        if (pt.y > maxY) maxY = pt.y;
    }
    double pad = std::max(eps * 10.0, std::max(maxX - minX, maxY - minY) * 0.1);
    if (pad <= 1e-9) pad = 1.0;
    minX -= pad; minY -= pad;
    maxX += pad; maxY += pad;
    int bp1 = gmsh::model::geo::addPoint(minX, minY, 0);
    int bp2 = gmsh::model::geo::addPoint(maxX, minY, 0);
    int bp3 = gmsh::model::geo::addPoint(maxX, maxY, 0);
    int bp4 = gmsh::model::geo::addPoint(minX, maxY, 0);
    int cl = gmsh::model::geo::addCurveLoop({
        gmsh::model::geo::addLine(bp1, bp2),
        gmsh::model::geo::addLine(bp2, bp3),
        gmsh::model::geo::addLine(bp3, bp4),
        gmsh::model::geo::addLine(bp4, bp1)
    });
    int surf = gmsh::model::geo::addPlaneSurface({cl});
    std::unordered_map<std::size_t, size_t> gmshTagToLocalIdx;
    std::vector<int> gmshPtTags(uniquePts.size());
    for (size_t i = 0; i < uniquePts.size(); ++i) {
        int tag = gmsh::model::geo::addPoint(uniquePts[i].x, uniquePts[i].y, 0);
        gmshPtTags[i] = tag;
        gmshTagToLocalIdx[tag] = i; 
    }
    std::vector<int> gmshLineTags;
    std::vector<bool> pointUsedInLine(uniquePts.size(), false);
    for (const EdgeKey& edge : CDTEdges) {
        gmshLineTags.push_back(gmsh::model::geo::addLine(gmshPtTags[edge.first], gmshPtTags[edge.second]));
        pointUsedInLine[edge.first] = true;
        pointUsedInLine[edge.second] = true;
    }
    std::vector<int> isolatedPtTags;
    for (size_t i = 0; i < uniquePts.size(); ++i) {
        if (!pointUsedInLine[i]) {
            isolatedPtTags.push_back(gmshPtTags[i]);
        }
    }
    gmsh::model::geo::synchronize();
    if (!gmshLineTags.empty()) {
        gmsh::model::mesh::embed(1, gmshLineTags, 2, surf);
    }
    if (!isolatedPtTags.empty()) {
        gmsh::model::mesh::embed(0, isolatedPtTags, 2, surf);
    }
    if (constrained){
        gmsh::option::setNumber("Mesh.MeshSizeMin", 1e22);
        gmsh::option::setNumber("Mesh.MeshSizeMax", 1e22);
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
        gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
        gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    }
    gmsh::option::setNumber("Mesh.Algorithm", 5); 
    gmsh::model::mesh::generate(2);
    std::vector<std::size_t> outNodeTags;
    std::vector<double> outNodeCoords, outNodeParametricCoords;
    gmsh::model::mesh::getNodes(outNodeTags, outNodeCoords, outNodeParametricCoords, -1, -1);
    for (size_t i = 0; i < outNodeTags.size(); ++i) {
        std::size_t tag = outNodeTags[i];
        if (gmshTagToLocalIdx.find(tag) == gmshTagToLocalIdx.end()) {
            Vec2 newPt = { outNodeCoords[i * 3], outNodeCoords[i * 3 + 1] };
            uniquePts.push_back(newPt);
            gmshTagToLocalIdx[tag] = uniquePts.size() - 1;
        }
    }
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags, nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2, surf);
    for (size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == 2) { 
            const auto& nTags = nodeTags[i];
            for (size_t t = 0; t < nTags.size() / 3; ++t) {
                std::size_t n1 = nTags[t * 3 + 0];
                std::size_t n2 = nTags[t * 3 + 1];
                std::size_t n3 = nTags[t * 3 + 2];
                if (gmshTagToLocalIdx.find(n1) == gmshTagToLocalIdx.end() ||
                    gmshTagToLocalIdx.find(n2) == gmshTagToLocalIdx.end() ||
                    gmshTagToLocalIdx.find(n3) == gmshTagToLocalIdx.end()) {
                    continue;
                }
                size_t idx1 = gmshTagToLocalIdx[n1];
                size_t idx2 = gmshTagToLocalIdx[n2];
                size_t idx3 = gmshTagToLocalIdx[n3];
                Vec2 centroid = {
                    (uniquePts[idx1].x + uniquePts[idx2].x + uniquePts[idx3].x) / 3.0,
                    (uniquePts[idx1].y + uniquePts[idx2].y + uniquePts[idx3].y) / 3.0
                };
                if (isPointInsidePolygon(centroid, boundary2D)) {
                    Triangle newTri;
                    newTri.v = {
                        nodeGrid.getOrAdd(frame.to3D(uniquePts[idx1])),
                        nodeGrid.getOrAdd(frame.to3D(uniquePts[idx2])),
                        nodeGrid.getOrAdd(frame.to3D(uniquePts[idx3]))
                    };
                    outTriangles.push_back(newTri);
                }
            }
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