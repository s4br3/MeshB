#include "triangulation.hpp"
#include "CDTUtils.h"
#include "bvh.hpp"
#include "bvh_collisions.hpp"
#include "math_utils.hpp"
#include "mesh_clean.hpp"
#include "mesh_types.hpp"
#include <unordered_map>
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
void intersect2DAllPoints(
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
            if (std::abs(t) <= eps) {
                addIdUnique(grid.getOrAdd(A));
            } else if (std::abs(t - 1.0) <= eps) {
                addIdUnique(grid.getOrAdd(B));
            }
            else if (std::abs(u) <= eps) {
                addIdUnique(grid.getOrAdd(C));
            } else if (std::abs(u - 1.0) <= eps) {
                addIdUnique(grid.getOrAdd(D));
            }
            else {
                addIdUnique(grid.getOrAdd( {A.x + t * x12, A.y + t * y12} ));
            }
        }
        return;
    }
    if (std::abs(x13 * y12 - y13 * x12) > eps) return;
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
    std::vector<BBox> boxes;
    std::vector<Vec3> centres;
    std::vector<std::vector<CDT::VertInd>> pointsOnSegments(segs.size());
    std::vector<std::pair<size_t, size_t>> stack;
    for (size_t i = 0; i < segs.size(); ++i) {
        CDT::V2d<double> A = initialPts[segs[i].first];
        CDT::V2d<double> B = initialPts[segs[i].second];
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
                    CDT::V2d<double> A = initialPts[segs[prim1Id].first], B = initialPts[segs[prim1Id].second];
                    CDT::V2d<double> C = initialPts[segs[prim2Id].first], D = initialPts[segs[prim2Id].second];
                    std::vector<CDT::VertInd> outs;
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
            parametrized.push_back( {t, id} );
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
void triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    const double eps, std::vector<Triangle>& outTriangles)
{
    if (polygonSegments.empty() && cuts.empty()) return;
    std::vector<CDT::V2d<double>> initialPts;
    std::vector<std::pair<size_t, size_t>> segs;
    for (const std::pair<Vec3, Vec3>& seg : polygonSegments) {
        size_t idx1 = nodeGrid.getOrAdd(seg.first);
        size_t idx2 = nodeGrid.getOrAdd(seg.second);
        if (idx1 != idx2) segs.push_back( {idx1, idx2} );
    }
    std::unordered_set<std::pair<size_t, size_t>> cutEdge;
    for (const std::pair<Vec3, Vec3>& cut : cuts) {
        size_t idx1 = nodeGrid.getOrAdd(cut.first);
        size_t idx2 = nodeGrid.getOrAdd(cut.second);
        if (idx1 != idx2) {
            cutEdge.insert(makeEdgeKey(idx1, idx2));
        }
    }
    for (const std::pair<size_t, size_t>& key : cutEdge) {
        segs.push_back(key);
    }
    for (const Vec3& vec : nodeGrid.getUniquePoints()) initialPts.push_back(frame.to2D(vec));
    std::vector<CDT::V2d<double>> uniquePts;
    std::vector<CDT::Edge> CDTEdges;
    buildSubdividedEdges(initialPts, segs, uniquePts, CDTEdges, eps);
    std::unordered_set<std::pair<CDT::VertInd, CDT::VertInd>> cleanEdges;
    std::vector<CDT::Edge> deduplicatedCDTEdges;
    for (const CDT::Edge& edge : CDTEdges) {
        CDT::VertInd v1 = edge.v1();
        CDT::VertInd v2 = edge.v2();
        if (v1 == v2) continue;
        if (cleanEdges.insert(makeEdgeKey(v1, v2)).second) {
            deduplicatedCDTEdges.push_back(CDT::Edge {v1, v2} );
        }
    }
    CDTEdges = std::move(deduplicatedCDTEdges);
    CDT_Triangulation cdt;
    cdt.insertVertices(uniquePts);
    cdt.insertEdges(CDTEdges);
    cdt.eraseOuterTriangles();
    std::vector<size_t> localToGlobal(cdt.vertices.size());
    for (size_t i = 0; i < cdt.vertices.size(); ++i) {
        localToGlobal[i] = nodeGrid.getOrAdd(frame.to3D(cdt.vertices[i]));
    }
    
    outTriangles.reserve(outTriangles.size() + cdt.triangles.size());
    for (const CDT::Triangle& tri : cdt.triangles) {
        Triangle newTri;
        newTri.v = {
            localToGlobal[tri.vertices[0]],
            localToGlobal[tri.vertices[1]],
            localToGlobal[tri.vertices[2]]
        };
        outTriangles.push_back(newTri);
    }
}
void cutTriangles(
    const std::vector<Triangle>& tris, const std::vector<Vec3>& nodes, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    const double eps, std::vector<Triangle>& outTriangles)
{
    std::unordered_map<std::pair<size_t, size_t>, size_t> edgeCounts;
    std::unordered_map<std::pair<size_t, size_t>, std::pair<Vec3, Vec3>> edgeGeom;
    for (const Triangle& tri : tris) {
        for (int i = 0; i < 3; ++i) {
            size_t u = tri.v[i];
            size_t v = tri.v[(i + 1) % 3];
            std::pair<size_t, size_t> key = makeEdgeKey(u, v);
            edgeCounts[key]++;
            if (edgeCounts[key] == 1) {
                edgeGeom[key] = {nodes[u], nodes[v]};
            }
        }
    }
    PolyLine segs;
    for (const auto& [key, count] : edgeCounts) {
        segs.push_back(edgeGeom[key]);
    }
    triangulate(segs, cuts, frame, nodeGrid, eps, outTriangles);
}
void dfsEdges(
    size_t curr,
    const std::unordered_map<size_t, std::vector<size_t>>& adjacencies,
    std::set<EdgeKey>& visitedEdges,
    std::vector<size_t>& pathStack,
    std::vector<std::vector<size_t>>& allCycles)
{
    auto it = adjacencies.find(curr);
    if (it == adjacencies.end()) return;
    for (size_t neighbour : it->second) {
        EdgeKey key = makeEdgeKey(curr, neighbour);
        if (visitedEdges.count(key)) continue;
        visitedEdges.insert(key);
        auto pathIt = std::find(pathStack.begin(), pathStack.end(), neighbour);
        if (pathIt != pathStack.end()) {
            std::vector<size_t> cycle(pathIt, pathStack.end());
            allCycles.push_back(cycle);
        } else {
            pathStack.push_back(neighbour);
            dfsEdges(neighbour, adjacencies, visitedEdges, pathStack, allCycles);
            pathStack.pop_back();
        }
    }
}
void findCycles(const PolyLine& edges, SpatialGrid3D& nodeGrid, std::vector<PolyLine>& outLoops) {
    if (edges.empty()) return;
    std::unordered_map<size_t, std::vector<size_t>> adjacencyMap;
    for (const std::pair<Vec3, Vec3>& seg : edges) {
        size_t u = nodeGrid.getOrAdd(seg.first);
        size_t v = nodeGrid.getOrAdd(seg.second);
        if (u == v) continue;
        adjacencyMap[u].push_back(v);
        adjacencyMap[v].push_back(u);
    }
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<size_t> toRemove;
        for (const auto& [node, neighbours] : adjacencyMap) {
            if (neighbours.size() <= 1) {
                toRemove.push_back(node);
            }
        }
        for (size_t node : toRemove) {
            auto itNode = adjacencyMap.find(node);
            if (itNode == adjacencyMap.end()) continue;
            for (size_t nb : itNode->second) {
                auto itNb = adjacencyMap.find(nb);
                if (itNb != adjacencyMap.end()) {
                    auto& vec = itNb->second;
                    vec.erase(std::remove(vec.begin(), vec.end(), node), vec.end());
                }
            }
            adjacencyMap.erase(itNode);
            changed = true;
        }
    }
    if (adjacencyMap.empty()) return;
    std::set<EdgeKey> visitedEdges;
    std::vector<std::vector<size_t>> allCycles;
    for (const auto& [startNode, _] : adjacencyMap) {
        std::vector<size_t> pathStack = {startNode};
        dfsEdges(startNode, adjacencyMap, visitedEdges, pathStack, allCycles);
    }
    
    const auto& uniquePoints = nodeGrid.getUniquePoints();
    for (const auto& cycleIndices : allCycles) {
        if (cycleIndices.size() < 3) continue;
        PolyLine loopPoly;
        for (size_t i = 0; i < cycleIndices.size(); ++i) {
            size_t u = cycleIndices[i];
            size_t v = cycleIndices[(i + 1) % cycleIndices.size()];
            loopPoly.push_back({uniquePoints[u], uniquePoints[v]});
        }
        outLoops.push_back(loopPoly);
    }
}
bool isCentroidInHole(
    const Vec3 centre,
    const std::vector<PolyLine>& allHoles, 
    const ProjectionFrame& frame) 
{
    CDT::V2d<double> pt = frame.to2D(centre);
    bool inside = false;
    for (const PolyLine& hole : allHoles) {
        for (const auto& edge : hole) {
            CDT::V2d<double> v1 = frame.to2D(edge.first);
            CDT::V2d<double> v2 = frame.to2D(edge.second);
            if (((v1.y > pt.y) != (v2.y > pt.y)) &&
                (pt.x < (v2.x - v1.x) * (pt.y - v1.y) / (v2.y - v1.y) + v1.x)) 
            {
                inside = !inside;
            }
        }
    }
    return inside;
}
MeshData cutMesh(
    const MeshData& meshData,
    const std::unordered_map<size_t, PolyLine>& NCCuts,
    const std::unordered_map<size_t, std::vector<PolyLine>>& CCuts,
    const double eps, bool removeTouchingSurfaces)
{
    size_t tagOffset = 0;
    for (size_t tag : meshData.tags) {
        tagOffset = std::max(tagOffset, tag + 1);
    }
    struct IterationResult {
        std::vector<TriVerts> mainTris;
        std::vector<TriVerts> extraTris;
    };
    const size_t numTris = meshData.triangles.size();
    std::vector<IterationResult> results(numTris);
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < numTris; ++i) {
        auto itNC = NCCuts.find(i);
        auto itC = CCuts.find(i);
        bool isModified = (itNC != NCCuts.end()) || (itC != CCuts.end());
        if (!isModified) {
            results[i].mainTris.push_back(meshData.triangles[i].getVertices(meshData.nodes));
            continue;
        }
        SpatialGrid3D localGrid(eps);
        Triangle t;
        t.v[0] = localGrid.getOrAdd(meshData.nodes[meshData.triangles[i].v[0]]);
        t.v[1] = localGrid.getOrAdd(meshData.nodes[meshData.triangles[i].v[1]]);
        t.v[2] = localGrid.getOrAdd(meshData.nodes[meshData.triangles[i].v[2]]);
        std::vector<Triangle> subTris;
        subTris.push_back(t);
        ProjectionFrame frame = computeSharedFrame(meshData.normals[i], meshData.centres[i]);
        std::vector<PolyLine> coplanarCuts;
        PolyLine noncoplanarCuts;
        if (itC != CCuts.end()) {
            coplanarCuts = itC->second;
        }
        if (itNC != NCCuts.end()) {
            noncoplanarCuts = itNC->second;
        }
        std::vector<PolyLine> holesNC;
        findCycles(noncoplanarCuts, localGrid, holesNC);
        std::vector<Triangle> holeTriangles;
        for (const PolyLine& hole : coplanarCuts) {
            if (hole.size() == 3) {
                Triangle holeTri {
                    localGrid.getOrAdd(hole[0].first),
                    localGrid.getOrAdd(hole[1].first),
                    localGrid.getOrAdd(hole[2].first)
                };
                holeTriangles.push_back(holeTri);
            } else {
                triangulate(hole, {}, frame, localGrid, eps, holeTriangles);
            }
        }
        subTris.insert(subTris.end(), holeTriangles.begin(), holeTriangles.end());
        std::vector<Triangle> cutTrisOut;
        cutTriangles(subTris, localGrid.getUniquePoints(), noncoplanarCuts, frame, localGrid, eps, cutTrisOut);
        if (removeTouchingSurfaces) {
            std::vector<PolyLine> holes = coplanarCuts;
            holes.insert(holes.end(), holesNC.begin(), holesNC.end());
            subTris.clear();
            for (const Triangle& tri : cutTrisOut) {
                if (!isCentroidInHole(tri.centre(localGrid.getUniquePoints()), holes, frame)) {
                    subTris.push_back(tri);
                }
            }
        } else {
            subTris = std::move(cutTrisOut);
        }
        if (!subTris.empty()) {
            const auto& localPoints = localGrid.getUniquePoints();
            results[i].mainTris.push_back(subTris[0].getVertices(localPoints));
            for (size_t k = 1; k < subTris.size(); ++k) {
                results[i].extraTris.push_back(subTris[k].getVertices(localPoints));
            }
        }
    }
    SpatialGrid3D globalGrid(eps);
    MeshData outData;
    for (size_t i = 0; i < numTris; ++i) {
        for (const TriVerts& tri : results[i].mainTris) {
            outData.triangles.emplace_back(
                globalGrid.getOrAdd(tri[0]),
                globalGrid.getOrAdd(tri[1]),
                globalGrid.getOrAdd(tri[2])
            );
            outData.tags.push_back(meshData.tags[i]);
        }
    }
    size_t extraCounter = 0;
    for (size_t i = 0; i < numTris; ++i) {
        for (const TriVerts& tri : results[i].extraTris) {
            outData.triangles.emplace_back(
                globalGrid.getOrAdd(tri[0]),
                globalGrid.getOrAdd(tri[1]),
                globalGrid.getOrAdd(tri[2])
            );
            outData.tags.push_back(tagOffset + extraCounter);
            extraCounter++;
        }
    }
    outData.nodes = std::move(globalGrid.getUniquePoints());
    cleanMesh(outData, eps);
    return outData;
}