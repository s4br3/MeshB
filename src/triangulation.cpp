#include "triangulation.hpp"
#include "CDTUtils.h"
#include "math_utils.hpp"
#include "mesh_clean.hpp"
#include "mesh_types.hpp"
#include <fstream>
#include <iostream>
#include <unordered_map>


void debugDumpToObj(const std::vector<Triangle>& tris, 
                    const std::vector<Vec3>& nodes, 
                    const std::string& filename, 
                    double eps) 
{
    std::ofstream out(filename);
    if (!out.is_open()) return;
    SpatialGrid3D debugGrid(eps);
    std::vector<size_t> objIndices;
    objIndices.reserve(tris.size() * 3);
    for (const Triangle& tri : tris) {
        for (int j = 0; j < 3; ++j) {
            objIndices.push_back(debugGrid.getOrAdd(nodes[tri.v[j]]));
        }
    }
    const std::vector<Vec3>& uniquePts = debugGrid.getUniquePoints();
    for (const Vec3& p : uniquePts) {
        out << "v " << p.x << " " << p.y << " " << p.z << "\n";
    }
    for (size_t i = 0; i < tris.size(); ++i) {
        out << "f " << (objIndices[i * 3 + 0] + 1) << " "
                    << (objIndices[i * 3 + 1] + 1) << " "
                    << (objIndices[i * 3 + 2] + 1) << "\n";
    }
}
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
    for (size_t i = 0; i < segs.size(); ++i) {
        CDT::V2d<double> A = initialPts[segs[i].first];
        CDT::V2d<double> B = initialPts[segs[i].second];
        pointsOnSegments[i].push_back(grid.getOrAdd(A));
        pointsOnSegments[i].push_back(grid.getOrAdd(B));
    }
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
std::vector<Triangle> triangulate(
    const PolyLine& polygonSegments, const PolyLine& cuts,
    const ProjectionFrame& frame, SpatialGrid3D& nodeGrid,
    const double eps)
{
    if (polygonSegments.empty() && cuts.empty()) return {};
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
    CDT::RemoveDuplicatesAndRemapEdges(uniquePts, CDTEdges);
    std::set<std::pair<CDT::VertInd, CDT::VertInd>> cleanEdges;
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
    const double eps)
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

    return triangulate(segs, cuts, frame, nodeGrid, eps);
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
std::vector<PolyLine> findCycles(const PolyLine& edges, SpatialGrid3D& nodeGrid) {
    if (edges.empty()) return {};
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
    if (adjacencyMap.empty()) return {};
    std::set<EdgeKey> visitedEdges;
    std::vector<std::vector<size_t>> allCycles;
    for (const auto& [startNode, _] : adjacencyMap) {
        std::vector<size_t> pathStack = {startNode};
        dfsEdges(startNode, adjacencyMap, visitedEdges, pathStack, allCycles);
    }
    std::vector<PolyLine> loops;
    const auto& uniquePoints = nodeGrid.getUniquePoints();
    for (const auto& cycleIndices : allCycles) {
        if (cycleIndices.size() < 3) continue;
        PolyLine loopPoly;
        for (size_t i = 0; i < cycleIndices.size(); ++i) {
            size_t u = cycleIndices[i];
            size_t v = cycleIndices[(i + 1) % cycleIndices.size()];
            loopPoly.push_back({uniquePoints[u], uniquePoints[v]});
        }
        loops.push_back(loopPoly);
    }
    return loops;
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
    const double eps, bool removeTouchingSurfaces, bool meshA)
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
        std::vector<PolyLine> holesNC = findCycles(noncoplanarCuts, nodeGrid);
        for (const PolyLine& hole : coplanarCuts) {
            if (hole.size() == 3) {
                Triangle holeTri {nodeGrid.getOrAdd(hole[0].first), nodeGrid.getOrAdd(hole[1].first), nodeGrid.getOrAdd(hole[2].first)};
                holeTriangles.push_back(holeTri);
                modifiedTris.push_back(holeTri);
            }
            else {
            std::vector<Triangle> currHoleTris = triangulate(hole, {}, frame, nodeGrid, eps);
            holeTriangles.insert(holeTriangles.end(), currHoleTris.begin(), currHoleTris.end());
            modifiedTris.insert(modifiedTris.end(), currHoleTris.begin(), currHoleTris.end());
            }
        }
        subTris.insert(subTris.end(), holeTriangles.begin(), holeTriangles.end());
        subTris = cutTriangles(subTris, nodeGrid.getUniquePoints(), noncoplanarCuts, frame, nodeGrid, eps);
        if(removeTouchingSurfaces){
            std::vector<PolyLine> holes = coplanarCuts;
            holes.insert(holes.end(),holesNC.begin(), holesNC.end());
            std::vector<Triangle> goodTris;
            for (const Triangle& tri : subTris){
                if (!isCentroidInHole(tri.centre(nodeGrid.getUniquePoints()), holes, frame)){
                    goodTris.push_back(tri);
                }
            }
            subTris = std::move(goodTris);
        }
        if(!subTris.empty()){
            outData.triangles.push_back(subTris[0]);
            outData.tags.push_back(meshData.tags[i]);
            extraTriangles.insert(extraTriangles.end(), subTris.begin() + 1, subTris.end());
            std::string s = meshA?"A":"B";
            std::string filename = "tri" + s + std::to_string(meshData.tags[i]) + ".obj";
            debugDumpToObj(modifiedTris, nodeGrid.getUniquePoints(), filename, eps);
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