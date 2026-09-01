#include "cdt.hpp"
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <optional>
#include <array>
static TriangleCDT makeCCWTriangle(Vec2 p1, Vec2 p2, Vec2 p3) {
    if (cross2D(p2 - p1, p3 - p1) < 0) {
        std::swap(p2, p3);
    }
    return TriangleCDT({p1, p2, p3});
}
static size_t getNeighbourAcrossEdge(const TriangleCDT& t, const Edge& e, double eps) {
    bool has1 = e.first.equals(t.p1, eps) || e.second.equals(t.p1, eps);
    bool has2 = e.first.equals(t.p2, eps) || e.second.equals(t.p2, eps);
    bool has3 = e.first.equals(t.p3, eps) || e.second.equals(t.p3, eps);
    if (has2 && has3) return t.n1;
    if (has3 && has1) return t.n2;
    if (has1 && has2) return t.n3;
    return NO_NEIGHBOUR;
}
static void setNeighbourAcrossEdge(TriangleCDT& t, const Edge& e, size_t neighbourIdx, double eps) {
    bool has1 = e.first.equals(t.p1, eps) || e.second.equals(t.p1, eps);
    bool has2 = e.first.equals(t.p2, eps) || e.second.equals(t.p2, eps);
    bool has3 = e.first.equals(t.p3, eps) || e.second.equals(t.p3, eps);
    if (has2 && has3) t.n1 = neighbourIdx;
    else if (has3 && has1) t.n2 = neighbourIdx;
    else if (has1 && has2) t.n3 = neighbourIdx;
}
static bool sameEdge(const Edge& e1, const Edge& e2, double eps) {
    return (e1.first.equals(e2.first, eps) && e1.second.equals(e2.second, eps)) ||
           (e1.first.equals(e2.second, eps) && e1.second.equals(e2.first, eps));
}
static bool isConstraintEdge(const Edge& e, const std::vector<Edge>& constraints, double eps) {
    for (const auto& c : constraints) {
        if (sameEdge(e, c, eps)) return true;
    }
    return false;
}
size_t hashEdgeSpatial(const Edge& e, double eps) {
    auto hashPoint = [eps](const Vec2& p) -> size_t {
        int64_t cx = static_cast<int64_t>(std::floor(p.x / eps));
        int64_t cy = static_cast<int64_t>(std::floor(p.y / eps));
        size_t h = std::hash<int64_t>{}(cx);
        return hashCombine(h, cy);
    };
    size_t h1 = hashPoint(e.first);
    size_t h2 = hashPoint(e.second);
    EdgeKey key = makeEdgeKey(h1, h2);
    return hashCombine(key.first, key.second);
}
bool pointInTriangle(const Vec2& p, const TriangleCDT& t, const double eps) {
    double c1 = cross2D(t.p2 - t.p1, p - t.p1);
    double c2 = cross2D(t.p3 - t.p2, p - t.p2);
    double c3 = cross2D(t.p1 - t.p3, p - t.p3);
    bool has_neg = (c1 < -eps) || (c2 < -eps) || (c3 < -eps);
    bool has_pos = (c1 > eps) || (c2 > eps) || (c3 > eps);
    return !(has_neg && has_pos);
}
static bool pointOnSegmentRobust(const Vec2& p, const Edge& e, double eps) {
    double minX = std::min(e.first.x, e.second.x) - eps;
    double maxX = std::max(e.first.x, e.second.x) + eps;
    double minY = std::min(e.first.y, e.second.y) - eps;
    double maxY = std::max(e.first.y, e.second.y) + eps;
    if (p.x < minX || p.x > maxX || p.y < minY || p.y > maxY) return false;
    
    Vec2 ab = e.second - e.first;
    Vec2 ap = p - e.first;
    double area = std::fabs(cross2D(ab, ap));
    double len = std::sqrt(dot(ab, ab));
    if (len <= eps) return p.equals(e.first, eps);
    return (area / len) <= eps;
}
double comparePoints(const Vec2& a, const Vec2& b, double eps) {
    if (std::abs(a.x - b.x) < eps) return a.y - b.y;
    return a.x - b.x;
}
bool doesEdgeExist(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Edge& e, const double eps) {
    for (size_t idx : activeTris) {
        if (allTris[idx].hasVertex(e.first, eps) && allTris[idx].hasVertex(e.second, eps)) return true;
    }
    return false;
}
bool areSegmentsCrossing(const Edge& e1, const Edge& e2, double eps) {
    Vec2 v12 = e1.second - e1.first;
    Vec2 v34 = e2.second - e2.first;
    Vec2 v13 = e2.first - e1.first;
    double denom = cross2D(v12, v34);
    if (std::fabs(denom) <= eps) return false; 
    double t = cross2D(v13, v34) / denom;
    double u = cross2D(v13, v12) / denom;
    return (t > eps && t < 1 - eps) && (u > eps && u < 1 - eps);
}
std::vector<Edge> findIntersectingEdges(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Edge& e, const double eps) {
    const Vec2& A = e.first;
    const Vec2& B = e.second;
    size_t startTriIdx = NO_NEIGHBOUR;
    for (size_t idx : activeTris) {
        if (allTris[idx].hasVertex(A, eps)) {
            startTriIdx = idx;
            break;
        }
    }
    if (startTriIdx == NO_NEIGHBOUR) return {};
    size_t currentTriIdx = NO_NEIGHBOUR;
    Edge exitEdge;
    bool foundStart = false;
    size_t searchTri = startTriIdx;
    std::unordered_set<size_t> visitedAroundA;
    while (searchTri != NO_NEIGHBOUR && activeTris.count(searchTri) && visitedAroundA.count(searchTri) == 0) {
        visitedAroundA.insert(searchTri);
        const TriangleCDT& t = allTris[searchTri];
        Vec2 v1, v2;
        if (t.p1.equals(A, eps)) { v1 = t.p2; v2 = t.p3; }
        else if (t.p2.equals(A, eps)) { v1 = t.p3; v2 = t.p1; }
        else { v1 = t.p1; v2 = t.p2; }
        Edge oppEdge = {v1, v2};
        if (areSegmentsCrossing(oppEdge, e, eps)) {
            currentTriIdx = searchTri;
            exitEdge = oppEdge;
            foundStart = true;
            break;
        }
        size_t next1 = getNeighbourAcrossEdge(t, Edge{A, v1}, eps);
        size_t next2 = getNeighbourAcrossEdge(t, Edge{A, v2}, eps);
        if (next1 != NO_NEIGHBOUR && activeTris.count(next1) && visitedAroundA.count(next1) == 0) {
            searchTri = next1;
        } else if (next2 != NO_NEIGHBOUR && activeTris.count(next2) && visitedAroundA.count(next2) == 0) {
            searchTri = next2;
        } else {
            break;
        }
    }
    if (!foundStart) {
        for (size_t idx : activeTris) {
            const TriangleCDT& t = allTris[idx];
            if (!t.hasVertex(A, eps)) continue;
            
            Vec2 v1, v2;
            if (t.p1.equals(A, eps)) { v1 = t.p2; v2 = t.p3; }
            else if (t.p2.equals(A, eps)) { v1 = t.p3; v2 = t.p1; }
            else { v1 = t.p1; v2 = t.p2; }
            
            Edge oppEdge = {v1, v2};
            if (areSegmentsCrossing(oppEdge, e, eps)) {
                currentTriIdx = idx;
                exitEdge = oppEdge;
                foundStart = true;
                break;
            }
        }
    }
    if (!foundStart) return {};
    std::vector<Edge> intersectingEdges;
    std::unordered_set<size_t> visitedWalk;
    while (currentTriIdx != NO_NEIGHBOUR && activeTris.count(currentTriIdx) && visitedWalk.count(currentTriIdx) == 0) {
        visitedWalk.insert(currentTriIdx);
        intersectingEdges.push_back(exitEdge);
        size_t nextTriIdx = getNeighbourAcrossEdge(allTris[currentTriIdx], exitEdge, eps);
        if (nextTriIdx == NO_NEIGHBOUR || !activeTris.count(nextTriIdx)) {
            break;
        }
        const TriangleCDT& nextTri = allTris[nextTriIdx];
        if (nextTri.hasVertex(B, eps)) {
            break;
        }
        TriangleEdges nEdges = nextTri.getEdges();
        Edge nextExitEdge;
        bool foundNextExit = false;
        for (int i = 0; i < 3; ++i) {
            if (sameEdge(nEdges[i], exitEdge, eps)) continue;
            if (areSegmentsCrossing(nEdges[i], e, eps)) {
                nextExitEdge = nEdges[i];
                foundNextExit = true;
                break;
            }
        }
        if (!foundNextExit) break;
        currentTriIdx = nextTriIdx;
        exitEdge = nextExitEdge;
    }
    return intersectingEdges;
}
size_t findTriangleWithPoint(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Vec2& p, const double eps) {
    for (size_t idx : activeTris) {
        if (pointInTriangle(p, allTris[idx], eps)) return idx;
    }
    for (size_t idx : activeTris) {
        if (pointInTriangle(p, allTris[idx], eps * 100.0)) return idx;
    }
    return NO_NEIGHBOUR;
}
std::optional<std::pair<size_t, size_t>> getTrianglesForEdge(
    const std::vector<TriangleCDT>& allTris, 
    const std::unordered_set<size_t>& activeTris, 
    const Edge& e, 
    const double eps) 
{
    std::vector<size_t> matches;
    for (size_t idx : activeTris) {
        if (allTris[idx].hasVertex(e.first, eps) && allTris[idx].hasVertex(e.second, eps)) {
            matches.push_back(idx);
            if (matches.size() == 2) return std::make_pair(matches[0], matches[1]);
        }
    }
    return std::nullopt;
}
bool isInCircumcircle(const TriangleCDT& t, const Vec2& d, const double eps) {
    Vec2 a = t.p1 - d;
    Vec2 b = t.p2 - d;
    Vec2 c = t.p3 - d;
    double sa = a.x * a.x + a.y * a.y;
    double sb = b.x * b.x + b.y * b.y;
    double sc = c.x * c.x + c.y * c.y;
    if (cross2D(t.p2 - t.p1, t.p3 - t.p1) < 0) {
        std::swap(b, c);
        std::swap(sb, sc);
    }
    double det = sa * (b.x * c.y - c.x * b.y) -
                 sb * (a.x * c.y - c.x * a.y) +
                 sc * (a.x * b.y - b.x * a.y);
    return det > eps; 
}
bool checkIfConvexQuadrilateral(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d, const double eps) {
    auto crossSign = [](const Vec2& p, const Vec2& lineStart, const Vec2& lineEnd) {
        return cross2D(lineEnd - lineStart, p - lineStart);
    };
    double s1 = crossSign(a, b, d);
    double s2 = crossSign(c, b, d);
    double s3 = crossSign(b, a, c);
    double s4 = crossSign(d, a, c);

    return (s1 * s2 < -eps) && (s3 * s4 < -eps);
}
std::array<Vec2, 4> getQuadrilateral(const std::vector<TriangleCDT>& allTris, const std::pair<size_t, size_t>& e, double eps) {
    const TriangleCDT& t1 = allTris[e.first];
    const TriangleCDT& t2 = allTris[e.second];
    auto nextCCW = [&eps](const TriangleCDT& t, const Vec2& p) -> Vec2 {
        if (t.p1.equals(p, eps)) return t.p2;
        if (t.p2.equals(p, eps)) return t.p3;
        if (t.p3.equals(p, eps)) return t.p1;
        throw std::runtime_error("Point not in triangle.");
    };
    Vec2 q0 = t1.p1;
    if (t2.hasVertex(q0, eps)) q0 = t1.p2;
    if (t2.hasVertex(q0, eps)) q0 = t1.p3;
    Vec2 q1 = nextCCW(t1, q0);
    Vec2 q3 = nextCCW(t1, q1);
    Vec2 q2 = t2.p1;
    if (t1.hasVertex(q2, eps)) q2 = t2.p2;
    if (t1.hasVertex(q2, eps)) q2 = t2.p3;
    return {q0, q1, q2, q3};
}
static void flipEdgeAndMaintainNeighbours(
    std::vector<TriangleCDT>& allTris,
    std::unordered_set<size_t>& activeTris,
    size_t t1Idx,
    size_t t2Idx,
    const std::array<Vec2, 4>& quad,
    double eps)
{
    size_t n01 = getNeighbourAcrossEdge(allTris[t1Idx], Edge{quad[0], quad[1]}, eps);
    size_t n12 = getNeighbourAcrossEdge(allTris[t2Idx], Edge{quad[1], quad[2]}, eps);
    size_t n23 = getNeighbourAcrossEdge(allTris[t2Idx], Edge{quad[2], quad[3]}, eps);
    size_t n30 = getNeighbourAcrossEdge(allTris[t1Idx], Edge{quad[3], quad[0]}, eps);
    activeTris.erase(t1Idx);
    activeTris.erase(t2Idx);
    size_t newIdx1 = allTris.size();
    allTris.push_back(makeCCWTriangle(quad[0], quad[1], quad[2]));
    activeTris.insert(newIdx1);
    size_t newIdx2 = allTris.size();
    allTris.push_back(makeCCWTriangle(quad[0], quad[2], quad[3]));
    activeTris.insert(newIdx2);
    setNeighbourAcrossEdge(allTris[newIdx1], Edge{quad[0], quad[2]}, newIdx2, eps);
    setNeighbourAcrossEdge(allTris[newIdx2], Edge{quad[0], quad[2]}, newIdx1, eps);
    if (n01 != NO_NEIGHBOUR && activeTris.count(n01)) {
        setNeighbourAcrossEdge(allTris[newIdx1], Edge{quad[0], quad[1]}, n01, eps);
        setNeighbourAcrossEdge(allTris[n01], Edge{quad[0], quad[1]}, newIdx1, eps);
    }
    if (n12 != NO_NEIGHBOUR && activeTris.count(n12)) {
        setNeighbourAcrossEdge(allTris[newIdx1], Edge{quad[1], quad[2]}, n12, eps);
        setNeighbourAcrossEdge(allTris[n12], Edge{quad[1], quad[2]}, newIdx1, eps);
    }
    if (n23 != NO_NEIGHBOUR && activeTris.count(n23)) {
        setNeighbourAcrossEdge(allTris[newIdx2], Edge{quad[2], quad[3]}, n23, eps);
        setNeighbourAcrossEdge(allTris[n23], Edge{quad[2], quad[3]}, newIdx2, eps);
    }
    if (n30 != NO_NEIGHBOUR && activeTris.count(n30)) {
        setNeighbourAcrossEdge(allTris[newIdx2], Edge{quad[3], quad[0]}, n30, eps);
        setNeighbourAcrossEdge(allTris[n30], Edge{quad[3], quad[0]}, newIdx2, eps);
    }
}
static std::array<size_t, 3> splitTriangle(
    std::vector<TriangleCDT>& allTris,
    std::unordered_set<size_t>& activeTris,
    size_t tIdx,
    const Vec2& pt,
    double eps) 
{
    TriangleCDT old = allTris[tIdx];
    activeTris.erase(tIdx);
    size_t t1 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, old.p1, old.p2));
    size_t t2 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, old.p2, old.p3));
    size_t t3 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, old.p3, old.p1));
    activeTris.insert(t1);
    activeTris.insert(t2);
    activeTris.insert(t3);
    setNeighbourAcrossEdge(allTris[t1], Edge{pt, old.p1}, t3, eps);
    setNeighbourAcrossEdge(allTris[t1], Edge{pt, old.p2}, t2, eps);
    setNeighbourAcrossEdge(allTris[t2], Edge{pt, old.p2}, t1, eps);
    setNeighbourAcrossEdge(allTris[t2], Edge{pt, old.p3}, t3, eps);
    setNeighbourAcrossEdge(allTris[t3], Edge{pt, old.p3}, t2, eps);
    setNeighbourAcrossEdge(allTris[t3], Edge{pt, old.p1}, t1, eps);
    if (old.n3 != NO_NEIGHBOUR && activeTris.count(old.n3)) {
        setNeighbourAcrossEdge(allTris[t1], Edge{old.p1, old.p2}, old.n3, eps);
        setNeighbourAcrossEdge(allTris[old.n3], Edge{old.p1, old.p2}, t1, eps);
    }
    if (old.n1 != NO_NEIGHBOUR && activeTris.count(old.n1)) {
        setNeighbourAcrossEdge(allTris[t2], Edge{old.p2, old.p3}, old.n1, eps);
        setNeighbourAcrossEdge(allTris[old.n1], Edge{old.p2, old.p3}, t2, eps);
    }
    if (old.n2 != NO_NEIGHBOUR && activeTris.count(old.n2)) {
        setNeighbourAcrossEdge(allTris[t3], Edge{old.p3, old.p1}, old.n2, eps);
        setNeighbourAcrossEdge(allTris[old.n2], Edge{old.p3, old.p1}, t3, eps);
    }
    return {t1, t2, t3};
}
static std::vector<size_t> splitEdge(
    std::vector<TriangleCDT>& allTris,
    std::unordered_set<size_t>& activeTris,
    size_t t1Idx,
    const Edge& e,
    const Vec2& pt,
    double eps)
{
    TriangleCDT t1 = allTris[t1Idx];
    Vec2 A = e.first;
    Vec2 B = e.second;
    Vec2 C;
    if (!t1.p1.equals(A, eps) && !t1.p1.equals(B, eps)) C = t1.p1;
    else if (!t1.p2.equals(A, eps) && !t1.p2.equals(B, eps)) C = t1.p2;
    else C = t1.p3;
    size_t t2Idx = getNeighbourAcrossEdge(t1, e, eps);
    if (t2Idx == NO_NEIGHBOUR || activeTris.count(t2Idx) == 0) {
        for (size_t idx : activeTris) {
            if (idx != t1Idx && allTris[idx].hasVertex(A, eps) && allTris[idx].hasVertex(B, eps)) {
                t2Idx = idx;
                break;
            }
        }
    }
    activeTris.erase(t1Idx);
    std::vector<size_t> created;
    if (t2Idx != NO_NEIGHBOUR && activeTris.count(t2Idx)) {
        TriangleCDT t2 = allTris[t2Idx];
        activeTris.erase(t2Idx);
        Vec2 D;
        if (!t2.p1.equals(A, eps) && !t2.p1.equals(B, eps)) D = t2.p1;
        else if (!t2.p2.equals(A, eps) && !t2.p2.equals(B, eps)) D = t2.p2;
        else D = t2.p3;
        size_t nAC = getNeighbourAcrossEdge(t1, Edge{A, C}, eps);
        size_t nCB = getNeighbourAcrossEdge(t1, Edge{C, B}, eps);
        size_t nBD = getNeighbourAcrossEdge(t2, Edge{B, D}, eps);
        size_t nDA = getNeighbourAcrossEdge(t2, Edge{D, A}, eps);
        size_t idx1 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, A, C));
        size_t idx2 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, C, B));
        size_t idx3 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, B, D));
        size_t idx4 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, D, A));
        activeTris.insert(idx1); activeTris.insert(idx2);
        activeTris.insert(idx3); activeTris.insert(idx4);
        setNeighbourAcrossEdge(allTris[idx1], Edge{pt, C}, idx2, eps);
        setNeighbourAcrossEdge(allTris[idx2], Edge{pt, C}, idx1, eps);
        setNeighbourAcrossEdge(allTris[idx2], Edge{pt, B}, idx3, eps);
        setNeighbourAcrossEdge(allTris[idx3], Edge{pt, B}, idx2, eps);
        setNeighbourAcrossEdge(allTris[idx3], Edge{pt, D}, idx4, eps);
        setNeighbourAcrossEdge(allTris[idx4], Edge{pt, D}, idx3, eps);
        setNeighbourAcrossEdge(allTris[idx4], Edge{pt, A}, idx1, eps);
        setNeighbourAcrossEdge(allTris[idx1], Edge{pt, A}, idx4, eps);
        if (nAC != NO_NEIGHBOUR && activeTris.count(nAC)) {
            setNeighbourAcrossEdge(allTris[idx1], Edge{A, C}, nAC, eps);
            setNeighbourAcrossEdge(allTris[nAC], Edge{A, C}, idx1, eps);
        }
        if (nCB != NO_NEIGHBOUR && activeTris.count(nCB)) {
            setNeighbourAcrossEdge(allTris[idx2], Edge{C, B}, nCB, eps);
            setNeighbourAcrossEdge(allTris[nCB], Edge{C, B}, idx2, eps);
        }
        if (nBD != NO_NEIGHBOUR && activeTris.count(nBD)) {
            setNeighbourAcrossEdge(allTris[idx3], Edge{B, D}, nBD, eps);
            setNeighbourAcrossEdge(allTris[nBD], Edge{B, D}, idx3, eps);
        }
        if (nDA != NO_NEIGHBOUR && activeTris.count(nDA)) {
            setNeighbourAcrossEdge(allTris[idx4], Edge{D, A}, nDA, eps);
            setNeighbourAcrossEdge(allTris[nDA], Edge{D, A}, idx4, eps);
        }
        created = {idx1, idx2, idx3, idx4};
    } else {
        size_t nAC = getNeighbourAcrossEdge(t1, Edge{A, C}, eps);
        size_t nCB = getNeighbourAcrossEdge(t1, Edge{C, B}, eps);
        size_t idx1 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, A, C));
        size_t idx2 = allTris.size(); allTris.emplace_back(makeCCWTriangle(pt, C, B));
        activeTris.insert(idx1); activeTris.insert(idx2);
        setNeighbourAcrossEdge(allTris[idx1], Edge{pt, C}, idx2, eps);
        setNeighbourAcrossEdge(allTris[idx2], Edge{pt, C}, idx1, eps);
        if (nAC != NO_NEIGHBOUR && activeTris.count(nAC)) {
            setNeighbourAcrossEdge(allTris[idx1], Edge{A, C}, nAC, eps);
            setNeighbourAcrossEdge(allTris[nAC], Edge{A, C}, idx1, eps);
        }
        if (nCB != NO_NEIGHBOUR && activeTris.count(nCB)) {
            setNeighbourAcrossEdge(allTris[idx2], Edge{C, B}, nCB, eps);
            setNeighbourAcrossEdge(allTris[nCB], Edge{C, B}, idx2, eps);
        }
        created = {idx1, idx2};
    }
    return created;
}
std::vector<TriangleCDT> calculateCDT(const std::vector<Vec2>& points, const std::vector<EdgeKey>& edges, double eps) {
    std::vector<TriangleCDT> allTris;
    std::unordered_set<size_t> activeTris;
    if (points.size() < 3) return allTris;
    std::vector<Edge> constraintList;
    constraintList.reserve(edges.size());
    for (const auto& ek : edges) {
        if (ek.first != ek.second && ek.first < points.size() && ek.second < points.size()) {
            constraintList.push_back({points[ek.first], points[ek.second]});
        }
    }
    double minX = points[0].x, minY = points[0].y;
    double maxX = minX, maxY = minY;
    for (const Vec2& p : points) {
        if (p.x < minX) minX = p.x;
        if (p.y < minY) minY = p.y;
        if (p.x > maxX) maxX = p.x;
        if (p.y > maxY) maxY = p.y;
    }
    double dMax = std::max(maxX - minX, maxY - minY);
    if (dMax <= eps) dMax = 1.0;
    double midX = (minX + maxX) / 2.0;
    double midY = (minY + maxY) / 2.0;
    Vec2 p1 = {midX - 100 * dMax, midY - 100 * dMax};
    Vec2 p2 = {midX + 100 * dMax, midY - 100 * dMax}; 
    Vec2 p3 = {midX, midY + 100 * dMax};
    allTris.emplace_back(std::array<Vec2, 3>{p1, p2, p3});
    activeTris.insert(0);
    for (const Vec2& pt : points) {
        size_t tIdx = findTriangleWithPoint(allTris, activeTris, pt, eps);
        if (tIdx == NO_NEIGHBOUR) continue;
        const TriangleCDT& targetTri = allTris[tIdx];
        if (targetTri.p1.equals(pt, eps) || targetTri.p2.equals(pt, eps) || targetTri.p3.equals(pt, eps)) continue;
        Edge edgeOn;
        bool isOnEdge = false;
        TriangleEdges tEdges = targetTri.getEdges();
        for (int i = 0; i < 3; ++i) {
            if (pointOnSegmentRobust(pt, tEdges[i], eps)) {
                edgeOn = tEdges[i];
                isOnEdge = true;
                break;
            }
        }
        std::vector<size_t> newTris;
        if (isOnEdge) {
            newTris = splitEdge(allTris, activeTris, tIdx, edgeOn, pt, eps);
        } else {
            std::array<size_t, 3> splitRes = splitTriangle(allTris, activeTris, tIdx, pt, eps);
            newTris.assign(splitRes.begin(), splitRes.end());
        }
        std::vector<size_t> stack = newTris;
        while (!stack.empty()) {
            size_t currTriIdx = stack.back();
            stack.pop_back();
            if (!activeTris.count(currTriIdx)) continue;
            TriangleCDT& currTri = allTris[currTriIdx];
            Vec2 vA, vB;
            if (currTri.p1.equals(pt, eps)) { vA = currTri.p2; vB = currTri.p3; }
            else if (currTri.p2.equals(pt, eps)) { vA = currTri.p1; vB = currTri.p3; }
            else if (currTri.p3.equals(pt, eps)) { vA = currTri.p1; vB = currTri.p2; }
            else continue; 
            size_t oppNeighbour = getNeighbourAcrossEdge(currTri, Edge{vA, vB}, eps);
            if (oppNeighbour != NO_NEIGHBOUR && activeTris.count(oppNeighbour)) {
                if (isInCircumcircle(allTris[oppNeighbour], pt, eps)) {
                    std::pair<size_t, size_t> adj = {currTriIdx, oppNeighbour};
                    std::array<Vec2, 4> quad;
                    try { 
                        quad = getQuadrilateral(allTris, adj, eps); 
                    } catch (...) { 
                        continue; 
                    }
                    if (checkIfConvexQuadrilateral(quad[0], quad[1], quad[2], quad[3], eps)) {
                        flipEdgeAndMaintainNeighbours(allTris, activeTris, currTriIdx, oppNeighbour, quad, eps);
                        size_t new1 = allTris.size() - 2;
                        size_t new2 = allTris.size() - 1;
                        stack.push_back(new1);
                        stack.push_back(new2);
                    }
                }
            }
        }
    }
    for (const Edge& constraint : constraintList) {
        if (doesEdgeExist(allTris, activeTris, constraint, eps)) continue;
        std::vector<Edge> intersectingEdges = findIntersectingEdges(allTris, activeTris, constraint, eps);
        int maxFlips = static_cast<int>(intersectingEdges.size() * intersectingEdges.size()) + 20;
        while (!intersectingEdges.empty() && maxFlips-- > 0) {
            Edge toFlip = intersectingEdges.front();
            intersectingEdges.erase(intersectingEdges.begin());
            std::optional<std::pair<size_t, size_t>> adj = getTrianglesForEdge(allTris, activeTris, toFlip, eps);
            if (!adj) continue;
            std::array<Vec2, 4> quad;
            try { quad = getQuadrilateral(allTris, *adj, eps); } catch (...) { continue; }
            if (checkIfConvexQuadrilateral(quad[0], quad[1], quad[2], quad[3], eps)) {
                Edge newDiag = {quad[0], quad[2]};
                if (areSegmentsCrossing(newDiag, constraint, eps)) {
                    intersectingEdges.push_back(toFlip); 
                } else {
                    flipEdgeAndMaintainNeighbours(allTris, activeTris, adj->first, adj->second, quad, eps);
                }
            } else {
                intersectingEdges.push_back(toFlip); 
            }
        }
    }
    bool flipped = true;
    int max_opt_flips = static_cast<int>(activeTris.size() * activeTris.size()) + 20;
    while (flipped && max_opt_flips-- > 0) {
        flipped = false;
        std::vector<size_t> currentActive(activeTris.begin(), activeTris.end());
        for (size_t tIdx : currentActive) {
            if (activeTris.count(tIdx) == 0) continue;
            TriangleEdges tEdges = allTris[tIdx].getEdges();
            for (int i = 0; i < 3; ++i) {
                Edge e = tEdges[i];
                if (isConstraintEdge(e, constraintList, eps)) continue;
                std::optional<std::pair<size_t, size_t>> adj = getTrianglesForEdge(allTris, activeTris, e, eps);
                if (!adj) continue;
                std::array<Vec2, 4> quad;
                try { quad = getQuadrilateral(allTris, *adj, eps); } catch (...) { continue; }
                if (checkIfConvexQuadrilateral(quad[0], quad[1], quad[2], quad[3], eps)) {
                    TriangleCDT t1_test({quad[0], quad[1], quad[3]});
                    Vec2 d = quad[2];
                    if (isInCircumcircle(t1_test, d, eps)) {
                        flipEdgeAndMaintainNeighbours(allTris, activeTris, adj->first, adj->second, quad, eps);
                        flipped = true;
                        break; 
                    }
                }
            }
            if (flipped) break;
        }
    }
    std::vector<TriangleCDT> result;
    for (size_t tIdx : activeTris) {
        const TriangleCDT& t = allTris[tIdx];
        if (!t.hasVertex(p1, eps) && !t.hasVertex(p2, eps) && !t.hasVertex(p3, eps)) {
            result.push_back(t);
        }
    }
    return result;
}
std::vector<TriangleCDT> calculateDelaunay(std::vector<Vec2>& points, double eps)
{
    if (points.size() < 3) return {};
    std::vector<Vec2> steinerPoints;
    if (points.size() >= 4) {
        Vec2 centroid{0.0, 0.0};
        for (const auto& p : points) {
            centroid.x += p.x;
            centroid.y += p.y;
        }
        centroid.x /= static_cast<double>(points.size());
        centroid.y /= static_cast<double>(points.size());
        points.push_back(centroid);
    }
    std::vector<TriangleCDT> initialTris = calculateCDT(points, {}, eps);
    for (const auto& tri : initialTris) {
        double d1 = dot(tri.p1 - tri.p2, tri.p1 - tri.p2);
        double d2 = dot(tri.p2 - tri.p3, tri.p2 - tri.p3);
        double d3 = dot(tri.p3 - tri.p1, tri.p3 - tri.p1);
        double maxEdgeSq = std::max({d1, d2, d3});
        if (maxEdgeSq > 100.0 * eps * eps) {
            Vec2 triCentroid = (tri.p1 + tri.p2 + tri.p3) / 3.0;
            steinerPoints.push_back(triCentroid);
        }
    }
    points.insert(points.end(), steinerPoints.begin(), steinerPoints.end());
    return calculateCDT(points, {}, eps);
}