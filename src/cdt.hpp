#pragma once
#include "vector.hpp"
#include "geom_2d.hpp"
#include <unordered_set>
#include <optional>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
using TriangleEdges = std::array<Edge, 3>;
constexpr size_t NO_NEIGHBOUR = std::numeric_limits<size_t>::max();
struct TriangleCDT {
    size_t n1 = NO_NEIGHBOUR;
    size_t n2 = NO_NEIGHBOUR;
    size_t n3 = NO_NEIGHBOUR;
    Vec2 p1;
    Vec2 p2;
    Vec2 p3;
    TriangleCDT(std::array<Vec2, 3> vertices) {
        p1 = vertices[0];
        p2 = vertices[1];
        p3 = vertices[2];
    }
    bool hasVertex(const Vec2& p, double eps) const {
        return p.equals(p1, eps) || p.equals(p2, eps) || p.equals(p3, eps);
    }
    const std::array<Vec2, 3> getVertices() const {
        return {p1, p2, p3};
    }
    const TriangleEdges getEdges() const {
        return {Edge{p1, p2}, Edge{p2, p3}, Edge{p3, p1}};
    }
    size_t getNeighbour(const Vec2& opposite, double eps) const {
        if (opposite.equals(p1, eps)) return n1;
        if (opposite.equals(p2, eps)) return n2;
        if (opposite.equals(p3, eps)) return n3;
        return NO_NEIGHBOUR;
    }
    void replaceNeighbour(size_t oldNeighbour, size_t newNeighbour) {
        if (n1 == oldNeighbour) n1 = newNeighbour;
        else if (n2 == oldNeighbour) n2 = newNeighbour;
        else if (n3 == oldNeighbour) n3 = newNeighbour;
    }
};

/**
* @brief Performs Constrained Delaunay Triangulation on a set of 2D points and constraint edges.
* @param[in] points Vector of unique 2D points.
* @param[in] edges Vector of constraint edges that must appear in the final triangulation.
* @param[in] eps Geometric tolerance.
* @return A vector of triangulated faces.
*/
std::vector<TriangleCDT> calculateCDT(const std::vector<Vec2>& points, const std::vector<EdgeKey>& edges, double eps);

/**
* @brief Checks if a 2D point lies within the boundaries of a triangle.
* @param[in] p The 2D point to test.
* @param[in] tri The triangle to test against.
* @param[in] eps Tolerance value.
* @return True if the point is inside or on the edge of the triangle.
*/
bool pointInTriangle(const Vec2& p, const TriangleCDT& tri, double eps);

/**
* @brief Compares two 2D points for sorting purposes.
* @param[in] a First point.
* @param[in] b Second point.
* @param[in] eps Tolerance value.
* @return Negative if a < b, positive if a > b, zero if equivalent.
*/
double comparePoints(const Vec2& a, const Vec2& b, double eps);

/**
* @brief Determines if two line segments cross each other.
* @param[in] e1 First edge segment.
* @param[in] e2 Second edge segment.
* @param[in] eps Tolerance value.
* @return True if the segments cross strictly within their endpoints.
*/
bool areSegmentsCrossing(const Edge& e1, const Edge& e2, double eps);

/**
* @brief Checks if an edge exists within a set of triangles.
* @param[in] allTris Vector of all triangles.
* @param[in] activeTris Set of indices of active triangles.
* @param[in] e The edge to search for.
* @param[in] eps Tolerance value.
* @return True if the edge is found, false otherwise.
*/
bool doesEdgeExist(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Edge& e, double eps);

/**
* @brief Finds all edges in the triangulation that intersect a given segment.
* @param[in] allTris Vector of all triangles.
* @param[in] activeTris Set of indices of active triangles.
* @param[in] e The segment to test for intersections.
* @param[in] eps Tolerance value.
* @return Vector of existing edges that cross the input segment.
*/
std::vector<Edge> findIntersectingEdges(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Edge& e, double eps);

/**
* @brief Finds the specific triangle index that contains a given point.
* @param[in] allTris Vector of all triangles.
* @param[in] activeTris Set of indices of active triangles.
* @param[in] p The 2D query point.
* @param[in] eps Tolerance value.
* @return The index of the containing triangle, or NO_NEIGHBOUR if none found.
*/
size_t findTriangleWithPoint(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Vec2& p, double eps);

/**
* @brief Retrieves the two adjacent triangle indices that share a specific edge.
* @param[in] allTris Vector of all triangles.
* @param[in] activeTris Set of indices of active triangles.
* @param[in] e The shared edge.
* @param[in] eps Tolerance value.
* @return Optional pair of triangle indices sharing the edge.
*/
std::optional<std::pair<size_t, size_t>> getTrianglesForEdge(const std::vector<TriangleCDT>& allTris, const std::unordered_set<size_t>& activeTris, const Edge& e, double eps);

/**
* @brief Checks if a point is strictly inside the circumcircle of a triangle.
* @param[in] t The triangle defining the circumcircle.
* @param[in] d The 2D query point.
* @return True if the point is inside the circumcircle.
*/
bool isInCircumcircle(const TriangleCDT& t, const Vec2& d);

/**
* @brief Validates if a quadrilateral formed by 4 points is strictly convex.
* @param[in] a First vertex.
* @param[in] b Second vertex.
* @param[in] c Third vertex.
* @param[in] d Fourth vertex.
* @return True if strictly convex.
*/
bool checkIfConvexQuadrilateral(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d);

/**
* @brief Extracts the 4 external vertices of a quadrilateral formed by two adjacent triangles.
* @param[in] allTris Vector of all triangles.
* @param[in] e Pair of adjacent triangle indices.
* @param[in] eps Tolerance value.
* @return Array of 4 points forming the outer boundary.
*/
std::array<Vec2, 4> getQuadrilateral(const std::vector<TriangleCDT>& allTris, const std::pair<size_t, size_t>& e, double eps);

/**
* @brief Computes a robust integer-based spatial hash for an edge.
* @param[in] e The edge to hash.
* @param[in] eps The cell size used for spatial discretization.
* @return Combined spatial hash value.
*/
size_t hashEdgeSpatial(const Edge& e, double eps);