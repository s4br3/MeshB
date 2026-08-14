#pragma once
#include "vector.hpp"
#include "bvh.hpp"
#include <cstdint>
#include <unordered_map>
#include <CDT.h>

/**
* @brief Formatted stream output operator for a 3D vector.
* @param[in,out] os - Output stream.
* @param[in] v - Vec3 to print.
* @return Reference to the output stream.
*/
/**
* @brief Combines a hash seed with a 64-bit integer value.
* @param[in] seed - Existing hash value.
* @param[in] v - Value to combine into hash.
* @return Combined hash value.
*/
size_t hashCombine(size_t seed, int64_t v);

/**
* @brief Computes adaptive dynamic geometric tolerance (epsilon) based on mesh bounding boxes.
* @param[in] boxA - Bounding box of mesh A.
* @param[in] boxB - Bounding box of mesh B.
* @return Computed absolute numerical tolerance.
*/
double computeMeshesEpsilon(const BBox& boxA, const BBox& boxB);

/**
* @brief Computes geometric tolerance for this single mesh (relevant to forced
* triangulation on loading meshes)
* @param[in] nodes - List of vector coordinates of vertices in the mesh
*/
double computeMeshEpsilon(const std::vector<Vec3>& nodes);
/**
* @brief Flattens a 2D nested vector into a 1D vector.
* @tparam T Type of elements contained in the vector.
* @param[in] v - Nested 2D vector.
* @return Flattened 1D vector containing all elements sequentially.
*/
template <class T>
std::vector<T> flattenVector(const std::vector<std::vector<T>>& v) {
    std::vector<T> out;
    size_t total = 0;
    for (const auto& row : v) total += row.size();
    out.reserve(total);
    for (const auto& row : v)
        for (const auto& x : row)
            out.push_back(x);
    return out;
}

/**
* @brief Evaluates the sign of a numerical value with geometric tolerance.
* @tparam T Numeric scalar type.
* @param[in] num - Value to check.
* @param[in] eps - Absolute tolerance threshold.
* @return -1 if negative beyond tolerance, +1 if positive beyond tolerance, 0 if within tolerance.
*/
template <class T>
int sign(const T& num, double eps) {
    if (num < -eps) return -1;
    if (num >  eps) return  1;
    return 0;
}

/**
* @brief Identifies the unique vertex index that lies on one side of a plane compared to the other two.
* @tparam T Signed distance scalar type.
* @param[in] dist - Array of signed distances for 3 triangle vertices.
* @param[in] eps - Geometric tolerance threshold.
* @return Index (0, 1, or 2) of isolated vertex, or -1 if no unique split exists.
*/
template <class T>
int uniqueSignIndex(const std::array<T,3>& dist, double eps) {
    int s0 = sign(dist[0], eps);
    int s1 = sign(dist[1], eps);
    int s2 = sign(dist[2], eps);
    if (s0 != s1 && s1 == s2) return 0;
    if (s1 != s0 && s0 == s2) return 1;
    if (s2 != s0 && s0 == s1) return 2;
    if (s0 != 0 && s1 != 0 && s0 != s1) return (s2 == 0) ? 2 : -1;
    if (s1 != 0 && s2 != 0 && s1 != s2) return (s0 == 0) ? 0 : -1;
    if (s2 != 0 && s0 != 0 && s2 != s0) return (s1 == 0) ? 1 : -1;
    return -1;
}

/**
* @class SpatialGrid3D
* @brief Spatial hashing grid for 3D point deduplication and fast proximity lookup.
* @details Divides 3D space into uniform cubic cells of size `eps` to query and merge near-duplicate coordinates.
*/
class SpatialGrid3D {
private:
    struct Key {
        int64_t x, y, z;
        bool operator==(const Key& o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<int64_t>{}(k.x);
            h = hashCombine(h, k.y);
            h = hashCombine(h, k.z);
            return h;
        }
    };
    double cellSize;
    double sqrEps;
    std::unordered_map<Key, std::vector<size_t>, KeyHash> grid;
    std::vector<Vec3> points;

public:
    /**
    * @brief Constructs a 3D spatial grid given a tolerance distance.
    * @param[in] eps - Cell size and maximum distance for duplicate point matching.
    */
    explicit SpatialGrid3D(double eps) : cellSize(eps), sqrEps(eps * eps) {}

    /**
    * @brief Queries an existing point within distance tolerance or adds a new point to the grid.
    * @param[in] p - 3D Point coordinate.
    * @return Index of existing or newly inserted point in the internal point list.
    * @details To avoid missing grid-misaligned degeneracy, this checks the current grid cube and its 26 neighbours
    */
    size_t getOrAdd(const Vec3& p) {
        int64_t cx = static_cast<int64_t>(std::floor(p[0] / cellSize));
        int64_t cy = static_cast<int64_t>(std::floor(p[1] / cellSize));
        int64_t cz = static_cast<int64_t>(std::floor(p[2] / cellSize));
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    Key neighborKey{cx + dx, cy + dy, cz + dz};
                    auto it = grid.find(neighborKey);
                    if (it != grid.end()) {
                        for (size_t idx : it->second) {
                            Vec3 diff = points[idx] - p;
                            double sqrDist = dot(diff, diff);
                            if (sqrDist < sqrEps) {
                                return idx;
                            }
                        }
                    }
                }
            }
        }
        size_t newIdx = points.size();
        points.push_back(p);
        grid[{cx, cy, cz}].push_back(newIdx);
        return newIdx;
    }

    /**
    * @brief Gets the list of unique points stored in the grid.
    * @return Reference to vector of unique 3D point positions.
    */
    const std::vector<Vec3>& getUniquePoints() const { return points; }

    /**
    * @brief Clears all grid cells and stored points.
    */
    void clear() { grid.clear(); points.clear(); }
};

/**
* @class SpatialGrid2D
* @brief Spatial hashing grid for 2D point deduplication on local coordinate projection planes.
*/
class SpatialGrid2D {
private:
    struct Key {
        int64_t x, y;
        bool operator==(const Key& o) const {
            return x == o.x && y == o.y;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<int64_t>{}(k.x);
            return hashCombine(h, k.y);
        }
    };
    double cellSize;
    double sqrEps;
    std::unordered_map<Key, std::vector<CDT::VertInd>, KeyHash> grid;
    std::vector<CDT::V2d<double>> points;

public:
    /**
    * @brief Constructs a 2D spatial grid.
    * @param[in] eps - Cell size and deduplication threshold.
    */
    explicit SpatialGrid2D(double eps) : cellSize(eps), sqrEps(eps * eps) {}

    /**
    * @brief Queries or inserts a 2D point into the grid.
    * @param[in] p - 2D Point coordinate.
    * @return Index of the existing or newly registered point.
    * @details To avoid missing grid-misaligned degeneracy, this checks the current grid square and its 8 neighbours
    */
    CDT::VertInd getOrAdd(const CDT::V2d<double>& p) {
        int64_t cx = static_cast<int64_t>(std::floor(p.x / cellSize));
        int64_t cy = static_cast<int64_t>(std::floor(p.y / cellSize));
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                Key neighborKey{cx + dx, cy + dy};
                auto it = grid.find(neighborKey);
                if (it != grid.end()) {
                    for (CDT::VertInd idx : it->second) {
                        double dx_p = points[idx].x - p.x;
                        double dy_p = points[idx].y - p.y;
                        if ((dx_p * dx_p + dy_p * dy_p) < sqrEps) {
                            return idx;
                        }
                    }
                }
            }
        }
        CDT::VertInd newIdx = static_cast<CDT::VertInd>(points.size());
        points.push_back(p);
        grid[{cx, cy}].push_back(newIdx);
        return newIdx;
    }

    /**
    * @brief Returns unique 2D point positions stored in the spatial grid.
    * @return Reference to array of 2D points.
    */
    const std::vector<CDT::V2d<double>>& getUniquePoints() const { return points; }

    /**
    * @brief Resets and clears the 2D grid structure.
    */
    void clear() { grid.clear(); points.clear(); }
};