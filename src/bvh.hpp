#pragma once
#include "vector.hpp"
#include <numeric>
#include <limits>
#include <vector>
#include <algorithm>

/**
* @brief Positive infinity floating-point constant.
*/
constexpr double inf = std::numeric_limits<double>::infinity();

/**
* @struct BBox
* @brief Represents a 3D Axis-Aligned Bounding Box (AABB)
*/
struct BBox {
    Vec3 min = Vec3{inf,  inf, inf};
    Vec3 max = Vec3{-inf, -inf, -inf};
    BBox() = default;

    /**
    * @brief Constructs a bounding box initialized with a list of 3D points.
    * @param[in] xs - Initializer list of 3D point coordinates.
    */
    BBox(std::initializer_list<Vec3> xs) {
        for (auto it = xs.begin(); it != xs.end(); ++it) {
            extend(*it);
        }
    }

    /**
    * @brief Extends the bounding box boundaries to enclose a 3D point.
    * @param[in] p - 3D point coordinate.
    */
    void extend(const Vec3& p) {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
        max.z = std::max(max.z, p.z);
    }

    /**
    * @brief Extends the bounding box boundaries to enclose another bounding box.
    * @param[in] b - Bounding box to merge.
    */
    void extend(const BBox& b) {
        min.x = std::min(min.x, b.min.x);
        min.y = std::min(min.y, b.min.y);
        min.z = std::min(min.z, b.min.z);
        max.x = std::max(max.x, b.max.x);
        max.y = std::max(max.y, b.max.y);
        max.z = std::max(max.z, b.max.z);
    }

    /**
    * @brief Calculates the total exterior surface area of the bounding box.
    * @return Surface area.
    */
    double surfaceArea() const {
        double dx = max.x - min.x;
        double dy = max.y - min.y;
        double dz = max.z - min.z;
        return 2.0* (dx* dy + dy* dz + dz* dx); 
    }
};

/**
* @struct Node
* @brief Node element contained within a Bounding Volume Hierarchy (BVH) spatial tree.
*/
struct Node {
    BBox boundingBox;
    size_t firstId = 0; 
    size_t primCount = 0;
    BBox getBbox() const { return boundingBox; }
    bool isLeaf() const { return primCount > 0; }
};

/**
* @struct Bin
* @brief Spatial partition bin used during Surface Area Heuristic (SAH) BVH splits.
*/
struct Bin {
    BBox box;
    size_t count = 0;
};

/**
* @struct BVH
* @brief Bounding Volume Hierarchy tree for spatial ray tracing and collision acceleration.
*/
struct BVH {
    const int binCount = 16;
    std::vector<size_t> primIds;
    std::vector<Node> nodes;

    /**
    * @struct BuildTask
    * @brief Traversal task unit used for iterative BVH tree construction.
    */
    struct BuildTask {
        size_t nodeIdx;
        size_t start;
        size_t count;
    };

    /**
    * @brief Builds a BVH spatial tree over primitive bounding boxes using Surface Area Heuristic (SAH).
    * @param[in] bboxes - Primitive bounding boxes.
    * @param[in] centres - Primitive centroids.
    * @param[in] eps - Distance tolerance threshold.
    */
    BVH(const std::vector<BBox>& bboxes, const std::vector<Vec3>& centres, double eps) {
        if (bboxes.empty()) return;
        primIds.resize(bboxes.size());
        std::iota(primIds.begin(), primIds.end(), 0);
        nodes.reserve(bboxes.size()* 2 - 1);
        nodes.emplace_back();
        std::vector<BuildTask> stack;
        stack.reserve(64);
        stack.push_back({0, 0, bboxes.size()});
        while (!stack.empty()) {
            BuildTask task = stack.back();
            stack.pop_back();
            for (size_t i = task.start; i < task.start + task.count; i++) {
                nodes[task.nodeIdx].boundingBox.extend(bboxes[primIds[i]]);
            }
            if (task.count <= 2) {
                nodes[task.nodeIdx].firstId = task.start;
                nodes[task.nodeIdx].primCount = task.count;
                continue;
            }
            BBox centBounds;
            for (size_t i = task.start; i < task.start + task.count; ++i) {
                centBounds.extend(centres[primIds[i]]);
            }
            Vec3 diag = centBounds.max - centBounds.min;
            int axis = getLargestAxis(diag);
            size_t leftCount = 0;
            if (diag[axis] <= eps) {
                leftCount = task.count / 2;
            } else {
                std::vector<Bin> bins(binCount);
                for (size_t i = task.start; i < task.start + task.count; ++i) {
                    size_t pId = primIds[i];
                    double pos = centres[pId][axis];
                    double ratio = (pos - centBounds.min[axis]) / diag[axis];
                    int binIdx = static_cast<int>(binCount* ratio);
                    binIdx = std::clamp(binIdx, 0, binCount - 1);
                    bins[binIdx].count++;
                    bins[binIdx].box.extend(bboxes[pId]);
                }
                double minCost = inf;
                int bestSplitBin = -1;
                for (int split = 0; split < binCount - 1; ++split) {
                    BBox leftBox, rightBox;
                    size_t lCount = 0, rCount = 0;
                    for (int i = 0; i <= split; ++i) {
                        if (bins[i].count > 0) {
                            leftBox.extend(bins[i].box);
                            lCount += bins[i].count;
                        }
                    }
                    for (int i = split + 1; i < binCount; ++i) {
                        if (bins[i].count > 0) {
                            rightBox.extend(bins[i].box);
                            rCount += bins[i].count;
                        }
                    }
                    double cost = leftBox.surfaceArea()* lCount + rightBox.surfaceArea()* rCount;
                    if (cost < minCost) {
                        minCost = cost;
                        bestSplitBin = split;
                    }
                }
                auto splitIt = std::partition(
                    primIds.begin() + task.start, 
                    primIds.begin() + task.start + task.count,
                    [&](size_t pId) {
                        double pos = centres[pId][axis];
                        double ratio = (pos - centBounds.min[axis]) / diag[axis];
                        int binIdx = static_cast<int>(binCount* ratio);
                        binIdx = std::clamp(binIdx, 0, binCount - 1);
                        return binIdx <= bestSplitBin;
                    }
                );
                leftCount = std::distance(primIds.begin() + task.start, splitIt);
                if (leftCount == 0 || leftCount == task.count) {
                    leftCount = task.count / 2;
                    std::nth_element(
                            primIds.begin() + task.start,
                            primIds.begin() + task.start + leftCount,
                            primIds.begin() + task.start + task.count,
                            [&](size_t a, size_t b) {
                                return centres[a][axis] < centres[b][axis];
                            }
                        );
                }
            }
            size_t leftChildIdx = nodes.size();
            nodes.emplace_back();
            nodes.emplace_back();
            nodes[task.nodeIdx].firstId = leftChildIdx;
            nodes[task.nodeIdx].primCount = 0;
            stack.push_back({leftChildIdx + 1, task.start + leftCount, task.count - leftCount});
            stack.push_back({leftChildIdx, task.start, leftCount});
        }
    }
};