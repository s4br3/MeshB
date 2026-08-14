#include "math_utils.hpp"
#include "bvh.hpp"

size_t hashCombine(size_t seed, int64_t v) {
    return seed ^ (std::hash<int64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}
double computeMeshesEpsilon(const BBox& boxA, const BBox& boxB) {
    BBox combined = boxA;
    combined.extend(boxB);
    Vec3 diag = combined.max - combined.min;
    int axis = getLargestAxis(diag);
    return std::max(diag[axis] * 1e-7, 1e-12);
}
double computeMeshEpsilon(const std::vector<Vec3>& nodes) {
    if (nodes.empty()) return 1e-6;
    BBox b;
    for (const auto& p : nodes) {
        b.extend(p);
    }
    double dx = b.max.x - b.min.x;
    double dy = b.max.y - b.min.y;
    double dz = b.max.z - b.min.z;
    double bboxDiag = std::sqrt(dx * dx + dy * dy + dz * dz);
    return std::max(bboxDiag * 1e-6, 1e-9);
}