#include "math_utils.hpp"

size_t hashCombine(size_t seed, int64_t v) {
    return seed ^ (std::hash<int64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}
double computeMeshEpsilon(const BBox& boxA, const BBox& boxB) {
    BBox combined = boxA;
    combined.extend(boxB);
    Vec3 diag = combined.max - combined.min;
    int axis = getLargestAxis(diag);
    return std::max(diag[axis] * 1e-7, 1e-12);
}