#include "math_utils.hpp"

size_t hashCombine(size_t seed, int64_t v) {
    return seed ^ (std::hash<int64_t>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}
double computeMeshEpsilon(const BBox& boxA, const BBox& boxB) {
    BBox combined = boxA;
    combined.extend(boxB);
    Vec3 diag = combined.max - combined.min;
    double max_dim = std::max({diag[0], diag[1], diag[2]});
    return std::max(max_dim * 1e-7, 1e-12);
}