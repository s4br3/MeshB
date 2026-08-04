#include "mesh_types.hpp"
MeshData filterMesh(const MeshData& mesh, const std::vector<bool>& remove);
void clearUnusedNodes(MeshData& mesh);
size_t findRoot(std::vector<size_t>& parent, size_t i);
void unionNodes(std::vector<size_t>& parent, size_t i, size_t j);
void collapseSlivers(MeshData& mesh, double eps);
void removeDegenerateTriangles(MeshData& mesh);
void cleanMesh(MeshData& mesh, double eps);
void invertWinding(MeshData& mesh);
MeshData combineMeshes(const std::vector<MeshData>& meshes, double eps);
