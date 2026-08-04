#pragma once
#include "mesh_types.hpp"
#include <source/geometry/mesh/triangle_mesh.hpp>
MeshData extractMeshData(const bem::TriangleMesh<3>& mesh);
void rebuildMesh(bem::TriangleMesh<3>& mesh, const MeshData& newMesh, double eps);