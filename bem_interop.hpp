#pragma once
#include "mesh_types.hpp"
#include <source/geometry/mesh/triangle_mesh.hpp>

/**
* @brief Converts an openBEM `bem::TriangleMesh<3>` structure into an internal `MeshData` representation.
* @param[in] mesh Source openBEM mesh object.
* @return Converted internal `MeshData` struct containing nodes, triangles, centres, normals, and tags.
*/
MeshData extractMeshData(const bem::TriangleMesh<3>& mesh);

/**
* @brief Updates and reconstructs an openBEM `bem::TriangleMesh<3>` instance from modified internal `MeshData`.
* @param[in,out] mesh - Destination openBEM triangle mesh object to update.
* @param[in] newMesh - Source mesh data struct.
* @param[in] eps - Distance tolerance threshold for vertex unification and cleanup.
*/
void rebuildMesh(bem::TriangleMesh<3>& mesh, const MeshData& newMesh, double eps);