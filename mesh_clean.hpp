#pragma once
#include "mesh_types.hpp"

/**
* @brief Filters a mesh by removing elements specified by a selection mask.
* @param[in] mesh - Source mesh data to filter.
* @param[in] remove - Boolean flag array where true indicates the triangle at that index should be removed.
* @return Filtered MeshData object containing remaining valid triangles and nodes.
*/
MeshData filterMesh(const MeshData& mesh, const std::vector<bool>& remove);

/**
* @brief Removes vertices from the mesh that are no longer referenced by any triangle elements.
* @param[in,out] mesh - Target MeshData object to prune in-place.
*/
void clearUnusedNodes(MeshData& mesh);

/**
* @brief Finds the representative root index of a node set in a Disjoint-Set Union (DSU) structure.
* @param[in,out] parent - Vector tracking node parent relationships.
* @param[in] i - Node index to query.
* @return Root node index representing the set.
*/
size_t findRoot(std::vector<size_t>& parent, size_t i);

/**
* @brief Unifies two subsets containing the specified node indices in a Disjoint-Set Union (DSU) structure.
* @param[in,out] parent - Vector tracking node parent relationships.
* @param[in] i - First node index.
* @param[in] j - Second node index.
*/
void unionNodes(std::vector<size_t>& parent, size_t i, size_t j);

/**
* @brief Collapses short edges and removes thin sliver triangles within a spatial tolerance distance.
* @param[in,out] mesh - Mesh object to modify in-place.
* @param[in] eps - Spatial distance threshold below which short edges are collapsed.
*/
void collapseSlivers(MeshData& mesh, double eps);

/**
* @brief Identifies and removes zero-area or duplicate-vertex degenerate triangles from the mesh.
* @param[in,out] mesh - Target MeshData object to clean in-place.
*/
void removeDegenerateTriangles(MeshData& mesh);

/**
* @brief Performs full topological cleanup on a mesh (collapses slivers, removes degeneracies, and purges unreferenced nodes).
* @param[in,out] mesh - Target MeshData object to process in-place.
* @param[in] eps - Geometric tolerance threshold for cleaning operations.
*/
void cleanMesh(MeshData& mesh, double eps);

/**
* @brief Reverses the vertex ordering/winding of all triangles to flip face normals.
* @param[in,out] mesh - Target MeshData object whose normal orientations will be inverted.
*/
void invertWinding(MeshData& mesh);

/**
* @brief Combines a list of separate mesh structures into a single unified mesh and welds nearby vertices.
* @param[in] meshes - Array of input MeshData instances to merge.
* @param[in] eps - Spatial tolerance distance threshold for welding co-located nodes.
* @return Unified combined MeshData structure.
*/
MeshData combineMeshes(const std::vector<MeshData>& meshes, double eps);