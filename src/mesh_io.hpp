#pragma once
#include "mesh_types.hpp"
#include "math_utils.hpp"
#include <string>

/**
* @brief Find extension of filename.
* @param[in] filename - The file to find the extension of.
* @return File format extension.
*/
std::string getExtension(const std::string& filename);

/**
* @brief Force triangulation of this polygon and add it to the mesh.
* @param[in] polyIndices - Vector of indices pointing to the current polygon's vertices.
* @param[in, out] mesh - Mesh to add the triangulated polygon to.
* @param [in] nodeGrid - Spatial node matching structure.
* @param[in] eps - Tolerance value.
* @param[in] tag - Triangle tag (primarily relevant to msh files).
*/
void addPolygonToMesh(const std::vector<size_t>& polyIndices, MeshData& mesh,
                      SpatialGrid3D& nodeGrid, double eps, size_t tag = 0);

/**
* @brief General loader dispatcher for mesh files; selects parser based on extension (.obj, .stl, .msh).
* @param[in] filename - Path to input mesh file.
* @return MeshData structure containing nodes, triangles, and metadata.
* @throws std::invalid_argument if file extension is unsupported.
* @throws std::runtime_error if file fails to open.
*/
MeshData loadMesh(const std::string& filename);

/**
* @brief General saver dispatcher for mesh files; selects writer based on extension (.obj, .stl, .msh).
* @param[in] mesh - Mesh structure to write.
* @param[in] filename - Target output file path.
* @throws std::invalid_argument if file extension is unsupported.
* @throws std::runtime_error if file fails to open for writing.
*/
void saveMesh(const MeshData& mesh, const std::string& filename);

/**
* @brief Parses Wavefront OBJ format 3D mesh files.
* @param[in] filename - Path to OBJ file.
* @return Loaded MeshData object.
* @throws std::runtime_error if file opening fails.
*/
MeshData loadOBJ(const std::string &filename);

/**
* @brief Writes mesh data to Wavefront OBJ file format (1-based indexing).
* @param[in] mesh - Mesh data to save.
* @param[in] filename - Target output file path.
* @throws std::runtime_error if file opening fails.
*/
void saveOBJ(const MeshData& mesh, const std::string& filename);

/**
* @brief Parses ASCII STL format 3D surface mesh files.
* @param[in] filename - Path to STL file.
* @return Loaded MeshData object.
* @throws std::runtime_error if file opening fails.
*/
MeshData loadSTL(const std::string &filename);

/**
* @brief Writes mesh data to ASCII STL file format with facet normals.
* @param[in] mesh - Mesh data to save.
* @param[in] filename - Target output file path.
* @throws std::runtime_error if file opening fails.
*/
void saveSTL(const MeshData& mesh, const std::string& filename);

/**
* @brief Parses Gmsh MSH (v2.2 format) mesh files.
* @param[in] filename - Path to MSH file.
* @return Loaded MeshData object.
* @throws std::runtime_error if file opening fails.
*/
MeshData loadMSH(const std::string &filename);

/**
* @brief Writes mesh data to Gmsh MSH (v2.2 format) file.
* @param[in] mesh - Mesh data to save.
* @param[in] filename - Target output file path.
* @throws std::runtime_error if file opening fails.
*/
void saveMSH(const MeshData& mesh, const std::string& filename);