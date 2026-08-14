#pragma once
#include <string>
#include <source/geometry/mesh/triangle_mesh.hpp>

/**
 * @brief Reads a .msh file directly into an openBEM TriangleMesh<3>.
 * @param[in] filename Path to the .msh file.
 * @return Fully initialized bem::TriangleMesh<3>.
 */
bem::TriangleMesh<3> loadMSH(const std::string& filename);

/**
 * @brief Saves an openBEM TriangleMesh<3> directly to a .msh file.
 * @param[in] mesh Source openBEM TriangleMesh<3>.
 * @param[in] filename Output .msh file path.
 */
void saveMSH(const bem::TriangleMesh<3>& mesh, const std::string& filename, const size_t tag = -1);