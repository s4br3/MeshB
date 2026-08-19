#include "mesh_io.hpp"
#include "triangulation.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <iostream>
std::string getExtension(const std::string& filename) {
    size_t pos = filename.find_last_of('.');
    if (pos == std::string::npos) return "";
    std::string ext = filename.substr(pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), 
                   [](unsigned char c){ return std::tolower(c); });
    return ext;
}
void addPolygonToMesh(const std::vector<size_t>& polyIndices,
                      MeshData& mesh,
                      SpatialGrid3D& nodeGrid,
                      double eps,
                      size_t tag) {
    if (polyIndices.size() < 3) return;
    if (polyIndices.size() == 3) {
        mesh.triangles.emplace_back(std::array<size_t, 3>{
            polyIndices[0], polyIndices[1], polyIndices[2]
        });
        mesh.tags.push_back(tag);
        return;
    }
    PolyLine outerSegments;
    outerSegments.reserve(polyIndices.size());
    Vec3 normal = {0, 0, 0};
    Vec3 center = {0, 0, 0};
    for (size_t i = 0; i < polyIndices.size(); ++i) {
        size_t idxA = polyIndices[i];
        size_t idxB = polyIndices[(i + 1) % polyIndices.size()];
        const Vec3& pA = mesh.nodes[idxA];
        const Vec3& pB = mesh.nodes[idxB];
        outerSegments.emplace_back(pA, pB);
        center = center + pA;
    }
    center = center * (1.0 / polyIndices.size());
    for (size_t i = 0; i < polyIndices.size(); ++i) {
        const Vec3& pCurr = mesh.nodes[polyIndices[i]];
        const Vec3& pNext = mesh.nodes[polyIndices[(i + 1) % polyIndices.size()]];
        normal.x += (pCurr.y - pNext.y) * (pCurr.z + pNext.z);
        normal.y += (pCurr.z - pNext.z) * (pCurr.x + pNext.x);
        normal.z += (pCurr.x - pNext.x) * (pCurr.y + pNext.y);
    }
    normal = normal.normalize();
    ProjectionFrame frame = computeSharedFrame(normal, center);
    std::vector<Triangle> tris;
    triangulate(outerSegments, {}, frame, nodeGrid, eps, tris);
    for (const auto& tri : tris) {
        mesh.triangles.push_back(tri);
        mesh.tags.push_back(tag);
    }
}
MeshData loadMesh(const std::string& filename) {
    std::string ext = getExtension(filename);
    if (ext == "msh") return loadMSH(filename);
    if (ext == "obj") return loadOBJ(filename);
    if (ext == "stl") return loadSTL(filename);
    throw std::invalid_argument("Unsupported file extension: " + ext);
}
void saveMesh(const MeshData& mesh, const std::string& filename) {
    std::string ext = getExtension(filename);
    if (ext == "msh") saveMSH(mesh, filename);
    else if (ext == "obj") saveOBJ(mesh, filename);
    else if (ext == "stl") saveSTL(mesh, filename);
    else throw std::invalid_argument("Unsupported file extension: " + ext);
}
MeshData loadOBJ(const std::string &filename) {
    PolygonSoup soup;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open OBJ file: " + filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if (prefix == "v") {
            double x, y, z;
            iss >> x >> y >> z;
            soup.nodes.push_back({x, y, z});
        } 
        else if (prefix == "f") {
            std::vector<size_t> vIndices;
            std::string vertexStr;
            while (iss >> vertexStr) {
                size_t slashPos = vertexStr.find('/');
                std::string vIdxStr = (slashPos != std::string::npos) ? 
                                      vertexStr.substr(0, slashPos) : vertexStr;
                long long idx = std::stoll(vIdxStr);
                if (idx > 0) {
                    vIndices.push_back(static_cast<size_t>(idx - 1));
                } else if (idx < 0) {
                    vIndices.push_back(static_cast<size_t>(soup.nodes.size() + idx));
                }
            }
            if (vIndices.size() >= 3) {
                soup.polygons.push_back(vIndices);
            }
        }
    }
    double eps = computeMeshEpsilon(soup.nodes);
    MeshData mesh;
    mesh.nodes = soup.nodes;
    SpatialGrid3D nodeGrid(eps);
    std::vector<size_t> indexMap(mesh.nodes.size());
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        indexMap[i] = nodeGrid.getOrAdd(mesh.nodes[i]);
    }
    for (const auto& poly : soup.polygons) {
        if (poly.size() < 3) continue;
        if (poly.size() == 3){
            mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[1]], indexMap[poly[2]]});
            continue;
        }
        if (poly.size() == 4) {
            double d1 = (mesh.nodes[poly[0]] - mesh.nodes[poly[2]]).length2();
            double d2 = (mesh.nodes[poly[1]] - mesh.nodes[poly[3]]).length2();
            if (d1 < d2) {
                mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[1]], indexMap[poly[2]]});
                mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[2]], indexMap[poly[3]]});
            } else {
                mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[1]], indexMap[poly[3]]});
                mesh.triangles.push_back({indexMap[poly[1]], indexMap[poly[2]], indexMap[poly[3]]});
            }
            continue;
        }
        addPolygonToMesh(poly, mesh, nodeGrid, eps, 0); 
    }
    mesh.nodes = nodeGrid.getUniquePoints();
    recomputeMeshData(mesh);
    return mesh;
}
void saveOBJ(const MeshData& mesh, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open file for writing: " + filename);
    for (const auto& node : mesh.nodes) {
        file << "v " << node.x << " " << node.y << " " << node.z << "\n";
    }
    for (const auto& tri : mesh.triangles) {
        // OBJ uses 1-based indexing
        file << "f " << (tri.v[0] + 1) << " " 
                     << (tri.v[1] + 1) << " " 
                     << (tri.v[2] + 1) << "\n";
    }
}
MeshData loadSTL(const std::string &filename) {
    MeshData mesh;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open STL file: " + filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        if (token == "vertex") {
            double x, y, z;
            iss >> x >> y >> z;
            mesh.nodes.push_back({x, y, z});
            if (mesh.nodes.size() % 3 == 0) {
                size_t idx = mesh.nodes.size();
                mesh.triangles.emplace_back(idx - 3, idx - 2, idx - 1);
                mesh.tags.push_back(0);
            }
        }
    }
    recomputeMeshData(mesh);
    return mesh;
}
void saveSTL(const MeshData& mesh, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open file for writing: " + filename);
    file << "solid mesh\n";
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        Vec3 n = mesh.normals.empty() ? mesh.triangles[i].normal(mesh.nodes) : mesh.normals[i];
        file << "  facet normal " << n.x << " " << n.y << " " << n.z << "\n";
        file << "    outer loop\n";
        for (int j = 0; j < 3; ++j) {
            Vec3 v = mesh.nodes[mesh.triangles[i].v[j]];
            file << "      vertex " << v.x << " " << v.y << " " << v.z << "\n";
        }
        file << "    endloop\n";
        file << "  endfacet\n";
    }
    file << "endsolid mesh\n";
}
MeshData loadMSH(const std::string &filename) {
    PolygonSoup soup;
    std::vector<int> polyTags;
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open MSH file: " + filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("$Nodes") != std::string::npos) {
            size_t numNodes;
            file >> numNodes;
            soup.nodes.resize(numNodes);
            for (size_t i = 0; i < numNodes; ++i) {
                size_t id;
                double x, y, z;
                file >> id >> x >> y >> z;
                soup.nodes[id - 1] = {x, y, z};
            }
        } 
        else if (line.find("$Elements") != std::string::npos) {
            size_t numElements;
            file >> numElements;
            for (size_t i = 0; i < numElements; ++i) {
                size_t id, type, numTags;
                file >> id >> type >> numTags;
                std::vector<int> tags(numTags);
                for (size_t j = 0; j < numTags; ++j) {
                    file >> tags[j];
                }
                if (type == 2 || type == 3) {
                    size_t numVerts = (type == 2) ? 3 : 4;
                    std::vector<size_t> poly(numVerts);
                    for (size_t k = 0; k < numVerts; ++k) {
                        size_t v;
                        file >> v;
                        poly[k] = v - 1;
                    }
                    soup.polygons.push_back(poly);
                    polyTags.push_back(tags.empty() ? 0 : tags[0]);
                } else {
                    std::string skip;
                    std::getline(file, skip); 
                }
            }
        }
    }
    double eps = computeMeshEpsilon(soup.nodes);
    MeshData mesh;
    mesh.nodes = soup.nodes;
    SpatialGrid3D nodeGrid(eps);
    std::vector<size_t> indexMap(mesh.nodes.size());
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        indexMap[i] = nodeGrid.getOrAdd(mesh.nodes[i]);
    }
    for (const auto& poly : soup.polygons) {
        if (poly.size() < 3) continue;
        if (poly.size() == 3){
            mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[1]], indexMap[poly[2]]});
            continue;
        }
        if (poly.size() == 4) {
            double d1 = (mesh.nodes[poly[0]] - mesh.nodes[poly[2]]).length2();
            double d2 = (mesh.nodes[poly[1]] - mesh.nodes[poly[3]]).length2();
            if (d1 < d2) {
                mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[1]], indexMap[poly[2]]});
                mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[2]], indexMap[poly[3]]});
            } else {
                mesh.triangles.push_back({indexMap[poly[0]], indexMap[poly[1]], indexMap[poly[3]]});
                mesh.triangles.push_back({indexMap[poly[1]], indexMap[poly[2]], indexMap[poly[3]]});
            }
            continue;
        }
        addPolygonToMesh(poly, mesh, nodeGrid, eps, 0); 
    }
    mesh.nodes = nodeGrid.getUniquePoints();
    recomputeMeshData(mesh);
    return mesh;
}
void saveMSH(const MeshData& mesh, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open file for writing: " + filename);
    file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
    file << "$Nodes\n" << mesh.nodes.size() << "\n";
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        file << (i + 1) << " " << mesh.nodes[i].x << " " 
             << mesh.nodes[i].y << " " << mesh.nodes[i].z << "\n";
    }
    file << "$EndNodes\n";
    file << "$Elements\n" << mesh.triangles.size() << "\n";
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        size_t tag = (i < mesh.tags.size()) ? mesh.tags[i] : 0;
        file << (i + 1) << " 2 2 " << tag << " " << tag << " " 
             << (mesh.triangles[i].v[0] + 1) << " " 
             << (mesh.triangles[i].v[1] + 1) << " " 
             << (mesh.triangles[i].v[2] + 1) << "\n";
    }
    file << "$EndElements\n";
}