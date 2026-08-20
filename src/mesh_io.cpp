#include "mesh_io.hpp"
#include "triangulation.hpp"
#include "geom_2d.hpp"
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
                   [](unsigned char c){return std::tolower(c); });
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
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open MSH file: " + filename);
    MeshData mesh;
    std::vector<size_t> tagToNodeIndex; 
    std::string line;
    bool gridInitialized = false;
    double eps = 0.0;
    SpatialGrid3D nodeGrid(1e-5);
    while (std::getline(file, line)) {
        if (line.find("$Nodes") != std::string::npos) {
            size_t numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
            file >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
            tagToNodeIndex.resize(maxNodeTag + 1);
            for (size_t b = 0; b < numEntityBlocks; ++b) {
                int entityDim, entityTag, parametric;
                size_t numNodesInBlock;
                file >> entityDim >> entityTag >> parametric >> numNodesInBlock;
                std::vector<size_t> nodeTags(numNodesInBlock);
                for (size_t i = 0; i < numNodesInBlock; ++i) {
                    file >> nodeTags[i];
                }
                for (size_t i = 0; i < numNodesInBlock; ++i) {
                    double x, y, z;
                    file >> x >> y >> z;
                    if (parametric) {
                        double u, v;
                        if (entityDim >= 1) file >> u;
                        if (entityDim >= 2) file >> v;
                    }
                    tagToNodeIndex[nodeTags[i]] = mesh.nodes.size();
                    mesh.nodes.push_back({x, y, z});
                }
            }
        }
        else if (line.find("$Elements") != std::string::npos) {
            if (!gridInitialized && !mesh.nodes.empty()) {
                eps = computeMeshEpsilon(mesh.nodes);
                nodeGrid = SpatialGrid3D(eps);
                for (const auto& node : mesh.nodes) {
                    nodeGrid.getOrAdd(node);
                }
                gridInitialized = true;
            }
            size_t numEntityBlocks, numElements, minElemTag, maxElemTag;
            file >> numEntityBlocks >> numElements >> minElemTag >> maxElemTag;
            for (size_t b = 0; b < numEntityBlocks; ++b) {
                int entityDim, entityTag, elementType;
                size_t numElementsInBlock;
                file >> entityDim >> entityTag >> elementType >> numElementsInBlock;
                for (size_t i = 0; i < numElementsInBlock; ++i) {
                    size_t elemTag;
                    file >> elemTag;
                    if (elementType == 2) {
                        size_t n1, n2, n3;
                        file >> n1 >> n2 >> n3;
                        mesh.triangles.push_back({tagToNodeIndex[n1], tagToNodeIndex[n2], tagToNodeIndex[n3]});
                        mesh.tags.push_back(entityTag);
                    } 
                    else if (elementType == 3) {
                        size_t n1, n2, n3, n4;
                        file >> n1 >> n2 >> n3 >> n4;
                        size_t idx0 = tagToNodeIndex[n1];
                        size_t idx1 = tagToNodeIndex[n2];
                        size_t idx2 = tagToNodeIndex[n3];
                        size_t idx3 = tagToNodeIndex[n4];
                        double d1 = (mesh.nodes[idx0] - mesh.nodes[idx2]).length2();
                        double d2 = (mesh.nodes[idx1] - mesh.nodes[idx3]).length2();
                        if (d1 < d2) {
                            mesh.triangles.push_back({idx0, idx1, idx2});
                            mesh.triangles.push_back({idx0, idx2, idx3});
                        } else {
                            mesh.triangles.push_back({idx0, idx1, idx3});
                            mesh.triangles.push_back({idx1, idx2, idx3});
                        }
                        mesh.tags.push_back(entityTag);
                        mesh.tags.push_back(entityTag);
                    } 
                    else {
                        std::string lineRest;
                        std::getline(file, lineRest);
                        std::istringstream iss(lineRest);
                        std::vector<size_t> polyIndices;
                        size_t nodeTag;
                        while (iss >> nodeTag) {
                            if (nodeTag < tagToNodeIndex.size()) {
                                polyIndices.push_back(tagToNodeIndex[nodeTag]);
                            }
                        }
                        if (polyIndices.size() >= 3) {
                            addPolygonToMesh(polyIndices, mesh, nodeGrid, eps, entityTag);
                        }
                    }
                }
            }
        }
    }
    recomputeMeshData(mesh);
    return mesh;
}
void saveMSH(const MeshData& mesh, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open file for writing: " + filename);
    int defaultTag = 1;
    if (filename.find('B') != std::string::npos) {
        defaultTag = 2;
    }
    size_t maxTag = 0;
    for (size_t tag : mesh.tags) {
        if (tag > maxTag) maxTag = tag;
    }
    maxTag = std::max(maxTag, static_cast<size_t>(defaultTag));
    std::vector<std::vector<size_t>> tagToTriIndices(maxTag + 1);
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        size_t tag = (i < mesh.tags.size() && mesh.tags[i] != 0) ? mesh.tags[i] : defaultTag;
        tagToTriIndices[tag].push_back(i);
    }
    size_t activeEntities = 0;
    for (const auto& indices : tagToTriIndices) {
        if (!indices.empty()) activeEntities++;
    }
    double minX = std::numeric_limits<double>::max(), minY = minX, minZ = minX;
    double maxX = std::numeric_limits<double>::lowest(), maxY = maxX, maxZ = maxX;
    for (const auto& node : mesh.nodes) {
        minX = std::min(minX, node.x); minY = std::min(minY, node.y); minZ = std::min(minZ, node.z);
        maxX = std::max(maxX, node.x); maxY = std::max(maxY, node.y); maxZ = std::max(maxZ, node.z);
    }
    if (mesh.nodes.empty()) {minX = minY = minZ = maxX = maxY = maxZ = 0.0; }
    file << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n";
    file << "$PhysicalNames\n" << activeEntities << "\n";
    for (size_t tag = 0; tag <= maxTag; ++tag) {
        if (!tagToTriIndices[tag].empty()) {
            file << "2 " << tag << " \"Surface_" << tag << "\"\n";
        }
    }
    file << "$EndPhysicalNames\n";
    file << "$Entities\n0 0 " << activeEntities << " 0\n";
    for (size_t tag = 0; tag <= maxTag; ++tag) {
        if (!tagToTriIndices[tag].empty()) {
            file << tag << " " << minX << " " << minY << " " << minZ << " "
                 << maxX << " " << maxY << " " << maxZ << " 1 " << tag << "\n";
        }
    }
    file << "$EndEntities\n";
    size_t numNodes = mesh.nodes.size();
    file << "$Nodes\n1 " << numNodes << " 1 " << (numNodes > 0 ? numNodes : 1) << "\n";
    if (numNodes > 0) {
        file << "2 1 0 " << numNodes << "\n";
        for (size_t i = 0; i < numNodes; ++i) file << (i + 1) << "\n";
        for (size_t i = 0; i < numNodes; ++i) {
            file << mesh.nodes[i].x << " " << mesh.nodes[i].y << " " << mesh.nodes[i].z << "\n";
        }
    }
    file << "$EndNodes\n";
    size_t numElements = mesh.triangles.size();
    file << "$Elements\n" << activeEntities << " " << numElements 
         << " 1 " << (numElements > 0 ? numElements : 1) << "\n";
    size_t elemId = 1;
    for (size_t tag = 0; tag <= maxTag; ++tag) {
        const auto& triIndices = tagToTriIndices[tag];
        if (triIndices.empty()) continue;
        
        file << "2 " << tag << " 2 " << triIndices.size() << "\n";
        for (size_t idx : triIndices) {
            file << elemId++ << " " 
                 << (mesh.triangles[idx].v[0] + 1) << " "
                 << (mesh.triangles[idx].v[1] + 1) << " "
                 << (mesh.triangles[idx].v[2] + 1) << "\n";
        }
    }
    file << "$EndElements\n";
}