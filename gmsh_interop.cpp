#include "gmsh_interop.hpp"
#include <gmsh.h>
#include <unordered_map>

bem::TriangleMesh<3> loadMSH(const std::string& filename) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(filename);
    std::vector<size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords, -1, -1);
    size_t numVerts = nodeTags.size();
    std::unordered_map<size_t, size_t> tagToIndex;
    tagToIndex.reserve(numVerts);
    Eigen::Matrix<double, 3, Eigen::Dynamic> verts(3, numVerts);
    for (size_t i = 0; i < numVerts; ++i) {
        tagToIndex[nodeTags[i]] = i;
        verts(0, i) = coords[3 * i + 0];
        verts(1, i) = coords[3 * i + 1];
        verts(2, i) = coords[3 * i + 2];
    }
    std::vector<int> elementTypes;
    std::vector<std::vector<size_t>> elementTags;
    std::vector<std::vector<size_t>> nodeTagsPerElement;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, -1, -1);
    size_t totalTriangles = 0;
    for (size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == 2) {
            totalTriangles += elementTags[i].size();
        }
    }
    Eigen::Matrix<size_t, 3, Eigen::Dynamic> elems(3, totalTriangles);
    Eigen::Matrix<size_t, 1, Eigen::Dynamic> elem_tags(1, totalTriangles);
    size_t currentElemIdx = 0;
    for (size_t i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] == 2) {
            const auto& nodeTagsList = nodeTagsPerElement[i];
            const auto& tagsList = elementTags[i];
            size_t numTris = tagsList.size();
            for (size_t j = 0; j < numTris; ++j) {
                elems(0, currentElemIdx) = tagToIndex.at(nodeTagsList[3 * j + 0]);
                elems(1, currentElemIdx) = tagToIndex.at(nodeTagsList[3 * j + 1]);
                elems(2, currentElemIdx) = tagToIndex.at(nodeTagsList[3 * j + 2]);
                elem_tags(0, currentElemIdx) = tagsList[j];
                currentElemIdx++;
            }
        }
    }
    gmsh::finalize();
    bem::TriangleMesh<3> bemMesh;
    bemMesh.set_data(verts, elems, elem_tags);
    return bemMesh;
}


void saveMSH(const bem::TriangleMesh<3>& mesh, const std::string& filename) {
    const auto& verts = mesh.verts();
    const auto& elems = mesh.elems();
    const auto& tags = mesh.elem_tags();
    size_t numVerts = verts.cols();
    size_t numElems = elems.cols();
    gmsh::initialize();
    gmsh::model::add("ExportedMesh");
    int surfaceTag = gmsh::model::addDiscreteEntity(2);
    std::vector<double> flatCoords;
    flatCoords.reserve(numVerts * 3);
    std::vector<size_t> nodeTags;
    nodeTags.reserve(numVerts);
    for (size_t i = 0; i < numVerts; ++i) {
        nodeTags.push_back(i + 1);
        flatCoords.push_back(verts(0, i));
        flatCoords.push_back(verts(1, i));
        flatCoords.push_back(verts(2, i));
    }
    gmsh::model::mesh::addNodes(2, surfaceTag, nodeTags, flatCoords);
    int elementType = 2;
    std::vector<size_t> elementNodeTags;
    elementNodeTags.reserve(numElems * 3);
    std::vector<size_t> elementTags;
    elementTags.reserve(numElems);
    for (size_t i = 0; i < numElems; ++i) {
        elementNodeTags.push_back(elems(0, i) + 1);
        elementNodeTags.push_back(elems(1, i) + 1);
        elementNodeTags.push_back(elems(2, i) + 1);
        if (tags.cols() == numElems) {
            elementTags.push_back(tags(0, i));
        }
    }
    gmsh::model::mesh::addElementsByType(surfaceTag, elementType, elementTags, elementNodeTags);
    gmsh::write(filename);
    gmsh::finalize();
}