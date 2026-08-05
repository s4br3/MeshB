#include "mesh_types.hpp"
std::ostream& operator<<(std::ostream& os, const std::array<Vec3, 3>& vertices){
    return os << "[" << vertices[0] << ", " << vertices[1] << ", " << vertices[2] << "]";
}

std::ostream& operator<<(std::ostream& os, const MeshData& mesh) {
    os << "MeshData{\n";
    os << "  nodes (" << mesh.nodes.size() << "):\n";
    for (const auto& n : mesh.nodes) {
        os << "    " << n << "\n";
    }
    os << "  triangles (" << mesh.triangles.size() << "):\n";
    for (const auto& t : mesh.triangles) {
        os << "    " << t.getVertices(mesh.nodes) << "\n";
    }
    os << "  tags (" << mesh.tags.size() << "): ";
    for (size_t i = 0; i < mesh.tags.size(); ++i) {
        if (i) os << ", ";
        os << mesh.tags[i];
    }
    os << "\n";
    os << "  centres (" << mesh.centres.size() << "):\n";
    for (const auto& c : mesh.centres) os << "    " << c << "\n";
    os << "  normals (" << mesh.normals.size() << "):\n";
    for (const auto& nn : mesh.normals) os << "    " << nn << "\n";

    os << "}\n";
    return os;
}

void recomputeMeshData(MeshData &mesh){
    int maxTag = 0;
    int tagged = mesh.tags.size();
    for (int i = 0; i < mesh.triangles.size(); i++){
        if (i < tagged){
            maxTag = maxTag > mesh.tags[i] ? maxTag : mesh.tags[i];
        }
        else{
            mesh.tags.push_back(maxTag + i - tagged);
        }
        mesh.centres.push_back(mesh.triangles[i].centre(mesh.nodes));
        mesh.normals.push_back(mesh.triangles[i].normal(mesh.nodes));
    }
}