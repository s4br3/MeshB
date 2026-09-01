#include "mesh_types.hpp"
std::ostream& operator<<(std::ostream& os, const std::array<Vec3, 3>& vertices) {
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
std::ostream& operator<<(std::ostream& os, const Connection& conn) {
    os << "Intersections A -> B:\n";
    for (const auto& [a, map] : conn.intersectionsAB) {
        os << "  Triangle A [" << a << "]:\n";
        for (const auto& [b, vecs] : map) {
            os << "    -> Triangle B [" << b << "]: ";
            for (size_t i = 0; i < vecs.size(); ++i) {
                os << vecs[i] << (i + 1 < vecs.size() ? ", " : "");
            }
            os << "\n";
        }
    }
    os << "\n";
    os << "Intersections B -> A:\n";
    for (const auto& [b, map] : conn.intersectionsBA) {
        os << "  Triangle B [" << b << "]:\n";
        for (const auto& [a, vecs] : map) {
            os << "    -> Triangle A [" << a << "]: ";
            for (size_t i = 0; i < vecs.size(); ++i) {
                os << vecs[i] << (i + 1 < vecs.size() ? ", " : "");
            }
            os << "\n";
        }
    }
    os << "\n";

    os << "Overlaps A -> B:\n";
    for (const auto& [a, list] : conn.overlapsAB) {
        os << "  Triangle A [" << a << "] overlaps with B tags: [";
        for (size_t i = 0; i < list.size(); ++i) {
            os << list[i] << (i + 1 < list.size() ? ", " : "");
        }
        os << "]\n";
    }
    os << "\n";
    os << "Overlaps B -> A:\n";
    for (const auto& [b, list] : conn.overlapsBA) {
        os << "  Triangle B [" << b << "] overlaps with A tags: [";
        for (size_t i = 0; i < list.size(); ++i) {
            os << list[i] << (i + 1 < list.size() ? ", " : "");
        }
        os << "]\n";
    }
    os << "\n";
    os << "Overlap Percentages A -> B:\n";
    for (const auto& [a, map] : conn.overlapPercentAB) {
        os << "  Triangle A [" << a << "]:\n";
        for (const auto& [b, percent] : map) {
            os << "    -> Triangle B [" << b << "]: " << percent << "%\n";
        }
    }
    os << "\n";
    os << "Overlap Percentages B -> A:\n";
    for (const auto& [b, map] : conn.overlapPercentBA) {
        os << "  Triangle B [" << b << "]:\n";
        for (const auto& [a, percent] : map) {
            os << "    -> Triangle A [" << a << "]: " << percent << "%\n";
        }
    }
    return os;
}
void recomputeMeshData(MeshData &mesh) {
    mesh.centres.clear();
    mesh.normals.clear();
    const size_t nTri = mesh.triangles.size();
    mesh.centres.reserve(nTri);
    mesh.normals.reserve(nTri);
    mesh.tags.reserve(nTri > mesh.tags.size() ? nTri - mesh.tags.size() : 0);
    size_t maxTag = 0;
    for (size_t t : mesh.tags) {
        if (t > maxTag) maxTag = t;
    }
    const size_t tagged0 = mesh.tags.size();
    while (mesh.tags.size() < nTri) {
        mesh.tags.push_back(maxTag + mesh.tags.size() - tagged0);
    }
    for (size_t i = 0; i < nTri; ++i) {
        mesh.centres.push_back(mesh.triangles[i].centre(mesh.nodes));
        mesh.normals.push_back(mesh.triangles[i].normal(mesh.nodes).normalize());
    }
}
