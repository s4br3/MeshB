#include "mesh_clean.hpp"
#include <numeric>

MeshData filterMesh(const MeshData& mesh, const std::vector<bool>& remove) {
    const size_t len = remove.size();
    MeshData out;
    out.nodes = mesh.nodes;
    out.triangles.reserve(len);
    out.tags.reserve(len);
    out.centres.reserve(len);
    out.normals.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        if (!remove[i]) {
            out.triangles.push_back(mesh.triangles[i]);
            if (!mesh.tags.empty()) out.tags.push_back(mesh.tags[i]);
            if (i < mesh.centres.size()) out.centres.push_back(mesh.centres[i]);
            if (i < mesh.normals.size()) out.normals.push_back(mesh.normals[i]);
        }
    }
    return out;
}

void clearUnusedNodes(MeshData& mesh) {
    const size_t n = mesh.nodes.size();
    if (n == 0) return;
    std::vector<uint8_t> used(n, 0);
    for (const Triangle& t : mesh.triangles) {
        used[t.v[0]] = 1;
        used[t.v[1]] = 1;
        used[t.v[2]] = 1;
    }
    const size_t SENT = static_cast<size_t>(-1);
    std::vector<size_t> remap(n, SENT);
    size_t newCount = 0;
    for (size_t old = 0; old < n; ++old) {
        if (used[old]) remap[old] = newCount++;
    }
    std::vector<Vec3> newNodes;
    for (size_t old = 0; old < n; ++old) {
        if (remap[old] != SENT) {
            newNodes.push_back(mesh.nodes[old]);
        }
    }
    for (Triangle& t : mesh.triangles) {
        t.v[0] = remap[t.v[0]];
        t.v[1] = remap[t.v[1]];
        t.v[2] = remap[t.v[2]];
    }
    mesh.nodes = std::move(newNodes);
}

size_t findRoot(std::vector<size_t>& parent, size_t i) {
    if (parent[i] == i)
        return i;
    return parent[i] = findRoot(parent, parent[i]); 
}

void unionNodes(std::vector<size_t>& parent, size_t i, size_t j) {
    size_t rootI = findRoot(parent, i);
    size_t rootJ = findRoot(parent, j);
    if (rootI != rootJ) {
        if (rootI < rootJ) parent[rootJ] = rootI;
        else parent[rootI] = rootJ;
    }
}

void collapseSlivers(MeshData& mesh, double eps) {
    const double eps2 = eps * eps;
    const size_t nNodes = mesh.nodes.size();
    if (nNodes == 0) return;
    std::vector<size_t> parent(nNodes);
    std::iota(parent.begin(), parent.end(), 0);
    for (const Triangle& t : mesh.triangles) {
        for (int i = 0; i < 3; ++i) {
            size_t v0 = t.v[i];
            size_t v1 = t.v[(i + 1) % 3];
            Vec3 diff = mesh.nodes[v1] - mesh.nodes[v0];
            if (BVH::dot(diff, diff) < eps2) {
                unionNodes(parent, v0, v1);
            }
        }
    }
    for (Triangle& t : mesh.triangles) {
        t.v[0] = findRoot(parent, t.v[0]);
        t.v[1] = findRoot(parent, t.v[1]);
        t.v[2] = findRoot(parent, t.v[2]);
    }
}

void removeDegenerateTriangles(MeshData& mesh) {
    std::vector<Triangle> goodTriangles;
    std::vector<size_t> goodTags;
    std::vector<Vec3> goodCentres;
    std::vector<Vec3> goodNormals;

    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const Triangle& t = mesh.triangles[i];
        if (t.v[0] != t.v[1] && t.v[1] != t.v[2] && t.v[2] != t.v[0]) {
            goodTriangles.push_back(t);
            goodTags.push_back(mesh.tags[i]);
            goodCentres.push_back(mesh.centres[i]);
            goodNormals.push_back(mesh.normals[i]);
        }
    }

    mesh.triangles = std::move(goodTriangles);
    mesh.tags = std::move(goodTags);
    mesh.centres = std::move(goodCentres);
    mesh.normals = std::move(goodNormals);
}

void cleanMesh(MeshData& mesh, double eps) {
    collapseSlivers(mesh, eps);
    removeDegenerateTriangles(mesh);
}

void invertWinding(MeshData& mesh) {
    for (Triangle& t : mesh.triangles) {
        std::swap(t.v[0], t.v[1]);
    }
    for (Vec3& n : mesh.normals) {
        n = n * (-1.0); 
    }
}

MeshData combineMeshes(const std::vector<MeshData>& meshes, double eps) {
    MeshData combined;
    SpatialGrid3D grid(eps);
    size_t tagOffset = 0;
    for (const auto& mesh : meshes) {
        std::vector<size_t> remap(mesh.nodes.size());
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            remap[i] = grid.getOrAdd(mesh.nodes[i]);
        }
        for (size_t i = 0; i < mesh.triangles.size(); ++i) {
            Triangle t = mesh.triangles[i];
            t.v[0] = remap[t.v[0]];
            t.v[1] = remap[t.v[1]];
            t.v[2] = remap[t.v[2]];
            combined.triangles.push_back(t);
            if (i < mesh.centres.size()) combined.centres.push_back(mesh.centres[i]);
            if (i < mesh.normals.size()) combined.normals.push_back(mesh.normals[i]);
        }
        size_t currentMeshMaxTag = 0;
        for (size_t tag : mesh.tags) {
            combined.tags.push_back(tag + tagOffset);
            currentMeshMaxTag = std::max(currentMeshMaxTag, tag);
        }
        if (!mesh.tags.empty()) {
            tagOffset += currentMeshMaxTag + 1;
        }
    }
    combined.nodes = grid.getUniquePoints();
    return combined;
}