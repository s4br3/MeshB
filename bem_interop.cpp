#include "bem_interop.hpp"

MeshData extractMeshData(const bem::TriangleMesh<3>& mesh) {
    MeshData data;
    const auto& verts = mesh.verts();
    const auto& elems = mesh.elems();
    const auto& tags = mesh.elem_tags();
    size_t num_verts = verts.cols();
    data.nodes.reserve(num_verts);
    for (size_t i = 0; i < num_verts; ++i) {
        data.nodes.push_back(Vec3(verts(0, i), verts(1, i), verts(2, i)));
    }
    size_t count = elems.cols();
    data.triangles.reserve(count);
    data.centres.reserve(count);
    data.normals.reserve(count);
    data.tags.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        Triangle tri = {
            (elems(0, i)),
            (elems(1, i)),
            (elems(2, i))
        };
        data.triangles.push_back(tri);
        data.centres.push_back(tri.centre(data.nodes));
        data.normals.push_back(tri.normal(data.nodes));
        data.tags.push_back(tags[i]);
    }
    return data;
}

void rebuildMesh(bem::TriangleMesh<3>& mesh, const MeshData& newMesh, double eps) {
    Eigen::Matrix<double, 3, Eigen::Dynamic> verts(3, newMesh.nodes.size());
    Eigen::Matrix<size_t, 3, Eigen::Dynamic> elems(3, newMesh.triangles.size());
    Eigen::Matrix<size_t, 1, Eigen::Dynamic> elem_tags(1, newMesh.tags.size());
    for (size_t i = 0; i < newMesh.nodes.size(); ++i) {
        verts(0, i) = newMesh.nodes[i][0];
        verts(1, i) = newMesh.nodes[i][1];
        verts(2, i) = newMesh.nodes[i][2];
    }
    for (size_t i = 0; i < newMesh.triangles.size(); ++i) {
        elems(0, i) = newMesh.triangles[i].v[0];
        elems(1, i) = newMesh.triangles[i].v[1];
        elems(2, i) = newMesh.triangles[i].v[2];
        elem_tags(0, i) = newMesh.tags[i];
    }
    mesh.set_data(verts, elems, elem_tags);
}