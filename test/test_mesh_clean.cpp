#include "mesh_clean.hpp"
#include "mesh_types.hpp"
#include <cassert>
#include <iostream>

constexpr double EPS = 1e-6;

void test_dsu_operations() {
    std::vector<size_t> parent(5);
    for (size_t i = 0; i < 5; ++i) parent[i] = i;

    unionNodes(parent, 1, 2);
    unionNodes(parent, 3, 4);
    unionNodes(parent, 2, 4);

    assert(findRoot(parent, 1) == findRoot(parent, 3));
    assert(findRoot(parent, 0) == 0);
}

void test_filter_and_clear() {
    MeshData mesh;
    mesh.nodes = {{0,0,0}, {1,0,0}, {0,1,0}, {10,10,10}};
    mesh.triangles = {{0,1,2}, {0,1,3}};
    mesh.tags = {1, 2};
    mesh.centres = {{0.3,0.3,0}, {3.6,3.6,3.3}};
    mesh.normals = {{0,0,1}, {0,0,1}};

    // Filter second triangle
    std::vector<bool> remove = {false, true};
    filterMesh(mesh, remove);

    assert(mesh.triangles.size() == 1);
    assert(mesh.tags.size() == 1 && mesh.tags[0] == 1);

    // Clear unused node {10,10,10} (index 3)
    clearUnusedNodes(mesh);
    assert(mesh.nodes.size() == 3);
    assert(mesh.triangles[0].v[0] == 0 && mesh.triangles[0].v[1] == 1 && mesh.triangles[0].v[2] == 2);
}

void test_sliver_and_degeneracy() {
    MeshData mesh;
    mesh.nodes = {{0,0,0}, {1e-8,0,0}, {0,1,0}, {1,0,0}};
    mesh.triangles = {
        {0, 1, 2},
        {0, 2, 3},
        {2, 2, 3}
    };
    mesh.tags = {10, 20, 30};

    collapseSlivers(mesh, EPS);
    assert(mesh.triangles[0].v[0] == mesh.triangles[0].v[1]);

    removeDegenerateTriangles(mesh);
    assert(mesh.triangles.size() == 1);
    assert(mesh.tags[0] == 20);

    // Test invert winding
    invertWinding(mesh);
    assert(mesh.triangles[0].v[0] == 2 && mesh.triangles[0].v[1] == 0);
}

void test_mesh_combine() {
    MeshData m1;
    m1.nodes = {{0,0,0}, {1,0,0}, {0,1,0}};
    m1.triangles = {{0,1,2}};
    m1.tags = {0};

    MeshData m2;
    m2.nodes = {{1,0,0}, {1,1,0}, {0,1,0}}; // Duplicates nodes {1,0,0} and {0,1,0}
    m2.triangles = {{0,1,2}};
    m2.tags = {0};

    MeshData combined = combineMeshes({m1, m2}, EPS);
    assert(combined.triangles.size() == 2);
    assert(combined.nodes.size() == 4); // Welded to 4 unique nodes
}

int main() {
    test_dsu_operations();
    test_filter_and_clear();
    test_sliver_and_degeneracy();
    test_mesh_combine();
    std::cout << "test_mesh_clean passed successfully.\n";
    return 0;
}