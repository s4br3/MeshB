# MeshB

**MeshB** is a fast boolean operations solver using floating-point arithmetic designed for electromagnetic flow modeling (and general mesh applications).

---

## Prerequisites & Dependencies
* **CMake**: 3.14 or newer
* **Compiler**: C++20-compatible
* **Tools & Libraries**: Git, OpenMP, and Gmsh (`sudo apt install libgmsh-dev`)
* **Network**: An active internet connection during the first configuration (OpenBEM is automatically fetched from GitHub)

---

## Build Instructions

Run the following commands from the repository root:

```bash
cmake -S . -B build
cmake --build build
```

Run the generated executable:
```bash
./build/main
```

---

## Testing
From root, run ```./build/(executable)``` where the executable is any of
```
test_bvh
test_classify
test_geom_2d
test_geom_3d
test_math_utils
test_mesh_clean
test_raycast
test_triangulation
test_vector
```

---

## Usage & Performance

### Inputs
Provide your mesh file inputs, the output folder destination, and the desired boolean operation.

### Programmatic API
You can call functions directly from `boolean_ops.hpp` or individual headers:
* Compatible with **Shashwat Sharma's OpenBEM** `TriangleMesh<3>` struct.
* Compatible with MeshB's native `MeshData` struct (holds `vector<Vec3<double>>` nodes and `vector<Triangle<array<size_t, 3>>>` triangles).

### Performance Note on Triangulation
MeshB requires a triangle mesh and forces triangulation via an internal subsystem:
* **Standard Meshes (mostly tris/quads, < 1,000 N-gons with > 4 vertices)**: Executes in **30 seconds or less**.
* **Heavy N-Gon Meshes (lots of N-gons with > 4 vertices)**: The current unoptimized triangulation subsystem can take **a few hours**.

### Notes
For some reason, the app is unable to successfully verify mesh conformity when the meshes have coplanar overlaps, though it does work for non-coplanar intersections. I feel like I had definitely seen it successfully work on coplanar intersections previously, but a lot has changed in the codebase so maybe that's why: maybe I'm misremembering. I don't quite know.
Running a self-boolean results in the filter function not properly knowing whether a face should be kept or not, so output ends up quite badly (I suppose that degeneracy is why they don't let you run self-booleans normally).

---

## Contributing & Support
* **Contributions**: Pull requests welcome. Feel free to send one over.
* **Questions / Issues**: If you run into any difficulties or have questions, please post them in the **Issues** section of the repository.

---

## Credits & Acknowledgments
* **S4BRE**: Creator and lead developer.
* **Shashwat Sharma**: Supervisor and project advisor.
* **Gmsh**: Excellent API used to cut triangles, and a solid 3D model visualization tool.
* **Blender**: Wonderful 3D model visualization tool.
* **My Mom**: Great audience, asked lots of questions, and was really supportive. She helped wrangle complex ideas.
