MeshB: A boolean operations solver that hopes to match up to predicate-based counterparts with floating-point arithmetic, built for researchers trying to model electromagnetic flow (though it can work for others just as well).

Usage: Your mesh file inputs, the folder where the mesh(es) should be outputted and the operation to be carried out.
If your mesh has lots of (>4) n-gons, please wait a few hours lol. This system needs a triangle mesh, so it forces triangulation through a bulky subsystem that I have yet to optimise. If your mesh is made almost exclusively of tris or quads (less than a thousand (>4) n-gons) then you have no worries. It'll work in 30 seconds tops.

Additionally, any program can directly call the functions from boolean_ops.hpp, whether on Shaswat Sharma's OpenBEM's TriangleMesh<3> struct or my own MeshData which holds vector(Vec3(double)) nodes and vector (Triangle(array(size_t, 3))) triangles, where the triangles are made up of indices of the associated nodes. You can also call individual functions directly from any header if you want to operate the pipeline your own way.

Build: Gmsh must be installed on your system. Besides this, the only library used is Shashwat Sharma's OpenBEM finite-element method solver (for TriangleMesh<3> compatibility), which is fetched by Cmake.

Contributing: Unlikely anybody will care about this library, but if you happen to be an odd manner of human, you can just send in pull requests, and I'll see what I can do.

Credits:

Shashwat Sharma: Supervisor for the project and great advisor for when I had no clue which way to go.

Blender: Wonderful 3D model visualisation tool. I ask that you donate the cost of a coffee once in a while. I have found much value from their work and have done accordingly.

Gmsh: Solid 3D model visualisation tool. Sorry, I'm not as emotionally attached to this software, as useful as it has been through this process. There's also no foundation you can give money to to support the development, which left me feeling a bit miffed (I did find great use from this app). The API, on the other hand, is superb: it's how this system cuts triangles.

My mom: Excellent audience when I'm explaining ideas I'm trying to wrangle for this program. Asked great questions as well, and was very supportive.

S4BRE: Myself. This project was heaps of fun and I am stoked to have brought it this far.

If you have any questions or difficulties with the software, feel free to post as such in the issues section of this repo. Enjoy your nicely triangulated meshes y'all.
