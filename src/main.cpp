#include "boolean_ops.hpp"
#include "mesh_io.hpp"
#include <iostream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

static std::string joinName(const fs::path& folder, const std::string& name) {
    fs::path p = folder;
    p /= name;
    return p.string();
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <meshA_path> <meshB_path> <outFolder> <operation> [removeTouchingSurfaces(0/1)]\n";
        std::cerr << "Operations available: union, intersect, difference, combine\n";
        return 1;
    }

    std::string meshAPath  = argv[1];
    std::string meshBPath  = argv[2];
    std::string outFolder = argv[3];
    std::string operation = argv[4];
    bool removeTouching = false;
    if (argc >= 6) {
        removeTouching = (std::string(argv[5]) == "1" || std::string(argv[5]) == "true");
    }
    if (!fs::exists(outFolder)) {
        fs::create_directories(outFolder);
    }
    std::cout << "Loading meshes...\n";
    MeshData meshA = loadMesh(meshAPath);
    MeshData meshB = loadMesh(meshBPath);
    std::cout << "Performing operation: " << operation << "...\n";
    if (operation == "union") {
        MeshData unionMesh = meshUnion(meshA, meshB);
        std::string outFile = "union_result."+getExtension(meshAPath);
        saveMesh(unionMesh, joinName(outFolder, outFile));
        std::cout << "Saved union results to output folder\n";
    } 
    else if (operation == "intersect") {
        MeshData intersectMesh = meshIntersect(meshA, meshB);
        std::string outFile = "intersect_result."+getExtension(meshAPath);
        saveMesh(intersectMesh, joinName(outFolder, outFile));
        std::cout << "Saved intersect results to output folder\n";
    } 
    else if (operation == "difference") {
        meshDifference(meshA, meshB);
        std::string outFileA = "difference_result_A."+getExtension(meshAPath);
        std::string outFileB = "difference_result_B."+getExtension(meshBPath);
        saveMesh(meshA, joinName(outFolder, outFileA));
        saveMesh(meshB, joinName(outFolder, outFileB));
        std::cout << "Saved difference results to output folder.\n";
    } 
    else if (operation == "combine") {
        meshCombine(meshA, meshB, removeTouching);
        std::string outFileA = "combine_result_A."+getExtension(meshAPath);
        std::string outFileB = "combine_result_B."+getExtension(meshBPath);
        saveMesh(meshA, joinName(outFolder, outFileA));
        saveMesh(meshB, joinName(outFolder, outFileB));
        std::cout << "Saved combine results to output folder.\n";
    } 
    else {
        std::cerr << "Error: Unknown operation '" << operation << "'\n";
        return 1;
    }

    std::cout << "Done!\n";
    return 0;
}