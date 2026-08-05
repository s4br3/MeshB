#include "boolean_ops.hpp"
#include <iostream>

int main(int argc, char** argv) {
    std::vector<Vec3> vertexPool = {
        Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0),
        Vec3(0, 0, -1), Vec3(1, 0, -1), Vec3(0, 1, 1)
    };
    Triangle t1 = {0, 1, 2};
    Triangle t2 = {3, 4, 5};

    std::cout << "Initialization successful!" << std::endl;
    return 0;
}