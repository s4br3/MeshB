#include <initializer_list>
#include <stdexcept>
struct Vector {
    double x = 0.0, y = 0.0, z = 0.0;
    int dim = 2;
    static Vector from2D(double x_, double y_) {
        Vector v; v.x = x_; v.y = y_; v.z = 0.0; v.dim = 2; return v;
    }
    static Vector from3D(double x_, double y_, double z_) {
        Vector v; v.x = x_; v.y = y_; v.z = z_; v.dim = 3; return v;
    }
    Vector() = default;
    Vector(std::initializer_list<double> xs) {
        auto it = xs.begin();
        int n = static_cast<int>(xs.size());
        if (n == 2) {
            x = *it++; y = *it++; z = 0.0; dim = 2;
        } else if (n == 3) {
            x = *it++; y = *it++; z = *it++; dim = 3;
        } else {
            throw std::invalid_argument("Vector must be constructed with {x,y} or {x,y,z}");
        }
    }
    double operator[](size_t i) const {
        switch (i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Vector index out of range");
        }
    }
    Vector operator+(const Vector& rhs) const {
        Vector r;
        r.x = x + rhs.x;
        r.y = y + rhs.y;
        r.z = z + rhs.z;
        r.dim = std::max(dim, rhs.dim);
        return r;
    }

    Vector operator-(const Vector& rhs) const {
        Vector r;
        r.x = x - rhs.x;
        r.y = y - rhs.y;
        r.z = z - rhs.z;
        r.dim = std::max(dim, rhs.dim);
        return r;
    }
    double dot(const Vector& rhs) const {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }
    Vector cross(const Vector& rhs) const {
        return Vector::from3D(
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        );
    }
    double cross2D(const Vector& rhs) const {
        return x * rhs.y - y * rhs.x;
    }
};
