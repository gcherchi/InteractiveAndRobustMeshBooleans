//
// Created by Michele on 20/06/24.
//

#ifndef GITBOOLEANS_INTERSECT_POINT_RATIONALS_H
#define GITBOOLEANS_INTERSECT_POINT_RATIONALS_H

#include <numerics.h>

class IntersectionPointRationals {

private:
    bigrational x,y,z;
    uint tri_id;
    uint patch_id;

public:

    // Costruttore di default
    IntersectionPointRationals() : x(), y(), z(), tri_id(-1), patch_id(-1) {}

    IntersectionPointRationals(bigrational x, bigrational y, bigrational z, uint tri_id, uint patch_id){
        this->x = x;
        this->y = y;
        this->z = z;
        this->tri_id = tri_id;
        this->patch_id = patch_id;
    }
    IntersectionPointRationals(bigrational x, bigrational y, bigrational z, uint tri_id){
        this->x = x;
        this->y = y;
        this->z = z;
        this->tri_id = tri_id;
    }

    bigrational getX() const {
        return x;
    }

    bigrational getY() const {
        return y;
    }

    bigrational getZ() const {
        return z;
    }

    uint getTriId() const {
        return tri_id;
    }

    uint getPatchId() const {
        return patch_id;
    }

    void setX(bigrational x) {
        this->x = x;
    }

    void setY(bigrational y) {
        this->y = y;
    }

    void setZ(bigrational z) {
        this->z = z;
    }

    void setTriId(uint tri_id) {
        this->tri_id = tri_id;
    }

    void setPatchId(uint patch_id) {
        this->patch_id = patch_id;
    }

    bool operator==(const IntersectionPointRationals &p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator!=(const IntersectionPointRationals &p) const {
        return x != p.x || y != p.y || z != p.z;
    }

    bool operator<(const IntersectionPointRationals &p) const {
        if(x < p.x) return true;
        if(x > p.x) return false;
        if(y < p.y) return true;
        if(y > p.y) return false;
        return z < p.z;
    }

    bool operator>(const IntersectionPointRationals &p) const {
        if(x > p.x) return true;
        if(x < p.x) return false;
        if(y > p.y) return true;
        if(y < p.y) return false;
        return z > p.z;
    }

    bool operator<=(const IntersectionPointRationals &p) const {
        return *this < p || *this == p;
    }

    bool operator>=(const IntersectionPointRationals &p) const {
        return *this > p || *this == p;
    }

    bool lessThanX(const IntersectionPointRationals &p) const {
        return x < p.x;
    }

    bool lessThanY(const IntersectionPointRationals &p) const {
        return y < p.y;
    }

    bool lessThanZ(const IntersectionPointRationals &p) const {
        return z < p.z;
    }

    bool greaterThanX(const IntersectionPointRationals &p) const {
        return x > p.x;
    }

    bool greaterThanY(const IntersectionPointRationals &p) const {
        return y > p.y;
    }

    bool greaterThanZ(const IntersectionPointRationals &p) const {
        return z > p.z;
    }


    static unsigned int findPatchIdByTriId(const std::vector<IntersectionPointRationals>& elements, unsigned int tri_id) {
        for (const auto& element : elements) {
            if (element.tri_id == tri_id) {
                return element.patch_id;
            }
        }
        //return an error message

        assert(true && "Error in findPatchIdByTriId");
    }

    bigrational& operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range("Index out of bounds for IntersectionPointRationals");
    }

    const bigrational& operator[](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range("Index out of bounds for IntersectionPointRationals");
    }

    void printIntersectionPoint(){
        std::cout << "##################################################################" << std::endl;
        std::cout << "Intersection point on triangle with t_id: " << tri_id << std::endl;
        std::cout << "Coords x y z rational: " << x << " --- " << y << " --- " << z << std::endl;
        std::cout << "Coords x y z double: " << x.get_d() << " --- " << y.get_d() << " --- " << z.get_d() << std::endl;
        std::cout << "##################################################################" << std::endl;
    }


};

#endif //GITBOOLEANS_INTERSECT_POINT_RATIONALS_H
