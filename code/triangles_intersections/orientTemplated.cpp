#include "orientTemplated.h"

template <typename T>
inline int orient3dT( const T * pa, const T * pb, const T * pc, const T * pd) {

    using ElementType = typename std::remove_extent<T>::type;

    if constexpr (std::is_same<ElementType,double>::value) {
        double result = cinolib::orient3d(pa, pb, pc, pd);
        if (result < 0) return -1;
        else if (result > 0) return 1;
        else return 0;
    }
    else if constexpr (std::is_same<ElementType, genericPoint>::value) {
        return -1 * genericPoint::orient3D(*pa, *pb, *pc, *pd);
    }

    else {
        T adx = pa[0] - pd[0];
        T bdx = pb[0] - pd[0];
        T cdx = pc[0] - pd[0];
        T ady = pa[1] - pd[1];
        T bdy = pb[1] - pd[1];
        T cdy = pc[1] - pd[1];
        T adz = pa[2] - pd[2];
        T bdz = pb[2] - pd[2];
        T cdz = pc[2] - pd[2];

        T res =  adx * (bdy * cdz - bdz * cdy)
             + bdx * (cdy * adz - cdz * ady)
             + cdx * (ady * bdz - adz * bdy);

        if constexpr (std::is_same<ElementType, bigrational>::value) {
            if (res == bigrational()) return 0;
            else if (res > bigrational()) return 1;
            else return -1;
        }
        else {
            if (res < 0) return -1;
            else if (res > 0) return 1;
            else return 0;
        }
    }
}

template <typename T>
inline int orient2dT( const T * pa, const T * pb, const T * pc) {

    using ElementType = typename std::remove_extent<T>::type;

    if constexpr (std::is_same<ElementType, double>::value) {
        //Double
        double result = cinolib::orient2d(pa, pb, pc);
        if (result < 0) return -1;
        else if (result > 0) return 1;
        else return 0;
    } else if constexpr (std::is_same<ElementType, genericPoint>::value) {
        return genericPoint::orient2D(*pa, *pb, *pc);
    } else {
        T acx = pa[0] - pc[0];
        T bcx = pb[0] - pc[0];
        T acy = pa[1] - pc[1];
        T bcy = pb[1] - pc[1];

        T res =  acx * bcy - acy * bcx;

        if constexpr (std::is_same<ElementType, bigrational>::value) {
            if (res == bigrational()) return 0;
            else if (res > bigrational()) return 1;
            else return -1;
        }
        else {
            if (res < 0) return -1;
            else if (res > 0) return 1;
            else return 0;
        }
    }
}