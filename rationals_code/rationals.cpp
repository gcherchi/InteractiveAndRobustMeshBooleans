/********************************************************************************
*  This file is part of CinoLib                                                 *
*  Copyright(C) 2023: Marco Livesu                                              *
*                                                                               *
*  The MIT License                                                              *
*                                                                               *
*  Permission is hereby granted, free of charge, to any person obtaining a      *
*  copy of this software and associated documentation files (the "Software"),   *
*  to deal in the Software without restriction, including without limitation    *
*  the rights to use, copy, modify, merge, publish, distribute, sublicense,     *
*  and/or sell copies of the Software, and to permit persons to whom the        *
*  Software is furnished to do so, subject to the following conditions:         *
*                                                                               *
*  The above copyright notice and this permission notice shall be included in   *
*  all copies or substantial portions of the Software.                          *
*                                                                               *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR   *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
*  FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      *
*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS *
*  IN THE SOFTWARE.                                                             *
*                                                                               *
*  Author(s):                                                                   *
*                                                                               *
*     Marco Livesu (marco.livesu@gmail.com)                                     *
*     http://pers.ge.imati.cnr.it/livesu/                                       *
*                                                                               *
*     Italian National Research Council (CNR)                                   *
*     Institute for Applied Mathematics and Information Technologies (IMATI)    *
*     Via de Marini, 6                                                          *
*     16149 Genoa,                                                              *
*     Italy                                                                     *
*********************************************************************************/
#include <rationals.h>
#include <cinolib/predicates.h>
#include <numerics.h>

namespace cinolib
{

CINO_INLINE
bool rationals_are_working()
{
    CGAL_Q a[3] = { 0,0,0 };
    CGAL_Q b[3] = { 1,2,3 };
    CGAL_Q c[3] = { 9,3,1 };
    CGAL_Q d[3];
    midpoint(a,b,c,d);

    double aa[3] = { 0,0,0 };
    double bb[3] = { 1,2,3 };
    double cc[3] = { 9,3,1 };
    double dd[3] = { (aa[0]+bb[0]+cc[0])/3,
                     (aa[1]+bb[1]+cc[1])/3,
                     (aa[2]+bb[2]+cc[2])/3 };

    if(orient3d(a,b,c,d)==0 && orient3d(aa,bb,cc,dd)!=0)
    {
        //std::cout << "CGAL::Lazy_exact_nt<CGAL::Gmpq> CHECK ... OK" << std::endl;
        return true;
    }
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

template <class T>
CINO_INLINE
    /*
T orient3d(const T * pa,
           const T * pb,
           const T * pc,
           const T * pd)
{

    std::cout << "NEW VALUES :::::::::" << std::endl;

    std::cout << "pa: " << pa[0] <<  " " << pa[1] << " " << pa[2] << std::endl;
    std::cout << "pb: " << pb[0] <<  " " << pb[1] << " " << pb[2] << std::endl;
    std::cout << "pc: " << pc[0] <<  " " << pc[1] << " " << pc[2] << std::endl;
    std::cout << "pd: " << pd[0] <<  " " << pd[1] << " " << pd[2] << std::endl;

    T adx = pa[0] - pd[0];
    T bdx = pb[0] - pd[0];
    T cdx = pc[0] - pd[0];
    T ady = pa[1] - pd[1];
    T bdy = pb[1] - pd[1];
    T cdy = pc[1] - pd[1];
    T adz = pa[2] - pd[2];
    T bdz = pb[2] - pd[2];
    T cdz = pc[2] - pd[2];

    std::cout << "adx: " << adx << std::endl;
    std::cout << "bdx: " << bdx << std::endl;
    std::cout << "cdx: " << cdx << std::endl;
    std::cout << "ady: " << ady << std::endl;
    std::cout << "bdy: " << bdy << std::endl;
    std::cout << "cdy: " << cdy << std::endl;
    std::cout << "adz: " << adz << std::endl;
    std::cout << "bdz: " << bdz << std::endl;
    std::cout << "cdz: " << cdz << std::endl;


    return adx * ((bdy * cdz) - (bdz * cdy))
         + bdx * ((cdy * adz) - (cdz * ady))
         + cdx * ((ady * bdz) - (adz * bdy));
}*/
    T orient3d(const T* pa, const T* pb, const T* pc, const T* pd)
{
    std::cout << "Il tipo della variabile è: " << typeid(decltype(pa)).name() << std::endl;

    std::array<bigrational,3 > s0_tmp = {bigrational(3074457343958179937,153576,1),
                          bigrational(3572095111,235216,1),
                          bigrational(9320677326,108497,1)
    };

    std::array<bigrational,3 > s1_tmp = {bigrational(116913819,89206,-1),
                         bigrational(2305843008117596711,17776,1),
                         bigrational()
    };

    std::array<bigrational,3 > s2_tmp = {bigrational(5504579610,1,-1),
                         bigrational(3074457344089985827,40257,-1),
                         bigrational()
    };

    std::array<bigrational,3 > s3_tmp = {bigrational(18446744070889277111,918766,-1),
                         bigrational(5870338011,958118,-1),
                         bigrational(18446744071096704124,682923,1)
    };
/*
    pa = s0_tmp.data();
    pb = s1_tmp.data();
    pc = s2_tmp.data();
    pd = s3_tmp.data();*/

    std::cout << "pa: " << pa[0] << " " << pa[1] << " " << pa[2] << std::endl;
    std::cout << "pb: " << pb[0] << " " << pb[1] << " " << pb[2] << std::endl;
    std::cout << "pc: " << pc[0] << " " << pc[1] << " " << pc[2] << std::endl;
    std::cout << "pd: " << pd[0] << " " << pd[1] << " " << pd[2] << std::endl;

    // Calcolo delle differenze tra i punti
    T adx = pa[0] - pd[0];
    T bdx = pb[0] - pd[0];
    T cdx = pc[0] - pd[0];
    T ady = pa[1] - pd[1];
    T bdy = pb[1] - pd[1];
    T cdy = pc[1] - pd[1];
    T adz = pa[2] - pd[2];
    T bdz = pb[2] - pd[2];
    T cdz = pc[2] - pd[2];


    // Stampa delle differenze
    std::cout << "adx: " << adx << std::endl;
    std::cout << "bdx: " << bdx << std::endl;
    std::cout << "cdx: " << cdx << std::endl;
    std::cout << "ady: " << ady << std::endl;
    std::cout << "bdy: " << bdy << std::endl;
    std::cout << "cdy: " << cdy << std::endl;
    std::cout << "adz: " << adz << std::endl;
    std::cout << "bdz: " << bdz << std::endl;
    std::cout << "cdz: " << cdz << std::endl;

    // Prima moltiplicazione e sottrazione
    T bdy_cdz = bdy * cdz;
    T bdz_cdy = bdz * cdy;
    std::cout << "bdy * cdz = " << bdy_cdz << std::endl;
    std::cout << "bdz * cdy = " << bdz_cdy << std::endl;
    T term1 = bdy_cdz - bdz_cdy;
    std::cout << "(bdy * cdz) - (bdz * cdy) = " << term1 << std::endl;

    // Moltiplicazione esterna alla prima parentesi
    T adx_term1 = adx * term1;
    std::cout << "adx * ((bdy * cdz) - (bdz * cdy)) = " << adx_term1 << std::endl;

    // Seconda moltiplicazione e sottrazione
    T cdy_adz = cdy * adz;
    T cdz_ady = cdz * ady;
    std::cout << "cdy * adz = " << cdy_adz << std::endl;
    std::cout << "cdz * ady = " << cdz_ady << std::endl;
    T term2 = cdy_adz - cdz_ady;
    std::cout << "(cdy * adz) - (cdz * ady) = " << term2 << std::endl;

    // Moltiplicazione esterna alla seconda parentesi
    T bdx_term2 = bdx * term2;
    std::cout << "bdx * ((cdy * adz) - (cdz * ady)) = " << bdx_term2 << std::endl;

    // Terza moltiplicazione e sottrazione
    T ady_bdz = ady * bdz;
    T adz_bdy = adz * bdy;
    std::cout << "ady * bdz = " << ady_bdz << std::endl;
    std::cout << "adz * bdy = " << adz_bdy << std::endl;
    T term3 = ady_bdz - adz_bdy;
    std::cout << "(ady * bdz) - (adz * bdy) = " << term3 << std::endl;

    // Moltiplicazione esterna alla terza parentesi
    T cdx_term3 = cdx * term3;
    std::cout << "cdx * ((ady * bdz) - (adz * bdy)) = " << cdx_term3 << std::endl;

    // Somma finale
    T result = adx_term1 + bdx_term2 + cdx_term3;
    std::cout << "Final result = " << result << std::endl;

    return result;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

template <class T>
CINO_INLINE
T orient2d(const T * pa,
           const T * pb,
           const T * pc)
{
    T acx = pa[0] - pc[0];
    T bcx = pb[0] - pc[0];
    T acy = pa[1] - pc[1];
    T bcy = pb[1] - pc[1];

    return acx * bcy - acy * bcx;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void midpoint(const CGAL_Q * pa,
              const CGAL_Q * pb,
              const CGAL_Q * pc,
                    CGAL_Q * res) // res = (pa + pb + pc)/3
{
    res[0] = (pa[0] + pb[0] + pc[0])/3;
    res[1] = (pa[1] + pb[1] + pc[1])/3;
    res[2] = (pa[2] + pb[2] + pc[2])/3;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void midpoint(const CGAL_Q * pa,
              const CGAL_Q * pb,
                    CGAL_Q * res) // res = (pa + pb)/2
{
    res[0] = (pa[0] + pb[0])/2;
    res[1] = (pa[1] + pb[1])/2;
    res[2] = (pa[2] + pb[2])/2;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void line_intersection2d(const CGAL_Q * pa,
                         const CGAL_Q * pb,
                         const CGAL_Q * pc,
                         const CGAL_Q * pd,
                               CGAL_Q * res) // intersection point of **NON PARALLEL** lines (pa,pb) and (pc,pd)
{
    // https://en.wikipedia.org/wiki/Line–line_intersection
    CGAL_Q det_ab   = pa[0]*pb[1] - pa[1]*pb[0];
    CGAL_Q det_cd   = pc[0]*pd[1] - pc[1]*pd[0];
    CGAL_Q det_ab_x = pa[0] - pb[0];
    CGAL_Q det_cd_x = pc[0] - pd[0];
    CGAL_Q det_ab_y = pa[1] - pb[1];
    CGAL_Q det_cd_y = pc[1] - pd[1];

    CGAL_Q den = det_ab_x*det_cd_y - det_ab_y*det_cd_x;
    assert(den!=0);

    res[0] = (det_ab*det_cd_x - det_ab_x*det_cd)/den;
    res[1] = (det_ab*det_cd_y - det_ab_y*det_cd)/den;

    // sanity checks
    assert(orient2d(pa,pb,res)==0);
    assert(orient2d(pc,pd,res)==0);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
CGAL_Q sqrd_distance2d(const CGAL_Q * pa,
                       const CGAL_Q * pb)
{
    return (pa[0]-pb[0])*(pa[0]-pb[0]) +
           (pa[1]-pb[1])*(pa[1]-pb[1]);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void copy(const CGAL_Q * pa, CGAL_Q * pb)
{
    pb[0] = pa[0];
    pb[1] = pa[1];
    pb[2] = pa[2];
}

}





