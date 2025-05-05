/****************************************************************************
* Indirect predicates for geometric constructions					        *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         * 
* IMATI-GE / CNR                                                            * 
*                                                                           *
* Authors: Marco Attene                                                     * 
* Copyright(C) 2019: IMATI-GE / CNR                                         * 
* All rights reserved.                                                      * 
*                                                                           *
* This program is free software; you can redistribute it and/or modify      * 
* it under the terms of the GNU Lesser General Public License as published  * 
* by the Free Software Foundation; either version 3 of the License, or (at  * 
* your option) any later version.                                           * 
*                                                                           *
* This program is distributed in the hope that it will be useful, but       * 
* WITHOUT ANY WARRANTY; without even the implied warranty of                * 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  * 
* General Public License for more details.                                  * 
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  * 
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/ 

/* This code was generated automatically. Do not edit unless you exactly   */
/* know what you are doing!                                                */

#include "implicit_point.h"

#pragma intrinsic(fabs)

// Uncomment the following to activate overflow/underflow checks
#define CHECK_FOR_XYZERFLOWS

inline int dotProductSign2D_filtered(double px, double py, double rx, double ry, double qx, double qy)
{
   const double lx = px - qx;
   const double ly = py - qy;
   const double gx = rx - qx;
   const double gy = ry - qy;
   const double dx = lx * gx;
   const double dy = ly * gy;
   const double d = dx + dy;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(lx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ly)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(gx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(gy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= 8.881784197001253e-16;
   if (d > epsilon) return IP_Sign::POSITIVE;
   if (-d > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int dotProductSign2D_interval(interval_number px, interval_number py, interval_number rx, interval_number ry, interval_number qx, interval_number qy)
{
   setFPUModeToRoundUP();
   const interval_number lx(px - qx);
   const interval_number ly(py - qy);
   const interval_number gx(rx - qx);
   const interval_number gy(ry - qy);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number d(dx + dy);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign2D_bigfloat(bigfloat px, bigfloat py, bigfloat rx, bigfloat ry, bigfloat qx, bigfloat qy)
{
   const bigfloat lx(px - qx);
   const bigfloat ly(py - qy);
   const bigfloat gx(rx - qx);
   const bigfloat gy(ry - qy);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat d(dx + dy);
   return sgn(d);
}

inline int dotProductSign2D_exact(double px, double py, double rx, double ry, double qx, double qy)
{
   double lx[2];
   expansionObject::two_Diff(px, qx, lx);
   double ly[2];
   expansionObject::two_Diff(py, qy, ly);
   double gx[2];
   expansionObject::two_Diff(rx, qx, gx);
   double gy[2];
   expansionObject::two_Diff(ry, qy, gy);
   double dx[8];
   int dx_len = expansionObject::Gen_Product(2, lx, 2, gx, dx);
   double dy[8];
   int dy_len = expansionObject::Gen_Product(2, ly, 2, gy, dy);
   double d[16];
   int d_len = expansionObject::Gen_Sum(dx_len, dx, dy_len, dy, d);

   double return_value = d[d_len - 1];

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int dotProductSign2D(double px, double py, double rx, double ry, double qx, double qy)
{
   int ret;
   ret = dotProductSign2D_filtered(px, py, rx, ry, qx, qy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = dotProductSign2D_interval(px, py, rx, ry, qx, qy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign2D_exact(px, py, rx, ry, qx, qy);
}

inline int dotProductSign3D_filtered(double px, double py, double pz, double rx, double ry, double rz, double qx, double qy, double qz)
{
   const double lx = px - qx;
   const double ly = py - qy;
   const double lz = pz - qz;
   const double gx = rx - qx;
   const double gy = ry - qy;
   const double gz = rz - qz;
   const double dx = lx * gx;
   const double dy = ly * gy;
   const double dz = lz * gz;
   const double d1 = dx + dy;
   const double d = d1 + dz;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(lx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ly)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(lz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(gx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(gy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(gz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= 1.443289932012704e-15;
   if (d > epsilon) return IP_Sign::POSITIVE;
   if (-d > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int dotProductSign3D_interval(interval_number px, interval_number py, interval_number pz, interval_number rx, interval_number ry, interval_number rz, interval_number qx, interval_number qy, interval_number qz)
{
   setFPUModeToRoundUP();
   const interval_number lx(px - qx);
   const interval_number ly(py - qy);
   const interval_number lz(pz - qz);
   const interval_number gx(rx - qx);
   const interval_number gy(ry - qy);
   const interval_number gz(rz - qz);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number dz(lz * gz);
   const interval_number d1(dx + dy);
   const interval_number d(d1 + dz);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign3D_bigfloat(bigfloat px, bigfloat py, bigfloat pz, bigfloat rx, bigfloat ry, bigfloat rz, bigfloat qx, bigfloat qy, bigfloat qz)
{
   const bigfloat lx(px - qx);
   const bigfloat ly(py - qy);
   const bigfloat lz(pz - qz);
   const bigfloat gx(rx - qx);
   const bigfloat gy(ry - qy);
   const bigfloat gz(rz - qz);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat dz(lz * gz);
   const bigfloat d1(dx + dy);
   const bigfloat d(d1 + dz);
   return sgn(d);
}

inline int dotProductSign3D_exact(double px, double py, double pz, double rx, double ry, double rz, double qx, double qy, double qz)
{
   double lx[2];
   expansionObject::two_Diff(px, qx, lx);
   double ly[2];
   expansionObject::two_Diff(py, qy, ly);
   double lz[2];
   expansionObject::two_Diff(pz, qz, lz);
   double gx[2];
   expansionObject::two_Diff(rx, qx, gx);
   double gy[2];
   expansionObject::two_Diff(ry, qy, gy);
   double gz[2];
   expansionObject::two_Diff(rz, qz, gz);
   double dx[8];
   int dx_len = expansionObject::Gen_Product(2, lx, 2, gx, dx);
   double dy[8];
   int dy_len = expansionObject::Gen_Product(2, ly, 2, gy, dy);
   double dz[8];
   int dz_len = expansionObject::Gen_Product(2, lz, 2, gz, dz);
   double d1[16];
   int d1_len = expansionObject::Gen_Sum(dx_len, dx, dy_len, dy, d1);
   double d[24];
   int d_len = expansionObject::Gen_Sum(d1_len, d1, dz_len, dz, d);

   double return_value = d[d_len - 1];

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int dotProductSign3D(double px, double py, double pz, double rx, double ry, double rz, double qx, double qy, double qz)
{
   int ret;
   ret = dotProductSign3D_filtered(px, py, pz, rx, ry, rz, qx, qy, qz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = dotProductSign3D_interval(px, py, pz, rx, ry, rz, qx, qy, qz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign3D_exact(px, py, pz, rx, ry, rz, qx, qy, qz);
}

inline int incircle_filtered(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   const double adx = pax - pdx;
   const double ady = pay - pdy;
   const double bdx = pbx - pdx;
   const double bdy = pby - pdy;
   const double cdx = pcx - pdx;
   const double cdy = pcy - pdy;
   const double abdeta = adx * bdy;
   const double abdetb = bdx * ady;
   const double abdet = abdeta - abdetb;
   const double bcdeta = bdx * cdy;
   const double bcdetb = cdx * bdy;
   const double bcdet = bcdeta - bcdetb;
   const double cadeta = cdx * ady;
   const double cadetb = adx * cdy;
   const double cadet = cadeta - cadetb;
   const double alifta = adx * adx;
   const double aliftb = ady * ady;
   const double alift = alifta + aliftb;
   const double blifta = bdx * bdx;
   const double bliftb = bdy * bdy;
   const double blift = blifta + bliftb;
   const double clifta = cdx * cdx;
   const double cliftb = cdy * cdy;
   const double clift = clifta + cliftb;
   const double la = alift * bcdet;
   const double lb = blift * cadet;
   const double lc = clift * abdet;
   const double lab = la + lb;
   const double L = lab + lc;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(adx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ady)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.376676550535195e-14;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int incircle_interval(interval_number pax, interval_number pay, interval_number pbx, interval_number pby, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   setFPUModeToRoundUP();
   const interval_number adx(pax - pdx);
   const interval_number ady(pay - pdy);
   const interval_number bdx(pbx - pdx);
   const interval_number bdy(pby - pdy);
   const interval_number cdx(pcx - pdx);
   const interval_number cdy(pcy - pdy);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number alift(alifta + aliftb);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number blift(blifta + bliftb);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number clift(clifta + cliftb);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab(la + lb);
   const interval_number L(lab + lc);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incircle_bigfloat(bigfloat pax, bigfloat pay, bigfloat pbx, bigfloat pby, bigfloat pcx, bigfloat pcy, bigfloat pdx, bigfloat pdy)
{
   const bigfloat adx(pax - pdx);
   const bigfloat ady(pay - pdy);
   const bigfloat bdx(pbx - pdx);
   const bigfloat bdy(pby - pdy);
   const bigfloat cdx(pcx - pdx);
   const bigfloat cdy(pcy - pdy);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat alift(alifta + aliftb);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat blift(blifta + bliftb);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat clift(clifta + cliftb);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab(la + lb);
   const bigfloat L(lab + lc);
   return sgn(L);
}

inline int incircle_exact(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   double adx[2];
   expansionObject::two_Diff(pax, pdx, adx);
   double ady[2];
   expansionObject::two_Diff(pay, pdy, ady);
   double bdx[2];
   expansionObject::two_Diff(pbx, pdx, bdx);
   double bdy[2];
   expansionObject::two_Diff(pby, pdy, bdy);
   double cdx[2];
   expansionObject::two_Diff(pcx, pdx, cdx);
   double cdy[2];
   expansionObject::two_Diff(pcy, pdy, cdy);
   double abdeta[8];
   int abdeta_len = expansionObject::Gen_Product(2, adx, 2, bdy, abdeta);
   double abdetb[8];
   int abdetb_len = expansionObject::Gen_Product(2, bdx, 2, ady, abdetb);
   double abdet[16];
   int abdet_len = expansionObject::Gen_Diff(abdeta_len, abdeta, abdetb_len, abdetb, abdet);
   double bcdeta[8];
   int bcdeta_len = expansionObject::Gen_Product(2, bdx, 2, cdy, bcdeta);
   double bcdetb[8];
   int bcdetb_len = expansionObject::Gen_Product(2, cdx, 2, bdy, bcdetb);
   double bcdet[16];
   int bcdet_len = expansionObject::Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
   double cadeta[8];
   int cadeta_len = expansionObject::Gen_Product(2, cdx, 2, ady, cadeta);
   double cadetb[8];
   int cadetb_len = expansionObject::Gen_Product(2, adx, 2, cdy, cadetb);
   double cadet[16];
   int cadet_len = expansionObject::Gen_Diff(cadeta_len, cadeta, cadetb_len, cadetb, cadet);
   double alifta[8];
   int alifta_len = expansionObject::Gen_Product(2, adx, 2, adx, alifta);
   double aliftb[8];
   int aliftb_len = expansionObject::Gen_Product(2, ady, 2, ady, aliftb);
   double alift[16];
   int alift_len = expansionObject::Gen_Sum(alifta_len, alifta, aliftb_len, aliftb, alift);
   double blifta[8];
   int blifta_len = expansionObject::Gen_Product(2, bdx, 2, bdx, blifta);
   double bliftb[8];
   int bliftb_len = expansionObject::Gen_Product(2, bdy, 2, bdy, bliftb);
   double blift[16];
   int blift_len = expansionObject::Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
   double clifta[8];
   int clifta_len = expansionObject::Gen_Product(2, cdx, 2, cdx, clifta);
   double cliftb[8];
   int cliftb_len = expansionObject::Gen_Product(2, cdy, 2, cdy, cliftb);
   double clift[16];
   int clift_len = expansionObject::Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
   double la_p[128], *la = la_p;
   int la_len = expansionObject::Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 128);
   double lb_p[128], *lb = lb_p;
   int lb_len = expansionObject::Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 128);
   double lc_p[128], *lc = lc_p;
   int lc_len = expansionObject::Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 128);
   double lab_p[128], *lab = lab_p;
   int lab_len = expansionObject::Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 128);
   double L_p[128], *L = L_p;
   int L_len = expansionObject::Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 128);

   double return_value = L[L_len - 1];
   if (L_p != L) FreeDoubles(L);
   if (lab_p != lab) FreeDoubles(lab);
   if (lc_p != lc) FreeDoubles(lc);
   if (lb_p != lb) FreeDoubles(lb);
   if (la_p != la) FreeDoubles(la);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int incircle(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_filtered(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_interval(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_exact(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
}

inline int inGabrielSphere_filtered(double qx, double qy, double qz, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   const double bax = bx - ax;
   const double bay = by - ay;
   const double baz = bz - az;
   const double cax = cx - ax;
   const double cay = cy - ay;
   const double caz = cz - az;
   const double qax = qx - ax;
   const double qay = qy - ay;
   const double qaz = qz - az;
   const double cx1 = bay * caz;
   const double cx2 = baz * cay;
   const double crossbcx = cx1 - cx2;
   const double cy1 = baz * cax;
   const double cy2 = bax * caz;
   const double crossbcy = cy1 - cy2;
   const double cz1 = bax * cay;
   const double cz2 = bay * cax;
   const double crossbcz = cz1 - cz2;
   const double ba2x = bax * bax;
   const double ba2y = bay * bay;
   const double ba2z = baz * baz;
   const double ba2t = ba2x + ba2y;
   const double ba2 = ba2t + ba2z;
   const double ca2x = cax * cax;
   const double ca2y = cay * cay;
   const double ca2z = caz * caz;
   const double ca2t = ca2x + ca2y;
   const double ca2 = ca2t + ca2z;
   const double calx = cax * ba2;
   const double caly = cay * ba2;
   const double calz = caz * ba2;
   const double balx = bax * ca2;
   const double baly = bay * ca2;
   const double balz = baz * ca2;
   const double abcx = calx - balx;
   const double abcy = caly - baly;
   const double abcz = calz - balz;
   const double kx1 = abcy * crossbcz;
   const double kx2 = abcz * crossbcy;
   const double ccax = kx1 - kx2;
   const double ky1 = abcz * crossbcx;
   const double ky2 = abcx * crossbcz;
   const double ccay = ky1 - ky2;
   const double kz1 = abcx * crossbcy;
   const double kz2 = abcy * crossbcx;
   const double ccaz = kz1 - kz2;
   const double cr2x = crossbcx * crossbcx;
   const double cr2y = crossbcy * crossbcy;
   const double cr2z = crossbcz * crossbcz;
   const double cr2t = cr2x + cr2y;
   const double c2 = cr2t + cr2z;
   const double c22 = 2 * c2;
   const double qa1x = qax * c22;
   const double qa1y = qay * c22;
   const double qa1z = qaz * c22;
   const double qa2x = qa1x - ccax;
   const double qa2y = qa1y - ccay;
   const double qa2z = qa1z - ccaz;
   const double r1x = qa2x * qa2x;
   const double r1y = qa2y * qa2y;
   const double r1z = qa2z * qa2z;
   const double r1t = r1x + r1y;
   const double r1 = r1t + r1z;
   const double r2x = ccax * ccax;
   const double r2y = ccay * ccay;
   const double r2z = ccaz * ccaz;
   const double r2t = r2x + r2y;
   const double r2 = r2t + r2z;
   const double ret = r1 - r2;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(bax)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bay)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(baz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cax)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cay)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(caz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qax)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qay)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qaz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.297895207448169e-11;
   if (ret > epsilon) return IP_Sign::POSITIVE;
   if (-ret > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int inGabrielSphere_interval(interval_number qx, interval_number qy, interval_number qz, interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   setFPUModeToRoundUP();
   const interval_number bax(bx - ax);
   const interval_number bay(by - ay);
   const interval_number baz(bz - az);
   const interval_number cax(cx - ax);
   const interval_number cay(cy - ay);
   const interval_number caz(cz - az);
   const interval_number qax(qx - ax);
   const interval_number qay(qy - ay);
   const interval_number qaz(qz - az);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_bigfloat(bigfloat qx, bigfloat qy, bigfloat qz, bigfloat ax, bigfloat ay, bigfloat az, bigfloat bx, bigfloat by, bigfloat bz, bigfloat cx, bigfloat cy, bigfloat cz)
{
   const bigfloat bax(bx - ax);
   const bigfloat bay(by - ay);
   const bigfloat baz(bz - az);
   const bigfloat cax(cx - ax);
   const bigfloat cay(cy - ay);
   const bigfloat caz(cz - az);
   const bigfloat qax(qx - ax);
   const bigfloat qay(qy - ay);
   const bigfloat qaz(qz - az);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere(double qx, double qy, double qz, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = inGabrielSphere_filtered(qx, qy, qz, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = inGabrielSphere_interval(qx, qy, qz, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_bigfloat(qx, qy, qz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline int inSphere_filtered(double pax, double pay, double paz, double pbx, double pby, double pbz, double pcx, double pcy, double pcz, double pdx, double pdy, double pdz, double pex, double pey, double pez)
{
   const double aex = pax - pex;
   const double aey = pay - pey;
   const double aez = paz - pez;
   const double bex = pbx - pex;
   const double bey = pby - pey;
   const double bez = pbz - pez;
   const double cex = pcx - pex;
   const double cey = pcy - pey;
   const double cez = pcz - pez;
   const double dex = pdx - pex;
   const double dey = pdy - pey;
   const double dez = pdz - pez;
   const double aexbey = aex * bey;
   const double bexaey = bex * aey;
   const double ab = aexbey - bexaey;
   const double bexcey = bex * cey;
   const double cexbey = cex * bey;
   const double bc = bexcey - cexbey;
   const double cexdey = cex * dey;
   const double dexcey = dex * cey;
   const double cd = cexdey - dexcey;
   const double dexaey = dex * aey;
   const double aexdey = aex * dey;
   const double da = dexaey - aexdey;
   const double aexcey = aex * cey;
   const double cexaey = cex * aey;
   const double ac = aexcey - cexaey;
   const double bexdey = bex * dey;
   const double dexbey = dex * bey;
   const double bd = bexdey - dexbey;
   const double abc1 = aez * bc;
   const double abc2 = bez * ac;
   const double abc3 = cez * ab;
   const double abc4 = abc1 + abc3;
   const double abc = abc4 - abc2;
   const double bcd1 = bez * cd;
   const double bcd2 = cez * bd;
   const double bcd3 = dez * bc;
   const double bcd4 = bcd1 + bcd3;
   const double bcd = bcd4 - bcd2;
   const double cda1 = cez * da;
   const double cda2 = dez * ac;
   const double cda3 = aez * cd;
   const double cda4 = cda1 + cda3;
   const double cda = cda4 + cda2;
   const double dab1 = dez * ab;
   const double dab2 = aez * bd;
   const double dab3 = bez * da;
   const double dab4 = dab1 + dab3;
   const double dab = dab4 + dab2;
   const double al1 = aex * aex;
   const double al2 = aey * aey;
   const double al3 = aez * aez;
   const double al4 = al1 + al2;
   const double alift = al4 + al3;
   const double bl1 = bex * bex;
   const double bl2 = bey * bey;
   const double bl3 = bez * bez;
   const double bl4 = bl1 + bl2;
   const double blift = bl4 + bl3;
   const double cl1 = cex * cex;
   const double cl2 = cey * cey;
   const double cl3 = cez * cez;
   const double cl4 = cl1 + cl2;
   const double clift = cl4 + cl3;
   const double dl1 = dex * dex;
   const double dl2 = dey * dey;
   const double dl3 = dez * dez;
   const double dl4 = dl1 + dl2;
   const double dlift = dl4 + dl3;
   const double ds1 = dlift * abc;
   const double ds2 = clift * dab;
   const double dl = ds2 - ds1;
   const double dr1 = blift * cda;
   const double dr2 = alift * bcd;
   const double dr = dr2 - dr1;
   const double det = dl + dr;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(aex)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(aey)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(aez)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bex)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bey)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bez)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cex)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cey)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cez)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dex)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dey)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dez)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= 1.145750161413163e-13;
   if (det > epsilon) return IP_Sign::POSITIVE;
   if (-det > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int inSphere_interval(interval_number pax, interval_number pay, interval_number paz, interval_number pbx, interval_number pby, interval_number pbz, interval_number pcx, interval_number pcy, interval_number pcz, interval_number pdx, interval_number pdy, interval_number pdz, interval_number pex, interval_number pey, interval_number pez)
{
   setFPUModeToRoundUP();
   const interval_number aex(pax - pex);
   const interval_number aey(pay - pey);
   const interval_number aez(paz - pez);
   const interval_number bex(pbx - pex);
   const interval_number bey(pby - pey);
   const interval_number bez(pbz - pez);
   const interval_number cex(pcx - pex);
   const interval_number cey(pcy - pey);
   const interval_number cez(pcz - pez);
   const interval_number dex(pdx - pex);
   const interval_number dey(pdy - pey);
   const interval_number dez(pdz - pez);
   const interval_number aexbey(aex * bey);
   const interval_number bexaey(bex * aey);
   const interval_number ab(aexbey - bexaey);
   const interval_number bexcey(bex * cey);
   const interval_number cexbey(cex * bey);
   const interval_number bc(bexcey - cexbey);
   const interval_number cexdey(cex * dey);
   const interval_number dexcey(dex * cey);
   const interval_number cd(cexdey - dexcey);
   const interval_number dexaey(dex * aey);
   const interval_number aexdey(aex * dey);
   const interval_number da(dexaey - aexdey);
   const interval_number aexcey(aex * cey);
   const interval_number cexaey(cex * aey);
   const interval_number ac(aexcey - cexaey);
   const interval_number bexdey(bex * dey);
   const interval_number dexbey(dex * bey);
   const interval_number bd(bexdey - dexbey);
   const interval_number abc1(aez * bc);
   const interval_number abc2(bez * ac);
   const interval_number abc3(cez * ab);
   const interval_number abc4(abc1 + abc3);
   const interval_number abc(abc4 - abc2);
   const interval_number bcd1(bez * cd);
   const interval_number bcd2(cez * bd);
   const interval_number bcd3(dez * bc);
   const interval_number bcd4(bcd1 + bcd3);
   const interval_number bcd(bcd4 - bcd2);
   const interval_number cda1(cez * da);
   const interval_number cda2(dez * ac);
   const interval_number cda3(aez * cd);
   const interval_number cda4(cda1 + cda3);
   const interval_number cda(cda4 + cda2);
   const interval_number dab1(dez * ab);
   const interval_number dab2(aez * bd);
   const interval_number dab3(bez * da);
   const interval_number dab4(dab1 + dab3);
   const interval_number dab(dab4 + dab2);
   const interval_number al1(aex * aex);
   const interval_number al2(aey * aey);
   const interval_number al3(aez * aez);
   const interval_number al4(al1 + al2);
   const interval_number alift(al4 + al3);
   const interval_number bl1(bex * bex);
   const interval_number bl2(bey * bey);
   const interval_number bl3(bez * bez);
   const interval_number bl4(bl1 + bl2);
   const interval_number blift(bl4 + bl3);
   const interval_number cl1(cex * cex);
   const interval_number cl2(cey * cey);
   const interval_number cl3(cez * cez);
   const interval_number cl4(cl1 + cl2);
   const interval_number clift(cl4 + cl3);
   const interval_number dl1(dex * dex);
   const interval_number dl2(dey * dey);
   const interval_number dl3(dez * dez);
   const interval_number dl4(dl1 + dl2);
   const interval_number dlift(dl4 + dl3);
   const interval_number ds1(dlift * abc);
   const interval_number ds2(clift * dab);
   const interval_number dl(ds2 - ds1);
   const interval_number dr1(blift * cda);
   const interval_number dr2(alift * bcd);
   const interval_number dr(dr2 - dr1);
   const interval_number det(dl + dr);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int inSphere_bigfloat(bigfloat pax, bigfloat pay, bigfloat paz, bigfloat pbx, bigfloat pby, bigfloat pbz, bigfloat pcx, bigfloat pcy, bigfloat pcz, bigfloat pdx, bigfloat pdy, bigfloat pdz, bigfloat pex, bigfloat pey, bigfloat pez)
{
   const bigfloat aex(pax - pex);
   const bigfloat aey(pay - pey);
   const bigfloat aez(paz - pez);
   const bigfloat bex(pbx - pex);
   const bigfloat bey(pby - pey);
   const bigfloat bez(pbz - pez);
   const bigfloat cex(pcx - pex);
   const bigfloat cey(pcy - pey);
   const bigfloat cez(pcz - pez);
   const bigfloat dex(pdx - pex);
   const bigfloat dey(pdy - pey);
   const bigfloat dez(pdz - pez);
   const bigfloat aexbey(aex * bey);
   const bigfloat bexaey(bex * aey);
   const bigfloat ab(aexbey - bexaey);
   const bigfloat bexcey(bex * cey);
   const bigfloat cexbey(cex * bey);
   const bigfloat bc(bexcey - cexbey);
   const bigfloat cexdey(cex * dey);
   const bigfloat dexcey(dex * cey);
   const bigfloat cd(cexdey - dexcey);
   const bigfloat dexaey(dex * aey);
   const bigfloat aexdey(aex * dey);
   const bigfloat da(dexaey - aexdey);
   const bigfloat aexcey(aex * cey);
   const bigfloat cexaey(cex * aey);
   const bigfloat ac(aexcey - cexaey);
   const bigfloat bexdey(bex * dey);
   const bigfloat dexbey(dex * bey);
   const bigfloat bd(bexdey - dexbey);
   const bigfloat abc1(aez * bc);
   const bigfloat abc2(bez * ac);
   const bigfloat abc3(cez * ab);
   const bigfloat abc4(abc1 + abc3);
   const bigfloat abc(abc4 - abc2);
   const bigfloat bcd1(bez * cd);
   const bigfloat bcd2(cez * bd);
   const bigfloat bcd3(dez * bc);
   const bigfloat bcd4(bcd1 + bcd3);
   const bigfloat bcd(bcd4 - bcd2);
   const bigfloat cda1(cez * da);
   const bigfloat cda2(dez * ac);
   const bigfloat cda3(aez * cd);
   const bigfloat cda4(cda1 + cda3);
   const bigfloat cda(cda4 + cda2);
   const bigfloat dab1(dez * ab);
   const bigfloat dab2(aez * bd);
   const bigfloat dab3(bez * da);
   const bigfloat dab4(dab1 + dab3);
   const bigfloat dab(dab4 + dab2);
   const bigfloat al1(aex * aex);
   const bigfloat al2(aey * aey);
   const bigfloat al3(aez * aez);
   const bigfloat al4(al1 + al2);
   const bigfloat alift(al4 + al3);
   const bigfloat bl1(bex * bex);
   const bigfloat bl2(bey * bey);
   const bigfloat bl3(bez * bez);
   const bigfloat bl4(bl1 + bl2);
   const bigfloat blift(bl4 + bl3);
   const bigfloat cl1(cex * cex);
   const bigfloat cl2(cey * cey);
   const bigfloat cl3(cez * cez);
   const bigfloat cl4(cl1 + cl2);
   const bigfloat clift(cl4 + cl3);
   const bigfloat dl1(dex * dex);
   const bigfloat dl2(dey * dey);
   const bigfloat dl3(dez * dez);
   const bigfloat dl4(dl1 + dl2);
   const bigfloat dlift(dl4 + dl3);
   const bigfloat ds1(dlift * abc);
   const bigfloat ds2(clift * dab);
   const bigfloat dl(ds2 - ds1);
   const bigfloat dr1(blift * cda);
   const bigfloat dr2(alift * bcd);
   const bigfloat dr(dr2 - dr1);
   const bigfloat det(dl + dr);
   return sgn(det);
}

inline int inSphere(double pax, double pay, double paz, double pbx, double pby, double pbz, double pcx, double pcy, double pcz, double pdx, double pdy, double pdz, double pex, double pey, double pez)
{
   int ret;
   ret = inSphere_filtered(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = inSphere_interval(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_bigfloat(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}

inline int dotProductSign2D_EEI_interval(const genericPoint& q, interval_number px, interval_number py, interval_number rx, interval_number ry)
{
   interval_number lqx, lqy, dq;
   if (
   !q.getIntervalLambda(lqx, lqy, dq)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pxq(px * dq);
   const interval_number pyq(py * dq);
   const interval_number rxq(rx * dq);
   const interval_number ryq(ry * dq);
   const interval_number lx(pxq - lqx);
   const interval_number ly(pyq - lqy);
   const interval_number gx(rxq - lqx);
   const interval_number gy(ryq - lqy);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number d(dx + dy);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign2D_EEI_bigfloat(const genericPoint& q, bigfloat px, bigfloat py, bigfloat rx, bigfloat ry)
{
   bigfloat lqx, lqy, dq;
   q.getBigfloatLambda(lqx, lqy, dq);
   const bigfloat pxq(px * dq);
   const bigfloat pyq(py * dq);
   const bigfloat rxq(rx * dq);
   const bigfloat ryq(ry * dq);
   const bigfloat lx(pxq - lqx);
   const bigfloat ly(pyq - lqy);
   const bigfloat gx(rxq - lqx);
   const bigfloat gy(ryq - lqy);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat d(dx + dy);
   return sgn(d);
}

inline int dotProductSign2D_EEI_exact(const genericPoint& q, double px, double py, double rx, double ry)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double lqx_p[128], *lqx = lqx_p, lqy_p[128], *lqy = lqy_p, dq_p[128], *dq = dq_p;
 int lqx_len = 128, lqy_len = 128, dq_len = 128;
 q.getExactLambda(&lqx, lqx_len, &lqy, lqy_len, &dq, dq_len);
 if ((dq[dq_len - 1] != 0))
 {
   double pxq_p[128], *pxq = pxq_p;
   int pxq_len = expansionObject::Gen_Scale_With_PreAlloc(dq_len, dq, px, &pxq, 128);
   double pyq_p[128], *pyq = pyq_p;
   int pyq_len = expansionObject::Gen_Scale_With_PreAlloc(dq_len, dq, py, &pyq, 128);
   double rxq_p[128], *rxq = rxq_p;
   int rxq_len = expansionObject::Gen_Scale_With_PreAlloc(dq_len, dq, rx, &rxq, 128);
   double ryq_p[128], *ryq = ryq_p;
   int ryq_len = expansionObject::Gen_Scale_With_PreAlloc(dq_len, dq, ry, &ryq, 128);
   double lx_p[128], *lx = lx_p;
   int lx_len = expansionObject::Gen_Diff_With_PreAlloc(pxq_len, pxq, lqx_len, lqx, &lx, 128);
   double ly_p[128], *ly = ly_p;
   int ly_len = expansionObject::Gen_Diff_With_PreAlloc(pyq_len, pyq, lqy_len, lqy, &ly, 128);
   double gx_p[128], *gx = gx_p;
   int gx_len = expansionObject::Gen_Diff_With_PreAlloc(rxq_len, rxq, lqx_len, lqx, &gx, 128);
   double gy_p[128], *gy = gy_p;
   int gy_len = expansionObject::Gen_Diff_With_PreAlloc(ryq_len, ryq, lqy_len, lqy, &gy, 128);
   double dx_p[128], *dx = dx_p;
   int dx_len = expansionObject::Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 128);
   double dy_p[128], *dy = dy_p;
   int dy_len = expansionObject::Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 128);
   double d_p[128], *d = d_p;
   int d_len = expansionObject::Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 128);

   return_value = d[d_len - 1];
   if (d_p != d) FreeDoubles(d);
   if (dy_p != dy) FreeDoubles(dy);
   if (dx_p != dx) FreeDoubles(dx);
   if (gy_p != gy) FreeDoubles(gy);
   if (gx_p != gx) FreeDoubles(gx);
   if (ly_p != ly) FreeDoubles(ly);
   if (lx_p != lx) FreeDoubles(lx);
   if (ryq_p != ryq) FreeDoubles(ryq);
   if (rxq_p != rxq) FreeDoubles(rxq);
   if (pyq_p != pyq) FreeDoubles(pyq);
   if (pxq_p != pxq) FreeDoubles(pxq);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = dotProductSign2D_EEI_bigfloat(q, px, py, rx, ry);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int dotProductSign2D_EEI(const genericPoint& q, double px, double py, double rx, double ry)
{
   int ret;
   ret = dotProductSign2D_EEI_interval(q, px, py, rx, ry);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign2D_EEI_exact(q, px, py, rx, ry);
}

inline int dotProductSign2D_IEE_interval(const genericPoint& p, interval_number rx, interval_number ry, interval_number qx, interval_number qy)
{
   interval_number lpx, lpy, dp;
   if (
   !p.getIntervalLambda(lpx, lpy, dp)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number qxd(qx * dp);
   const interval_number qyd(qy * dp);
   const interval_number lx(lpx - qxd);
   const interval_number ly(lpy - qyd);
   const interval_number gx(rx - qx);
   const interval_number gy(ry - qy);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number d(dx + dy);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign2D_IEE_bigfloat(const genericPoint& p, bigfloat rx, bigfloat ry, bigfloat qx, bigfloat qy)
{
   bigfloat lpx, lpy, dp;
   p.getBigfloatLambda(lpx, lpy, dp);
   const bigfloat qxd(qx * dp);
   const bigfloat qyd(qy * dp);
   const bigfloat lx(lpx - qxd);
   const bigfloat ly(lpy - qyd);
   const bigfloat gx(rx - qx);
   const bigfloat gy(ry - qy);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat d(dx + dy);
   return sgn(d);
}

inline int dotProductSign2D_IEE_exact(const genericPoint& p, double rx, double ry, double qx, double qy)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double lpx_p[128], *lpx = lpx_p, lpy_p[128], *lpy = lpy_p, dp_p[128], *dp = dp_p;
 int lpx_len = 128, lpy_len = 128, dp_len = 128;
 p.getExactLambda(&lpx, lpx_len, &lpy, lpy_len, &dp, dp_len);
 if ((dp[dp_len - 1] != 0))
 {
   double qxd_p[128], *qxd = qxd_p;
   int qxd_len = expansionObject::Gen_Scale_With_PreAlloc(dp_len, dp, qx, &qxd, 128);
   double qyd_p[128], *qyd = qyd_p;
   int qyd_len = expansionObject::Gen_Scale_With_PreAlloc(dp_len, dp, qy, &qyd, 128);
   double lx_p[128], *lx = lx_p;
   int lx_len = expansionObject::Gen_Diff_With_PreAlloc(lpx_len, lpx, qxd_len, qxd, &lx, 128);
   double ly_p[128], *ly = ly_p;
   int ly_len = expansionObject::Gen_Diff_With_PreAlloc(lpy_len, lpy, qyd_len, qyd, &ly, 128);
   double gx[2];
   expansionObject::two_Diff(rx, qx, gx);
   double gy[2];
   expansionObject::two_Diff(ry, qy, gy);
   double dx_p[128], *dx = dx_p;
   int dx_len = expansionObject::Gen_Product_With_PreAlloc(lx_len, lx, 2, gx, &dx, 128);
   double dy_p[128], *dy = dy_p;
   int dy_len = expansionObject::Gen_Product_With_PreAlloc(ly_len, ly, 2, gy, &dy, 128);
   double d_p[128], *d = d_p;
   int d_len = expansionObject::Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 128);

   return_value = d[d_len - 1];
   if (d_p != d) FreeDoubles(d);
   if (dy_p != dy) FreeDoubles(dy);
   if (dx_p != dx) FreeDoubles(dx);
   if (ly_p != ly) FreeDoubles(ly);
   if (lx_p != lx) FreeDoubles(lx);
   if (qyd_p != qyd) FreeDoubles(qyd);
   if (qxd_p != qxd) FreeDoubles(qxd);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = dotProductSign2D_IEE_bigfloat(p, rx, ry, qx, qy);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int dotProductSign2D_IEE(const genericPoint& p, double rx, double ry, double qx, double qy)
{
   int ret;
   ret = dotProductSign2D_IEE_interval(p, rx, ry, qx, qy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign2D_IEE_exact(p, rx, ry, qx, qy);
}

inline int dotProductSign2D_IEI_interval(const genericPoint& p, const genericPoint& q, interval_number rx, interval_number ry)
{
   interval_number lpx, lpy, dp, lqx, lqy, dq;
   if (
   !p.getIntervalLambda(lpx, lpy, dp)
   || !q.getIntervalLambda(lqx, lqy, dq)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number dqp(dq * dp);
   const interval_number pxq(lpx * dqp);
   const interval_number pyq(lpy * dqp);
   const interval_number rxq(rx * dq);
   const interval_number ryq(ry * dq);
   const interval_number lqxd(lqx * dp);
   const interval_number lqyd(lqy * dp);
   const interval_number lx(pxq - lqxd);
   const interval_number ly(pyq - lqyd);
   const interval_number gx(rxq - lqx);
   const interval_number gy(ryq - lqy);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number d(dx + dy);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign2D_IEI_bigfloat(const genericPoint& p, const genericPoint& q, bigfloat rx, bigfloat ry)
{
   bigfloat lpx, lpy, dp, lqx, lqy, dq;
   p.getBigfloatLambda(lpx, lpy, dp);
   q.getBigfloatLambda(lqx, lqy, dq);
   const bigfloat dqp(dq * dp);
   const bigfloat pxq(lpx * dqp);
   const bigfloat pyq(lpy * dqp);
   const bigfloat rxq(rx * dq);
   const bigfloat ryq(ry * dq);
   const bigfloat lqxd(lqx * dp);
   const bigfloat lqyd(lqy * dp);
   const bigfloat lx(pxq - lqxd);
   const bigfloat ly(pyq - lqyd);
   const bigfloat gx(rxq - lqx);
   const bigfloat gy(ryq - lqy);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat d(dx + dy);
   return sgn(d);
}

inline int dotProductSign2D_IEI(const genericPoint& p, const genericPoint& q, double rx, double ry)
{
   int ret;
   ret = dotProductSign2D_IEI_interval(p, q, rx, ry);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign2D_IEI_bigfloat(p, q, rx, ry);
}

inline int dotProductSign2D_IIE_interval(const genericPoint& p, const genericPoint& r, interval_number qx, interval_number qy)
{
   interval_number lpx, lpy, dp, lrx, lry, dr;
   if (
   !p.getIntervalLambda(lpx, lpy, dp)
   || !r.getIntervalLambda(lrx, lry, dr)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number qxd(qx * dp);
   const interval_number qyd(qy * dp);
   const interval_number lx(lpx - qxd);
   const interval_number ly(lpy - qyd);
   const interval_number qxr(qx * dr);
   const interval_number qyr(qy * dr);
   const interval_number gx(lrx - qxr);
   const interval_number gy(lry - qyr);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number d(dx + dy);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign2D_IIE_bigfloat(const genericPoint& p, const genericPoint& r, bigfloat qx, bigfloat qy)
{
   bigfloat lpx, lpy, dp, lrx, lry, dr;
   p.getBigfloatLambda(lpx, lpy, dp);
   r.getBigfloatLambda(lrx, lry, dr);
   const bigfloat qxd(qx * dp);
   const bigfloat qyd(qy * dp);
   const bigfloat lx(lpx - qxd);
   const bigfloat ly(lpy - qyd);
   const bigfloat qxr(qx * dr);
   const bigfloat qyr(qy * dr);
   const bigfloat gx(lrx - qxr);
   const bigfloat gy(lry - qyr);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat d(dx + dy);
   return sgn(d);
}

inline int dotProductSign2D_IIE(const genericPoint& p, const genericPoint& r, double qx, double qy)
{
   int ret;
   ret = dotProductSign2D_IIE_interval(p, r, qx, qy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign2D_IIE_bigfloat(p, r, qx, qy);
}

inline int dotProductSign2D_III_interval(const genericPoint& p, const genericPoint& r, const genericPoint& q)
{
   interval_number lpx, lpy, dp, lrx, lry, dr, lqx, lqy, dq;
   if (
   !p.getIntervalLambda(lpx, lpy, dp)
   || !r.getIntervalLambda(lrx, lry, dr)
   || !q.getIntervalLambda(lqx, lqy, dq)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number qxd(lqx * dp);
   const interval_number qyd(lqy * dp);
   const interval_number lpxq(lpx * dq);
   const interval_number lpyq(lpy * dq);
   const interval_number lx(lpxq - qxd);
   const interval_number ly(lpyq - qyd);
   const interval_number qxr(lqx * dr);
   const interval_number qyr(lqy * dr);
   const interval_number lrxq(lrx * dq);
   const interval_number lryq(lry * dq);
   const interval_number gx(lrxq - qxr);
   const interval_number gy(lryq - qyr);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number d(dx + dy);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign2D_III_bigfloat(const genericPoint& p, const genericPoint& r, const genericPoint& q)
{
   bigfloat lpx, lpy, dp, lrx, lry, dr, lqx, lqy, dq;
   p.getBigfloatLambda(lpx, lpy, dp);
   r.getBigfloatLambda(lrx, lry, dr);
   q.getBigfloatLambda(lqx, lqy, dq);
   const bigfloat qxd(lqx * dp);
   const bigfloat qyd(lqy * dp);
   const bigfloat lpxq(lpx * dq);
   const bigfloat lpyq(lpy * dq);
   const bigfloat lx(lpxq - qxd);
   const bigfloat ly(lpyq - qyd);
   const bigfloat qxr(lqx * dr);
   const bigfloat qyr(lqy * dr);
   const bigfloat lrxq(lrx * dq);
   const bigfloat lryq(lry * dq);
   const bigfloat gx(lrxq - qxr);
   const bigfloat gy(lryq - qyr);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat d(dx + dy);
   return sgn(d);
}

inline int dotProductSign2D_III(const genericPoint& p, const genericPoint& r, const genericPoint& q)
{
   int ret;
   ret = dotProductSign2D_III_interval(p, r, q);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign2D_III_bigfloat(p, r, q);
}

inline int dotProductSign3D_EEI_interval(const genericPoint& q, interval_number px, interval_number py, interval_number pz, interval_number rx, interval_number ry, interval_number rz)
{
   interval_number lqx, lqy, lqz, dq;
   if (
   !q.getIntervalLambda(lqx, lqy, lqz, dq)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pxq(px * dq);
   const interval_number pyq(py * dq);
   const interval_number pzq(pz * dq);
   const interval_number rxq(rx * dq);
   const interval_number ryq(ry * dq);
   const interval_number rzq(rz * dq);
   const interval_number lx(pxq - lqx);
   const interval_number ly(pyq - lqy);
   const interval_number lz(pzq - lqz);
   const interval_number gx(rxq - lqx);
   const interval_number gy(ryq - lqy);
   const interval_number gz(rzq - lqz);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number dz(lz * gz);
   const interval_number d1(dx + dy);
   const interval_number d(d1 + dz);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign3D_EEI_bigfloat(const genericPoint& q, bigfloat px, bigfloat py, bigfloat pz, bigfloat rx, bigfloat ry, bigfloat rz)
{
   bigfloat lqx, lqy, lqz, dq;
   q.getBigfloatLambda(lqx, lqy, lqz, dq);
   const bigfloat pxq(px * dq);
   const bigfloat pyq(py * dq);
   const bigfloat pzq(pz * dq);
   const bigfloat rxq(rx * dq);
   const bigfloat ryq(ry * dq);
   const bigfloat rzq(rz * dq);
   const bigfloat lx(pxq - lqx);
   const bigfloat ly(pyq - lqy);
   const bigfloat lz(pzq - lqz);
   const bigfloat gx(rxq - lqx);
   const bigfloat gy(ryq - lqy);
   const bigfloat gz(rzq - lqz);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat dz(lz * gz);
   const bigfloat d1(dx + dy);
   const bigfloat d(d1 + dz);
   return sgn(d);
}

inline int dotProductSign3D_EEI(const genericPoint& q, double px, double py, double pz, double rx, double ry, double rz)
{
   int ret;
   ret = dotProductSign3D_EEI_interval(q, px, py, pz, rx, ry, rz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign3D_EEI_bigfloat(q, px, py, pz, rx, ry, rz);
}

inline int dotProductSign3D_IEE_interval(const genericPoint& p, interval_number rx, interval_number ry, interval_number rz, interval_number qx, interval_number qy, interval_number qz)
{
   interval_number lpx, lpy, lpz, dp;
   if (
   !p.getIntervalLambda(lpx, lpy, lpz, dp)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number qxd(qx * dp);
   const interval_number qyd(qy * dp);
   const interval_number qzd(qz * dp);
   const interval_number lx(lpx - qxd);
   const interval_number ly(lpy - qyd);
   const interval_number lz(lpz - qzd);
   const interval_number gx(rx - qx);
   const interval_number gy(ry - qy);
   const interval_number gz(rz - qz);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number dz(lz * gz);
   const interval_number d1(dx + dy);
   const interval_number d(d1 + dz);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign3D_IEE_bigfloat(const genericPoint& p, bigfloat rx, bigfloat ry, bigfloat rz, bigfloat qx, bigfloat qy, bigfloat qz)
{
   bigfloat lpx, lpy, lpz, dp;
   p.getBigfloatLambda(lpx, lpy, lpz, dp);
   const bigfloat qxd(qx * dp);
   const bigfloat qyd(qy * dp);
   const bigfloat qzd(qz * dp);
   const bigfloat lx(lpx - qxd);
   const bigfloat ly(lpy - qyd);
   const bigfloat lz(lpz - qzd);
   const bigfloat gx(rx - qx);
   const bigfloat gy(ry - qy);
   const bigfloat gz(rz - qz);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat dz(lz * gz);
   const bigfloat d1(dx + dy);
   const bigfloat d(d1 + dz);
   return sgn(d);
}

inline int dotProductSign3D_IEE_exact(const genericPoint& p, double rx, double ry, double rz, double qx, double qy, double qz)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double lpx_p[128], *lpx = lpx_p, lpy_p[128], *lpy = lpy_p, lpz_p[128], *lpz = lpz_p, dp_p[128], *dp = dp_p;
 int lpx_len = 128, lpy_len = 128, lpz_len = 128, dp_len = 128;
 p.getExactLambda(&lpx, lpx_len, &lpy, lpy_len, &lpz, lpz_len, &dp, dp_len);
 if ((dp[dp_len - 1] != 0))
 {
   double qxd_p[128], *qxd = qxd_p;
   int qxd_len = expansionObject::Gen_Scale_With_PreAlloc(dp_len, dp, qx, &qxd, 128);
   double qyd_p[128], *qyd = qyd_p;
   int qyd_len = expansionObject::Gen_Scale_With_PreAlloc(dp_len, dp, qy, &qyd, 128);
   double qzd_p[128], *qzd = qzd_p;
   int qzd_len = expansionObject::Gen_Scale_With_PreAlloc(dp_len, dp, qz, &qzd, 128);
   double lx_p[128], *lx = lx_p;
   int lx_len = expansionObject::Gen_Diff_With_PreAlloc(lpx_len, lpx, qxd_len, qxd, &lx, 128);
   double ly_p[128], *ly = ly_p;
   int ly_len = expansionObject::Gen_Diff_With_PreAlloc(lpy_len, lpy, qyd_len, qyd, &ly, 128);
   double lz_p[128], *lz = lz_p;
   int lz_len = expansionObject::Gen_Diff_With_PreAlloc(lpz_len, lpz, qzd_len, qzd, &lz, 128);
   double gx[2];
   expansionObject::two_Diff(rx, qx, gx);
   double gy[2];
   expansionObject::two_Diff(ry, qy, gy);
   double gz[2];
   expansionObject::two_Diff(rz, qz, gz);
   double dx_p[128], *dx = dx_p;
   int dx_len = expansionObject::Gen_Product_With_PreAlloc(lx_len, lx, 2, gx, &dx, 128);
   double dy_p[128], *dy = dy_p;
   int dy_len = expansionObject::Gen_Product_With_PreAlloc(ly_len, ly, 2, gy, &dy, 128);
   double dz_p[128], *dz = dz_p;
   int dz_len = expansionObject::Gen_Product_With_PreAlloc(lz_len, lz, 2, gz, &dz, 128);
   double d1_p[128], *d1 = d1_p;
   int d1_len = expansionObject::Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d1, 128);
   double d_p[128], *d = d_p;
   int d_len = expansionObject::Gen_Sum_With_PreAlloc(d1_len, d1, dz_len, dz, &d, 128);

   return_value = d[d_len - 1];
   if (d_p != d) FreeDoubles(d);
   if (d1_p != d1) FreeDoubles(d1);
   if (dz_p != dz) FreeDoubles(dz);
   if (dy_p != dy) FreeDoubles(dy);
   if (dx_p != dx) FreeDoubles(dx);
   if (lz_p != lz) FreeDoubles(lz);
   if (ly_p != ly) FreeDoubles(ly);
   if (lx_p != lx) FreeDoubles(lx);
   if (qzd_p != qzd) FreeDoubles(qzd);
   if (qyd_p != qyd) FreeDoubles(qyd);
   if (qxd_p != qxd) FreeDoubles(qxd);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = dotProductSign3D_IEE_bigfloat(p, rx, ry, rz, qx, qy, qz);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int dotProductSign3D_IEE(const genericPoint& p, double rx, double ry, double rz, double qx, double qy, double qz)
{
   int ret;
   ret = dotProductSign3D_IEE_interval(p, rx, ry, rz, qx, qy, qz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign3D_IEE_exact(p, rx, ry, rz, qx, qy, qz);
}

inline int dotProductSign3D_IEI_interval(const genericPoint& p, const genericPoint& q, interval_number rx, interval_number ry, interval_number rz)
{
   interval_number lpx, lpy, lpz, dp, lqx, lqy, lqz, dq;
   if (
   !p.getIntervalLambda(lpx, lpy, lpz, dp)
   || !q.getIntervalLambda(lqx, lqy, lqz, dq)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number dqp(dq * dp);
   const interval_number pxq(lpx * dqp);
   const interval_number pyq(lpy * dqp);
   const interval_number pzq(lpz * dqp);
   const interval_number rxq(rx * dq);
   const interval_number ryq(ry * dq);
   const interval_number rzq(rz * dq);
   const interval_number lqxd(lqx * dp);
   const interval_number lqyd(lqy * dp);
   const interval_number lqzd(lqz * dp);
   const interval_number lx(pxq - lqxd);
   const interval_number ly(pyq - lqyd);
   const interval_number lz(pzq - lqzd);
   const interval_number gx(rxq - lqx);
   const interval_number gy(ryq - lqy);
   const interval_number gz(rzq - lqz);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number dz(lz * gz);
   const interval_number d1(dx + dy);
   const interval_number d(d1 + dz);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign3D_IEI_bigfloat(const genericPoint& p, const genericPoint& q, bigfloat rx, bigfloat ry, bigfloat rz)
{
   bigfloat lpx, lpy, lpz, dp, lqx, lqy, lqz, dq;
   p.getBigfloatLambda(lpx, lpy, lpz, dp);
   q.getBigfloatLambda(lqx, lqy, lqz, dq);
   const bigfloat dqp(dq * dp);
   const bigfloat pxq(lpx * dqp);
   const bigfloat pyq(lpy * dqp);
   const bigfloat pzq(lpz * dqp);
   const bigfloat rxq(rx * dq);
   const bigfloat ryq(ry * dq);
   const bigfloat rzq(rz * dq);
   const bigfloat lqxd(lqx * dp);
   const bigfloat lqyd(lqy * dp);
   const bigfloat lqzd(lqz * dp);
   const bigfloat lx(pxq - lqxd);
   const bigfloat ly(pyq - lqyd);
   const bigfloat lz(pzq - lqzd);
   const bigfloat gx(rxq - lqx);
   const bigfloat gy(ryq - lqy);
   const bigfloat gz(rzq - lqz);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat dz(lz * gz);
   const bigfloat d1(dx + dy);
   const bigfloat d(d1 + dz);
   return sgn(d);
}

inline int dotProductSign3D_IEI(const genericPoint& p, const genericPoint& q, double rx, double ry, double rz)
{
   int ret;
   ret = dotProductSign3D_IEI_interval(p, q, rx, ry, rz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign3D_IEI_bigfloat(p, q, rx, ry, rz);
}

inline int dotProductSign3D_IIE_interval(const genericPoint& p, const genericPoint& r, interval_number qx, interval_number qy, interval_number qz)
{
   interval_number lpx, lpy, lpz, dp, lrx, lry, lrz, dr;
   if (
   !p.getIntervalLambda(lpx, lpy, lpz, dp)
   || !r.getIntervalLambda(lrx, lry, lrz, dr)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number qxd(qx * dp);
   const interval_number qyd(qy * dp);
   const interval_number qzd(qz * dp);
   const interval_number lx(lpx - qxd);
   const interval_number ly(lpy - qyd);
   const interval_number lz(lpz - qzd);
   const interval_number qxr(qx * dr);
   const interval_number qyr(qy * dr);
   const interval_number qzr(qz * dr);
   const interval_number gx(lrx - qxr);
   const interval_number gy(lry - qyr);
   const interval_number gz(lrz - qzr);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number dz(lz * gz);
   const interval_number d1(dx + dy);
   const interval_number d(d1 + dz);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign3D_IIE_bigfloat(const genericPoint& p, const genericPoint& r, bigfloat qx, bigfloat qy, bigfloat qz)
{
   bigfloat lpx, lpy, lpz, dp, lrx, lry, lrz, dr;
   p.getBigfloatLambda(lpx, lpy, lpz, dp);
   r.getBigfloatLambda(lrx, lry, lrz, dr);
   const bigfloat qxd(qx * dp);
   const bigfloat qyd(qy * dp);
   const bigfloat qzd(qz * dp);
   const bigfloat lx(lpx - qxd);
   const bigfloat ly(lpy - qyd);
   const bigfloat lz(lpz - qzd);
   const bigfloat qxr(qx * dr);
   const bigfloat qyr(qy * dr);
   const bigfloat qzr(qz * dr);
   const bigfloat gx(lrx - qxr);
   const bigfloat gy(lry - qyr);
   const bigfloat gz(lrz - qzr);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat dz(lz * gz);
   const bigfloat d1(dx + dy);
   const bigfloat d(d1 + dz);
   return sgn(d);
}

inline int dotProductSign3D_IIE(const genericPoint& p, const genericPoint& r, double qx, double qy, double qz)
{
   int ret;
   ret = dotProductSign3D_IIE_interval(p, r, qx, qy, qz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign3D_IIE_bigfloat(p, r, qx, qy, qz);
}

inline int dotProductSign3D_III_interval(const genericPoint& p, const genericPoint& r, const genericPoint& q)
{
   interval_number lpx, lpy, lpz, dp, lrx, lry, lrz, dr, lqx, lqy, lqz, dq;
   if (
   !p.getIntervalLambda(lpx, lpy, lpz, dp)
   || !r.getIntervalLambda(lrx, lry, lrz, dr)
   || !q.getIntervalLambda(lqx, lqy, lqz, dq)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number qxd(lqx * dp);
   const interval_number qyd(lqy * dp);
   const interval_number qzd(lqz * dp);
   const interval_number lpxq(lpx * dq);
   const interval_number lpyq(lpy * dq);
   const interval_number lpzq(lpz * dq);
   const interval_number lx(lpxq - qxd);
   const interval_number ly(lpyq - qyd);
   const interval_number lz(lpzq - qzd);
   const interval_number qxr(lqx * dr);
   const interval_number qyr(lqy * dr);
   const interval_number qzr(lqz * dr);
   const interval_number lrxq(lrx * dq);
   const interval_number lryq(lry * dq);
   const interval_number lrzq(lrz * dq);
   const interval_number gx(lrxq - qxr);
   const interval_number gy(lryq - qyr);
   const interval_number gz(lrzq - qzr);
   const interval_number dx(lx * gx);
   const interval_number dy(ly * gy);
   const interval_number dz(lz * gz);
   const interval_number d1(dx + dy);
   const interval_number d(d1 + dz);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

inline int dotProductSign3D_III_bigfloat(const genericPoint& p, const genericPoint& r, const genericPoint& q)
{
   bigfloat lpx, lpy, lpz, dp, lrx, lry, lrz, dr, lqx, lqy, lqz, dq;
   p.getBigfloatLambda(lpx, lpy, lpz, dp);
   r.getBigfloatLambda(lrx, lry, lrz, dr);
   q.getBigfloatLambda(lqx, lqy, lqz, dq);
   const bigfloat qxd(lqx * dp);
   const bigfloat qyd(lqy * dp);
   const bigfloat qzd(lqz * dp);
   const bigfloat lpxq(lpx * dq);
   const bigfloat lpyq(lpy * dq);
   const bigfloat lpzq(lpz * dq);
   const bigfloat lx(lpxq - qxd);
   const bigfloat ly(lpyq - qyd);
   const bigfloat lz(lpzq - qzd);
   const bigfloat qxr(lqx * dr);
   const bigfloat qyr(lqy * dr);
   const bigfloat qzr(lqz * dr);
   const bigfloat lrxq(lrx * dq);
   const bigfloat lryq(lry * dq);
   const bigfloat lrzq(lrz * dq);
   const bigfloat gx(lrxq - qxr);
   const bigfloat gy(lryq - qyr);
   const bigfloat gz(lrzq - qzr);
   const bigfloat dx(lx * gx);
   const bigfloat dy(ly * gy);
   const bigfloat dz(lz * gz);
   const bigfloat d1(dx + dy);
   const bigfloat d(d1 + dz);
   return sgn(d);
}

inline int dotProductSign3D_III(const genericPoint& p, const genericPoint& r, const genericPoint& q)
{
   int ret;
   ret = dotProductSign3D_III_interval(p, r, q);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return dotProductSign3D_III_bigfloat(p, r, q);
}

inline int incirclexy_indirect_IEEE_interval(const genericPoint& p1, interval_number pbx, interval_number pby, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pdxt(pdx * d1);
   const interval_number pdyt(pdy * d1);
   const interval_number adx(l1x - pdxt);
   const interval_number ady(l1y - pdyt);
   const interval_number bdx(pbx - pdx);
   const interval_number bdy(pby - pdy);
   const interval_number cdx(pcx - pdx);
   const interval_number cdy(pcy - pdy);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number alift(alifta + aliftb);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number blift(blifta + bliftb);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number clift(clifta + cliftb);
   const interval_number la(alift * bcdet);
   const interval_number lbt(blift * cadet);
   const interval_number lb(lbt * d1);
   const interval_number lct(clift * abdet);
   const interval_number lc(lct * d1);
   const interval_number lab(la + lb);
   const interval_number L(lab + lc);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incirclexy_indirect_IEEE_bigfloat(const genericPoint& p1, bigfloat pbx, bigfloat pby, bigfloat pcx, bigfloat pcy, bigfloat pdx, bigfloat pdy)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat pdxt(pdx * d1);
   const bigfloat pdyt(pdy * d1);
   const bigfloat adx(l1x - pdxt);
   const bigfloat ady(l1y - pdyt);
   const bigfloat bdx(pbx - pdx);
   const bigfloat bdy(pby - pdy);
   const bigfloat cdx(pcx - pdx);
   const bigfloat cdy(pcy - pdy);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat alift(alifta + aliftb);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat blift(blifta + bliftb);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat clift(clifta + cliftb);
   const bigfloat la(alift * bcdet);
   const bigfloat lbt(blift * cadet);
   const bigfloat lb(lbt * d1);
   const bigfloat lct(clift * abdet);
   const bigfloat lc(lct * d1);
   const bigfloat lab(la + lb);
   const bigfloat L(lab + lc);
   return sgn(L);
}

inline int incirclexy_indirect_IEEE(const genericPoint& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incirclexy_indirect_IEEE_interval(p1, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incirclexy_indirect_IEEE_bigfloat(p1, pbx, pby, pcx, pcy, pdx, pdy);
}

inline int incirclexy_indirect_IIEE_interval(const genericPoint& p1, const genericPoint& p2, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pdx1(pdx * d1);
   const interval_number pdy1(pdy * d1);
   const interval_number adx(l1x - pdx1);
   const interval_number ady(l1y - pdy1);
   const interval_number pdx2(pdx * d2);
   const interval_number pdy2(pdy * d2);
   const interval_number bdx(l2x - pdx2);
   const interval_number bdy(l2y - pdy2);
   const interval_number cdx(pcx - pdx);
   const interval_number cdy(pcy - pdy);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number aliftt(alifta + aliftb);
   const interval_number alift(aliftt * d2);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number blift(blifta + bliftb);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number cliftt(clifta + cliftb);
   const interval_number clift(cliftt * d2);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab(lc + lb);
   const interval_number lab2(lab * d1);
   const interval_number L(lab2 + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incirclexy_indirect_IIEE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat pcx, bigfloat pcy, bigfloat pdx, bigfloat pdy)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat pdx1(pdx * d1);
   const bigfloat pdy1(pdy * d1);
   const bigfloat adx(l1x - pdx1);
   const bigfloat ady(l1y - pdy1);
   const bigfloat pdx2(pdx * d2);
   const bigfloat pdy2(pdy * d2);
   const bigfloat bdx(l2x - pdx2);
   const bigfloat bdy(l2y - pdy2);
   const bigfloat cdx(pcx - pdx);
   const bigfloat cdy(pcy - pdy);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat aliftt(alifta + aliftb);
   const bigfloat alift(aliftt * d2);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat blift(blifta + bliftb);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat cliftt(clifta + cliftb);
   const bigfloat clift(cliftt * d2);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab(lc + lb);
   const bigfloat lab2(lab * d1);
   const bigfloat L(lab2 + la);
   return sgn(L);
}

inline int incirclexy_indirect_IIEE(const genericPoint& p1, const genericPoint& p2, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incirclexy_indirect_IIEE_interval(p1, p2, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incirclexy_indirect_IIEE_bigfloat(p1, p2, pcx, pcy, pdx, pdy);
}

inline int incirclexy_indirect_IIIE_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pdx1(pdx * d1);
   const interval_number pdy1(pdy * d1);
   const interval_number adx(l1x - pdx1);
   const interval_number ady(l1y - pdy1);
   const interval_number pdx2(pdx * d2);
   const interval_number pdy2(pdy * d2);
   const interval_number bdx(l2x - pdx2);
   const interval_number bdy(l2y - pdy2);
   const interval_number pdx3(pdx * d3);
   const interval_number pdy3(pdy * d3);
   const interval_number cdx(l3x - pdx3);
   const interval_number cdy(l3y - pdy3);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number aliftt(alifta + aliftb);
   const interval_number alift2(aliftt * d2);
   const interval_number alift(alift2 * d3);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number bliftt(blifta + bliftb);
   const interval_number blift(bliftt * d3);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number cliftt(clifta + cliftb);
   const interval_number clift(cliftt * d2);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab2(lc + lb);
   const interval_number lab(lab2 * d1);
   const interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incirclexy_indirect_IIIE_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, bigfloat pdx, bigfloat pdy)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   const bigfloat pdx1(pdx * d1);
   const bigfloat pdy1(pdy * d1);
   const bigfloat adx(l1x - pdx1);
   const bigfloat ady(l1y - pdy1);
   const bigfloat pdx2(pdx * d2);
   const bigfloat pdy2(pdy * d2);
   const bigfloat bdx(l2x - pdx2);
   const bigfloat bdy(l2y - pdy2);
   const bigfloat pdx3(pdx * d3);
   const bigfloat pdy3(pdy * d3);
   const bigfloat cdx(l3x - pdx3);
   const bigfloat cdy(l3y - pdy3);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat aliftt(alifta + aliftb);
   const bigfloat alift2(aliftt * d2);
   const bigfloat alift(alift2 * d3);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat bliftt(blifta + bliftb);
   const bigfloat blift(bliftt * d3);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat cliftt(clifta + cliftb);
   const bigfloat clift(cliftt * d2);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab2(lc + lb);
   const bigfloat lab(lab2 * d1);
   const bigfloat L(lab + la);
   return sgn(L);
}

inline int incirclexy_indirect_IIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, double pdx, double pdy)
{
   int ret;
   ret = incirclexy_indirect_IIIE_interval(p1, p2, p3, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incirclexy_indirect_IIIE_bigfloat(p1, p2, p3, pdx, pdy);
}

inline int incirclexy_indirect_IIII_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   || !p4.getIntervalLambda(l4x, l4y, l4z, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number l1xt(l1x * d4);
   const interval_number l1yt(l1y * d4);
   const interval_number l2xt(l2x * d4);
   const interval_number l2yt(l2y * d4);
   const interval_number l3xt(l3x * d4);
   const interval_number l3yt(l3y * d4);
   const interval_number l4x1(l4x * d1);
   const interval_number l4y1(l4y * d1);
   const interval_number adx(l1xt - l4x1);
   const interval_number ady(l1yt - l4y1);
   const interval_number l4x2(l4x * d2);
   const interval_number l4y2(l4y * d2);
   const interval_number bdx(l2xt - l4x2);
   const interval_number bdy(l2yt - l4y2);
   const interval_number l4x3(l4x * d3);
   const interval_number l4y3(l4y * d3);
   const interval_number cdx(l3xt - l4x3);
   const interval_number cdy(l3yt - l4y3);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number aliftt(alifta + aliftb);
   const interval_number alift2(aliftt * d2);
   const interval_number alift(alift2 * d3);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number bliftt(blifta + bliftb);
   const interval_number blift(bliftt * d3);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number cliftt(clifta + cliftb);
   const interval_number clift(cliftt * d2);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab2(lc + lb);
   const interval_number lab(lab2 * d1);
   const interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incirclexy_indirect_IIII_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   p4.getBigfloatLambda(l4x, l4y, l4z, d4);
   const bigfloat l1xt(l1x * d4);
   const bigfloat l1yt(l1y * d4);
   const bigfloat l2xt(l2x * d4);
   const bigfloat l2yt(l2y * d4);
   const bigfloat l3xt(l3x * d4);
   const bigfloat l3yt(l3y * d4);
   const bigfloat l4x1(l4x * d1);
   const bigfloat l4y1(l4y * d1);
   const bigfloat adx(l1xt - l4x1);
   const bigfloat ady(l1yt - l4y1);
   const bigfloat l4x2(l4x * d2);
   const bigfloat l4y2(l4y * d2);
   const bigfloat bdx(l2xt - l4x2);
   const bigfloat bdy(l2yt - l4y2);
   const bigfloat l4x3(l4x * d3);
   const bigfloat l4y3(l4y * d3);
   const bigfloat cdx(l3xt - l4x3);
   const bigfloat cdy(l3yt - l4y3);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat aliftt(alifta + aliftb);
   const bigfloat alift2(aliftt * d2);
   const bigfloat alift(alift2 * d3);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat bliftt(blifta + bliftb);
   const bigfloat blift(bliftt * d3);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat cliftt(clifta + cliftb);
   const bigfloat clift(cliftt * d2);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab2(lc + lb);
   const bigfloat lab(lab2 * d1);
   const bigfloat L(lab + la);
   return sgn(L);
}

inline int incirclexy_indirect_IIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   int ret;
   ret = incirclexy_indirect_IIII_interval(p1, p2, p3, p4);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incirclexy_indirect_IIII_bigfloat(p1, p2, p3, p4);
}

inline int incircle_indirect_IEEE_interval(const genericPoint& p1, interval_number pbx, interval_number pby, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pdxt(pdx * d1);
   const interval_number pdyt(pdy * d1);
   const interval_number adx(l1x - pdxt);
   const interval_number ady(l1y - pdyt);
   const interval_number bdx(pbx - pdx);
   const interval_number bdy(pby - pdy);
   const interval_number cdx(pcx - pdx);
   const interval_number cdy(pcy - pdy);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number alift(alifta + aliftb);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number blift(blifta + bliftb);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number clift(clifta + cliftb);
   const interval_number la(alift * bcdet);
   const interval_number lbt(blift * cadet);
   const interval_number lb(lbt * d1);
   const interval_number lct(clift * abdet);
   const interval_number lc(lct * d1);
   const interval_number lab(la + lb);
   const interval_number L(lab + lc);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incircle_indirect_IEEE_bigfloat(const genericPoint& p1, bigfloat pbx, bigfloat pby, bigfloat pcx, bigfloat pcy, bigfloat pdx, bigfloat pdy)
{
   bigfloat l1x, l1y, d1;
   p1.getBigfloatLambda(l1x, l1y, d1);
   const bigfloat pdxt(pdx * d1);
   const bigfloat pdyt(pdy * d1);
   const bigfloat adx(l1x - pdxt);
   const bigfloat ady(l1y - pdyt);
   const bigfloat bdx(pbx - pdx);
   const bigfloat bdy(pby - pdy);
   const bigfloat cdx(pcx - pdx);
   const bigfloat cdy(pcy - pdy);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat alift(alifta + aliftb);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat blift(blifta + bliftb);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat clift(clifta + cliftb);
   const bigfloat la(alift * bcdet);
   const bigfloat lbt(blift * cadet);
   const bigfloat lb(lbt * d1);
   const bigfloat lct(clift * abdet);
   const bigfloat lc(lct * d1);
   const bigfloat lab(la + lb);
   const bigfloat L(lab + lc);
   return sgn(L);
}

inline int incircle_indirect_IEEE(const genericPoint& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_IEEE_interval(p1, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_IEEE_bigfloat(p1, pbx, pby, pcx, pcy, pdx, pdy);
}

inline int incircle_indirect_IIEE_interval(const genericPoint& p1, const genericPoint& p2, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pdx1(pdx * d1);
   const interval_number pdy1(pdy * d1);
   const interval_number adx(l1x - pdx1);
   const interval_number ady(l1y - pdy1);
   const interval_number pdx2(pdx * d2);
   const interval_number pdy2(pdy * d2);
   const interval_number bdx(l2x - pdx2);
   const interval_number bdy(l2y - pdy2);
   const interval_number cdx(pcx - pdx);
   const interval_number cdy(pcy - pdy);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number aliftt(alifta + aliftb);
   const interval_number alift(aliftt * d2);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number blift(blifta + bliftb);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number cliftt(clifta + cliftb);
   const interval_number clift(cliftt * d2);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab(lc + lb);
   const interval_number lab2(lab * d1);
   const interval_number L(lab2 + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incircle_indirect_IIEE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat pcx, bigfloat pcy, bigfloat pdx, bigfloat pdy)
{
   bigfloat l1x, l1y, d1, l2x, l2y, d2;
   p1.getBigfloatLambda(l1x, l1y, d1);
   p2.getBigfloatLambda(l2x, l2y, d2);
   const bigfloat pdx1(pdx * d1);
   const bigfloat pdy1(pdy * d1);
   const bigfloat adx(l1x - pdx1);
   const bigfloat ady(l1y - pdy1);
   const bigfloat pdx2(pdx * d2);
   const bigfloat pdy2(pdy * d2);
   const bigfloat bdx(l2x - pdx2);
   const bigfloat bdy(l2y - pdy2);
   const bigfloat cdx(pcx - pdx);
   const bigfloat cdy(pcy - pdy);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat aliftt(alifta + aliftb);
   const bigfloat alift(aliftt * d2);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat blift(blifta + bliftb);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat cliftt(clifta + cliftb);
   const bigfloat clift(cliftt * d2);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab(lc + lb);
   const bigfloat lab2(lab * d1);
   const bigfloat L(lab2 + la);
   return sgn(L);
}

inline int incircle_indirect_IIEE(const genericPoint& p1, const genericPoint& p2, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_IIEE_interval(p1, p2, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_IIEE_bigfloat(p1, p2, pcx, pcy, pdx, pdy);
}

inline int incircle_indirect_IIIE_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   || !p3.getIntervalLambda(l3x, l3y, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pdx1(pdx * d1);
   const interval_number pdy1(pdy * d1);
   const interval_number adx(l1x - pdx1);
   const interval_number ady(l1y - pdy1);
   const interval_number pdx2(pdx * d2);
   const interval_number pdy2(pdy * d2);
   const interval_number bdx(l2x - pdx2);
   const interval_number bdy(l2y - pdy2);
   const interval_number pdx3(pdx * d3);
   const interval_number pdy3(pdy * d3);
   const interval_number cdx(l3x - pdx3);
   const interval_number cdy(l3y - pdy3);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number aliftt(alifta + aliftb);
   const interval_number alift2(aliftt * d2);
   const interval_number alift(alift2 * d3);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number bliftt(blifta + bliftb);
   const interval_number blift(bliftt * d3);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number cliftt(clifta + cliftb);
   const interval_number clift(cliftt * d2);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab2(lc + lb);
   const interval_number lab(lab2 * d1);
   const interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incircle_indirect_IIIE_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, bigfloat pdx, bigfloat pdy)
{
   bigfloat l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
   p1.getBigfloatLambda(l1x, l1y, d1);
   p2.getBigfloatLambda(l2x, l2y, d2);
   p3.getBigfloatLambda(l3x, l3y, d3);
   const bigfloat pdx1(pdx * d1);
   const bigfloat pdy1(pdy * d1);
   const bigfloat adx(l1x - pdx1);
   const bigfloat ady(l1y - pdy1);
   const bigfloat pdx2(pdx * d2);
   const bigfloat pdy2(pdy * d2);
   const bigfloat bdx(l2x - pdx2);
   const bigfloat bdy(l2y - pdy2);
   const bigfloat pdx3(pdx * d3);
   const bigfloat pdy3(pdy * d3);
   const bigfloat cdx(l3x - pdx3);
   const bigfloat cdy(l3y - pdy3);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat aliftt(alifta + aliftb);
   const bigfloat alift2(aliftt * d2);
   const bigfloat alift(alift2 * d3);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat bliftt(blifta + bliftb);
   const bigfloat blift(bliftt * d3);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat cliftt(clifta + cliftb);
   const bigfloat clift(cliftt * d2);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab2(lc + lb);
   const bigfloat lab(lab2 * d1);
   const bigfloat L(lab + la);
   return sgn(L);
}

inline int incircle_indirect_IIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_IIIE_interval(p1, p2, p3, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_IIIE_bigfloat(p1, p2, p3, pdx, pdy);
}

inline int incircle_indirect_IIII_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, l4x, l4y, d4;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   || !p3.getIntervalLambda(l3x, l3y, d3)
   || !p4.getIntervalLambda(l4x, l4y, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number l1xt(l1x * d4);
   const interval_number l1yt(l1y * d4);
   const interval_number l2xt(l2x * d4);
   const interval_number l2yt(l2y * d4);
   const interval_number l3xt(l3x * d4);
   const interval_number l3yt(l3y * d4);
   const interval_number l4x1(l4x * d1);
   const interval_number l4y1(l4y * d1);
   const interval_number adx(l1xt - l4x1);
   const interval_number ady(l1yt - l4y1);
   const interval_number l4x2(l4x * d2);
   const interval_number l4y2(l4y * d2);
   const interval_number bdx(l2xt - l4x2);
   const interval_number bdy(l2yt - l4y2);
   const interval_number l4x3(l4x * d3);
   const interval_number l4y3(l4y * d3);
   const interval_number cdx(l3xt - l4x3);
   const interval_number cdy(l3yt - l4y3);
   const interval_number abdeta(adx * bdy);
   const interval_number abdetb(bdx * ady);
   const interval_number abdet(abdeta - abdetb);
   const interval_number bcdeta(bdx * cdy);
   const interval_number bcdetb(cdx * bdy);
   const interval_number bcdet(bcdeta - bcdetb);
   const interval_number cadeta(cdx * ady);
   const interval_number cadetb(adx * cdy);
   const interval_number cadet(cadeta - cadetb);
   const interval_number alifta(adx * adx);
   const interval_number aliftb(ady * ady);
   const interval_number aliftt(alifta + aliftb);
   const interval_number alift2(aliftt * d2);
   const interval_number alift(alift2 * d3);
   const interval_number blifta(bdx * bdx);
   const interval_number bliftb(bdy * bdy);
   const interval_number bliftt(blifta + bliftb);
   const interval_number blift(bliftt * d3);
   const interval_number clifta(cdx * cdx);
   const interval_number cliftb(cdy * cdy);
   const interval_number cliftt(clifta + cliftb);
   const interval_number clift(cliftt * d2);
   const interval_number la(alift * bcdet);
   const interval_number lb(blift * cadet);
   const interval_number lc(clift * abdet);
   const interval_number lab2(lc + lb);
   const interval_number lab(lab2 * d1);
   const interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int incircle_indirect_IIII_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   bigfloat l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, l4x, l4y, d4;
   p1.getBigfloatLambda(l1x, l1y, d1);
   p2.getBigfloatLambda(l2x, l2y, d2);
   p3.getBigfloatLambda(l3x, l3y, d3);
   p4.getBigfloatLambda(l4x, l4y, d4);
   const bigfloat l1xt(l1x * d4);
   const bigfloat l1yt(l1y * d4);
   const bigfloat l2xt(l2x * d4);
   const bigfloat l2yt(l2y * d4);
   const bigfloat l3xt(l3x * d4);
   const bigfloat l3yt(l3y * d4);
   const bigfloat l4x1(l4x * d1);
   const bigfloat l4y1(l4y * d1);
   const bigfloat adx(l1xt - l4x1);
   const bigfloat ady(l1yt - l4y1);
   const bigfloat l4x2(l4x * d2);
   const bigfloat l4y2(l4y * d2);
   const bigfloat bdx(l2xt - l4x2);
   const bigfloat bdy(l2yt - l4y2);
   const bigfloat l4x3(l4x * d3);
   const bigfloat l4y3(l4y * d3);
   const bigfloat cdx(l3xt - l4x3);
   const bigfloat cdy(l3yt - l4y3);
   const bigfloat abdeta(adx * bdy);
   const bigfloat abdetb(bdx * ady);
   const bigfloat abdet(abdeta - abdetb);
   const bigfloat bcdeta(bdx * cdy);
   const bigfloat bcdetb(cdx * bdy);
   const bigfloat bcdet(bcdeta - bcdetb);
   const bigfloat cadeta(cdx * ady);
   const bigfloat cadetb(adx * cdy);
   const bigfloat cadet(cadeta - cadetb);
   const bigfloat alifta(adx * adx);
   const bigfloat aliftb(ady * ady);
   const bigfloat aliftt(alifta + aliftb);
   const bigfloat alift2(aliftt * d2);
   const bigfloat alift(alift2 * d3);
   const bigfloat blifta(bdx * bdx);
   const bigfloat bliftb(bdy * bdy);
   const bigfloat bliftt(blifta + bliftb);
   const bigfloat blift(bliftt * d3);
   const bigfloat clifta(cdx * cdx);
   const bigfloat cliftb(cdy * cdy);
   const bigfloat cliftt(clifta + cliftb);
   const bigfloat clift(cliftt * d2);
   const bigfloat la(alift * bcdet);
   const bigfloat lb(blift * cadet);
   const bigfloat lc(clift * abdet);
   const bigfloat lab2(lc + lb);
   const bigfloat lab(lab2 * d1);
   const bigfloat L(lab + la);
   return sgn(L);
}

inline int incircle_indirect_IIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   int ret;
   ret = incircle_indirect_IIII_interval(p1, p2, p3, p4);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_IIII_bigfloat(p1, p2, p3, p4);
}

inline int inGabrielSphere_EIEE_interval(const genericPoint& a, interval_number qx, interval_number qy, interval_number qz, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l2x, l2y, l2z, d2;
   if (
   !a.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number bxd2(bx * d2);
   const interval_number byd2(by * d2);
   const interval_number bzd2(bz * d2);
   const interval_number cxd2(cx * d2);
   const interval_number cyd2(cy * d2);
   const interval_number czd2(cz * d2);
   const interval_number qxd2(qx * d2);
   const interval_number qyd2(qy * d2);
   const interval_number qzd2(qz * d2);
   const interval_number bax(bxd2 - l2x);
   const interval_number bay(byd2 - l2y);
   const interval_number baz(bzd2 - l2z);
   const interval_number cax(cxd2 - l2x);
   const interval_number cay(cyd2 - l2y);
   const interval_number caz(czd2 - l2z);
   const interval_number qax(qxd2 - l2x);
   const interval_number qay(qyd2 - l2y);
   const interval_number qaz(qzd2 - l2z);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_EIEE_bigfloat(const genericPoint& a, bigfloat qx, bigfloat qy, bigfloat qz, bigfloat bx, bigfloat by, bigfloat bz, bigfloat cx, bigfloat cy, bigfloat cz)
{
   bigfloat l2x, l2y, l2z, d2;
   a.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat bxd2(bx * d2);
   const bigfloat byd2(by * d2);
   const bigfloat bzd2(bz * d2);
   const bigfloat cxd2(cx * d2);
   const bigfloat cyd2(cy * d2);
   const bigfloat czd2(cz * d2);
   const bigfloat qxd2(qx * d2);
   const bigfloat qyd2(qy * d2);
   const bigfloat qzd2(qz * d2);
   const bigfloat bax(bxd2 - l2x);
   const bigfloat bay(byd2 - l2y);
   const bigfloat baz(bzd2 - l2z);
   const bigfloat cax(cxd2 - l2x);
   const bigfloat cay(cyd2 - l2y);
   const bigfloat caz(czd2 - l2z);
   const bigfloat qax(qxd2 - l2x);
   const bigfloat qay(qyd2 - l2y);
   const bigfloat qaz(qzd2 - l2z);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_EIEE(const genericPoint& a, double qx, double qy, double qz, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = inGabrielSphere_EIEE_interval(a, qx, qy, qz, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_EIEE_bigfloat(a, qx, qy, qz, bx, by, bz, cx, cy, cz);
}

inline int inGabrielSphere_EIIE_interval(const genericPoint& a, const genericPoint& b, interval_number qx, interval_number qy, interval_number qz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !a.getIntervalLambda(l2x, l2y, l2z, d2)
   || !b.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d23(d2 * d3);
   const interval_number bxd2(l3x * d2);
   const interval_number byd2(l3y * d2);
   const interval_number bzd2(l3z * d2);
   const interval_number cxd2(cx * d23);
   const interval_number cyd2(cy * d23);
   const interval_number czd2(cz * d23);
   const interval_number qxd2(qx * d23);
   const interval_number qyd2(qy * d23);
   const interval_number qzd2(qz * d23);
   const interval_number l2x3(l2x * d3);
   const interval_number l2y3(l2y * d3);
   const interval_number l2z3(l2z * d3);
   const interval_number bax(bxd2 - l2x3);
   const interval_number bay(byd2 - l2y3);
   const interval_number baz(bzd2 - l2z3);
   const interval_number cax(cxd2 - l2x3);
   const interval_number cay(cyd2 - l2y3);
   const interval_number caz(czd2 - l2z3);
   const interval_number qax(qxd2 - l2x3);
   const interval_number qay(qyd2 - l2y3);
   const interval_number qaz(qzd2 - l2z3);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_EIIE_bigfloat(const genericPoint& a, const genericPoint& b, bigfloat qx, bigfloat qy, bigfloat qz, bigfloat cx, bigfloat cy, bigfloat cz)
{
   bigfloat l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   a.getBigfloatLambda(l2x, l2y, l2z, d2);
   b.getBigfloatLambda(l3x, l3y, l3z, d3);
   const bigfloat d23(d2 * d3);
   const bigfloat bxd2(l3x * d2);
   const bigfloat byd2(l3y * d2);
   const bigfloat bzd2(l3z * d2);
   const bigfloat cxd2(cx * d23);
   const bigfloat cyd2(cy * d23);
   const bigfloat czd2(cz * d23);
   const bigfloat qxd2(qx * d23);
   const bigfloat qyd2(qy * d23);
   const bigfloat qzd2(qz * d23);
   const bigfloat l2x3(l2x * d3);
   const bigfloat l2y3(l2y * d3);
   const bigfloat l2z3(l2z * d3);
   const bigfloat bax(bxd2 - l2x3);
   const bigfloat bay(byd2 - l2y3);
   const bigfloat baz(bzd2 - l2z3);
   const bigfloat cax(cxd2 - l2x3);
   const bigfloat cay(cyd2 - l2y3);
   const bigfloat caz(czd2 - l2z3);
   const bigfloat qax(qxd2 - l2x3);
   const bigfloat qay(qyd2 - l2y3);
   const bigfloat qaz(qzd2 - l2z3);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_EIIE(const genericPoint& a, const genericPoint& b, double qx, double qy, double qz, double cx, double cy, double cz)
{
   int ret;
   ret = inGabrielSphere_EIIE_interval(a, b, qx, qy, qz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_EIIE_bigfloat(a, b, qx, qy, qz, cx, cy, cz);
}

inline int inGabrielSphere_EIII_interval(const genericPoint& a, const genericPoint& b, const genericPoint& c, interval_number qx, interval_number qy, interval_number qz)
{
   interval_number l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   if (
   !a.getIntervalLambda(l2x, l2y, l2z, d2)
   || !b.getIntervalLambda(l3x, l3y, l3z, d3)
   || !c.getIntervalLambda(l4x, l4y, l4z, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d23(d2 * d3);
   const interval_number d24(d2 * d4);
   const interval_number d34(d3 * d4);
   const interval_number d234(d23 * d4);
   const interval_number bxd2(l3x * d24);
   const interval_number byd2(l3y * d24);
   const interval_number bzd2(l3z * d24);
   const interval_number cxd2(l4x * d23);
   const interval_number cyd2(l4y * d23);
   const interval_number czd2(l4z * d23);
   const interval_number qxd2(qx * d234);
   const interval_number qyd2(qy * d234);
   const interval_number qzd2(qz * d234);
   const interval_number l2xm(l2x * d34);
   const interval_number l2ym(l2y * d34);
   const interval_number l2zm(l2z * d34);
   const interval_number bax(bxd2 - l2xm);
   const interval_number bay(byd2 - l2ym);
   const interval_number baz(bzd2 - l2zm);
   const interval_number cax(cxd2 - l2xm);
   const interval_number cay(cyd2 - l2ym);
   const interval_number caz(czd2 - l2zm);
   const interval_number qax(qxd2 - l2xm);
   const interval_number qay(qyd2 - l2ym);
   const interval_number qaz(qzd2 - l2zm);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_EIII_bigfloat(const genericPoint& a, const genericPoint& b, const genericPoint& c, bigfloat qx, bigfloat qy, bigfloat qz)
{
   bigfloat l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   a.getBigfloatLambda(l2x, l2y, l2z, d2);
   b.getBigfloatLambda(l3x, l3y, l3z, d3);
   c.getBigfloatLambda(l4x, l4y, l4z, d4);
   const bigfloat d23(d2 * d3);
   const bigfloat d24(d2 * d4);
   const bigfloat d34(d3 * d4);
   const bigfloat d234(d23 * d4);
   const bigfloat bxd2(l3x * d24);
   const bigfloat byd2(l3y * d24);
   const bigfloat bzd2(l3z * d24);
   const bigfloat cxd2(l4x * d23);
   const bigfloat cyd2(l4y * d23);
   const bigfloat czd2(l4z * d23);
   const bigfloat qxd2(qx * d234);
   const bigfloat qyd2(qy * d234);
   const bigfloat qzd2(qz * d234);
   const bigfloat l2xm(l2x * d34);
   const bigfloat l2ym(l2y * d34);
   const bigfloat l2zm(l2z * d34);
   const bigfloat bax(bxd2 - l2xm);
   const bigfloat bay(byd2 - l2ym);
   const bigfloat baz(bzd2 - l2zm);
   const bigfloat cax(cxd2 - l2xm);
   const bigfloat cay(cyd2 - l2ym);
   const bigfloat caz(czd2 - l2zm);
   const bigfloat qax(qxd2 - l2xm);
   const bigfloat qay(qyd2 - l2ym);
   const bigfloat qaz(qzd2 - l2zm);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_EIII(const genericPoint& a, const genericPoint& b, const genericPoint& c, double qx, double qy, double qz)
{
   int ret;
   ret = inGabrielSphere_EIII_interval(a, b, c, qx, qy, qz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_EIII_bigfloat(a, b, c, qx, qy, qz);
}

inline int inGabrielSphere_IEEE_interval(const genericPoint& q, interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !q.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number baxs(bx - ax);
   const interval_number bays(by - ay);
   const interval_number bazs(bz - az);
   const interval_number caxs(cx - ax);
   const interval_number cays(cy - ay);
   const interval_number cazs(cz - az);
   const interval_number bax(baxs * d1);
   const interval_number bay(bays * d1);
   const interval_number baz(bazs * d1);
   const interval_number cax(caxs * d1);
   const interval_number cay(cays * d1);
   const interval_number caz(cazs * d1);
   const interval_number axs(ax * d1);
   const interval_number ays(ay * d1);
   const interval_number azs(az * d1);
   const interval_number qax(l1x - axs);
   const interval_number qay(l1y - ays);
   const interval_number qaz(l1z - azs);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_IEEE_bigfloat(const genericPoint& q, bigfloat ax, bigfloat ay, bigfloat az, bigfloat bx, bigfloat by, bigfloat bz, bigfloat cx, bigfloat cy, bigfloat cz)
{
   bigfloat l1x, l1y, l1z, d1;
   q.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat baxs(bx - ax);
   const bigfloat bays(by - ay);
   const bigfloat bazs(bz - az);
   const bigfloat caxs(cx - ax);
   const bigfloat cays(cy - ay);
   const bigfloat cazs(cz - az);
   const bigfloat bax(baxs * d1);
   const bigfloat bay(bays * d1);
   const bigfloat baz(bazs * d1);
   const bigfloat cax(caxs * d1);
   const bigfloat cay(cays * d1);
   const bigfloat caz(cazs * d1);
   const bigfloat axs(ax * d1);
   const bigfloat ays(ay * d1);
   const bigfloat azs(az * d1);
   const bigfloat qax(l1x - axs);
   const bigfloat qay(l1y - ays);
   const bigfloat qaz(l1z - azs);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_IEEE(const genericPoint& q, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = inGabrielSphere_IEEE_interval(q, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_IEEE_bigfloat(q, ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline int inGabrielSphere_IIEE_interval(const genericPoint& q, const genericPoint& a, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !q.getIntervalLambda(l1x, l1y, l1z, d1)
   || !a.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number bxd2(bx * d2);
   const interval_number byd2(by * d2);
   const interval_number bzd2(bz * d2);
   const interval_number cxd2(cx * d2);
   const interval_number cyd2(cy * d2);
   const interval_number czd2(cz * d2);
   const interval_number l1x2(l1x * d2);
   const interval_number l1y2(l1y * d2);
   const interval_number l1z2(l1z * d2);
   const interval_number baxs(bxd2 - l2x);
   const interval_number bays(byd2 - l2y);
   const interval_number bazs(bzd2 - l2z);
   const interval_number caxs(cxd2 - l2x);
   const interval_number cays(cyd2 - l2y);
   const interval_number cazs(czd2 - l2z);
   const interval_number bax(baxs * d1);
   const interval_number bay(bays * d1);
   const interval_number baz(bazs * d1);
   const interval_number cax(caxs * d1);
   const interval_number cay(cays * d1);
   const interval_number caz(cazs * d1);
   const interval_number l2x1(l2x * d1);
   const interval_number l2y1(l2y * d1);
   const interval_number l2z1(l2z * d1);
   const interval_number qax(l1x2 - l2x1);
   const interval_number qay(l1y2 - l2y1);
   const interval_number qaz(l1z2 - l2z1);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_IIEE_bigfloat(const genericPoint& q, const genericPoint& a, bigfloat bx, bigfloat by, bigfloat bz, bigfloat cx, bigfloat cy, bigfloat cz)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   q.getBigfloatLambda(l1x, l1y, l1z, d1);
   a.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat bxd2(bx * d2);
   const bigfloat byd2(by * d2);
   const bigfloat bzd2(bz * d2);
   const bigfloat cxd2(cx * d2);
   const bigfloat cyd2(cy * d2);
   const bigfloat czd2(cz * d2);
   const bigfloat l1x2(l1x * d2);
   const bigfloat l1y2(l1y * d2);
   const bigfloat l1z2(l1z * d2);
   const bigfloat baxs(bxd2 - l2x);
   const bigfloat bays(byd2 - l2y);
   const bigfloat bazs(bzd2 - l2z);
   const bigfloat caxs(cxd2 - l2x);
   const bigfloat cays(cyd2 - l2y);
   const bigfloat cazs(czd2 - l2z);
   const bigfloat bax(baxs * d1);
   const bigfloat bay(bays * d1);
   const bigfloat baz(bazs * d1);
   const bigfloat cax(caxs * d1);
   const bigfloat cay(cays * d1);
   const bigfloat caz(cazs * d1);
   const bigfloat l2x1(l2x * d1);
   const bigfloat l2y1(l2y * d1);
   const bigfloat l2z1(l2z * d1);
   const bigfloat qax(l1x2 - l2x1);
   const bigfloat qay(l1y2 - l2y1);
   const bigfloat qaz(l1z2 - l2z1);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_IIEE(const genericPoint& q, const genericPoint& a, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = inGabrielSphere_IIEE_interval(q, a, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_IIEE_bigfloat(q, a, bx, by, bz, cx, cy, cz);
}

inline int inGabrielSphere_IIIE_interval(const genericPoint& q, const genericPoint& a, const genericPoint& b, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !q.getIntervalLambda(l1x, l1y, l1z, d1)
   || !a.getIntervalLambda(l2x, l2y, l2z, d2)
   || !b.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d12(d1 * d2);
   const interval_number d23(d2 * d3);
   const interval_number d13(d1 * d3);
   const interval_number d123(d12 * d3);
   const interval_number bxd2(l3x * d12);
   const interval_number byd2(l3y * d12);
   const interval_number bzd2(l3z * d12);
   const interval_number cxd2(cx * d123);
   const interval_number cyd2(cy * d123);
   const interval_number czd2(cz * d123);
   const interval_number qxd2(l1x * d23);
   const interval_number qyd2(l1y * d23);
   const interval_number qzd2(l1z * d23);
   const interval_number l2xm(l2x * d13);
   const interval_number l2ym(l2y * d13);
   const interval_number l2zm(l2z * d13);
   const interval_number bax(bxd2 - l2xm);
   const interval_number bay(byd2 - l2ym);
   const interval_number baz(bzd2 - l2zm);
   const interval_number cax(cxd2 - l2xm);
   const interval_number cay(cyd2 - l2ym);
   const interval_number caz(czd2 - l2zm);
   const interval_number qax(qxd2 - l2xm);
   const interval_number qay(qyd2 - l2ym);
   const interval_number qaz(qzd2 - l2zm);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_IIIE_bigfloat(const genericPoint& q, const genericPoint& a, const genericPoint& b, bigfloat cx, bigfloat cy, bigfloat cz)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   q.getBigfloatLambda(l1x, l1y, l1z, d1);
   a.getBigfloatLambda(l2x, l2y, l2z, d2);
   b.getBigfloatLambda(l3x, l3y, l3z, d3);
   const bigfloat d12(d1 * d2);
   const bigfloat d23(d2 * d3);
   const bigfloat d13(d1 * d3);
   const bigfloat d123(d12 * d3);
   const bigfloat bxd2(l3x * d12);
   const bigfloat byd2(l3y * d12);
   const bigfloat bzd2(l3z * d12);
   const bigfloat cxd2(cx * d123);
   const bigfloat cyd2(cy * d123);
   const bigfloat czd2(cz * d123);
   const bigfloat qxd2(l1x * d23);
   const bigfloat qyd2(l1y * d23);
   const bigfloat qzd2(l1z * d23);
   const bigfloat l2xm(l2x * d13);
   const bigfloat l2ym(l2y * d13);
   const bigfloat l2zm(l2z * d13);
   const bigfloat bax(bxd2 - l2xm);
   const bigfloat bay(byd2 - l2ym);
   const bigfloat baz(bzd2 - l2zm);
   const bigfloat cax(cxd2 - l2xm);
   const bigfloat cay(cyd2 - l2ym);
   const bigfloat caz(czd2 - l2zm);
   const bigfloat qax(qxd2 - l2xm);
   const bigfloat qay(qyd2 - l2ym);
   const bigfloat qaz(qzd2 - l2zm);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_IIIE(const genericPoint& q, const genericPoint& a, const genericPoint& b, double cx, double cy, double cz)
{
   int ret;
   ret = inGabrielSphere_IIIE_interval(q, a, b, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_IIIE_bigfloat(q, a, b, cx, cy, cz);
}

inline int inGabrielSphere_IIII_interval(const genericPoint& q, const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   if (
   !q.getIntervalLambda(l1x, l1y, l1z, d1)
   || !a.getIntervalLambda(l2x, l2y, l2z, d2)
   || !b.getIntervalLambda(l3x, l3y, l3z, d3)
   || !c.getIntervalLambda(l4x, l4y, l4z, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d23(d2 * d3);
   const interval_number d24(d2 * d4);
   const interval_number d34(d3 * d4);
   const interval_number d123(d1 * d23);
   const interval_number d124(d1 * d24);
   const interval_number d134(d1 * d34);
   const interval_number d234(d23 * d4);
   const interval_number bxd2(l3x * d124);
   const interval_number byd2(l3y * d124);
   const interval_number bzd2(l3z * d124);
   const interval_number cxd2(l4x * d123);
   const interval_number cyd2(l4y * d123);
   const interval_number czd2(l4z * d123);
   const interval_number qxd2(l1x * d234);
   const interval_number qyd2(l1y * d234);
   const interval_number qzd2(l1z * d234);
   const interval_number l2xm(l2x * d134);
   const interval_number l2ym(l2y * d134);
   const interval_number l2zm(l2z * d134);
   const interval_number bax(bxd2 - l2xm);
   const interval_number bay(byd2 - l2ym);
   const interval_number baz(bzd2 - l2zm);
   const interval_number cax(cxd2 - l2xm);
   const interval_number cay(cyd2 - l2ym);
   const interval_number caz(czd2 - l2zm);
   const interval_number qax(qxd2 - l2xm);
   const interval_number qay(qyd2 - l2ym);
   const interval_number qaz(qzd2 - l2zm);
   const interval_number cx1(bay * caz);
   const interval_number cx2(baz * cay);
   const interval_number crossbcx(cx1 - cx2);
   const interval_number cy1(baz * cax);
   const interval_number cy2(bax * caz);
   const interval_number crossbcy(cy1 - cy2);
   const interval_number cz1(bax * cay);
   const interval_number cz2(bay * cax);
   const interval_number crossbcz(cz1 - cz2);
   const interval_number ba2x(bax * bax);
   const interval_number ba2y(bay * bay);
   const interval_number ba2z(baz * baz);
   const interval_number ba2t(ba2x + ba2y);
   const interval_number ba2(ba2t + ba2z);
   const interval_number ca2x(cax * cax);
   const interval_number ca2y(cay * cay);
   const interval_number ca2z(caz * caz);
   const interval_number ca2t(ca2x + ca2y);
   const interval_number ca2(ca2t + ca2z);
   const interval_number calx(cax * ba2);
   const interval_number caly(cay * ba2);
   const interval_number calz(caz * ba2);
   const interval_number balx(bax * ca2);
   const interval_number baly(bay * ca2);
   const interval_number balz(baz * ca2);
   const interval_number abcx(calx - balx);
   const interval_number abcy(caly - baly);
   const interval_number abcz(calz - balz);
   const interval_number kx1(abcy * crossbcz);
   const interval_number kx2(abcz * crossbcy);
   const interval_number ccax(kx1 - kx2);
   const interval_number ky1(abcz * crossbcx);
   const interval_number ky2(abcx * crossbcz);
   const interval_number ccay(ky1 - ky2);
   const interval_number kz1(abcx * crossbcy);
   const interval_number kz2(abcy * crossbcx);
   const interval_number ccaz(kz1 - kz2);
   const interval_number cr2x(crossbcx * crossbcx);
   const interval_number cr2y(crossbcy * crossbcy);
   const interval_number cr2z(crossbcz * crossbcz);
   const interval_number cr2t(cr2x + cr2y);
   const interval_number c2(cr2t + cr2z);
   const interval_number c22(c2 * 2);
   const interval_number qa1x(qax * c22);
   const interval_number qa1y(qay * c22);
   const interval_number qa1z(qaz * c22);
   const interval_number qa2x(qa1x - ccax);
   const interval_number qa2y(qa1y - ccay);
   const interval_number qa2z(qa1z - ccaz);
   const interval_number r1x(qa2x * qa2x);
   const interval_number r1y(qa2y * qa2y);
   const interval_number r1z(qa2z * qa2z);
   const interval_number r1t(r1x + r1y);
   const interval_number r1(r1t + r1z);
   const interval_number r2x(ccax * ccax);
   const interval_number r2y(ccay * ccay);
   const interval_number r2z(ccaz * ccaz);
   const interval_number r2t(r2x + r2y);
   const interval_number r2(r2t + r2z);
   const interval_number ret(r1 - r2);
   setFPUModeToRoundNEAR();

   if (!ret.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ret.sign();
}

inline int inGabrielSphere_IIII_bigfloat(const genericPoint& q, const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   q.getBigfloatLambda(l1x, l1y, l1z, d1);
   a.getBigfloatLambda(l2x, l2y, l2z, d2);
   b.getBigfloatLambda(l3x, l3y, l3z, d3);
   c.getBigfloatLambda(l4x, l4y, l4z, d4);
   const bigfloat d23(d2 * d3);
   const bigfloat d24(d2 * d4);
   const bigfloat d34(d3 * d4);
   const bigfloat d123(d1 * d23);
   const bigfloat d124(d1 * d24);
   const bigfloat d134(d1 * d34);
   const bigfloat d234(d23 * d4);
   const bigfloat bxd2(l3x * d124);
   const bigfloat byd2(l3y * d124);
   const bigfloat bzd2(l3z * d124);
   const bigfloat cxd2(l4x * d123);
   const bigfloat cyd2(l4y * d123);
   const bigfloat czd2(l4z * d123);
   const bigfloat qxd2(l1x * d234);
   const bigfloat qyd2(l1y * d234);
   const bigfloat qzd2(l1z * d234);
   const bigfloat l2xm(l2x * d134);
   const bigfloat l2ym(l2y * d134);
   const bigfloat l2zm(l2z * d134);
   const bigfloat bax(bxd2 - l2xm);
   const bigfloat bay(byd2 - l2ym);
   const bigfloat baz(bzd2 - l2zm);
   const bigfloat cax(cxd2 - l2xm);
   const bigfloat cay(cyd2 - l2ym);
   const bigfloat caz(czd2 - l2zm);
   const bigfloat qax(qxd2 - l2xm);
   const bigfloat qay(qyd2 - l2ym);
   const bigfloat qaz(qzd2 - l2zm);
   const bigfloat cx1(bay * caz);
   const bigfloat cx2(baz * cay);
   const bigfloat crossbcx(cx1 - cx2);
   const bigfloat cy1(baz * cax);
   const bigfloat cy2(bax * caz);
   const bigfloat crossbcy(cy1 - cy2);
   const bigfloat cz1(bax * cay);
   const bigfloat cz2(bay * cax);
   const bigfloat crossbcz(cz1 - cz2);
   const bigfloat ba2x(bax * bax);
   const bigfloat ba2y(bay * bay);
   const bigfloat ba2z(baz * baz);
   const bigfloat ba2t(ba2x + ba2y);
   const bigfloat ba2(ba2t + ba2z);
   const bigfloat ca2x(cax * cax);
   const bigfloat ca2y(cay * cay);
   const bigfloat ca2z(caz * caz);
   const bigfloat ca2t(ca2x + ca2y);
   const bigfloat ca2(ca2t + ca2z);
   const bigfloat calx(cax * ba2);
   const bigfloat caly(cay * ba2);
   const bigfloat calz(caz * ba2);
   const bigfloat balx(bax * ca2);
   const bigfloat baly(bay * ca2);
   const bigfloat balz(baz * ca2);
   const bigfloat abcx(calx - balx);
   const bigfloat abcy(caly - baly);
   const bigfloat abcz(calz - balz);
   const bigfloat kx1(abcy * crossbcz);
   const bigfloat kx2(abcz * crossbcy);
   const bigfloat ccax(kx1 - kx2);
   const bigfloat ky1(abcz * crossbcx);
   const bigfloat ky2(abcx * crossbcz);
   const bigfloat ccay(ky1 - ky2);
   const bigfloat kz1(abcx * crossbcy);
   const bigfloat kz2(abcy * crossbcx);
   const bigfloat ccaz(kz1 - kz2);
   const bigfloat cr2x(crossbcx * crossbcx);
   const bigfloat cr2y(crossbcy * crossbcy);
   const bigfloat cr2z(crossbcz * crossbcz);
   const bigfloat cr2t(cr2x + cr2y);
   const bigfloat c2(cr2t + cr2z);
   const bigfloat c22(c2 * 2);
   const bigfloat qa1x(qax * c22);
   const bigfloat qa1y(qay * c22);
   const bigfloat qa1z(qaz * c22);
   const bigfloat qa2x(qa1x - ccax);
   const bigfloat qa2y(qa1y - ccay);
   const bigfloat qa2z(qa1z - ccaz);
   const bigfloat r1x(qa2x * qa2x);
   const bigfloat r1y(qa2y * qa2y);
   const bigfloat r1z(qa2z * qa2z);
   const bigfloat r1t(r1x + r1y);
   const bigfloat r1(r1t + r1z);
   const bigfloat r2x(ccax * ccax);
   const bigfloat r2y(ccay * ccay);
   const bigfloat r2z(ccaz * ccaz);
   const bigfloat r2t(r2x + r2y);
   const bigfloat r2(r2t + r2z);
   const bigfloat ret(r1 - r2);
   return sgn(ret);
}

inline int inGabrielSphere_IIII(const genericPoint& q, const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
   int ret;
   ret = inGabrielSphere_IIII_interval(q, a, b, c);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inGabrielSphere_IIII_bigfloat(q, a, b, c);
}

inline int inSphere_IEEEE_interval(const genericPoint& p1, interval_number pbx, interval_number pby, interval_number pbz, interval_number pcx, interval_number pcy, interval_number pcz, interval_number pdx, interval_number pdy, interval_number pdz, interval_number pex, interval_number pey, interval_number pez)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pexd(pex * d1);
   const interval_number peyd(pey * d1);
   const interval_number pezd(pez * d1);
   const interval_number aex(l1x - pexd);
   const interval_number aey(l1y - peyd);
   const interval_number aez(l1z - pezd);
   const interval_number bex(pbx - pex);
   const interval_number bey(pby - pey);
   const interval_number bez(pbz - pez);
   const interval_number cex(pcx - pex);
   const interval_number cey(pcy - pey);
   const interval_number cez(pcz - pez);
   const interval_number dex(pdx - pex);
   const interval_number dey(pdy - pey);
   const interval_number dez(pdz - pez);
   const interval_number aexbey(aex * bey);
   const interval_number bexaey(bex * aey);
   const interval_number ab(aexbey - bexaey);
   const interval_number bexcey(bex * cey);
   const interval_number cexbey(cex * bey);
   const interval_number bc(bexcey - cexbey);
   const interval_number cexdey(cex * dey);
   const interval_number dexcey(dex * cey);
   const interval_number cd(cexdey - dexcey);
   const interval_number dexaey(dex * aey);
   const interval_number aexdey(aex * dey);
   const interval_number da(dexaey - aexdey);
   const interval_number aexcey(aex * cey);
   const interval_number cexaey(cex * aey);
   const interval_number ac(aexcey - cexaey);
   const interval_number bexdey(bex * dey);
   const interval_number dexbey(dex * bey);
   const interval_number bd(bexdey - dexbey);
   const interval_number abc1(aez * bc);
   const interval_number abc2(bez * ac);
   const interval_number abc3(cez * ab);
   const interval_number abc4(abc1 + abc3);
   const interval_number abc(abc4 - abc2);
   const interval_number bcd1(bez * cd);
   const interval_number bcd2(cez * bd);
   const interval_number bcd3(dez * bc);
   const interval_number bcd4(bcd1 + bcd3);
   const interval_number bcd(bcd4 - bcd2);
   const interval_number cda1(cez * da);
   const interval_number cda2(dez * ac);
   const interval_number cda3(aez * cd);
   const interval_number cda4(cda1 + cda3);
   const interval_number cda(cda4 + cda2);
   const interval_number dab1(dez * ab);
   const interval_number dab2(aez * bd);
   const interval_number dab3(bez * da);
   const interval_number dab4(dab1 + dab3);
   const interval_number dab(dab4 + dab2);
   const interval_number al1(aex * aex);
   const interval_number al2(aey * aey);
   const interval_number al3(aez * aez);
   const interval_number al4(al1 + al2);
   const interval_number alift(al4 + al3);
   const interval_number bl1(bex * bex);
   const interval_number bl2(bey * bey);
   const interval_number bl3(bez * bez);
   const interval_number bl4(bl1 + bl2);
   const interval_number blift(bl4 + bl3);
   const interval_number cl1(cex * cex);
   const interval_number cl2(cey * cey);
   const interval_number cl3(cez * cez);
   const interval_number cl4(cl1 + cl2);
   const interval_number clift(cl4 + cl3);
   const interval_number dl1(dex * dex);
   const interval_number dl2(dey * dey);
   const interval_number dl3(dez * dez);
   const interval_number dl4(dl1 + dl2);
   const interval_number dlift(dl4 + dl3);
   const interval_number ds1(dlift * abc);
   const interval_number ds2(clift * dab);
   const interval_number dlp(ds2 - ds1);
   const interval_number dl(dlp * d1);
   const interval_number dr1p(blift * cda);
   const interval_number dr1(dr1p * d1);
   const interval_number dr2(alift * bcd);
   const interval_number dr(dr2 - dr1);
   const interval_number det(dl + dr);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int inSphere_IEEEE_bigfloat(const genericPoint& p1, bigfloat pbx, bigfloat pby, bigfloat pbz, bigfloat pcx, bigfloat pcy, bigfloat pcz, bigfloat pdx, bigfloat pdy, bigfloat pdz, bigfloat pex, bigfloat pey, bigfloat pez)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat pexd(pex * d1);
   const bigfloat peyd(pey * d1);
   const bigfloat pezd(pez * d1);
   const bigfloat aex(l1x - pexd);
   const bigfloat aey(l1y - peyd);
   const bigfloat aez(l1z - pezd);
   const bigfloat bex(pbx - pex);
   const bigfloat bey(pby - pey);
   const bigfloat bez(pbz - pez);
   const bigfloat cex(pcx - pex);
   const bigfloat cey(pcy - pey);
   const bigfloat cez(pcz - pez);
   const bigfloat dex(pdx - pex);
   const bigfloat dey(pdy - pey);
   const bigfloat dez(pdz - pez);
   const bigfloat aexbey(aex * bey);
   const bigfloat bexaey(bex * aey);
   const bigfloat ab(aexbey - bexaey);
   const bigfloat bexcey(bex * cey);
   const bigfloat cexbey(cex * bey);
   const bigfloat bc(bexcey - cexbey);
   const bigfloat cexdey(cex * dey);
   const bigfloat dexcey(dex * cey);
   const bigfloat cd(cexdey - dexcey);
   const bigfloat dexaey(dex * aey);
   const bigfloat aexdey(aex * dey);
   const bigfloat da(dexaey - aexdey);
   const bigfloat aexcey(aex * cey);
   const bigfloat cexaey(cex * aey);
   const bigfloat ac(aexcey - cexaey);
   const bigfloat bexdey(bex * dey);
   const bigfloat dexbey(dex * bey);
   const bigfloat bd(bexdey - dexbey);
   const bigfloat abc1(aez * bc);
   const bigfloat abc2(bez * ac);
   const bigfloat abc3(cez * ab);
   const bigfloat abc4(abc1 + abc3);
   const bigfloat abc(abc4 - abc2);
   const bigfloat bcd1(bez * cd);
   const bigfloat bcd2(cez * bd);
   const bigfloat bcd3(dez * bc);
   const bigfloat bcd4(bcd1 + bcd3);
   const bigfloat bcd(bcd4 - bcd2);
   const bigfloat cda1(cez * da);
   const bigfloat cda2(dez * ac);
   const bigfloat cda3(aez * cd);
   const bigfloat cda4(cda1 + cda3);
   const bigfloat cda(cda4 + cda2);
   const bigfloat dab1(dez * ab);
   const bigfloat dab2(aez * bd);
   const bigfloat dab3(bez * da);
   const bigfloat dab4(dab1 + dab3);
   const bigfloat dab(dab4 + dab2);
   const bigfloat al1(aex * aex);
   const bigfloat al2(aey * aey);
   const bigfloat al3(aez * aez);
   const bigfloat al4(al1 + al2);
   const bigfloat alift(al4 + al3);
   const bigfloat bl1(bex * bex);
   const bigfloat bl2(bey * bey);
   const bigfloat bl3(bez * bez);
   const bigfloat bl4(bl1 + bl2);
   const bigfloat blift(bl4 + bl3);
   const bigfloat cl1(cex * cex);
   const bigfloat cl2(cey * cey);
   const bigfloat cl3(cez * cez);
   const bigfloat cl4(cl1 + cl2);
   const bigfloat clift(cl4 + cl3);
   const bigfloat dl1(dex * dex);
   const bigfloat dl2(dey * dey);
   const bigfloat dl3(dez * dez);
   const bigfloat dl4(dl1 + dl2);
   const bigfloat dlift(dl4 + dl3);
   const bigfloat ds1(dlift * abc);
   const bigfloat ds2(clift * dab);
   const bigfloat dlp(ds2 - ds1);
   const bigfloat dl(dlp * d1);
   const bigfloat dr1p(blift * cda);
   const bigfloat dr1(dr1p * d1);
   const bigfloat dr2(alift * bcd);
   const bigfloat dr(dr2 - dr1);
   const bigfloat det(dl + dr);
   return sgn(det);
}

inline int inSphere_IEEEE(const genericPoint& p1, double pbx, double pby, double pbz, double pcx, double pcy, double pcz, double pdx, double pdy, double pdz, double pex, double pey, double pez)
{
   int ret;
   ret = inSphere_IEEEE_interval(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_IEEEE_bigfloat(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}

inline int inSphere_IIEEE_interval(const genericPoint& p1, const genericPoint& p2, interval_number pcx, interval_number pcy, interval_number pcz, interval_number pdx, interval_number pdy, interval_number pdz, interval_number pex, interval_number pey, interval_number pez)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pexd(pex * d1);
   const interval_number peyd(pey * d1);
   const interval_number pezd(pez * d1);
   const interval_number aex(l1x - pexd);
   const interval_number aey(l1y - peyd);
   const interval_number aez(l1z - pezd);
   const interval_number pexd2(pex * d2);
   const interval_number peyd2(pey * d2);
   const interval_number pezd2(pez * d2);
   const interval_number bex(l2x - pexd2);
   const interval_number bey(l2y - peyd2);
   const interval_number bez(l2z - pezd2);
   const interval_number cex(pcx - pex);
   const interval_number cey(pcy - pey);
   const interval_number cez(pcz - pez);
   const interval_number dex(pdx - pex);
   const interval_number dey(pdy - pey);
   const interval_number dez(pdz - pez);
   const interval_number aexbey(aex * bey);
   const interval_number bexaey(bex * aey);
   const interval_number ab(aexbey - bexaey);
   const interval_number bexcey(bex * cey);
   const interval_number cexbey(cex * bey);
   const interval_number bc(bexcey - cexbey);
   const interval_number cexdey(cex * dey);
   const interval_number dexcey(dex * cey);
   const interval_number cd(cexdey - dexcey);
   const interval_number dexaey(dex * aey);
   const interval_number aexdey(aex * dey);
   const interval_number da(dexaey - aexdey);
   const interval_number aexcey(aex * cey);
   const interval_number cexaey(cex * aey);
   const interval_number ac(aexcey - cexaey);
   const interval_number bexdey(bex * dey);
   const interval_number dexbey(dex * bey);
   const interval_number bd(bexdey - dexbey);
   const interval_number abc1(aez * bc);
   const interval_number abc2(bez * ac);
   const interval_number abc3(cez * ab);
   const interval_number abc4(abc1 + abc3);
   const interval_number abc(abc4 - abc2);
   const interval_number bcd1(bez * cd);
   const interval_number bcd2(cez * bd);
   const interval_number bcd3(dez * bc);
   const interval_number bcd4(bcd1 + bcd3);
   const interval_number bcd(bcd4 - bcd2);
   const interval_number cda1(cez * da);
   const interval_number cda2(dez * ac);
   const interval_number cda3(aez * cd);
   const interval_number cda4(cda1 + cda3);
   const interval_number cda(cda4 + cda2);
   const interval_number dab1(dez * ab);
   const interval_number dab2(aez * bd);
   const interval_number dab3(bez * da);
   const interval_number dab4(dab1 + dab3);
   const interval_number dab(dab4 + dab2);
   const interval_number al1(aex * aex);
   const interval_number al2(aey * aey);
   const interval_number al3(aez * aez);
   const interval_number al4(al1 + al2);
   const interval_number alift(al4 + al3);
   const interval_number bl1(bex * bex);
   const interval_number bl2(bey * bey);
   const interval_number bl3(bez * bez);
   const interval_number bl4(bl1 + bl2);
   const interval_number blift(bl4 + bl3);
   const interval_number cl1(cex * cex);
   const interval_number cl2(cey * cey);
   const interval_number cl3(cez * cez);
   const interval_number cl4(cl1 + cl2);
   const interval_number clift(cl4 + cl3);
   const interval_number dl1(dex * dex);
   const interval_number dl2(dey * dey);
   const interval_number dl3(dez * dez);
   const interval_number dl4(dl1 + dl2);
   const interval_number dlift(dl4 + dl3);
   const interval_number ds1(dlift * abc);
   const interval_number ds2(clift * dab);
   const interval_number dl(ds2 - ds1);
   const interval_number dll(dl * d1);
   const interval_number dlll(dll * d2);
   const interval_number dr1(blift * cda);
   const interval_number dr12(dr1 * d1);
   const interval_number dr2(alift * bcd);
   const interval_number dr22(dr2 * d2);
   const interval_number dr(dr22 - dr12);
   const interval_number det(dlll + dr);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int inSphere_IIEEE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat pcx, bigfloat pcy, bigfloat pcz, bigfloat pdx, bigfloat pdy, bigfloat pdz, bigfloat pex, bigfloat pey, bigfloat pez)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat pexd(pex * d1);
   const bigfloat peyd(pey * d1);
   const bigfloat pezd(pez * d1);
   const bigfloat aex(l1x - pexd);
   const bigfloat aey(l1y - peyd);
   const bigfloat aez(l1z - pezd);
   const bigfloat pexd2(pex * d2);
   const bigfloat peyd2(pey * d2);
   const bigfloat pezd2(pez * d2);
   const bigfloat bex(l2x - pexd2);
   const bigfloat bey(l2y - peyd2);
   const bigfloat bez(l2z - pezd2);
   const bigfloat cex(pcx - pex);
   const bigfloat cey(pcy - pey);
   const bigfloat cez(pcz - pez);
   const bigfloat dex(pdx - pex);
   const bigfloat dey(pdy - pey);
   const bigfloat dez(pdz - pez);
   const bigfloat aexbey(aex * bey);
   const bigfloat bexaey(bex * aey);
   const bigfloat ab(aexbey - bexaey);
   const bigfloat bexcey(bex * cey);
   const bigfloat cexbey(cex * bey);
   const bigfloat bc(bexcey - cexbey);
   const bigfloat cexdey(cex * dey);
   const bigfloat dexcey(dex * cey);
   const bigfloat cd(cexdey - dexcey);
   const bigfloat dexaey(dex * aey);
   const bigfloat aexdey(aex * dey);
   const bigfloat da(dexaey - aexdey);
   const bigfloat aexcey(aex * cey);
   const bigfloat cexaey(cex * aey);
   const bigfloat ac(aexcey - cexaey);
   const bigfloat bexdey(bex * dey);
   const bigfloat dexbey(dex * bey);
   const bigfloat bd(bexdey - dexbey);
   const bigfloat abc1(aez * bc);
   const bigfloat abc2(bez * ac);
   const bigfloat abc3(cez * ab);
   const bigfloat abc4(abc1 + abc3);
   const bigfloat abc(abc4 - abc2);
   const bigfloat bcd1(bez * cd);
   const bigfloat bcd2(cez * bd);
   const bigfloat bcd3(dez * bc);
   const bigfloat bcd4(bcd1 + bcd3);
   const bigfloat bcd(bcd4 - bcd2);
   const bigfloat cda1(cez * da);
   const bigfloat cda2(dez * ac);
   const bigfloat cda3(aez * cd);
   const bigfloat cda4(cda1 + cda3);
   const bigfloat cda(cda4 + cda2);
   const bigfloat dab1(dez * ab);
   const bigfloat dab2(aez * bd);
   const bigfloat dab3(bez * da);
   const bigfloat dab4(dab1 + dab3);
   const bigfloat dab(dab4 + dab2);
   const bigfloat al1(aex * aex);
   const bigfloat al2(aey * aey);
   const bigfloat al3(aez * aez);
   const bigfloat al4(al1 + al2);
   const bigfloat alift(al4 + al3);
   const bigfloat bl1(bex * bex);
   const bigfloat bl2(bey * bey);
   const bigfloat bl3(bez * bez);
   const bigfloat bl4(bl1 + bl2);
   const bigfloat blift(bl4 + bl3);
   const bigfloat cl1(cex * cex);
   const bigfloat cl2(cey * cey);
   const bigfloat cl3(cez * cez);
   const bigfloat cl4(cl1 + cl2);
   const bigfloat clift(cl4 + cl3);
   const bigfloat dl1(dex * dex);
   const bigfloat dl2(dey * dey);
   const bigfloat dl3(dez * dez);
   const bigfloat dl4(dl1 + dl2);
   const bigfloat dlift(dl4 + dl3);
   const bigfloat ds1(dlift * abc);
   const bigfloat ds2(clift * dab);
   const bigfloat dl(ds2 - ds1);
   const bigfloat dll(dl * d1);
   const bigfloat dlll(dll * d2);
   const bigfloat dr1(blift * cda);
   const bigfloat dr12(dr1 * d1);
   const bigfloat dr2(alift * bcd);
   const bigfloat dr22(dr2 * d2);
   const bigfloat dr(dr22 - dr12);
   const bigfloat det(dlll + dr);
   return sgn(det);
}

inline int inSphere_IIEEE(const genericPoint& p1, const genericPoint& p2, double pcx, double pcy, double pcz, double pdx, double pdy, double pdz, double pex, double pey, double pez)
{
   int ret;
   ret = inSphere_IIEEE_interval(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_IIEEE_bigfloat(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez);
}

inline int inSphere_IIIEE_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, interval_number pdx, interval_number pdy, interval_number pdz, interval_number pex, interval_number pey, interval_number pez)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pexd(pex * d1);
   const interval_number peyd(pey * d1);
   const interval_number pezd(pez * d1);
   const interval_number aex(l1x - pexd);
   const interval_number aey(l1y - peyd);
   const interval_number aez(l1z - pezd);
   const interval_number pexd2(pex * d2);
   const interval_number peyd2(pey * d2);
   const interval_number pezd2(pez * d2);
   const interval_number bex(l2x - pexd2);
   const interval_number bey(l2y - peyd2);
   const interval_number bez(l2z - pezd2);
   const interval_number pexd3(pex * d3);
   const interval_number peyd3(pey * d3);
   const interval_number pezd3(pez * d3);
   const interval_number cex(l3x - pexd3);
   const interval_number cey(l3y - peyd3);
   const interval_number cez(l3z - pezd3);
   const interval_number dex(pdx - pex);
   const interval_number dey(pdy - pey);
   const interval_number dez(pdz - pez);
   const interval_number aexbey(aex * bey);
   const interval_number bexaey(bex * aey);
   const interval_number ab(aexbey - bexaey);
   const interval_number bexcey(bex * cey);
   const interval_number cexbey(cex * bey);
   const interval_number bc(bexcey - cexbey);
   const interval_number cexdey(cex * dey);
   const interval_number dexcey(dex * cey);
   const interval_number cd(cexdey - dexcey);
   const interval_number dexaey(dex * aey);
   const interval_number aexdey(aex * dey);
   const interval_number da(dexaey - aexdey);
   const interval_number aexcey(aex * cey);
   const interval_number cexaey(cex * aey);
   const interval_number ac(aexcey - cexaey);
   const interval_number bexdey(bex * dey);
   const interval_number dexbey(dex * bey);
   const interval_number bd(bexdey - dexbey);
   const interval_number abc1(aez * bc);
   const interval_number abc2(bez * ac);
   const interval_number abc3(cez * ab);
   const interval_number abc4(abc1 + abc3);
   const interval_number abc(abc4 - abc2);
   const interval_number bcd1(bez * cd);
   const interval_number bcd2(cez * bd);
   const interval_number bcd3(dez * bc);
   const interval_number bcd4(bcd1 + bcd3);
   const interval_number bcd(bcd4 - bcd2);
   const interval_number cda1(cez * da);
   const interval_number cda2(dez * ac);
   const interval_number cda3(aez * cd);
   const interval_number cda4(cda1 + cda3);
   const interval_number cda(cda4 + cda2);
   const interval_number dab1(dez * ab);
   const interval_number dab2(aez * bd);
   const interval_number dab3(bez * da);
   const interval_number dab4(dab1 + dab3);
   const interval_number dab(dab4 + dab2);
   const interval_number al1(aex * aex);
   const interval_number al2(aey * aey);
   const interval_number al3(aez * aez);
   const interval_number al4(al1 + al2);
   const interval_number alift(al4 + al3);
   const interval_number bl1(bex * bex);
   const interval_number bl2(bey * bey);
   const interval_number bl3(bez * bez);
   const interval_number bl4(bl1 + bl2);
   const interval_number blift(bl4 + bl3);
   const interval_number cl1(cex * cex);
   const interval_number cl2(cey * cey);
   const interval_number cl3(cez * cez);
   const interval_number cl4(cl1 + cl2);
   const interval_number clift(cl4 + cl3);
   const interval_number dl1(dex * dex);
   const interval_number dl2(dey * dey);
   const interval_number dl3(dez * dez);
   const interval_number dl4(dl1 + dl2);
   const interval_number dlift(dl4 + dl3);
   const interval_number ds1(dlift * abc);
   const interval_number ds1n(ds1 * d3);
   const interval_number ds2(clift * dab);
   const interval_number dl(ds2 - ds1n);
   const interval_number dlm(dl * d1);
   const interval_number dln(dlm * d2);
   const interval_number dr1(blift * cda);
   const interval_number dr1n(dr1 * d1);
   const interval_number dr2(alift * bcd);
   const interval_number dr2n(dr2 * d2);
   const interval_number dr(dr2n - dr1n);
   const interval_number drn(dr * d3);
   const interval_number det(dln + drn);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int inSphere_IIIEE_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, bigfloat pdx, bigfloat pdy, bigfloat pdz, bigfloat pex, bigfloat pey, bigfloat pez)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   const bigfloat pexd(pex * d1);
   const bigfloat peyd(pey * d1);
   const bigfloat pezd(pez * d1);
   const bigfloat aex(l1x - pexd);
   const bigfloat aey(l1y - peyd);
   const bigfloat aez(l1z - pezd);
   const bigfloat pexd2(pex * d2);
   const bigfloat peyd2(pey * d2);
   const bigfloat pezd2(pez * d2);
   const bigfloat bex(l2x - pexd2);
   const bigfloat bey(l2y - peyd2);
   const bigfloat bez(l2z - pezd2);
   const bigfloat pexd3(pex * d3);
   const bigfloat peyd3(pey * d3);
   const bigfloat pezd3(pez * d3);
   const bigfloat cex(l3x - pexd3);
   const bigfloat cey(l3y - peyd3);
   const bigfloat cez(l3z - pezd3);
   const bigfloat dex(pdx - pex);
   const bigfloat dey(pdy - pey);
   const bigfloat dez(pdz - pez);
   const bigfloat aexbey(aex * bey);
   const bigfloat bexaey(bex * aey);
   const bigfloat ab(aexbey - bexaey);
   const bigfloat bexcey(bex * cey);
   const bigfloat cexbey(cex * bey);
   const bigfloat bc(bexcey - cexbey);
   const bigfloat cexdey(cex * dey);
   const bigfloat dexcey(dex * cey);
   const bigfloat cd(cexdey - dexcey);
   const bigfloat dexaey(dex * aey);
   const bigfloat aexdey(aex * dey);
   const bigfloat da(dexaey - aexdey);
   const bigfloat aexcey(aex * cey);
   const bigfloat cexaey(cex * aey);
   const bigfloat ac(aexcey - cexaey);
   const bigfloat bexdey(bex * dey);
   const bigfloat dexbey(dex * bey);
   const bigfloat bd(bexdey - dexbey);
   const bigfloat abc1(aez * bc);
   const bigfloat abc2(bez * ac);
   const bigfloat abc3(cez * ab);
   const bigfloat abc4(abc1 + abc3);
   const bigfloat abc(abc4 - abc2);
   const bigfloat bcd1(bez * cd);
   const bigfloat bcd2(cez * bd);
   const bigfloat bcd3(dez * bc);
   const bigfloat bcd4(bcd1 + bcd3);
   const bigfloat bcd(bcd4 - bcd2);
   const bigfloat cda1(cez * da);
   const bigfloat cda2(dez * ac);
   const bigfloat cda3(aez * cd);
   const bigfloat cda4(cda1 + cda3);
   const bigfloat cda(cda4 + cda2);
   const bigfloat dab1(dez * ab);
   const bigfloat dab2(aez * bd);
   const bigfloat dab3(bez * da);
   const bigfloat dab4(dab1 + dab3);
   const bigfloat dab(dab4 + dab2);
   const bigfloat al1(aex * aex);
   const bigfloat al2(aey * aey);
   const bigfloat al3(aez * aez);
   const bigfloat al4(al1 + al2);
   const bigfloat alift(al4 + al3);
   const bigfloat bl1(bex * bex);
   const bigfloat bl2(bey * bey);
   const bigfloat bl3(bez * bez);
   const bigfloat bl4(bl1 + bl2);
   const bigfloat blift(bl4 + bl3);
   const bigfloat cl1(cex * cex);
   const bigfloat cl2(cey * cey);
   const bigfloat cl3(cez * cez);
   const bigfloat cl4(cl1 + cl2);
   const bigfloat clift(cl4 + cl3);
   const bigfloat dl1(dex * dex);
   const bigfloat dl2(dey * dey);
   const bigfloat dl3(dez * dez);
   const bigfloat dl4(dl1 + dl2);
   const bigfloat dlift(dl4 + dl3);
   const bigfloat ds1(dlift * abc);
   const bigfloat ds1n(ds1 * d3);
   const bigfloat ds2(clift * dab);
   const bigfloat dl(ds2 - ds1n);
   const bigfloat dlm(dl * d1);
   const bigfloat dln(dlm * d2);
   const bigfloat dr1(blift * cda);
   const bigfloat dr1n(dr1 * d1);
   const bigfloat dr2(alift * bcd);
   const bigfloat dr2n(dr2 * d2);
   const bigfloat dr(dr2n - dr1n);
   const bigfloat drn(dr * d3);
   const bigfloat det(dln + drn);
   return sgn(det);
}

inline int inSphere_IIIEE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, double pdx, double pdy, double pdz, double pex, double pey, double pez)
{
   int ret;
   ret = inSphere_IIIEE_interval(p1, p2, p3, pdx, pdy, pdz, pex, pey, pez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_IIIEE_bigfloat(p1, p2, p3, pdx, pdy, pdz, pex, pey, pez);
}

inline int inSphere_IIIIE_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, interval_number pex, interval_number pey, interval_number pez)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   || !p4.getIntervalLambda(l4x, l4y, l4z, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pexd(pex * d1);
   const interval_number peyd(pey * d1);
   const interval_number pezd(pez * d1);
   const interval_number aex(l1x - pexd);
   const interval_number aey(l1y - peyd);
   const interval_number aez(l1z - pezd);
   const interval_number pexd2(pex * d2);
   const interval_number peyd2(pey * d2);
   const interval_number pezd2(pez * d2);
   const interval_number bex(l2x - pexd2);
   const interval_number bey(l2y - peyd2);
   const interval_number bez(l2z - pezd2);
   const interval_number pexd3(pex * d3);
   const interval_number peyd3(pey * d3);
   const interval_number pezd3(pez * d3);
   const interval_number cex(l3x - pexd3);
   const interval_number cey(l3y - peyd3);
   const interval_number cez(l3z - pezd3);
   const interval_number pexd4(pex * d4);
   const interval_number peyd4(pey * d4);
   const interval_number pezd4(pez * d4);
   const interval_number dex(l4x - pexd4);
   const interval_number dey(l4y - peyd4);
   const interval_number dez(l4z - pezd4);
   const interval_number aexbey(aex * bey);
   const interval_number bexaey(bex * aey);
   const interval_number ab(aexbey - bexaey);
   const interval_number bexcey(bex * cey);
   const interval_number cexbey(cex * bey);
   const interval_number bc(bexcey - cexbey);
   const interval_number cexdey(cex * dey);
   const interval_number dexcey(dex * cey);
   const interval_number cd(cexdey - dexcey);
   const interval_number dexaey(dex * aey);
   const interval_number aexdey(aex * dey);
   const interval_number da(dexaey - aexdey);
   const interval_number aexcey(aex * cey);
   const interval_number cexaey(cex * aey);
   const interval_number ac(aexcey - cexaey);
   const interval_number bexdey(bex * dey);
   const interval_number dexbey(dex * bey);
   const interval_number bd(bexdey - dexbey);
   const interval_number abc1(aez * bc);
   const interval_number abc2(bez * ac);
   const interval_number abc3(cez * ab);
   const interval_number abc4(abc1 + abc3);
   const interval_number abc(abc4 - abc2);
   const interval_number bcd1(bez * cd);
   const interval_number bcd2(cez * bd);
   const interval_number bcd3(dez * bc);
   const interval_number bcd4(bcd1 + bcd3);
   const interval_number bcd(bcd4 - bcd2);
   const interval_number cda1(cez * da);
   const interval_number cda2(dez * ac);
   const interval_number cda3(aez * cd);
   const interval_number cda4(cda1 + cda3);
   const interval_number cda(cda4 + cda2);
   const interval_number dab1(dez * ab);
   const interval_number dab2(aez * bd);
   const interval_number dab3(bez * da);
   const interval_number dab4(dab1 + dab3);
   const interval_number dab(dab4 + dab2);
   const interval_number al1(aex * aex);
   const interval_number al2(aey * aey);
   const interval_number al3(aez * aez);
   const interval_number al4(al1 + al2);
   const interval_number alift(al4 + al3);
   const interval_number bl1(bex * bex);
   const interval_number bl2(bey * bey);
   const interval_number bl3(bez * bez);
   const interval_number bl4(bl1 + bl2);
   const interval_number blift(bl4 + bl3);
   const interval_number cl1(cex * cex);
   const interval_number cl2(cey * cey);
   const interval_number cl3(cez * cez);
   const interval_number cl4(cl1 + cl2);
   const interval_number clift(cl4 + cl3);
   const interval_number dl1(dex * dex);
   const interval_number dl2(dey * dey);
   const interval_number dl3(dez * dez);
   const interval_number dl4(dl1 + dl2);
   const interval_number dlift(dl4 + dl3);
   const interval_number ds1(dlift * abc);
   const interval_number ds12(ds1 * d3);
   const interval_number ds2(clift * dab);
   const interval_number ds22(ds2 * d4);
   const interval_number dl(ds22 - ds12);
   const interval_number dlx1(dl * d1);
   const interval_number dlx2(dlx1 * d2);
   const interval_number dr1(blift * cda);
   const interval_number dr12(dr1 * d1);
   const interval_number dr2(alift * bcd);
   const interval_number dr22(dr2 * d2);
   const interval_number dr(dr22 - dr12);
   const interval_number drx1(dr * d3);
   const interval_number drx2(drx1 * d4);
   const interval_number det(dlx2 + drx2);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int inSphere_IIIIE_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, bigfloat pex, bigfloat pey, bigfloat pez)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   p4.getBigfloatLambda(l4x, l4y, l4z, d4);
   const bigfloat pexd(pex * d1);
   const bigfloat peyd(pey * d1);
   const bigfloat pezd(pez * d1);
   const bigfloat aex(l1x - pexd);
   const bigfloat aey(l1y - peyd);
   const bigfloat aez(l1z - pezd);
   const bigfloat pexd2(pex * d2);
   const bigfloat peyd2(pey * d2);
   const bigfloat pezd2(pez * d2);
   const bigfloat bex(l2x - pexd2);
   const bigfloat bey(l2y - peyd2);
   const bigfloat bez(l2z - pezd2);
   const bigfloat pexd3(pex * d3);
   const bigfloat peyd3(pey * d3);
   const bigfloat pezd3(pez * d3);
   const bigfloat cex(l3x - pexd3);
   const bigfloat cey(l3y - peyd3);
   const bigfloat cez(l3z - pezd3);
   const bigfloat pexd4(pex * d4);
   const bigfloat peyd4(pey * d4);
   const bigfloat pezd4(pez * d4);
   const bigfloat dex(l4x - pexd4);
   const bigfloat dey(l4y - peyd4);
   const bigfloat dez(l4z - pezd4);
   const bigfloat aexbey(aex * bey);
   const bigfloat bexaey(bex * aey);
   const bigfloat ab(aexbey - bexaey);
   const bigfloat bexcey(bex * cey);
   const bigfloat cexbey(cex * bey);
   const bigfloat bc(bexcey - cexbey);
   const bigfloat cexdey(cex * dey);
   const bigfloat dexcey(dex * cey);
   const bigfloat cd(cexdey - dexcey);
   const bigfloat dexaey(dex * aey);
   const bigfloat aexdey(aex * dey);
   const bigfloat da(dexaey - aexdey);
   const bigfloat aexcey(aex * cey);
   const bigfloat cexaey(cex * aey);
   const bigfloat ac(aexcey - cexaey);
   const bigfloat bexdey(bex * dey);
   const bigfloat dexbey(dex * bey);
   const bigfloat bd(bexdey - dexbey);
   const bigfloat abc1(aez * bc);
   const bigfloat abc2(bez * ac);
   const bigfloat abc3(cez * ab);
   const bigfloat abc4(abc1 + abc3);
   const bigfloat abc(abc4 - abc2);
   const bigfloat bcd1(bez * cd);
   const bigfloat bcd2(cez * bd);
   const bigfloat bcd3(dez * bc);
   const bigfloat bcd4(bcd1 + bcd3);
   const bigfloat bcd(bcd4 - bcd2);
   const bigfloat cda1(cez * da);
   const bigfloat cda2(dez * ac);
   const bigfloat cda3(aez * cd);
   const bigfloat cda4(cda1 + cda3);
   const bigfloat cda(cda4 + cda2);
   const bigfloat dab1(dez * ab);
   const bigfloat dab2(aez * bd);
   const bigfloat dab3(bez * da);
   const bigfloat dab4(dab1 + dab3);
   const bigfloat dab(dab4 + dab2);
   const bigfloat al1(aex * aex);
   const bigfloat al2(aey * aey);
   const bigfloat al3(aez * aez);
   const bigfloat al4(al1 + al2);
   const bigfloat alift(al4 + al3);
   const bigfloat bl1(bex * bex);
   const bigfloat bl2(bey * bey);
   const bigfloat bl3(bez * bez);
   const bigfloat bl4(bl1 + bl2);
   const bigfloat blift(bl4 + bl3);
   const bigfloat cl1(cex * cex);
   const bigfloat cl2(cey * cey);
   const bigfloat cl3(cez * cez);
   const bigfloat cl4(cl1 + cl2);
   const bigfloat clift(cl4 + cl3);
   const bigfloat dl1(dex * dex);
   const bigfloat dl2(dey * dey);
   const bigfloat dl3(dez * dez);
   const bigfloat dl4(dl1 + dl2);
   const bigfloat dlift(dl4 + dl3);
   const bigfloat ds1(dlift * abc);
   const bigfloat ds12(ds1 * d3);
   const bigfloat ds2(clift * dab);
   const bigfloat ds22(ds2 * d4);
   const bigfloat dl(ds22 - ds12);
   const bigfloat dlx1(dl * d1);
   const bigfloat dlx2(dlx1 * d2);
   const bigfloat dr1(blift * cda);
   const bigfloat dr12(dr1 * d1);
   const bigfloat dr2(alift * bcd);
   const bigfloat dr22(dr2 * d2);
   const bigfloat dr(dr22 - dr12);
   const bigfloat drx1(dr * d3);
   const bigfloat drx2(drx1 * d4);
   const bigfloat det(dlx2 + drx2);
   return sgn(det);
}

inline int inSphere_IIIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, double pex, double pey, double pez)
{
   int ret;
   ret = inSphere_IIIIE_interval(p1, p2, p3, p4, pex, pey, pez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_IIIIE_bigfloat(p1, p2, p3, p4, pex, pey, pez);
}

inline int inSphere_IIIII_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const genericPoint& p5)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4, l5x, l5y, l5z, d5;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   || !p4.getIntervalLambda(l4x, l4y, l4z, d4)
   || !p5.getIntervalLambda(l5x, l5y, l5z, d5)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number pexd(l5x * d1);
   const interval_number peyd(l5y * d1);
   const interval_number pezd(l5z * d1);
   const interval_number ll1x(l1x * d5);
   const interval_number ll1y(l1y * d5);
   const interval_number ll1z(l1z * d5);
   const interval_number aex(ll1x - pexd);
   const interval_number aey(ll1y - peyd);
   const interval_number aez(ll1z - pezd);
   const interval_number pexd2(l5x * d2);
   const interval_number peyd2(l5y * d2);
   const interval_number pezd2(l5z * d2);
   const interval_number ll2x(l2x * d5);
   const interval_number ll2y(l2y * d5);
   const interval_number ll2z(l2z * d5);
   const interval_number bex(ll2x - pexd2);
   const interval_number bey(ll2y - peyd2);
   const interval_number bez(ll2z - pezd2);
   const interval_number pexd3(l5x * d3);
   const interval_number peyd3(l5y * d3);
   const interval_number pezd3(l5z * d3);
   const interval_number ll3x(l3x * d5);
   const interval_number ll3y(l3y * d5);
   const interval_number ll3z(l3z * d5);
   const interval_number cex(ll3x - pexd3);
   const interval_number cey(ll3y - peyd3);
   const interval_number cez(ll3z - pezd3);
   const interval_number pexd4(l5x * d4);
   const interval_number peyd4(l5y * d4);
   const interval_number pezd4(l5z * d4);
   const interval_number ll4x(l4x * d5);
   const interval_number ll4y(l4y * d5);
   const interval_number ll4z(l4z * d5);
   const interval_number dex(ll4x - pexd4);
   const interval_number dey(ll4y - peyd4);
   const interval_number dez(ll4z - pezd4);
   const interval_number aexbey(aex * bey);
   const interval_number bexaey(bex * aey);
   const interval_number ab(aexbey - bexaey);
   const interval_number bexcey(bex * cey);
   const interval_number cexbey(cex * bey);
   const interval_number bc(bexcey - cexbey);
   const interval_number cexdey(cex * dey);
   const interval_number dexcey(dex * cey);
   const interval_number cd(cexdey - dexcey);
   const interval_number dexaey(dex * aey);
   const interval_number aexdey(aex * dey);
   const interval_number da(dexaey - aexdey);
   const interval_number aexcey(aex * cey);
   const interval_number cexaey(cex * aey);
   const interval_number ac(aexcey - cexaey);
   const interval_number bexdey(bex * dey);
   const interval_number dexbey(dex * bey);
   const interval_number bd(bexdey - dexbey);
   const interval_number abc1(aez * bc);
   const interval_number abc2(bez * ac);
   const interval_number abc3(cez * ab);
   const interval_number abc4(abc1 + abc3);
   const interval_number abc(abc4 - abc2);
   const interval_number bcd1(bez * cd);
   const interval_number bcd2(cez * bd);
   const interval_number bcd3(dez * bc);
   const interval_number bcd4(bcd1 + bcd3);
   const interval_number bcd(bcd4 - bcd2);
   const interval_number cda1(cez * da);
   const interval_number cda2(dez * ac);
   const interval_number cda3(aez * cd);
   const interval_number cda4(cda1 + cda3);
   const interval_number cda(cda4 + cda2);
   const interval_number dab1(dez * ab);
   const interval_number dab2(aez * bd);
   const interval_number dab3(bez * da);
   const interval_number dab4(dab1 + dab3);
   const interval_number dab(dab4 + dab2);
   const interval_number al1(aex * aex);
   const interval_number al2(aey * aey);
   const interval_number al3(aez * aez);
   const interval_number al4(al1 + al2);
   const interval_number alift(al4 + al3);
   const interval_number bl1(bex * bex);
   const interval_number bl2(bey * bey);
   const interval_number bl3(bez * bez);
   const interval_number bl4(bl1 + bl2);
   const interval_number blift(bl4 + bl3);
   const interval_number cl1(cex * cex);
   const interval_number cl2(cey * cey);
   const interval_number cl3(cez * cez);
   const interval_number cl4(cl1 + cl2);
   const interval_number clift(cl4 + cl3);
   const interval_number dl1(dex * dex);
   const interval_number dl2(dey * dey);
   const interval_number dl3(dez * dez);
   const interval_number dl4(dl1 + dl2);
   const interval_number dlift(dl4 + dl3);
   const interval_number ds1(dlift * abc);
   const interval_number ds1n(ds1 * d3);
   const interval_number ds2(clift * dab);
   const interval_number ds2n(ds2 * d4);
   const interval_number dl(ds2n - ds1n);
   const interval_number dla(dl * d1);
   const interval_number dlb(dla * d2);
   const interval_number dr1(blift * cda);
   const interval_number dr1n(dr1 * d1);
   const interval_number dr2(alift * bcd);
   const interval_number dr2n(dr2 * d2);
   const interval_number dr(dr2n - dr1n);
   const interval_number dra(dr * d3);
   const interval_number drb(dra * d4);
   const interval_number det(dlb + drb);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int inSphere_IIIII_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const genericPoint& p5)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4, l5x, l5y, l5z, d5;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   p4.getBigfloatLambda(l4x, l4y, l4z, d4);
   p5.getBigfloatLambda(l5x, l5y, l5z, d5);
   const bigfloat pexd(l5x * d1);
   const bigfloat peyd(l5y * d1);
   const bigfloat pezd(l5z * d1);
   const bigfloat ll1x(l1x * d5);
   const bigfloat ll1y(l1y * d5);
   const bigfloat ll1z(l1z * d5);
   const bigfloat aex(ll1x - pexd);
   const bigfloat aey(ll1y - peyd);
   const bigfloat aez(ll1z - pezd);
   const bigfloat pexd2(l5x * d2);
   const bigfloat peyd2(l5y * d2);
   const bigfloat pezd2(l5z * d2);
   const bigfloat ll2x(l2x * d5);
   const bigfloat ll2y(l2y * d5);
   const bigfloat ll2z(l2z * d5);
   const bigfloat bex(ll2x - pexd2);
   const bigfloat bey(ll2y - peyd2);
   const bigfloat bez(ll2z - pezd2);
   const bigfloat pexd3(l5x * d3);
   const bigfloat peyd3(l5y * d3);
   const bigfloat pezd3(l5z * d3);
   const bigfloat ll3x(l3x * d5);
   const bigfloat ll3y(l3y * d5);
   const bigfloat ll3z(l3z * d5);
   const bigfloat cex(ll3x - pexd3);
   const bigfloat cey(ll3y - peyd3);
   const bigfloat cez(ll3z - pezd3);
   const bigfloat pexd4(l5x * d4);
   const bigfloat peyd4(l5y * d4);
   const bigfloat pezd4(l5z * d4);
   const bigfloat ll4x(l4x * d5);
   const bigfloat ll4y(l4y * d5);
   const bigfloat ll4z(l4z * d5);
   const bigfloat dex(ll4x - pexd4);
   const bigfloat dey(ll4y - peyd4);
   const bigfloat dez(ll4z - pezd4);
   const bigfloat aexbey(aex * bey);
   const bigfloat bexaey(bex * aey);
   const bigfloat ab(aexbey - bexaey);
   const bigfloat bexcey(bex * cey);
   const bigfloat cexbey(cex * bey);
   const bigfloat bc(bexcey - cexbey);
   const bigfloat cexdey(cex * dey);
   const bigfloat dexcey(dex * cey);
   const bigfloat cd(cexdey - dexcey);
   const bigfloat dexaey(dex * aey);
   const bigfloat aexdey(aex * dey);
   const bigfloat da(dexaey - aexdey);
   const bigfloat aexcey(aex * cey);
   const bigfloat cexaey(cex * aey);
   const bigfloat ac(aexcey - cexaey);
   const bigfloat bexdey(bex * dey);
   const bigfloat dexbey(dex * bey);
   const bigfloat bd(bexdey - dexbey);
   const bigfloat abc1(aez * bc);
   const bigfloat abc2(bez * ac);
   const bigfloat abc3(cez * ab);
   const bigfloat abc4(abc1 + abc3);
   const bigfloat abc(abc4 - abc2);
   const bigfloat bcd1(bez * cd);
   const bigfloat bcd2(cez * bd);
   const bigfloat bcd3(dez * bc);
   const bigfloat bcd4(bcd1 + bcd3);
   const bigfloat bcd(bcd4 - bcd2);
   const bigfloat cda1(cez * da);
   const bigfloat cda2(dez * ac);
   const bigfloat cda3(aez * cd);
   const bigfloat cda4(cda1 + cda3);
   const bigfloat cda(cda4 + cda2);
   const bigfloat dab1(dez * ab);
   const bigfloat dab2(aez * bd);
   const bigfloat dab3(bez * da);
   const bigfloat dab4(dab1 + dab3);
   const bigfloat dab(dab4 + dab2);
   const bigfloat al1(aex * aex);
   const bigfloat al2(aey * aey);
   const bigfloat al3(aez * aez);
   const bigfloat al4(al1 + al2);
   const bigfloat alift(al4 + al3);
   const bigfloat bl1(bex * bex);
   const bigfloat bl2(bey * bey);
   const bigfloat bl3(bez * bez);
   const bigfloat bl4(bl1 + bl2);
   const bigfloat blift(bl4 + bl3);
   const bigfloat cl1(cex * cex);
   const bigfloat cl2(cey * cey);
   const bigfloat cl3(cez * cez);
   const bigfloat cl4(cl1 + cl2);
   const bigfloat clift(cl4 + cl3);
   const bigfloat dl1(dex * dex);
   const bigfloat dl2(dey * dey);
   const bigfloat dl3(dez * dez);
   const bigfloat dl4(dl1 + dl2);
   const bigfloat dlift(dl4 + dl3);
   const bigfloat ds1(dlift * abc);
   const bigfloat ds1n(ds1 * d3);
   const bigfloat ds2(clift * dab);
   const bigfloat ds2n(ds2 * d4);
   const bigfloat dl(ds2n - ds1n);
   const bigfloat dla(dl * d1);
   const bigfloat dlb(dla * d2);
   const bigfloat dr1(blift * cda);
   const bigfloat dr1n(dr1 * d1);
   const bigfloat dr2(alift * bcd);
   const bigfloat dr2n(dr2 * d2);
   const bigfloat dr(dr2n - dr1n);
   const bigfloat dra(dr * d3);
   const bigfloat drb(dra * d4);
   const bigfloat det(dlb + drb);
   return sgn(det);
}

inline int inSphere_IIIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4, const genericPoint& p5)
{
   int ret;
   ret = inSphere_IIIII_interval(p1, p2, p3, p4, p5);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_IIIII_bigfloat(p1, p2, p3, p4, p5);
}

inline bool lambda2d_SSI_interval(interval_number ea1x, interval_number ea1y, interval_number ea2x, interval_number ea2y, interval_number eb1x, interval_number eb1y, interval_number eb2x, interval_number eb2y, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_det)
{
   setFPUModeToRoundUP();
   const interval_number t1a(ea1x * ea2y);
   const interval_number t1b(ea2x * ea1y);
   const interval_number t1(t1a - t1b);
   const interval_number tx2(eb1x - eb2x);
   const interval_number t3a(eb1x * eb2y);
   const interval_number t3b(eb2x * eb1y);
   const interval_number t3(t3a - t3b);
   const interval_number tx4(ea1x - ea2x);
   const interval_number ty2(eb1y - eb2y);
   const interval_number ty4(ea1y - ea2y);
   const interval_number lxa(t1 * tx2);
   const interval_number lxb(t3 * tx4);
   lambda_x = lxa - lxb;
   const interval_number lya(t1 * ty2);
   const interval_number lyb(t3 * ty4);
   lambda_y = lya - lyb;
   const interval_number deta(tx4 * ty2);
   const interval_number detb(tx2 * ty4);
   lambda_det = deta - detb;
   setFPUModeToRoundNEAR();

   return lambda_det.signIsReliable();
}

inline void lambda2d_SSI_bigfloat(bigfloat ea1x, bigfloat ea1y, bigfloat ea2x, bigfloat ea2y, bigfloat eb1x, bigfloat eb1y, bigfloat eb2x, bigfloat eb2y, bigfloat& lambda_x, bigfloat& lambda_y, bigfloat& lambda_det)
{
   const bigfloat t1a(ea1x * ea2y);
   const bigfloat t1b(ea2x * ea1y);
   const bigfloat t1(t1a - t1b);
   const bigfloat tx2(eb1x - eb2x);
   const bigfloat t3a(eb1x * eb2y);
   const bigfloat t3b(eb2x * eb1y);
   const bigfloat t3(t3a - t3b);
   const bigfloat tx4(ea1x - ea2x);
   const bigfloat ty2(eb1y - eb2y);
   const bigfloat ty4(ea1y - ea2y);
   const bigfloat lxa(t1 * tx2);
   const bigfloat lxb(t3 * tx4);
   lambda_x = lxa - lxb;
   const bigfloat lya(t1 * ty2);
   const bigfloat lyb(t3 * ty4);
   lambda_y = lya - lyb;
   const bigfloat deta(tx4 * ty2);
   const bigfloat detb(tx2 * ty4);
   lambda_det = deta - detb;
}

inline void lambda2d_SSI_exact(double ea1x, double ea1y, double ea2x, double ea2y, double eb1x, double eb1y, double eb2x, double eb2y, double **lambda_x, int& lambda_x_len, double **lambda_y, int& lambda_y_len, double **lambda_det, int& lambda_det_len)
{
   double t1a[2];
   expansionObject::Two_Prod(ea1x, ea2y, t1a);
   double t1b[2];
   expansionObject::Two_Prod(ea2x, ea1y, t1b);
   double t1[4];
   expansionObject::Two_Two_Diff(t1a, t1b, t1);
   double tx2[2];
   expansionObject::two_Diff(eb1x, eb2x, tx2);
   double t3a[2];
   expansionObject::Two_Prod(eb1x, eb2y, t3a);
   double t3b[2];
   expansionObject::Two_Prod(eb2x, eb1y, t3b);
   double t3[4];
   expansionObject::Two_Two_Diff(t3a, t3b, t3);
   double tx4[2];
   expansionObject::two_Diff(ea1x, ea2x, tx4);
   double ty2[2];
   expansionObject::two_Diff(eb1y, eb2y, ty2);
   double ty4[2];
   expansionObject::two_Diff(ea1y, ea2y, ty4);
   double lxa[16];
   int lxa_len = expansionObject::Gen_Product(4, t1, 2, tx2, lxa);
   double lxb[16];
   int lxb_len = expansionObject::Gen_Product(4, t3, 2, tx4, lxb);
   lambda_x_len = expansionObject::Gen_Diff(lxa_len, lxa, lxb_len, lxb, *lambda_x);
   double lya[16];
   int lya_len = expansionObject::Gen_Product(4, t1, 2, ty2, lya);
   double lyb[16];
   int lyb_len = expansionObject::Gen_Product(4, t3, 2, ty4, lyb);
   lambda_y_len = expansionObject::Gen_Diff(lya_len, lya, lyb_len, lyb, *lambda_y);
   double deta[8];
   int deta_len = expansionObject::Gen_Product(2, tx4, 2, ty2, deta);
   double detb[8];
   int detb_len = expansionObject::Gen_Product(2, tx2, 2, ty4, detb);
   lambda_det_len = expansionObject::Gen_Diff(deta_len, deta, detb_len, detb, *lambda_det);

}

inline bool lambda3d_BPT_interval(interval_number px, interval_number py, interval_number pz, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number u, interval_number v, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z, interval_number& lambda_d)
{
   setFPUModeToRoundUP();
   const interval_number dvx(px - rx);
   const interval_number dvy(py - ry);
   const interval_number dvz(pz - rz);
   const interval_number dux(qx - rx);
   const interval_number duy(qy - ry);
   const interval_number duz(qz - rz);
   const interval_number evx(dvx * v);
   const interval_number evy(dvy * v);
   const interval_number evz(dvz * v);
   const interval_number eux(dux * u);
   const interval_number euy(duy * u);
   const interval_number euz(duz * u);
   const interval_number rux(eux + rx);
   const interval_number ruy(euy + ry);
   const interval_number ruz(euz + rz);
   lambda_x = evx + rux;
   lambda_y = evy + ruy;
   lambda_z = evz + ruz;
   lambda_d = 1;
   setFPUModeToRoundNEAR();

   return true;
}

inline void lambda3d_BPT_bigfloat(bigfloat px, bigfloat py, bigfloat pz, bigfloat qx, bigfloat qy, bigfloat qz, bigfloat rx, bigfloat ry, bigfloat rz, bigfloat u, bigfloat v, bigfloat& lambda_x, bigfloat& lambda_y, bigfloat& lambda_z, bigfloat& lambda_d)
{
   const bigfloat dvx(px - rx);
   const bigfloat dvy(py - ry);
   const bigfloat dvz(pz - rz);
   const bigfloat dux(qx - rx);
   const bigfloat duy(qy - ry);
   const bigfloat duz(qz - rz);
   const bigfloat evx(dvx * v);
   const bigfloat evy(dvy * v);
   const bigfloat evz(dvz * v);
   const bigfloat eux(dux * u);
   const bigfloat euy(duy * u);
   const bigfloat euz(duz * u);
   const bigfloat rux(eux + rx);
   const bigfloat ruy(euy + ry);
   const bigfloat ruz(euz + rz);
   lambda_x = evx + rux;
   lambda_y = evy + ruy;
   lambda_z = evz + ruz;
   lambda_d = 1;
}

inline void lambda3d_BPT_exact(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double u, double v, double **lambda_x, int& lambda_x_len, double **lambda_y, int& lambda_y_len, double **lambda_z, int& lambda_z_len, double **lambda_d, int& lambda_d_len)
{
   double dvx[2];
   expansionObject::two_Diff(px, rx, dvx);
   double dvy[2];
   expansionObject::two_Diff(py, ry, dvy);
   double dvz[2];
   expansionObject::two_Diff(pz, rz, dvz);
   double dux[2];
   expansionObject::two_Diff(qx, rx, dux);
   double duy[2];
   expansionObject::two_Diff(qy, ry, duy);
   double duz[2];
   expansionObject::two_Diff(qz, rz, duz);
   double evx[4];
   expansionObject::Two_One_Prod(dvx, v, evx);
   double evy[4];
   expansionObject::Two_One_Prod(dvy, v, evy);
   double evz[4];
   expansionObject::Two_One_Prod(dvz, v, evz);
   double eux[4];
   expansionObject::Two_One_Prod(dux, u, eux);
   double euy[4];
   expansionObject::Two_One_Prod(duy, u, euy);
   double euz[4];
   expansionObject::Two_One_Prod(duz, u, euz);
   double rux[5];
   int rux_len = expansionObject::Gen_Sum(4, eux, 1, &rx, rux);
   double ruy[5];
   int ruy_len = expansionObject::Gen_Sum(4, euy, 1, &ry, ruy);
   double ruz[5];
   int ruz_len = expansionObject::Gen_Sum(4, euz, 1, &rz, ruz);
   lambda_x_len = expansionObject::Gen_Sum(4, evx, rux_len, rux, *lambda_x);
   lambda_y_len = expansionObject::Gen_Sum(4, evy, ruy_len, ruy, *lambda_y);
   lambda_z_len = expansionObject::Gen_Sum(4, evz, ruz_len, ruz, *lambda_z);
   (*lambda_d)[0] = 1;
   lambda_d_len = 1;

}

inline bool lambda3d_LNC_interval(interval_number px, interval_number py, interval_number pz, interval_number qx, interval_number qy, interval_number qz, interval_number t, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z, interval_number& lambda_d)
{
   setFPUModeToRoundUP();
   const interval_number vx(px - qx);
   const interval_number vy(py - qy);
   const interval_number vz(pz - qz);
   const interval_number vxt(vx * t);
   const interval_number vyt(vy * t);
   const interval_number vzt(vz * t);
   lambda_x = px - vxt;
   lambda_y = py - vyt;
   lambda_z = pz - vzt;
   lambda_d = 1;
   setFPUModeToRoundNEAR();

   return true;
}

inline void lambda3d_LNC_bigfloat(bigfloat px, bigfloat py, bigfloat pz, bigfloat qx, bigfloat qy, bigfloat qz, bigfloat t, bigfloat& lambda_x, bigfloat& lambda_y, bigfloat& lambda_z, bigfloat& lambda_d)
{
   const bigfloat vx(px - qx);
   const bigfloat vy(py - qy);
   const bigfloat vz(pz - qz);
   const bigfloat vxt(vx * t);
   const bigfloat vyt(vy * t);
   const bigfloat vzt(vz * t);
   lambda_x = px - vxt;
   lambda_y = py - vyt;
   lambda_z = pz - vzt;
   lambda_d = 1;
}

inline void lambda3d_LNC_exact(double px, double py, double pz, double qx, double qy, double qz, double t, double **lambda_x, int& lambda_x_len, double **lambda_y, int& lambda_y_len, double **lambda_z, int& lambda_z_len, double **lambda_d, int& lambda_d_len)
{
   double vx[2];
   expansionObject::two_Diff(px, qx, vx);
   double vy[2];
   expansionObject::two_Diff(py, qy, vy);
   double vz[2];
   expansionObject::two_Diff(pz, qz, vz);
   double vxt[4];
   expansionObject::Two_One_Prod(vx, t, vxt);
   double vyt[4];
   expansionObject::Two_One_Prod(vy, t, vyt);
   double vzt[4];
   expansionObject::Two_One_Prod(vz, t, vzt);
   lambda_x_len = expansionObject::Gen_Diff(1, &px, 4, vxt, *lambda_x);
   lambda_y_len = expansionObject::Gen_Diff(1, &py, 4, vyt, *lambda_y);
   lambda_z_len = expansionObject::Gen_Diff(1, &pz, 4, vzt, *lambda_z);
   (*lambda_d)[0] = 1;
   lambda_d_len = 1;

}

inline bool lambda3d_LPI_interval(interval_number px, interval_number py, interval_number pz, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number sx, interval_number sy, interval_number sz, interval_number tx, interval_number ty, interval_number tz, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z, interval_number& lambda_d)
{
   setFPUModeToRoundUP();
   const interval_number a11(px - qx);
   const interval_number a12(py - qy);
   const interval_number a13(pz - qz);
   const interval_number a21(sx - rx);
   const interval_number a22(sy - ry);
   const interval_number a23(sz - rz);
   const interval_number a31(tx - rx);
   const interval_number a32(ty - ry);
   const interval_number a33(tz - rz);
   const interval_number tv1(a22 * a33);
   const interval_number tv2(a23 * a32);
   const interval_number a2233(tv1 - tv2);
   const interval_number tv3(a21 * a33);
   const interval_number tv4(a23 * a31);
   const interval_number a2133(tv3 - tv4);
   const interval_number tv5(a21 * a32);
   const interval_number tv6(a22 * a31);
   const interval_number a2132(tv5 - tv6);
   const interval_number tv7(a11 * a2233);
   const interval_number tv8(a12 * a2133);
   const interval_number tv9(a13 * a2132);
   const interval_number tt1(tv7 - tv8);
   const interval_number ld(tt1 + tv9);
   const interval_number px_rx(px - rx);
   const interval_number py_ry(py - ry);
   const interval_number pz_rz(pz - rz);
   const interval_number tt2(py_ry * a2133);
   const interval_number tt3(px_rx * a2233);
   const interval_number tt4(pz_rz * a2132);
   const interval_number tt5(tt3 + tt4);
   const interval_number n(tt5 - tt2);
   const interval_number ax(a11 * n);
   const interval_number ay(a12 * n);
   const interval_number az(a13 * n);
   const interval_number dpx(ld * px);
   const interval_number dpy(ld * py);
   const interval_number dpz(ld * pz);
   lambda_x = dpx - ax;
   lambda_y = dpy - ay;
   lambda_z = dpz - az;
   lambda_d = tt1 + tv9;
   setFPUModeToRoundNEAR();

   return lambda_d.signIsReliable();
}

inline void lambda3d_LPI_bigfloat(bigfloat px, bigfloat py, bigfloat pz, bigfloat qx, bigfloat qy, bigfloat qz, bigfloat rx, bigfloat ry, bigfloat rz, bigfloat sx, bigfloat sy, bigfloat sz, bigfloat tx, bigfloat ty, bigfloat tz, bigfloat& lambda_x, bigfloat& lambda_y, bigfloat& lambda_z, bigfloat& lambda_d)
{
   const bigfloat a11(px - qx);
   const bigfloat a12(py - qy);
   const bigfloat a13(pz - qz);
   const bigfloat a21(sx - rx);
   const bigfloat a22(sy - ry);
   const bigfloat a23(sz - rz);
   const bigfloat a31(tx - rx);
   const bigfloat a32(ty - ry);
   const bigfloat a33(tz - rz);
   const bigfloat tv1(a22 * a33);
   const bigfloat tv2(a23 * a32);
   const bigfloat a2233(tv1 - tv2);
   const bigfloat tv3(a21 * a33);
   const bigfloat tv4(a23 * a31);
   const bigfloat a2133(tv3 - tv4);
   const bigfloat tv5(a21 * a32);
   const bigfloat tv6(a22 * a31);
   const bigfloat a2132(tv5 - tv6);
   const bigfloat tv7(a11 * a2233);
   const bigfloat tv8(a12 * a2133);
   const bigfloat tv9(a13 * a2132);
   const bigfloat tt1(tv7 - tv8);
   const bigfloat ld(tt1 + tv9);
   const bigfloat px_rx(px - rx);
   const bigfloat py_ry(py - ry);
   const bigfloat pz_rz(pz - rz);
   const bigfloat tt2(py_ry * a2133);
   const bigfloat tt3(px_rx * a2233);
   const bigfloat tt4(pz_rz * a2132);
   const bigfloat tt5(tt3 + tt4);
   const bigfloat n(tt5 - tt2);
   const bigfloat ax(a11 * n);
   const bigfloat ay(a12 * n);
   const bigfloat az(a13 * n);
   const bigfloat dpx(ld * px);
   const bigfloat dpy(ld * py);
   const bigfloat dpz(ld * pz);
   lambda_x = dpx - ax;
   lambda_y = dpy - ay;
   lambda_z = dpz - az;
   lambda_d = tt1 + tv9;
}

inline void lambda3d_LPI_exact(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz, double **lambda_x, int& lambda_x_len, double **lambda_y, int& lambda_y_len, double **lambda_z, int& lambda_z_len, double **lambda_d, int& lambda_d_len)
{
   double a11[2];
   expansionObject::two_Diff(px, qx, a11);
   double a12[2];
   expansionObject::two_Diff(py, qy, a12);
   double a13[2];
   expansionObject::two_Diff(pz, qz, a13);
   double a21[2];
   expansionObject::two_Diff(sx, rx, a21);
   double a22[2];
   expansionObject::two_Diff(sy, ry, a22);
   double a23[2];
   expansionObject::two_Diff(sz, rz, a23);
   double a31[2];
   expansionObject::two_Diff(tx, rx, a31);
   double a32[2];
   expansionObject::two_Diff(ty, ry, a32);
   double a33[2];
   expansionObject::two_Diff(tz, rz, a33);
   double tv1[8];
   int tv1_len = expansionObject::Gen_Product(2, a22, 2, a33, tv1);
   double tv2[8];
   int tv2_len = expansionObject::Gen_Product(2, a23, 2, a32, tv2);
   double a2233[16];
   int a2233_len = expansionObject::Gen_Diff(tv1_len, tv1, tv2_len, tv2, a2233);
   double tv3[8];
   int tv3_len = expansionObject::Gen_Product(2, a21, 2, a33, tv3);
   double tv4[8];
   int tv4_len = expansionObject::Gen_Product(2, a23, 2, a31, tv4);
   double a2133[16];
   int a2133_len = expansionObject::Gen_Diff(tv3_len, tv3, tv4_len, tv4, a2133);
   double tv5[8];
   int tv5_len = expansionObject::Gen_Product(2, a21, 2, a32, tv5);
   double tv6[8];
   int tv6_len = expansionObject::Gen_Product(2, a22, 2, a31, tv6);
   double a2132[16];
   int a2132_len = expansionObject::Gen_Diff(tv5_len, tv5, tv6_len, tv6, a2132);
   double tv7[64];
   int tv7_len = expansionObject::Gen_Product(2, a11, a2233_len, a2233, tv7);
   double tv8[64];
   int tv8_len = expansionObject::Gen_Product(2, a12, a2133_len, a2133, tv8);
   double tv9[64];
   int tv9_len = expansionObject::Gen_Product(2, a13, a2132_len, a2132, tv9);
   double tt1[128];
   int tt1_len = expansionObject::Gen_Diff(tv7_len, tv7, tv8_len, tv8, tt1);
   double ld_p[128], *ld = ld_p;
   int ld_len = expansionObject::Gen_Sum_With_PreAlloc(tt1_len, tt1, tv9_len, tv9, &ld, 128);
   double px_rx[2];
   expansionObject::two_Diff(px, rx, px_rx);
   double py_ry[2];
   expansionObject::two_Diff(py, ry, py_ry);
   double pz_rz[2];
   expansionObject::two_Diff(pz, rz, pz_rz);
   double tt2[64];
   int tt2_len = expansionObject::Gen_Product(2, py_ry, a2133_len, a2133, tt2);
   double tt3[64];
   int tt3_len = expansionObject::Gen_Product(2, px_rx, a2233_len, a2233, tt3);
   double tt4[64];
   int tt4_len = expansionObject::Gen_Product(2, pz_rz, a2132_len, a2132, tt4);
   double tt5[128];
   int tt5_len = expansionObject::Gen_Sum(tt3_len, tt3, tt4_len, tt4, tt5);
   double n_p[128], *n = n_p;
   int n_len = expansionObject::Gen_Diff_With_PreAlloc(tt5_len, tt5, tt2_len, tt2, &n, 128);
   double ax_p[128], *ax = ax_p;
   int ax_len = expansionObject::Gen_Product_With_PreAlloc(2, a11, n_len, n, &ax, 128);
   double ay_p[128], *ay = ay_p;
   int ay_len = expansionObject::Gen_Product_With_PreAlloc(2, a12, n_len, n, &ay, 128);
   double az_p[128], *az = az_p;
   int az_len = expansionObject::Gen_Product_With_PreAlloc(2, a13, n_len, n, &az, 128);
   double dpx_p[128], *dpx = dpx_p;
   int dpx_len = expansionObject::Gen_Scale_With_PreAlloc(ld_len, ld, px, &dpx, 128);
   double dpy_p[128], *dpy = dpy_p;
   int dpy_len = expansionObject::Gen_Scale_With_PreAlloc(ld_len, ld, py, &dpy, 128);
   double dpz_p[128], *dpz = dpz_p;
   int dpz_len = expansionObject::Gen_Scale_With_PreAlloc(ld_len, ld, pz, &dpz, 128);
   lambda_x_len = expansionObject::Gen_Diff_With_PreAlloc(dpx_len, dpx, ax_len, ax, lambda_x, lambda_x_len);
   lambda_y_len = expansionObject::Gen_Diff_With_PreAlloc(dpy_len, dpy, ay_len, ay, lambda_y, lambda_y_len);
   lambda_z_len = expansionObject::Gen_Diff_With_PreAlloc(dpz_len, dpz, az_len, az, lambda_z, lambda_z_len);
   lambda_d_len = expansionObject::Gen_Sum_With_PreAlloc(tt1_len, tt1, tv9_len, tv9, lambda_d, lambda_d_len);

   if (dpz_p != dpz) FreeDoubles(dpz);
   if (dpy_p != dpy) FreeDoubles(dpy);
   if (dpx_p != dpx) FreeDoubles(dpx);
   if (az_p != az) FreeDoubles(az);
   if (ay_p != ay) FreeDoubles(ay);
   if (ax_p != ax) FreeDoubles(ax);
   if (n_p != n) FreeDoubles(n);
   if (ld_p != ld) FreeDoubles(ld);
}

inline bool lambda3d_TPI_interval(interval_number ov1x, interval_number ov1y, interval_number ov1z, interval_number ov2x, interval_number ov2y, interval_number ov2z, interval_number ov3x, interval_number ov3y, interval_number ov3z, interval_number ow1x, interval_number ow1y, interval_number ow1z, interval_number ow2x, interval_number ow2y, interval_number ow2z, interval_number ow3x, interval_number ow3y, interval_number ow3z, interval_number ou1x, interval_number ou1y, interval_number ou1z, interval_number ou2x, interval_number ou2y, interval_number ou2z, interval_number ou3x, interval_number ou3y, interval_number ou3z, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z, interval_number& lambda_d)
{
   setFPUModeToRoundUP();
   const interval_number v3x(ov3x - ov2x);
   const interval_number v3y(ov3y - ov2y);
   const interval_number v3z(ov3z - ov2z);
   const interval_number v2x(ov2x - ov1x);
   const interval_number v2y(ov2y - ov1y);
   const interval_number v2z(ov2z - ov1z);
   const interval_number w3x(ow3x - ow2x);
   const interval_number w3y(ow3y - ow2y);
   const interval_number w3z(ow3z - ow2z);
   const interval_number w2x(ow2x - ow1x);
   const interval_number w2y(ow2y - ow1y);
   const interval_number w2z(ow2z - ow1z);
   const interval_number u3x(ou3x - ou2x);
   const interval_number u3y(ou3y - ou2y);
   const interval_number u3z(ou3z - ou2z);
   const interval_number u2x(ou2x - ou1x);
   const interval_number u2y(ou2y - ou1y);
   const interval_number u2z(ou2z - ou1z);
   const interval_number nvx1(v2y * v3z);
   const interval_number nvx2(v2z * v3y);
   const interval_number nvx(nvx1 - nvx2);
   const interval_number nvy1(v3x * v2z);
   const interval_number nvy2(v3z * v2x);
   const interval_number nvy(nvy1 - nvy2);
   const interval_number nvz1(v2x * v3y);
   const interval_number nvz2(v2y * v3x);
   const interval_number nvz(nvz1 - nvz2);
   const interval_number nwx1(w2y * w3z);
   const interval_number nwx2(w2z * w3y);
   const interval_number nwx(nwx1 - nwx2);
   const interval_number nwy1(w3x * w2z);
   const interval_number nwy2(w3z * w2x);
   const interval_number nwy(nwy1 - nwy2);
   const interval_number nwz1(w2x * w3y);
   const interval_number nwz2(w2y * w3x);
   const interval_number nwz(nwz1 - nwz2);
   const interval_number nux1(u2y * u3z);
   const interval_number nux2(u2z * u3y);
   const interval_number nux(nux1 - nux2);
   const interval_number nuy1(u3x * u2z);
   const interval_number nuy2(u3z * u2x);
   const interval_number nuy(nuy1 - nuy2);
   const interval_number nuz1(u2x * u3y);
   const interval_number nuz2(u2y * u3x);
   const interval_number nuz(nuz1 - nuz2);
   const interval_number nwyuz1(nwy * nuz);
   const interval_number nwyuz2(nwz * nuy);
   const interval_number nwyuz(nwyuz1 - nwyuz2);
   const interval_number nwxuz1(nwx * nuz);
   const interval_number nwxuz2(nwz * nux);
   const interval_number nwxuz(nwxuz1 - nwxuz2);
   const interval_number nwxuy1(nwx * nuy);
   const interval_number nwxuy2(nwy * nux);
   const interval_number nwxuy(nwxuy1 - nwxuy2);
   const interval_number nvyuz1(nvy * nuz);
   const interval_number nvyuz2(nvz * nuy);
   const interval_number nvyuz(nvyuz1 - nvyuz2);
   const interval_number nvxuz1(nvx * nuz);
   const interval_number nvxuz2(nvz * nux);
   const interval_number nvxuz(nvxuz1 - nvxuz2);
   const interval_number nvxuy1(nvx * nuy);
   const interval_number nvxuy2(nvy * nux);
   const interval_number nvxuy(nvxuy1 - nvxuy2);
   const interval_number nvywz1(nvy * nwz);
   const interval_number nvywz2(nvz * nwy);
   const interval_number nvywz(nvywz1 - nvywz2);
   const interval_number nvxwz1(nvx * nwz);
   const interval_number nvxwz2(nvz * nwx);
   const interval_number nvxwz(nvxwz1 - nvxwz2);
   const interval_number nvxwy1(nvx * nwy);
   const interval_number nvxwy2(nvy * nwx);
   const interval_number nvxwy(nvxwy1 - nvxwy2);
   const interval_number p1a(nvx * ov1x);
   const interval_number p1b(nvy * ov1y);
   const interval_number p1c(nvz * ov1z);
   const interval_number p1ab(p1a + p1b);
   const interval_number p1(p1ab + p1c);
   const interval_number p2a(nwx * ow1x);
   const interval_number p2b(nwy * ow1y);
   const interval_number p2c(nwz * ow1z);
   const interval_number p2ab(p2a + p2b);
   const interval_number p2(p2ab + p2c);
   const interval_number p3a(nux * ou1x);
   const interval_number p3b(nuy * ou1y);
   const interval_number p3c(nuz * ou1z);
   const interval_number p3ab(p3a + p3b);
   const interval_number p3(p3ab + p3c);
   const interval_number lxa(p1 * nwyuz);
   const interval_number lxb(p3 * nvywz);
   const interval_number lxc(p2 * nvyuz);
   const interval_number lxab(lxa + lxb);
   lambda_x = lxab - lxc;
   const interval_number lya(p2 * nvxuz);
   const interval_number lyb(p3 * nvxwz);
   const interval_number lyc(p1 * nwxuz);
   const interval_number lybc(lyc + lyb);
   lambda_y = lya - lybc;
   const interval_number lza(p3 * nvxwy);
   const interval_number lzb(p1 * nwxuy);
   const interval_number lzc(p2 * nvxuy);
   const interval_number lzab(lza + lzb);
   lambda_z = lzab - lzc;
   const interval_number da(nvx * nwyuz);
   const interval_number db(nvz * nwxuy);
   const interval_number dc(nvy * nwxuz);
   const interval_number dab(da + db);
   lambda_d = dab - dc;
   setFPUModeToRoundNEAR();

   return lambda_d.signIsReliable();
}

inline void lambda3d_TPI_bigfloat(bigfloat ov1x, bigfloat ov1y, bigfloat ov1z, bigfloat ov2x, bigfloat ov2y, bigfloat ov2z, bigfloat ov3x, bigfloat ov3y, bigfloat ov3z, bigfloat ow1x, bigfloat ow1y, bigfloat ow1z, bigfloat ow2x, bigfloat ow2y, bigfloat ow2z, bigfloat ow3x, bigfloat ow3y, bigfloat ow3z, bigfloat ou1x, bigfloat ou1y, bigfloat ou1z, bigfloat ou2x, bigfloat ou2y, bigfloat ou2z, bigfloat ou3x, bigfloat ou3y, bigfloat ou3z, bigfloat& lambda_x, bigfloat& lambda_y, bigfloat& lambda_z, bigfloat& lambda_d)
{
   const bigfloat v3x(ov3x - ov2x);
   const bigfloat v3y(ov3y - ov2y);
   const bigfloat v3z(ov3z - ov2z);
   const bigfloat v2x(ov2x - ov1x);
   const bigfloat v2y(ov2y - ov1y);
   const bigfloat v2z(ov2z - ov1z);
   const bigfloat w3x(ow3x - ow2x);
   const bigfloat w3y(ow3y - ow2y);
   const bigfloat w3z(ow3z - ow2z);
   const bigfloat w2x(ow2x - ow1x);
   const bigfloat w2y(ow2y - ow1y);
   const bigfloat w2z(ow2z - ow1z);
   const bigfloat u3x(ou3x - ou2x);
   const bigfloat u3y(ou3y - ou2y);
   const bigfloat u3z(ou3z - ou2z);
   const bigfloat u2x(ou2x - ou1x);
   const bigfloat u2y(ou2y - ou1y);
   const bigfloat u2z(ou2z - ou1z);
   const bigfloat nvx1(v2y * v3z);
   const bigfloat nvx2(v2z * v3y);
   const bigfloat nvx(nvx1 - nvx2);
   const bigfloat nvy1(v3x * v2z);
   const bigfloat nvy2(v3z * v2x);
   const bigfloat nvy(nvy1 - nvy2);
   const bigfloat nvz1(v2x * v3y);
   const bigfloat nvz2(v2y * v3x);
   const bigfloat nvz(nvz1 - nvz2);
   const bigfloat nwx1(w2y * w3z);
   const bigfloat nwx2(w2z * w3y);
   const bigfloat nwx(nwx1 - nwx2);
   const bigfloat nwy1(w3x * w2z);
   const bigfloat nwy2(w3z * w2x);
   const bigfloat nwy(nwy1 - nwy2);
   const bigfloat nwz1(w2x * w3y);
   const bigfloat nwz2(w2y * w3x);
   const bigfloat nwz(nwz1 - nwz2);
   const bigfloat nux1(u2y * u3z);
   const bigfloat nux2(u2z * u3y);
   const bigfloat nux(nux1 - nux2);
   const bigfloat nuy1(u3x * u2z);
   const bigfloat nuy2(u3z * u2x);
   const bigfloat nuy(nuy1 - nuy2);
   const bigfloat nuz1(u2x * u3y);
   const bigfloat nuz2(u2y * u3x);
   const bigfloat nuz(nuz1 - nuz2);
   const bigfloat nwyuz1(nwy * nuz);
   const bigfloat nwyuz2(nwz * nuy);
   const bigfloat nwyuz(nwyuz1 - nwyuz2);
   const bigfloat nwxuz1(nwx * nuz);
   const bigfloat nwxuz2(nwz * nux);
   const bigfloat nwxuz(nwxuz1 - nwxuz2);
   const bigfloat nwxuy1(nwx * nuy);
   const bigfloat nwxuy2(nwy * nux);
   const bigfloat nwxuy(nwxuy1 - nwxuy2);
   const bigfloat nvyuz1(nvy * nuz);
   const bigfloat nvyuz2(nvz * nuy);
   const bigfloat nvyuz(nvyuz1 - nvyuz2);
   const bigfloat nvxuz1(nvx * nuz);
   const bigfloat nvxuz2(nvz * nux);
   const bigfloat nvxuz(nvxuz1 - nvxuz2);
   const bigfloat nvxuy1(nvx * nuy);
   const bigfloat nvxuy2(nvy * nux);
   const bigfloat nvxuy(nvxuy1 - nvxuy2);
   const bigfloat nvywz1(nvy * nwz);
   const bigfloat nvywz2(nvz * nwy);
   const bigfloat nvywz(nvywz1 - nvywz2);
   const bigfloat nvxwz1(nvx * nwz);
   const bigfloat nvxwz2(nvz * nwx);
   const bigfloat nvxwz(nvxwz1 - nvxwz2);
   const bigfloat nvxwy1(nvx * nwy);
   const bigfloat nvxwy2(nvy * nwx);
   const bigfloat nvxwy(nvxwy1 - nvxwy2);
   const bigfloat p1a(nvx * ov1x);
   const bigfloat p1b(nvy * ov1y);
   const bigfloat p1c(nvz * ov1z);
   const bigfloat p1ab(p1a + p1b);
   const bigfloat p1(p1ab + p1c);
   const bigfloat p2a(nwx * ow1x);
   const bigfloat p2b(nwy * ow1y);
   const bigfloat p2c(nwz * ow1z);
   const bigfloat p2ab(p2a + p2b);
   const bigfloat p2(p2ab + p2c);
   const bigfloat p3a(nux * ou1x);
   const bigfloat p3b(nuy * ou1y);
   const bigfloat p3c(nuz * ou1z);
   const bigfloat p3ab(p3a + p3b);
   const bigfloat p3(p3ab + p3c);
   const bigfloat lxa(p1 * nwyuz);
   const bigfloat lxb(p3 * nvywz);
   const bigfloat lxc(p2 * nvyuz);
   const bigfloat lxab(lxa + lxb);
   lambda_x = lxab - lxc;
   const bigfloat lya(p2 * nvxuz);
   const bigfloat lyb(p3 * nvxwz);
   const bigfloat lyc(p1 * nwxuz);
   const bigfloat lybc(lyc + lyb);
   lambda_y = lya - lybc;
   const bigfloat lza(p3 * nvxwy);
   const bigfloat lzb(p1 * nwxuy);
   const bigfloat lzc(p2 * nvxuy);
   const bigfloat lzab(lza + lzb);
   lambda_z = lzab - lzc;
   const bigfloat da(nvx * nwyuz);
   const bigfloat db(nvz * nwxuy);
   const bigfloat dc(nvy * nwxuz);
   const bigfloat dab(da + db);
   lambda_d = dab - dc;
}

inline void lambda3d_TPI_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z, double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z, double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z, double ou3x, double ou3y, double ou3z, double **lambda_x, int& lambda_x_len, double **lambda_y, int& lambda_y_len, double **lambda_z, int& lambda_z_len, double **lambda_d, int& lambda_d_len)
{
   double v3x[2];
   expansionObject::two_Diff(ov3x, ov2x, v3x);
   double v3y[2];
   expansionObject::two_Diff(ov3y, ov2y, v3y);
   double v3z[2];
   expansionObject::two_Diff(ov3z, ov2z, v3z);
   double v2x[2];
   expansionObject::two_Diff(ov2x, ov1x, v2x);
   double v2y[2];
   expansionObject::two_Diff(ov2y, ov1y, v2y);
   double v2z[2];
   expansionObject::two_Diff(ov2z, ov1z, v2z);
   double w3x[2];
   expansionObject::two_Diff(ow3x, ow2x, w3x);
   double w3y[2];
   expansionObject::two_Diff(ow3y, ow2y, w3y);
   double w3z[2];
   expansionObject::two_Diff(ow3z, ow2z, w3z);
   double w2x[2];
   expansionObject::two_Diff(ow2x, ow1x, w2x);
   double w2y[2];
   expansionObject::two_Diff(ow2y, ow1y, w2y);
   double w2z[2];
   expansionObject::two_Diff(ow2z, ow1z, w2z);
   double u3x[2];
   expansionObject::two_Diff(ou3x, ou2x, u3x);
   double u3y[2];
   expansionObject::two_Diff(ou3y, ou2y, u3y);
   double u3z[2];
   expansionObject::two_Diff(ou3z, ou2z, u3z);
   double u2x[2];
   expansionObject::two_Diff(ou2x, ou1x, u2x);
   double u2y[2];
   expansionObject::two_Diff(ou2y, ou1y, u2y);
   double u2z[2];
   expansionObject::two_Diff(ou2z, ou1z, u2z);
   double nvx1[8];
   int nvx1_len = expansionObject::Gen_Product(2, v2y, 2, v3z, nvx1);
   double nvx2[8];
   int nvx2_len = expansionObject::Gen_Product(2, v2z, 2, v3y, nvx2);
   double nvx[16];
   int nvx_len = expansionObject::Gen_Diff(nvx1_len, nvx1, nvx2_len, nvx2, nvx);
   double nvy1[8];
   int nvy1_len = expansionObject::Gen_Product(2, v3x, 2, v2z, nvy1);
   double nvy2[8];
   int nvy2_len = expansionObject::Gen_Product(2, v3z, 2, v2x, nvy2);
   double nvy[16];
   int nvy_len = expansionObject::Gen_Diff(nvy1_len, nvy1, nvy2_len, nvy2, nvy);
   double nvz1[8];
   int nvz1_len = expansionObject::Gen_Product(2, v2x, 2, v3y, nvz1);
   double nvz2[8];
   int nvz2_len = expansionObject::Gen_Product(2, v2y, 2, v3x, nvz2);
   double nvz[16];
   int nvz_len = expansionObject::Gen_Diff(nvz1_len, nvz1, nvz2_len, nvz2, nvz);
   double nwx1[8];
   int nwx1_len = expansionObject::Gen_Product(2, w2y, 2, w3z, nwx1);
   double nwx2[8];
   int nwx2_len = expansionObject::Gen_Product(2, w2z, 2, w3y, nwx2);
   double nwx[16];
   int nwx_len = expansionObject::Gen_Diff(nwx1_len, nwx1, nwx2_len, nwx2, nwx);
   double nwy1[8];
   int nwy1_len = expansionObject::Gen_Product(2, w3x, 2, w2z, nwy1);
   double nwy2[8];
   int nwy2_len = expansionObject::Gen_Product(2, w3z, 2, w2x, nwy2);
   double nwy[16];
   int nwy_len = expansionObject::Gen_Diff(nwy1_len, nwy1, nwy2_len, nwy2, nwy);
   double nwz1[8];
   int nwz1_len = expansionObject::Gen_Product(2, w2x, 2, w3y, nwz1);
   double nwz2[8];
   int nwz2_len = expansionObject::Gen_Product(2, w2y, 2, w3x, nwz2);
   double nwz[16];
   int nwz_len = expansionObject::Gen_Diff(nwz1_len, nwz1, nwz2_len, nwz2, nwz);
   double nux1[8];
   int nux1_len = expansionObject::Gen_Product(2, u2y, 2, u3z, nux1);
   double nux2[8];
   int nux2_len = expansionObject::Gen_Product(2, u2z, 2, u3y, nux2);
   double nux[16];
   int nux_len = expansionObject::Gen_Diff(nux1_len, nux1, nux2_len, nux2, nux);
   double nuy1[8];
   int nuy1_len = expansionObject::Gen_Product(2, u3x, 2, u2z, nuy1);
   double nuy2[8];
   int nuy2_len = expansionObject::Gen_Product(2, u3z, 2, u2x, nuy2);
   double nuy[16];
   int nuy_len = expansionObject::Gen_Diff(nuy1_len, nuy1, nuy2_len, nuy2, nuy);
   double nuz1[8];
   int nuz1_len = expansionObject::Gen_Product(2, u2x, 2, u3y, nuz1);
   double nuz2[8];
   int nuz2_len = expansionObject::Gen_Product(2, u2y, 2, u3x, nuz2);
   double nuz[16];
   int nuz_len = expansionObject::Gen_Diff(nuz1_len, nuz1, nuz2_len, nuz2, nuz);
   double nwyuz1_p[16], *nwyuz1 = nwyuz1_p;
   int nwyuz1_len = expansionObject::Gen_Product_With_PreAlloc(nwy_len, nwy, nuz_len, nuz, &nwyuz1, 16);
   double nwyuz2_p[16], *nwyuz2 = nwyuz2_p;
   int nwyuz2_len = expansionObject::Gen_Product_With_PreAlloc(nwz_len, nwz, nuy_len, nuy, &nwyuz2, 16);
   double nwyuz_p[16], *nwyuz = nwyuz_p;
   int nwyuz_len = expansionObject::Gen_Diff_With_PreAlloc(nwyuz1_len, nwyuz1, nwyuz2_len, nwyuz2, &nwyuz, 16);
   double nwxuz1_p[16], *nwxuz1 = nwxuz1_p;
   int nwxuz1_len = expansionObject::Gen_Product_With_PreAlloc(nwx_len, nwx, nuz_len, nuz, &nwxuz1, 16);
   double nwxuz2_p[16], *nwxuz2 = nwxuz2_p;
   int nwxuz2_len = expansionObject::Gen_Product_With_PreAlloc(nwz_len, nwz, nux_len, nux, &nwxuz2, 16);
   double nwxuz_p[16], *nwxuz = nwxuz_p;
   int nwxuz_len = expansionObject::Gen_Diff_With_PreAlloc(nwxuz1_len, nwxuz1, nwxuz2_len, nwxuz2, &nwxuz, 16);
   double nwxuy1_p[16], *nwxuy1 = nwxuy1_p;
   int nwxuy1_len = expansionObject::Gen_Product_With_PreAlloc(nwx_len, nwx, nuy_len, nuy, &nwxuy1, 16);
   double nwxuy2_p[16], *nwxuy2 = nwxuy2_p;
   int nwxuy2_len = expansionObject::Gen_Product_With_PreAlloc(nwy_len, nwy, nux_len, nux, &nwxuy2, 16);
   double nwxuy_p[16], *nwxuy = nwxuy_p;
   int nwxuy_len = expansionObject::Gen_Diff_With_PreAlloc(nwxuy1_len, nwxuy1, nwxuy2_len, nwxuy2, &nwxuy, 16);
   double nvyuz1_p[16], *nvyuz1 = nvyuz1_p;
   int nvyuz1_len = expansionObject::Gen_Product_With_PreAlloc(nvy_len, nvy, nuz_len, nuz, &nvyuz1, 16);
   double nvyuz2_p[16], *nvyuz2 = nvyuz2_p;
   int nvyuz2_len = expansionObject::Gen_Product_With_PreAlloc(nvz_len, nvz, nuy_len, nuy, &nvyuz2, 16);
   double nvyuz_p[16], *nvyuz = nvyuz_p;
   int nvyuz_len = expansionObject::Gen_Diff_With_PreAlloc(nvyuz1_len, nvyuz1, nvyuz2_len, nvyuz2, &nvyuz, 16);
   double nvxuz1_p[16], *nvxuz1 = nvxuz1_p;
   int nvxuz1_len = expansionObject::Gen_Product_With_PreAlloc(nvx_len, nvx, nuz_len, nuz, &nvxuz1, 16);
   double nvxuz2_p[16], *nvxuz2 = nvxuz2_p;
   int nvxuz2_len = expansionObject::Gen_Product_With_PreAlloc(nvz_len, nvz, nux_len, nux, &nvxuz2, 16);
   double nvxuz_p[16], *nvxuz = nvxuz_p;
   int nvxuz_len = expansionObject::Gen_Diff_With_PreAlloc(nvxuz1_len, nvxuz1, nvxuz2_len, nvxuz2, &nvxuz, 16);
   double nvxuy1_p[16], *nvxuy1 = nvxuy1_p;
   int nvxuy1_len = expansionObject::Gen_Product_With_PreAlloc(nvx_len, nvx, nuy_len, nuy, &nvxuy1, 16);
   double nvxuy2_p[16], *nvxuy2 = nvxuy2_p;
   int nvxuy2_len = expansionObject::Gen_Product_With_PreAlloc(nvy_len, nvy, nux_len, nux, &nvxuy2, 16);
   double nvxuy_p[16], *nvxuy = nvxuy_p;
   int nvxuy_len = expansionObject::Gen_Diff_With_PreAlloc(nvxuy1_len, nvxuy1, nvxuy2_len, nvxuy2, &nvxuy, 16);
   double nvywz1_p[16], *nvywz1 = nvywz1_p;
   int nvywz1_len = expansionObject::Gen_Product_With_PreAlloc(nvy_len, nvy, nwz_len, nwz, &nvywz1, 16);
   double nvywz2_p[16], *nvywz2 = nvywz2_p;
   int nvywz2_len = expansionObject::Gen_Product_With_PreAlloc(nvz_len, nvz, nwy_len, nwy, &nvywz2, 16);
   double nvywz_p[16], *nvywz = nvywz_p;
   int nvywz_len = expansionObject::Gen_Diff_With_PreAlloc(nvywz1_len, nvywz1, nvywz2_len, nvywz2, &nvywz, 16);
   double nvxwz1_p[16], *nvxwz1 = nvxwz1_p;
   int nvxwz1_len = expansionObject::Gen_Product_With_PreAlloc(nvx_len, nvx, nwz_len, nwz, &nvxwz1, 16);
   double nvxwz2_p[16], *nvxwz2 = nvxwz2_p;
   int nvxwz2_len = expansionObject::Gen_Product_With_PreAlloc(nvz_len, nvz, nwx_len, nwx, &nvxwz2, 16);
   double nvxwz_p[16], *nvxwz = nvxwz_p;
   int nvxwz_len = expansionObject::Gen_Diff_With_PreAlloc(nvxwz1_len, nvxwz1, nvxwz2_len, nvxwz2, &nvxwz, 16);
   double nvxwy1_p[16], *nvxwy1 = nvxwy1_p;
   int nvxwy1_len = expansionObject::Gen_Product_With_PreAlloc(nvx_len, nvx, nwy_len, nwy, &nvxwy1, 16);
   double nvxwy2_p[16], *nvxwy2 = nvxwy2_p;
   int nvxwy2_len = expansionObject::Gen_Product_With_PreAlloc(nvy_len, nvy, nwx_len, nwx, &nvxwy2, 16);
   double nvxwy_p[16], *nvxwy = nvxwy_p;
   int nvxwy_len = expansionObject::Gen_Diff_With_PreAlloc(nvxwy1_len, nvxwy1, nvxwy2_len, nvxwy2, &nvxwy, 16);
   double p1a_p[16], *p1a = p1a_p;
   int p1a_len = expansionObject::Gen_Scale_With_PreAlloc(nvx_len, nvx, ov1x, &p1a, 16);
   double p1b_p[16], *p1b = p1b_p;
   int p1b_len = expansionObject::Gen_Scale_With_PreAlloc(nvy_len, nvy, ov1y, &p1b, 16);
   double p1c_p[16], *p1c = p1c_p;
   int p1c_len = expansionObject::Gen_Scale_With_PreAlloc(nvz_len, nvz, ov1z, &p1c, 16);
   double p1ab_p[16], *p1ab = p1ab_p;
   int p1ab_len = expansionObject::Gen_Sum_With_PreAlloc(p1a_len, p1a, p1b_len, p1b, &p1ab, 16);
   double p1_p[16], *p1 = p1_p;
   int p1_len = expansionObject::Gen_Sum_With_PreAlloc(p1ab_len, p1ab, p1c_len, p1c, &p1, 16);
   double p2a_p[16], *p2a = p2a_p;
   int p2a_len = expansionObject::Gen_Scale_With_PreAlloc(nwx_len, nwx, ow1x, &p2a, 16);
   double p2b_p[16], *p2b = p2b_p;
   int p2b_len = expansionObject::Gen_Scale_With_PreAlloc(nwy_len, nwy, ow1y, &p2b, 16);
   double p2c_p[16], *p2c = p2c_p;
   int p2c_len = expansionObject::Gen_Scale_With_PreAlloc(nwz_len, nwz, ow1z, &p2c, 16);
   double p2ab_p[16], *p2ab = p2ab_p;
   int p2ab_len = expansionObject::Gen_Sum_With_PreAlloc(p2a_len, p2a, p2b_len, p2b, &p2ab, 16);
   double p2_p[16], *p2 = p2_p;
   int p2_len = expansionObject::Gen_Sum_With_PreAlloc(p2ab_len, p2ab, p2c_len, p2c, &p2, 16);
   double p3a_p[16], *p3a = p3a_p;
   int p3a_len = expansionObject::Gen_Scale_With_PreAlloc(nux_len, nux, ou1x, &p3a, 16);
   double p3b_p[16], *p3b = p3b_p;
   int p3b_len = expansionObject::Gen_Scale_With_PreAlloc(nuy_len, nuy, ou1y, &p3b, 16);
   double p3c_p[16], *p3c = p3c_p;
   int p3c_len = expansionObject::Gen_Scale_With_PreAlloc(nuz_len, nuz, ou1z, &p3c, 16);
   double p3ab_p[16], *p3ab = p3ab_p;
   int p3ab_len = expansionObject::Gen_Sum_With_PreAlloc(p3a_len, p3a, p3b_len, p3b, &p3ab, 16);
   double p3_p[16], *p3 = p3_p;
   int p3_len = expansionObject::Gen_Sum_With_PreAlloc(p3ab_len, p3ab, p3c_len, p3c, &p3, 16);
   double lxa_p[16], *lxa = lxa_p;
   int lxa_len = expansionObject::Gen_Product_With_PreAlloc(p1_len, p1, nwyuz_len, nwyuz, &lxa, 16);
   double lxb_p[16], *lxb = lxb_p;
   int lxb_len = expansionObject::Gen_Product_With_PreAlloc(p3_len, p3, nvywz_len, nvywz, &lxb, 16);
   double lxc_p[16], *lxc = lxc_p;
   int lxc_len = expansionObject::Gen_Product_With_PreAlloc(p2_len, p2, nvyuz_len, nvyuz, &lxc, 16);
   double lxab_p[16], *lxab = lxab_p;
   int lxab_len = expansionObject::Gen_Sum_With_PreAlloc(lxa_len, lxa, lxb_len, lxb, &lxab, 16);
   lambda_x_len = expansionObject::Gen_Diff_With_PreAlloc(lxab_len, lxab, lxc_len, lxc, lambda_x, lambda_x_len);
   double lya_p[16], *lya = lya_p;
   int lya_len = expansionObject::Gen_Product_With_PreAlloc(p2_len, p2, nvxuz_len, nvxuz, &lya, 16);
   double lyb_p[16], *lyb = lyb_p;
   int lyb_len = expansionObject::Gen_Product_With_PreAlloc(p3_len, p3, nvxwz_len, nvxwz, &lyb, 16);
   double lyc_p[16], *lyc = lyc_p;
   int lyc_len = expansionObject::Gen_Product_With_PreAlloc(p1_len, p1, nwxuz_len, nwxuz, &lyc, 16);
   double lybc_p[16], *lybc = lybc_p;
   int lybc_len = expansionObject::Gen_Sum_With_PreAlloc(lyc_len, lyc, lyb_len, lyb, &lybc, 16);
   lambda_y_len = expansionObject::Gen_Diff_With_PreAlloc(lya_len, lya, lybc_len, lybc, lambda_y, lambda_y_len);
   double lza_p[16], *lza = lza_p;
   int lza_len = expansionObject::Gen_Product_With_PreAlloc(p3_len, p3, nvxwy_len, nvxwy, &lza, 16);
   double lzb_p[16], *lzb = lzb_p;
   int lzb_len = expansionObject::Gen_Product_With_PreAlloc(p1_len, p1, nwxuy_len, nwxuy, &lzb, 16);
   double lzc_p[16], *lzc = lzc_p;
   int lzc_len = expansionObject::Gen_Product_With_PreAlloc(p2_len, p2, nvxuy_len, nvxuy, &lzc, 16);
   double lzab_p[16], *lzab = lzab_p;
   int lzab_len = expansionObject::Gen_Sum_With_PreAlloc(lza_len, lza, lzb_len, lzb, &lzab, 16);
   lambda_z_len = expansionObject::Gen_Diff_With_PreAlloc(lzab_len, lzab, lzc_len, lzc, lambda_z, lambda_z_len);
   double da_p[16], *da = da_p;
   int da_len = expansionObject::Gen_Product_With_PreAlloc(nvx_len, nvx, nwyuz_len, nwyuz, &da, 16);
   double db_p[16], *db = db_p;
   int db_len = expansionObject::Gen_Product_With_PreAlloc(nvz_len, nvz, nwxuy_len, nwxuy, &db, 16);
   double dc_p[16], *dc = dc_p;
   int dc_len = expansionObject::Gen_Product_With_PreAlloc(nvy_len, nvy, nwxuz_len, nwxuz, &dc, 16);
   double dab_p[16], *dab = dab_p;
   int dab_len = expansionObject::Gen_Sum_With_PreAlloc(da_len, da, db_len, db, &dab, 16);
   lambda_d_len = expansionObject::Gen_Diff_With_PreAlloc(dab_len, dab, dc_len, dc, lambda_d, lambda_d_len);

   if (dab_p != dab) FreeDoubles(dab);
   if (dc_p != dc) FreeDoubles(dc);
   if (db_p != db) FreeDoubles(db);
   if (da_p != da) FreeDoubles(da);
   if (lzab_p != lzab) FreeDoubles(lzab);
   if (lzc_p != lzc) FreeDoubles(lzc);
   if (lzb_p != lzb) FreeDoubles(lzb);
   if (lza_p != lza) FreeDoubles(lza);
   if (lybc_p != lybc) FreeDoubles(lybc);
   if (lyc_p != lyc) FreeDoubles(lyc);
   if (lyb_p != lyb) FreeDoubles(lyb);
   if (lya_p != lya) FreeDoubles(lya);
   if (lxab_p != lxab) FreeDoubles(lxab);
   if (lxc_p != lxc) FreeDoubles(lxc);
   if (lxb_p != lxb) FreeDoubles(lxb);
   if (lxa_p != lxa) FreeDoubles(lxa);
   if (p3_p != p3) FreeDoubles(p3);
   if (p3ab_p != p3ab) FreeDoubles(p3ab);
   if (p3c_p != p3c) FreeDoubles(p3c);
   if (p3b_p != p3b) FreeDoubles(p3b);
   if (p3a_p != p3a) FreeDoubles(p3a);
   if (p2_p != p2) FreeDoubles(p2);
   if (p2ab_p != p2ab) FreeDoubles(p2ab);
   if (p2c_p != p2c) FreeDoubles(p2c);
   if (p2b_p != p2b) FreeDoubles(p2b);
   if (p2a_p != p2a) FreeDoubles(p2a);
   if (p1_p != p1) FreeDoubles(p1);
   if (p1ab_p != p1ab) FreeDoubles(p1ab);
   if (p1c_p != p1c) FreeDoubles(p1c);
   if (p1b_p != p1b) FreeDoubles(p1b);
   if (p1a_p != p1a) FreeDoubles(p1a);
   if (nvxwy_p != nvxwy) FreeDoubles(nvxwy);
   if (nvxwy2_p != nvxwy2) FreeDoubles(nvxwy2);
   if (nvxwy1_p != nvxwy1) FreeDoubles(nvxwy1);
   if (nvxwz_p != nvxwz) FreeDoubles(nvxwz);
   if (nvxwz2_p != nvxwz2) FreeDoubles(nvxwz2);
   if (nvxwz1_p != nvxwz1) FreeDoubles(nvxwz1);
   if (nvywz_p != nvywz) FreeDoubles(nvywz);
   if (nvywz2_p != nvywz2) FreeDoubles(nvywz2);
   if (nvywz1_p != nvywz1) FreeDoubles(nvywz1);
   if (nvxuy_p != nvxuy) FreeDoubles(nvxuy);
   if (nvxuy2_p != nvxuy2) FreeDoubles(nvxuy2);
   if (nvxuy1_p != nvxuy1) FreeDoubles(nvxuy1);
   if (nvxuz_p != nvxuz) FreeDoubles(nvxuz);
   if (nvxuz2_p != nvxuz2) FreeDoubles(nvxuz2);
   if (nvxuz1_p != nvxuz1) FreeDoubles(nvxuz1);
   if (nvyuz_p != nvyuz) FreeDoubles(nvyuz);
   if (nvyuz2_p != nvyuz2) FreeDoubles(nvyuz2);
   if (nvyuz1_p != nvyuz1) FreeDoubles(nvyuz1);
   if (nwxuy_p != nwxuy) FreeDoubles(nwxuy);
   if (nwxuy2_p != nwxuy2) FreeDoubles(nwxuy2);
   if (nwxuy1_p != nwxuy1) FreeDoubles(nwxuy1);
   if (nwxuz_p != nwxuz) FreeDoubles(nwxuz);
   if (nwxuz2_p != nwxuz2) FreeDoubles(nwxuz2);
   if (nwxuz1_p != nwxuz1) FreeDoubles(nwxuz1);
   if (nwyuz_p != nwyuz) FreeDoubles(nwyuz);
   if (nwyuz2_p != nwyuz2) FreeDoubles(nwyuz2);
   if (nwyuz1_p != nwyuz1) FreeDoubles(nwyuz1);
}

inline int lessThanOnX_IE_interval(const genericPoint& p1, interval_number bx)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number dbx(bx * d1);
   const interval_number kx(l1x - dbx);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return kx.sign();
}

inline int lessThanOnX_IE_bigfloat(const genericPoint& p1, bigfloat bx)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat dbx(bx * d1);
   const bigfloat kx(l1x - dbx);
   return sgn(kx);
}

inline int lessThanOnX_IE_exact(const genericPoint& p1, double bx)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double dbx_p[128], *dbx = dbx_p;
   int dbx_len = expansionObject::Gen_Scale_With_PreAlloc(d1_len, d1, bx, &dbx, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = expansionObject::Gen_Diff_With_PreAlloc(l1x_len, l1x, dbx_len, dbx, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) FreeDoubles(kx);
   if (dbx_p != dbx) FreeDoubles(dbx);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = lessThanOnX_IE_bigfloat(p1, bx);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int lessThanOnX_IE(const genericPoint& p1, double bx)
{
   int ret;
   ret = lessThanOnX_IE_interval(p1, bx);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_IE_exact(p1, bx);
}

inline int lessThanOnX_II_interval(const genericPoint& p1, const genericPoint& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number k1(d2 * l1x);
   const interval_number k2(d1 * l2x);
   const interval_number kx(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return kx.sign();
}

inline int lessThanOnX_II_bigfloat(const genericPoint& p1, const genericPoint& p2)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat k1(d2 * l1x);
   const bigfloat k2(d1 * l2x);
   const bigfloat kx(k1 - k2);
   return sgn(kx);
}

inline int lessThanOnX_II_exact(const genericPoint& p1, const genericPoint& p2)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128, l2y_len = 128, l2z_len = 128, d2_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 p2.getExactLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   double k1_p[128], *k1 = k1_p;
   int k1_len = expansionObject::Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &k2, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = expansionObject::Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) FreeDoubles(kx);
   if (k2_p != k2) FreeDoubles(k2);
   if (k1_p != k1) FreeDoubles(k1);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = lessThanOnX_II_bigfloat(p1, p2);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int lessThanOnX_II(const genericPoint& p1, const genericPoint& p2)
{
   int ret;
   ret = lessThanOnX_II_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_II_exact(p1, p2);
}

inline int lessThanOnY_IE_interval(const genericPoint& p1, interval_number by)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number dby(by * d1);
   const interval_number ky(l1y - dby);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ky.sign();
}

inline int lessThanOnY_IE_bigfloat(const genericPoint& p1, bigfloat by)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat dby(by * d1);
   const bigfloat ky(l1y - dby);
   return sgn(ky);
}

inline int lessThanOnY_IE_exact(const genericPoint& p1, double by)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double dby_p[128], *dby = dby_p;
   int dby_len = expansionObject::Gen_Scale_With_PreAlloc(d1_len, d1, by, &dby, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = expansionObject::Gen_Diff_With_PreAlloc(l1y_len, l1y, dby_len, dby, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) FreeDoubles(ky);
   if (dby_p != dby) FreeDoubles(dby);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = lessThanOnY_IE_bigfloat(p1, by);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int lessThanOnY_IE(const genericPoint& p1, double by)
{
   int ret;
   ret = lessThanOnY_IE_interval(p1, by);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_IE_exact(p1, by);
}

inline int lessThanOnY_II_interval(const genericPoint& p1, const genericPoint& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number k1(d2 * l1y);
   const interval_number k2(d1 * l2y);
   const interval_number ky(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return ky.sign();
}

inline int lessThanOnY_II_bigfloat(const genericPoint& p1, const genericPoint& p2)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat k1(d2 * l1y);
   const bigfloat k2(d1 * l2y);
   const bigfloat ky(k1 - k2);
   return sgn(ky);
}

inline int lessThanOnY_II_exact(const genericPoint& p1, const genericPoint& p2)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128, l2y_len = 128, l2z_len = 128, d2_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 p2.getExactLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   double k1_p[128], *k1 = k1_p;
   int k1_len = expansionObject::Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &k2, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = expansionObject::Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) FreeDoubles(ky);
   if (k2_p != k2) FreeDoubles(k2);
   if (k1_p != k1) FreeDoubles(k1);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = lessThanOnY_II_bigfloat(p1, p2);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int lessThanOnY_II(const genericPoint& p1, const genericPoint& p2)
{
   int ret;
   ret = lessThanOnY_II_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_II_exact(p1, p2);
}

inline int lessThanOnZ_IE_interval(const genericPoint& p1, interval_number bz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number dbz(bz * d1);
   const interval_number kz(l1z - dbz);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return kz.sign();
}

inline int lessThanOnZ_IE_bigfloat(const genericPoint& p1, bigfloat bz)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat dbz(bz * d1);
   const bigfloat kz(l1z - dbz);
   return sgn(kz);
}

inline int lessThanOnZ_IE_exact(const genericPoint& p1, double bz)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double dbz_p[128], *dbz = dbz_p;
   int dbz_len = expansionObject::Gen_Scale_With_PreAlloc(d1_len, d1, bz, &dbz, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = expansionObject::Gen_Diff_With_PreAlloc(l1z_len, l1z, dbz_len, dbz, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) FreeDoubles(kz);
   if (dbz_p != dbz) FreeDoubles(dbz);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = lessThanOnZ_IE_bigfloat(p1, bz);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int lessThanOnZ_IE(const genericPoint& p1, double bz)
{
   int ret;
   ret = lessThanOnZ_IE_interval(p1, bz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_IE_exact(p1, bz);
}

inline int lessThanOnZ_II_interval(const genericPoint& p1, const genericPoint& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number k1(d2 * l1z);
   const interval_number k2(d1 * l2z);
   const interval_number kz(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return kz.sign();
}

inline int lessThanOnZ_II_bigfloat(const genericPoint& p1, const genericPoint& p2)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat k1(d2 * l1z);
   const bigfloat k2(d1 * l2z);
   const bigfloat kz(k1 - k2);
   return sgn(kz);
}

inline int lessThanOnZ_II_exact(const genericPoint& p1, const genericPoint& p2)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128, l2y_len = 128, l2z_len = 128, d2_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 p2.getExactLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   double k1_p[128], *k1 = k1_p;
   int k1_len = expansionObject::Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &k2, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = expansionObject::Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) FreeDoubles(kz);
   if (k2_p != k2) FreeDoubles(k2);
   if (k1_p != k1) FreeDoubles(k1);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = lessThanOnZ_II_bigfloat(p1, p2);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int lessThanOnZ_II(const genericPoint& p1, const genericPoint& p2)
{
   int ret;
   ret = lessThanOnZ_II_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_II_exact(p1, p2);
}

inline int orient2dxy_indirect_IEE_interval(const genericPoint& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number t1x(p2y - p3y);
   const interval_number t1y(p3x - p2x);
   const interval_number e2(l1x * t1x);
   const interval_number e3(l1y * t1y);
   const interval_number e(e2 + e3);
   const interval_number pr1(p2x * p3y);
   const interval_number pr2(p2y * p3x);
   const interval_number pr(pr1 - pr2);
   const interval_number dpr(d1 * pr);
   const interval_number det(dpr + e);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int orient2dxy_indirect_IEE_bigfloat(const genericPoint& p1, bigfloat p2x, bigfloat p2y, bigfloat p3x, bigfloat p3y)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat t1x(p2y - p3y);
   const bigfloat t1y(p3x - p2x);
   const bigfloat e2(l1x * t1x);
   const bigfloat e3(l1y * t1y);
   const bigfloat e(e2 + e3);
   const bigfloat pr1(p2x * p3y);
   const bigfloat pr2(p2y * p3x);
   const bigfloat pr(pr1 - pr2);
   const bigfloat dpr(d1 * pr);
   const bigfloat det(dpr + e);
   return sgn(det);
}

inline int orient2dxy_indirect_IEE_exact(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double t1x[2];
   expansionObject::two_Diff(p2y, p3y, t1x);
   double t1y[2];
   expansionObject::two_Diff(p3x, p2x, t1y);
   double e2_p[128], *e2 = e2_p;
   int e2_len = expansionObject::Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e2, 128);
   double e3_p[128], *e3 = e3_p;
   int e3_len = expansionObject::Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e3, 128);
   double e_p[128], *e = e_p;
   int e_len = expansionObject::Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
   double pr1[2];
   expansionObject::Two_Prod(p2x, p3y, pr1);
   double pr2[2];
   expansionObject::Two_Prod(p2y, p3x, pr2);
   double pr[4];
   expansionObject::Two_Two_Diff(pr1, pr2, pr);
   double dpr_p[128], *dpr = dpr_p;
   int dpr_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
   double det_p[128], *det = det_p;
   int det_len = expansionObject::Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

   return_value = det[det_len - 1];
   if (det_p != det) FreeDoubles(det);
   if (dpr_p != dpr) FreeDoubles(dpr);
   if (e_p != e) FreeDoubles(e);
   if (e3_p != e3) FreeDoubles(e3);
   if (e2_p != e2) FreeDoubles(e2);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = orient2dxy_indirect_IEE_bigfloat(p1, p2x, p2y, p3x, p3y);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int orient2dxy_indirect_IEE(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2dxy_indirect_IEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dxy_indirect_IEE_exact(p1, p2x, p2y, p3x, p3y);
}

inline int orient2dxy_indirect_IIE_interval(const genericPoint& p1, const genericPoint& p2, interval_number op3x, interval_number op3y)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * op3y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * op3x);
   const interval_number ab(a - b);
   const interval_number cd(c - l1y);
   const interval_number ef(e - f);
   const interval_number gh(g - l1x);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2dxy_indirect_IIE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat op3x, bigfloat op3y)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * op3y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * op3x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - l1y);
   const bigfloat ef(e - f);
   const bigfloat gh(g - l1x);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2dxy_indirect_IIE(const genericPoint& p1, const genericPoint& p2, double op3x, double op3y)
{
   int ret;
   ret = orient2dxy_indirect_IIE_interval(p1, p2, op3x, op3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dxy_indirect_IIE_bigfloat(p1, p2, op3x, op3y);
}

inline int orient2dxy_indirect_III_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * l3y);
   const interval_number d(d3 * l1y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * l3x);
   const interval_number h(d3 * l1x);
   const interval_number ab(a - b);
   const interval_number cd(c - d);
   const interval_number ef(e - f);
   const interval_number gh(g - h);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2dxy_indirect_III_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * l3y);
   const bigfloat d(d3 * l1y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * l3x);
   const bigfloat h(d3 * l1x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - d);
   const bigfloat ef(e - f);
   const bigfloat gh(g - h);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2dxy_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   int ret;
   ret = orient2dxy_indirect_III_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dxy_indirect_III_bigfloat(p1, p2, p3);
}

inline int orient2dyz_indirect_IEE_interval(const genericPoint& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1z, l1x, l1y, d1;
   if (
   !p1.getIntervalLambda(l1z, l1x, l1y, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number t1x(p2y - p3y);
   const interval_number t1y(p3x - p2x);
   const interval_number e2(l1x * t1x);
   const interval_number e3(l1y * t1y);
   const interval_number e(e2 + e3);
   const interval_number pr1(p2x * p3y);
   const interval_number pr2(p2y * p3x);
   const interval_number pr(pr1 - pr2);
   const interval_number dpr(d1 * pr);
   const interval_number det(dpr + e);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int orient2dyz_indirect_IEE_bigfloat(const genericPoint& p1, bigfloat p2x, bigfloat p2y, bigfloat p3x, bigfloat p3y)
{
   bigfloat l1z, l1x, l1y, d1;
   p1.getBigfloatLambda(l1z, l1x, l1y, d1);
   const bigfloat t1x(p2y - p3y);
   const bigfloat t1y(p3x - p2x);
   const bigfloat e2(l1x * t1x);
   const bigfloat e3(l1y * t1y);
   const bigfloat e(e2 + e3);
   const bigfloat pr1(p2x * p3y);
   const bigfloat pr2(p2y * p3x);
   const bigfloat pr(pr1 - pr2);
   const bigfloat dpr(d1 * pr);
   const bigfloat det(dpr + e);
   return sgn(det);
}

inline int orient2dyz_indirect_IEE_exact(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1z_p[128], *l1z = l1z_p, l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, d1_p[128], *d1 = d1_p;
 int l1z_len = 128, l1x_len = 128, l1y_len = 128, d1_len = 128;
 p1.getExactLambda(&l1z, l1z_len, &l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double t1x[2];
   expansionObject::two_Diff(p2y, p3y, t1x);
   double t1y[2];
   expansionObject::two_Diff(p3x, p2x, t1y);
   double e2_p[128], *e2 = e2_p;
   int e2_len = expansionObject::Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e2, 128);
   double e3_p[128], *e3 = e3_p;
   int e3_len = expansionObject::Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e3, 128);
   double e_p[128], *e = e_p;
   int e_len = expansionObject::Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
   double pr1[2];
   expansionObject::Two_Prod(p2x, p3y, pr1);
   double pr2[2];
   expansionObject::Two_Prod(p2y, p3x, pr2);
   double pr[4];
   expansionObject::Two_Two_Diff(pr1, pr2, pr);
   double dpr_p[128], *dpr = dpr_p;
   int dpr_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
   double det_p[128], *det = det_p;
   int det_len = expansionObject::Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

   return_value = det[det_len - 1];
   if (det_p != det) FreeDoubles(det);
   if (dpr_p != dpr) FreeDoubles(dpr);
   if (e_p != e) FreeDoubles(e);
   if (e3_p != e3) FreeDoubles(e3);
   if (e2_p != e2) FreeDoubles(e2);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = orient2dyz_indirect_IEE_bigfloat(p1, p2x, p2y, p3x, p3y);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int orient2dyz_indirect_IEE(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2dyz_indirect_IEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dyz_indirect_IEE_exact(p1, p2x, p2y, p3x, p3y);
}

inline int orient2dyz_indirect_IIE_interval(const genericPoint& p1, const genericPoint& p2, interval_number op3x, interval_number op3y)
{
   interval_number l1z, l1x, l1y, d1, l2z, l2x, l2y, d2;
   if (
   !p1.getIntervalLambda(l1z, l1x, l1y, d1)
   || !p2.getIntervalLambda(l2z, l2x, l2y, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * op3y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * op3x);
   const interval_number ab(a - b);
   const interval_number cd(c - l1y);
   const interval_number ef(e - f);
   const interval_number gh(g - l1x);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2dyz_indirect_IIE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat op3x, bigfloat op3y)
{
   bigfloat l1z, l1x, l1y, d1, l2z, l2x, l2y, d2;
   p1.getBigfloatLambda(l1z, l1x, l1y, d1);
   p2.getBigfloatLambda(l2z, l2x, l2y, d2);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * op3y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * op3x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - l1y);
   const bigfloat ef(e - f);
   const bigfloat gh(g - l1x);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2dyz_indirect_IIE(const genericPoint& p1, const genericPoint& p2, double op3x, double op3y)
{
   int ret;
   ret = orient2dyz_indirect_IIE_interval(p1, p2, op3x, op3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dyz_indirect_IIE_bigfloat(p1, p2, op3x, op3y);
}

inline int orient2dyz_indirect_III_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   interval_number l1z, l1x, l1y, d1, l2z, l2x, l2y, d2, l3z, l3x, l3y, d3;
   if (
   !p1.getIntervalLambda(l1z, l1x, l1y, d1)
   || !p2.getIntervalLambda(l2z, l2x, l2y, d2)
   || !p3.getIntervalLambda(l3z, l3x, l3y, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * l3y);
   const interval_number d(d3 * l1y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * l3x);
   const interval_number h(d3 * l1x);
   const interval_number ab(a - b);
   const interval_number cd(c - d);
   const interval_number ef(e - f);
   const interval_number gh(g - h);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2dyz_indirect_III_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   bigfloat l1z, l1x, l1y, d1, l2z, l2x, l2y, d2, l3z, l3x, l3y, d3;
   p1.getBigfloatLambda(l1z, l1x, l1y, d1);
   p2.getBigfloatLambda(l2z, l2x, l2y, d2);
   p3.getBigfloatLambda(l3z, l3x, l3y, d3);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * l3y);
   const bigfloat d(d3 * l1y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * l3x);
   const bigfloat h(d3 * l1x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - d);
   const bigfloat ef(e - f);
   const bigfloat gh(g - h);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2dyz_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   int ret;
   ret = orient2dyz_indirect_III_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dyz_indirect_III_bigfloat(p1, p2, p3);
}

inline int orient2dzx_indirect_IEE_interval(const genericPoint& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1y, l1z, l1x, d1;
   if (
   !p1.getIntervalLambda(l1y, l1z, l1x, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number t1x(p2y - p3y);
   const interval_number t1y(p3x - p2x);
   const interval_number e2(l1x * t1x);
   const interval_number e3(l1y * t1y);
   const interval_number e(e2 + e3);
   const interval_number pr1(p2x * p3y);
   const interval_number pr2(p2y * p3x);
   const interval_number pr(pr1 - pr2);
   const interval_number dpr(d1 * pr);
   const interval_number det(dpr + e);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int orient2dzx_indirect_IEE_bigfloat(const genericPoint& p1, bigfloat p2x, bigfloat p2y, bigfloat p3x, bigfloat p3y)
{
   bigfloat l1y, l1z, l1x, d1;
   p1.getBigfloatLambda(l1y, l1z, l1x, d1);
   const bigfloat t1x(p2y - p3y);
   const bigfloat t1y(p3x - p2x);
   const bigfloat e2(l1x * t1x);
   const bigfloat e3(l1y * t1y);
   const bigfloat e(e2 + e3);
   const bigfloat pr1(p2x * p3y);
   const bigfloat pr2(p2y * p3x);
   const bigfloat pr(pr1 - pr2);
   const bigfloat dpr(d1 * pr);
   const bigfloat det(dpr + e);
   return sgn(det);
}

inline int orient2dzx_indirect_IEE_exact(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, l1x_p[128], *l1x = l1x_p, d1_p[128], *d1 = d1_p;
 int l1y_len = 128, l1z_len = 128, l1x_len = 128, d1_len = 128;
 p1.getExactLambda(&l1y, l1y_len, &l1z, l1z_len, &l1x, l1x_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double t1x[2];
   expansionObject::two_Diff(p2y, p3y, t1x);
   double t1y[2];
   expansionObject::two_Diff(p3x, p2x, t1y);
   double e2_p[128], *e2 = e2_p;
   int e2_len = expansionObject::Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e2, 128);
   double e3_p[128], *e3 = e3_p;
   int e3_len = expansionObject::Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e3, 128);
   double e_p[128], *e = e_p;
   int e_len = expansionObject::Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
   double pr1[2];
   expansionObject::Two_Prod(p2x, p3y, pr1);
   double pr2[2];
   expansionObject::Two_Prod(p2y, p3x, pr2);
   double pr[4];
   expansionObject::Two_Two_Diff(pr1, pr2, pr);
   double dpr_p[128], *dpr = dpr_p;
   int dpr_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
   double det_p[128], *det = det_p;
   int det_len = expansionObject::Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

   return_value = det[det_len - 1];
   if (det_p != det) FreeDoubles(det);
   if (dpr_p != dpr) FreeDoubles(dpr);
   if (e_p != e) FreeDoubles(e);
   if (e3_p != e3) FreeDoubles(e3);
   if (e2_p != e2) FreeDoubles(e2);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = orient2dzx_indirect_IEE_bigfloat(p1, p2x, p2y, p3x, p3y);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int orient2dzx_indirect_IEE(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2dzx_indirect_IEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dzx_indirect_IEE_exact(p1, p2x, p2y, p3x, p3y);
}

inline int orient2dzx_indirect_IIE_interval(const genericPoint& p1, const genericPoint& p2, interval_number op3x, interval_number op3y)
{
   interval_number l1y, l1z, l1x, d1, l2y, l2z, l2x, d2;
   if (
   !p1.getIntervalLambda(l1y, l1z, l1x, d1)
   || !p2.getIntervalLambda(l2y, l2z, l2x, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * op3y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * op3x);
   const interval_number ab(a - b);
   const interval_number cd(c - l1y);
   const interval_number ef(e - f);
   const interval_number gh(g - l1x);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2dzx_indirect_IIE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat op3x, bigfloat op3y)
{
   bigfloat l1y, l1z, l1x, d1, l2y, l2z, l2x, d2;
   p1.getBigfloatLambda(l1y, l1z, l1x, d1);
   p2.getBigfloatLambda(l2y, l2z, l2x, d2);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * op3y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * op3x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - l1y);
   const bigfloat ef(e - f);
   const bigfloat gh(g - l1x);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2dzx_indirect_IIE(const genericPoint& p1, const genericPoint& p2, double op3x, double op3y)
{
   int ret;
   ret = orient2dzx_indirect_IIE_interval(p1, p2, op3x, op3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dzx_indirect_IIE_bigfloat(p1, p2, op3x, op3y);
}

inline int orient2dzx_indirect_III_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   interval_number l1y, l1z, l1x, d1, l2y, l2z, l2x, d2, l3y, l3z, l3x, d3;
   if (
   !p1.getIntervalLambda(l1y, l1z, l1x, d1)
   || !p2.getIntervalLambda(l2y, l2z, l2x, d2)
   || !p3.getIntervalLambda(l3y, l3z, l3x, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * l3y);
   const interval_number d(d3 * l1y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * l3x);
   const interval_number h(d3 * l1x);
   const interval_number ab(a - b);
   const interval_number cd(c - d);
   const interval_number ef(e - f);
   const interval_number gh(g - h);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2dzx_indirect_III_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   bigfloat l1y, l1z, l1x, d1, l2y, l2z, l2x, d2, l3y, l3z, l3x, d3;
   p1.getBigfloatLambda(l1y, l1z, l1x, d1);
   p2.getBigfloatLambda(l2y, l2z, l2x, d2);
   p3.getBigfloatLambda(l3y, l3z, l3x, d3);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * l3y);
   const bigfloat d(d3 * l1y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * l3x);
   const bigfloat h(d3 * l1x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - d);
   const bigfloat ef(e - f);
   const bigfloat gh(g - h);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2dzx_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   int ret;
   ret = orient2dzx_indirect_III_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2dzx_indirect_III_bigfloat(p1, p2, p3);
}

inline int orient2d_indirect_IEE_interval(const genericPoint& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number t1x(p2y - p3y);
   const interval_number t1y(p3x - p2x);
   const interval_number e2(l1x * t1x);
   const interval_number e3(l1y * t1y);
   const interval_number e(e2 + e3);
   const interval_number pr1(p2x * p3y);
   const interval_number pr2(p2y * p3x);
   const interval_number pr(pr1 - pr2);
   const interval_number dpr(d1 * pr);
   const interval_number det(dpr + e);
   setFPUModeToRoundNEAR();

   if (!det.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return det.sign();
}

inline int orient2d_indirect_IEE_bigfloat(const genericPoint& p1, bigfloat p2x, bigfloat p2y, bigfloat p3x, bigfloat p3y)
{
   bigfloat l1x, l1y, d1;
   p1.getBigfloatLambda(l1x, l1y, d1);
   const bigfloat t1x(p2y - p3y);
   const bigfloat t1y(p3x - p2x);
   const bigfloat e2(l1x * t1x);
   const bigfloat e3(l1y * t1y);
   const bigfloat e(e2 + e3);
   const bigfloat pr1(p2x * p3y);
   const bigfloat pr2(p2y * p3x);
   const bigfloat pr(pr1 - pr2);
   const bigfloat dpr(d1 * pr);
   const bigfloat det(dpr + e);
   return sgn(det);
}

inline int orient2d_indirect_IEE_exact(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = NAN;
#ifdef CHECK_FOR_XYZERFLOWS
   feclearexcept(FE_ALL_EXCEPT);
#endif
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, d1_p[128], *d1 = d1_p;
 int l1x_len = 128, l1y_len = 128, d1_len = 128;
 p1.getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   double t1x[2];
   expansionObject::two_Diff(p2y, p3y, t1x);
   double t1y[2];
   expansionObject::two_Diff(p3x, p2x, t1y);
   double e2_p[128], *e2 = e2_p;
   int e2_len = expansionObject::Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e2, 128);
   double e3_p[128], *e3 = e3_p;
   int e3_len = expansionObject::Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e3, 128);
   double e_p[128], *e = e_p;
   int e_len = expansionObject::Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
   double pr1[2];
   expansionObject::Two_Prod(p2x, p3y, pr1);
   double pr2[2];
   expansionObject::Two_Prod(p2y, p3x, pr2);
   double pr[4];
   expansionObject::Two_Two_Diff(pr1, pr2, pr);
   double dpr_p[128], *dpr = dpr_p;
   int dpr_len = expansionObject::Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
   double det_p[128], *det = det_p;
   int det_len = expansionObject::Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

   return_value = det[det_len - 1];
   if (det_p != det) FreeDoubles(det);
   if (dpr_p != dpr) FreeDoubles(dpr);
   if (e_p != e) FreeDoubles(e);
   if (e3_p != e3) FreeDoubles(e3);
   if (e2_p != e2) FreeDoubles(e2);
 }

#ifdef CHECK_FOR_XYZERFLOWS
   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = orient2d_indirect_IEE_bigfloat(p1, p2x, p2y, p3x, p3y);
#endif


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int orient2d_indirect_IEE(const genericPoint& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2d_indirect_IEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_indirect_IEE_exact(p1, p2x, p2y, p3x, p3y);
}

inline int orient2d_indirect_IIE_interval(const genericPoint& p1, const genericPoint& p2, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * p3y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * p3x);
   const interval_number ab(a - b);
   const interval_number cd(c - l1y);
   const interval_number ef(e - f);
   const interval_number gh(g - l1x);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2d_indirect_IIE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat p3x, bigfloat p3y)
{
   bigfloat l1x, l1y, d1, l2x, l2y, d2;
   p1.getBigfloatLambda(l1x, l1y, d1);
   p2.getBigfloatLambda(l2x, l2y, d2);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * p3y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * p3x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - l1y);
   const bigfloat ef(e - f);
   const bigfloat gh(g - l1x);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2d_indirect_IIE(const genericPoint& p1, const genericPoint& p2, double p3x, double p3y)
{
   int ret;
   ret = orient2d_indirect_IIE_interval(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_indirect_IIE_bigfloat(p1, p2, p3x, p3y);
}

inline int orient2d_indirect_III_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   || !p3.getIntervalLambda(l3x, l3y, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number a(d1 * l2x);
   const interval_number b(d2 * l1x);
   const interval_number c(d1 * l3y);
   const interval_number d(d3 * l1y);
   const interval_number e(d1 * l2y);
   const interval_number f(d2 * l1y);
   const interval_number g(d1 * l3x);
   const interval_number h(d3 * l1x);
   const interval_number ab(a - b);
   const interval_number cd(c - d);
   const interval_number ef(e - f);
   const interval_number gh(g - h);
   const interval_number abcd(ab * cd);
   const interval_number efgh(ef * gh);
   const interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

inline int orient2d_indirect_III_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   bigfloat l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
   p1.getBigfloatLambda(l1x, l1y, d1);
   p2.getBigfloatLambda(l2x, l2y, d2);
   p3.getBigfloatLambda(l3x, l3y, d3);
   const bigfloat a(d1 * l2x);
   const bigfloat b(d2 * l1x);
   const bigfloat c(d1 * l3y);
   const bigfloat d(d3 * l1y);
   const bigfloat e(d1 * l2y);
   const bigfloat f(d2 * l1y);
   const bigfloat g(d1 * l3x);
   const bigfloat h(d3 * l1x);
   const bigfloat ab(a - b);
   const bigfloat cd(c - d);
   const bigfloat ef(e - f);
   const bigfloat gh(g - h);
   const bigfloat abcd(ab * cd);
   const bigfloat efgh(ef * gh);
   const bigfloat L(abcd - efgh);
   return sgn(L);
}

inline int orient2d_indirect_III(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3)
{
   int ret;
   ret = orient2d_indirect_III_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_indirect_III_bigfloat(p1, p2, p3);
}

inline int orient3d_indirect_IEEE_interval(const genericPoint& p1, interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number dcx(d1 * cx);
   const interval_number dcy(d1 * cy);
   const interval_number dcz(d1 * cz);
   const interval_number ix_cx(l1x - dcx);
   const interval_number iy_cy(l1y - dcy);
   const interval_number ax_cx(ax - cx);
   const interval_number ay_cy(ay - cy);
   const interval_number az_cz(az - cz);
   const interval_number iz_cz(l1z - dcz);
   const interval_number bx_cx(bx - cx);
   const interval_number by_cy(by - cy);
   const interval_number bz_cz(bz - cz);
   const interval_number tmc_a(ix_cx * ay_cy);
   const interval_number tmc_b(iy_cy * ax_cx);
   const interval_number m01(tmc_a - tmc_b);
   const interval_number tmi_a(ix_cx * az_cz);
   const interval_number tmi_b(iz_cz * ax_cx);
   const interval_number m02(tmi_a - tmi_b);
   const interval_number tma_a(iy_cy * az_cz);
   const interval_number tma_b(iz_cz * ay_cy);
   const interval_number m12(tma_a - tma_b);
   const interval_number mt1(m01 * bz_cz);
   const interval_number mt2(m02 * by_cy);
   const interval_number mt3(m12 * bx_cx);
   const interval_number mtt(mt2 - mt1);
   const interval_number m012(mtt - mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return m012.sign();
}

inline int orient3d_indirect_IEEE_bigfloat(const genericPoint& p1, bigfloat ax, bigfloat ay, bigfloat az, bigfloat bx, bigfloat by, bigfloat bz, bigfloat cx, bigfloat cy, bigfloat cz)
{
   bigfloat l1x, l1y, l1z, d1;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   const bigfloat dcx(d1 * cx);
   const bigfloat dcy(d1 * cy);
   const bigfloat dcz(d1 * cz);
   const bigfloat ix_cx(l1x - dcx);
   const bigfloat iy_cy(l1y - dcy);
   const bigfloat ax_cx(ax - cx);
   const bigfloat ay_cy(ay - cy);
   const bigfloat az_cz(az - cz);
   const bigfloat iz_cz(l1z - dcz);
   const bigfloat bx_cx(bx - cx);
   const bigfloat by_cy(by - cy);
   const bigfloat bz_cz(bz - cz);
   const bigfloat tmc_a(ix_cx * ay_cy);
   const bigfloat tmc_b(iy_cy * ax_cx);
   const bigfloat m01(tmc_a - tmc_b);
   const bigfloat tmi_a(ix_cx * az_cz);
   const bigfloat tmi_b(iz_cz * ax_cx);
   const bigfloat m02(tmi_a - tmi_b);
   const bigfloat tma_a(iy_cy * az_cz);
   const bigfloat tma_b(iz_cz * ay_cy);
   const bigfloat m12(tma_a - tma_b);
   const bigfloat mt1(m01 * bz_cz);
   const bigfloat mt2(m02 * by_cy);
   const bigfloat mt3(m12 * bx_cx);
   const bigfloat mtt(mt2 - mt1);
   const bigfloat m012(mtt - mt3);
   return sgn(m012);
}

inline int orient3d_indirect_IEEE(const genericPoint& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = orient3d_indirect_IEEE_interval(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_indirect_IEEE_bigfloat(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline int orient3d_indirect_IIEE_interval(const genericPoint& p1, const genericPoint& p2, interval_number p3x, interval_number p3y, interval_number p3z, interval_number p4x, interval_number p4y, interval_number p4z)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d1p4x(d1 * p4x);
   const interval_number d1p4y(d1 * p4y);
   const interval_number d1p4z(d1 * p4z);
   const interval_number d2p4x(d2 * p4x);
   const interval_number d2p4y(d2 * p4y);
   const interval_number d2p4z(d2 * p4z);
   const interval_number p1p4x(l1x - d1p4x);
   const interval_number p1p4y(l1y - d1p4y);
   const interval_number p1p4z(l1z - d1p4z);
   const interval_number p2p4x(l2x - d2p4x);
   const interval_number p2p4y(l2y - d2p4y);
   const interval_number p2p4z(l2z - d2p4z);
   const interval_number p3p4x(p3x - p4x);
   const interval_number p3p4y(p3y - p4y);
   const interval_number p3p4z(p3z - p4z);
   const interval_number tmc_a(p1p4x * p2p4y);
   const interval_number tmc_b(p1p4y * p2p4x);
   const interval_number m01(tmc_a - tmc_b);
   const interval_number tmi_a(p1p4x * p2p4z);
   const interval_number tmi_b(p1p4z * p2p4x);
   const interval_number m02(tmi_a - tmi_b);
   const interval_number tma_a(p1p4y * p2p4z);
   const interval_number tma_b(p1p4z * p2p4y);
   const interval_number m12(tma_a - tma_b);
   const interval_number mt1(m01 * p3p4z);
   const interval_number mt2(m02 * p3p4y);
   const interval_number mt3(m12 * p3p4x);
   const interval_number mtt(mt2 - mt1);
   const interval_number m012(mtt - mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return m012.sign();
}

inline int orient3d_indirect_IIEE_bigfloat(const genericPoint& p1, const genericPoint& p2, bigfloat p3x, bigfloat p3y, bigfloat p3z, bigfloat p4x, bigfloat p4y, bigfloat p4z)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   const bigfloat d1p4x(d1 * p4x);
   const bigfloat d1p4y(d1 * p4y);
   const bigfloat d1p4z(d1 * p4z);
   const bigfloat d2p4x(d2 * p4x);
   const bigfloat d2p4y(d2 * p4y);
   const bigfloat d2p4z(d2 * p4z);
   const bigfloat p1p4x(l1x - d1p4x);
   const bigfloat p1p4y(l1y - d1p4y);
   const bigfloat p1p4z(l1z - d1p4z);
   const bigfloat p2p4x(l2x - d2p4x);
   const bigfloat p2p4y(l2y - d2p4y);
   const bigfloat p2p4z(l2z - d2p4z);
   const bigfloat p3p4x(p3x - p4x);
   const bigfloat p3p4y(p3y - p4y);
   const bigfloat p3p4z(p3z - p4z);
   const bigfloat tmc_a(p1p4x * p2p4y);
   const bigfloat tmc_b(p1p4y * p2p4x);
   const bigfloat m01(tmc_a - tmc_b);
   const bigfloat tmi_a(p1p4x * p2p4z);
   const bigfloat tmi_b(p1p4z * p2p4x);
   const bigfloat m02(tmi_a - tmi_b);
   const bigfloat tma_a(p1p4y * p2p4z);
   const bigfloat tma_b(p1p4z * p2p4y);
   const bigfloat m12(tma_a - tma_b);
   const bigfloat mt1(m01 * p3p4z);
   const bigfloat mt2(m02 * p3p4y);
   const bigfloat mt3(m12 * p3p4x);
   const bigfloat mtt(mt2 - mt1);
   const bigfloat m012(mtt - mt3);
   return sgn(m012);
}

inline int orient3d_indirect_IIEE(const genericPoint& p1, const genericPoint& p2, double p3x, double p3y, double p3z, double p4x, double p4y, double p4z)
{
   int ret;
   ret = orient3d_indirect_IIEE_interval(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_indirect_IIEE_bigfloat(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
}

inline int orient3d_indirect_IIIE_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, interval_number p4x, interval_number p4y, interval_number p4z)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d1p4x(d1 * p4x);
   const interval_number d1p4y(d1 * p4y);
   const interval_number d1p4z(d1 * p4z);
   const interval_number d2p4x(d2 * p4x);
   const interval_number d2p4y(d2 * p4y);
   const interval_number d2p4z(d2 * p4z);
   const interval_number d3p4x(d3 * p4x);
   const interval_number d3p4y(d3 * p4y);
   const interval_number d3p4z(d3 * p4z);
   const interval_number p1p4x(l1x - d1p4x);
   const interval_number p1p4y(l1y - d1p4y);
   const interval_number p1p4z(l1z - d1p4z);
   const interval_number p2p4x(l2x - d2p4x);
   const interval_number p2p4y(l2y - d2p4y);
   const interval_number p2p4z(l2z - d2p4z);
   const interval_number p3p4x(l3x - d3p4x);
   const interval_number p3p4y(l3y - d3p4y);
   const interval_number p3p4z(l3z - d3p4z);
   const interval_number tmc_a(p1p4x * p2p4y);
   const interval_number tmc_b(p1p4y * p2p4x);
   const interval_number m01(tmc_a - tmc_b);
   const interval_number tmi_a(p1p4x * p2p4z);
   const interval_number tmi_b(p1p4z * p2p4x);
   const interval_number m02(tmi_a - tmi_b);
   const interval_number tma_a(p1p4y * p2p4z);
   const interval_number tma_b(p1p4z * p2p4y);
   const interval_number m12(tma_a - tma_b);
   const interval_number mt1(m01 * p3p4z);
   const interval_number mt2(m02 * p3p4y);
   const interval_number mt3(m12 * p3p4x);
   const interval_number mtt(mt2 - mt1);
   const interval_number m012(mtt - mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return m012.sign();
}

inline int orient3d_indirect_IIIE_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, bigfloat p4x, bigfloat p4y, bigfloat p4z)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   const bigfloat d1p4x(d1 * p4x);
   const bigfloat d1p4y(d1 * p4y);
   const bigfloat d1p4z(d1 * p4z);
   const bigfloat d2p4x(d2 * p4x);
   const bigfloat d2p4y(d2 * p4y);
   const bigfloat d2p4z(d2 * p4z);
   const bigfloat d3p4x(d3 * p4x);
   const bigfloat d3p4y(d3 * p4y);
   const bigfloat d3p4z(d3 * p4z);
   const bigfloat p1p4x(l1x - d1p4x);
   const bigfloat p1p4y(l1y - d1p4y);
   const bigfloat p1p4z(l1z - d1p4z);
   const bigfloat p2p4x(l2x - d2p4x);
   const bigfloat p2p4y(l2y - d2p4y);
   const bigfloat p2p4z(l2z - d2p4z);
   const bigfloat p3p4x(l3x - d3p4x);
   const bigfloat p3p4y(l3y - d3p4y);
   const bigfloat p3p4z(l3z - d3p4z);
   const bigfloat tmc_a(p1p4x * p2p4y);
   const bigfloat tmc_b(p1p4y * p2p4x);
   const bigfloat m01(tmc_a - tmc_b);
   const bigfloat tmi_a(p1p4x * p2p4z);
   const bigfloat tmi_b(p1p4z * p2p4x);
   const bigfloat m02(tmi_a - tmi_b);
   const bigfloat tma_a(p1p4y * p2p4z);
   const bigfloat tma_b(p1p4z * p2p4y);
   const bigfloat m12(tma_a - tma_b);
   const bigfloat mt1(m01 * p3p4z);
   const bigfloat mt2(m02 * p3p4y);
   const bigfloat mt3(m12 * p3p4x);
   const bigfloat mtt(mt2 - mt1);
   const bigfloat m012(mtt - mt3);
   return sgn(m012);
}

inline int orient3d_indirect_IIIE(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, double p4x, double p4y, double p4z)
{
   int ret;
   ret = orient3d_indirect_IIIE_interval(p1, p2, p3, p4x, p4y, p4z);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_indirect_IIIE_bigfloat(p1, p2, p3, p4x, p4y, p4z);
}

inline int orient3d_indirect_IIII_interval(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   || !p4.getIntervalLambda(l4x, l4y, l4z, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   const interval_number d1p4x(d1 * l4x);
   const interval_number d1p4y(d1 * l4y);
   const interval_number d1p4z(d1 * l4z);
   const interval_number d2p4x(d2 * l4x);
   const interval_number d2p4y(d2 * l4y);
   const interval_number d2p4z(d2 * l4z);
   const interval_number d3p4x(d3 * l4x);
   const interval_number d3p4y(d3 * l4y);
   const interval_number d3p4z(d3 * l4z);
   const interval_number d4l1x(d4 * l1x);
   const interval_number d4l1y(d4 * l1y);
   const interval_number d4l1z(d4 * l1z);
   const interval_number d4l2x(d4 * l2x);
   const interval_number d4l2y(d4 * l2y);
   const interval_number d4l2z(d4 * l2z);
   const interval_number d4l3x(d4 * l3x);
   const interval_number d4l3y(d4 * l3y);
   const interval_number d4l3z(d4 * l3z);
   const interval_number p1p4x(d4l1x - d1p4x);
   const interval_number p1p4y(d4l1y - d1p4y);
   const interval_number p1p4z(d4l1z - d1p4z);
   const interval_number p2p4x(d4l2x - d2p4x);
   const interval_number p2p4y(d4l2y - d2p4y);
   const interval_number p2p4z(d4l2z - d2p4z);
   const interval_number p3p4x(d4l3x - d3p4x);
   const interval_number p3p4y(d4l3y - d3p4y);
   const interval_number p3p4z(d4l3z - d3p4z);
   const interval_number tmc_a(p1p4x * p2p4y);
   const interval_number tmc_b(p1p4y * p2p4x);
   const interval_number m01(tmc_a - tmc_b);
   const interval_number tmi_a(p1p4x * p2p4z);
   const interval_number tmi_b(p1p4z * p2p4x);
   const interval_number m02(tmi_a - tmi_b);
   const interval_number tma_a(p1p4y * p2p4z);
   const interval_number tma_b(p1p4z * p2p4y);
   const interval_number m12(tma_a - tma_b);
   const interval_number mt1(m01 * p3p4z);
   const interval_number mt2(m02 * p3p4y);
   const interval_number mt3(m12 * p3p4x);
   const interval_number mtt(mt2 - mt1);
   const interval_number m012(mtt - mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return m012.sign();
}

inline int orient3d_indirect_IIII_bigfloat(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   bigfloat l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   p1.getBigfloatLambda(l1x, l1y, l1z, d1);
   p2.getBigfloatLambda(l2x, l2y, l2z, d2);
   p3.getBigfloatLambda(l3x, l3y, l3z, d3);
   p4.getBigfloatLambda(l4x, l4y, l4z, d4);
   const bigfloat d1p4x(d1 * l4x);
   const bigfloat d1p4y(d1 * l4y);
   const bigfloat d1p4z(d1 * l4z);
   const bigfloat d2p4x(d2 * l4x);
   const bigfloat d2p4y(d2 * l4y);
   const bigfloat d2p4z(d2 * l4z);
   const bigfloat d3p4x(d3 * l4x);
   const bigfloat d3p4y(d3 * l4y);
   const bigfloat d3p4z(d3 * l4z);
   const bigfloat d4l1x(d4 * l1x);
   const bigfloat d4l1y(d4 * l1y);
   const bigfloat d4l1z(d4 * l1z);
   const bigfloat d4l2x(d4 * l2x);
   const bigfloat d4l2y(d4 * l2y);
   const bigfloat d4l2z(d4 * l2z);
   const bigfloat d4l3x(d4 * l3x);
   const bigfloat d4l3y(d4 * l3y);
   const bigfloat d4l3z(d4 * l3z);
   const bigfloat p1p4x(d4l1x - d1p4x);
   const bigfloat p1p4y(d4l1y - d1p4y);
   const bigfloat p1p4z(d4l1z - d1p4z);
   const bigfloat p2p4x(d4l2x - d2p4x);
   const bigfloat p2p4y(d4l2y - d2p4y);
   const bigfloat p2p4z(d4l2z - d2p4z);
   const bigfloat p3p4x(d4l3x - d3p4x);
   const bigfloat p3p4y(d4l3y - d3p4y);
   const bigfloat p3p4z(d4l3z - d3p4z);
   const bigfloat tmc_a(p1p4x * p2p4y);
   const bigfloat tmc_b(p1p4y * p2p4x);
   const bigfloat m01(tmc_a - tmc_b);
   const bigfloat tmi_a(p1p4x * p2p4z);
   const bigfloat tmi_b(p1p4z * p2p4x);
   const bigfloat m02(tmi_a - tmi_b);
   const bigfloat tma_a(p1p4y * p2p4z);
   const bigfloat tma_b(p1p4z * p2p4y);
   const bigfloat m12(tma_a - tma_b);
   const bigfloat mt1(m01 * p3p4z);
   const bigfloat mt2(m02 * p3p4y);
   const bigfloat mt3(m12 * p3p4x);
   const bigfloat mtt(mt2 - mt1);
   const bigfloat m012(mtt - mt3);
   return sgn(m012);
}

inline int orient3d_indirect_IIII(const genericPoint& p1, const genericPoint& p2, const genericPoint& p3, const genericPoint& p4)
{
   int ret;
   ret = orient3d_indirect_IIII_interval(p1, p2, p3, p4);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_indirect_IIII_bigfloat(p1, p2, p3, p4);
}

