#include "implicit_point.h"

#pragma intrinsic(fabs)

// Uncomment the following to activate overflow/underflow checks
#define CHECK_FOR_XYZERFLOWS

inline int inSphere_LLLLE_filtered(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double ex, double ey, double ez)
{
   const double avx = a2x - a1x;
   const double avxt = avx * at;
   const double ax = a1x + avxt;
   const double avy = a2y - a1y;
   const double avyt = avy * at;
   const double ay = a1y + avyt;
   const double avz = a2z - a1z;
   const double avzt = avz * at;
   const double az = a1z + avzt;
   const double bvx = b2x - b1x;
   const double bvxt = bvx * bt;
   const double bx = b1x + bvxt;
   const double bvy = b2y - b1y;
   const double bvyt = bvy * bt;
   const double by = b1y + bvyt;
   const double bvz = b2z - b1z;
   const double bvzt = bvz * bt;
   const double bz = b1z + bvzt;
   const double cvx = c2x - c1x;
   const double cvxt = cvx * ct;
   const double cx = c1x + cvxt;
   const double cvy = c2y - c1y;
   const double cvyt = cvy * ct;
   const double cy = c1y + cvyt;
   const double cvz = c2z - c1z;
   const double cvzt = cvz * ct;
   const double cz = c1z + cvzt;
   const double dvx = d2x - d1x;
   const double dvxt = dvx * dt;
   const double dx = d1x + dvxt;
   const double dvy = d2y - d1y;
   const double dvyt = dvy * dt;
   const double dy = d1y + dvyt;
   const double dvz = d2z - d1z;
   const double dvzt = dvz * dt;
   const double dz = d1z + dvzt;
   const double aex = ax - ex;
   const double aey = ay - ey;
   const double aez = az - ez;
   const double bex = bx - ex;
   const double bey = by - ey;
   const double bez = bz - ez;
   const double cex = cx - ex;
   const double cey = cy - ey;
   const double cez = cz - ez;
   const double dex = dx - ex;
   const double dey = dy - ey;
   const double dez = dz - ez;
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
   if ((_tmp_fabs = fabs(a1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(at)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(b1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(b1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(b1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bt)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(c1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(c1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(c1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ct)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(d1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(d1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(d1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dt)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ex)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ey)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ez)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(avx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(avy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(avz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bvx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bvy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bvz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cvx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cvy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cvz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dvx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dvy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dvz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.581668295282733e-11;
   if (det > epsilon) return IP_Sign::POSITIVE;
   if (-det > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int inSphere_LLLLE_interval(interval_number a1x, interval_number a1y, interval_number a1z, interval_number a2x, interval_number a2y, interval_number a2z, interval_number at, interval_number b1x, interval_number b1y, interval_number b1z, interval_number b2x, interval_number b2y, interval_number b2z, interval_number bt, interval_number c1x, interval_number c1y, interval_number c1z, interval_number c2x, interval_number c2y, interval_number c2z, interval_number ct, interval_number d1x, interval_number d1y, interval_number d1z, interval_number d2x, interval_number d2y, interval_number d2z, interval_number dt, interval_number ex, interval_number ey, interval_number ez)
{
   setFPUModeToRoundUP();
   const interval_number avx(a2x - a1x);
   const interval_number avxt(avx * at);
   const interval_number ax(a1x + avxt);
   const interval_number avy(a2y - a1y);
   const interval_number avyt(avy * at);
   const interval_number ay(a1y + avyt);
   const interval_number avz(a2z - a1z);
   const interval_number avzt(avz * at);
   const interval_number az(a1z + avzt);
   const interval_number bvx(b2x - b1x);
   const interval_number bvxt(bvx * bt);
   const interval_number bx(b1x + bvxt);
   const interval_number bvy(b2y - b1y);
   const interval_number bvyt(bvy * bt);
   const interval_number by(b1y + bvyt);
   const interval_number bvz(b2z - b1z);
   const interval_number bvzt(bvz * bt);
   const interval_number bz(b1z + bvzt);
   const interval_number cvx(c2x - c1x);
   const interval_number cvxt(cvx * ct);
   const interval_number cx(c1x + cvxt);
   const interval_number cvy(c2y - c1y);
   const interval_number cvyt(cvy * ct);
   const interval_number cy(c1y + cvyt);
   const interval_number cvz(c2z - c1z);
   const interval_number cvzt(cvz * ct);
   const interval_number cz(c1z + cvzt);
   const interval_number dvx(d2x - d1x);
   const interval_number dvxt(dvx * dt);
   const interval_number dx(d1x + dvxt);
   const interval_number dvy(d2y - d1y);
   const interval_number dvyt(dvy * dt);
   const interval_number dy(d1y + dvyt);
   const interval_number dvz(d2z - d1z);
   const interval_number dvzt(dvz * dt);
   const interval_number dz(d1z + dvzt);
   const interval_number aex(ax - ex);
   const interval_number aey(ay - ey);
   const interval_number aez(az - ez);
   const interval_number bex(bx - ex);
   const interval_number bey(by - ey);
   const interval_number bez(bz - ez);
   const interval_number cex(cx - ex);
   const interval_number cey(cy - ey);
   const interval_number cez(cz - ez);
   const interval_number dex(dx - ex);
   const interval_number dey(dy - ey);
   const interval_number dez(dz - ez);
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

inline int inSphere_LLLLE_bigfloat(bigfloat a1x, bigfloat a1y, bigfloat a1z, bigfloat a2x, bigfloat a2y, bigfloat a2z, bigfloat at, bigfloat b1x, bigfloat b1y, bigfloat b1z, bigfloat b2x, bigfloat b2y, bigfloat b2z, bigfloat bt, bigfloat c1x, bigfloat c1y, bigfloat c1z, bigfloat c2x, bigfloat c2y, bigfloat c2z, bigfloat ct, bigfloat d1x, bigfloat d1y, bigfloat d1z, bigfloat d2x, bigfloat d2y, bigfloat d2z, bigfloat dt, bigfloat ex, bigfloat ey, bigfloat ez)
{
   const bigfloat avx(a2x - a1x);
   const bigfloat avxt(avx * at);
   const bigfloat ax(a1x + avxt);
   const bigfloat avy(a2y - a1y);
   const bigfloat avyt(avy * at);
   const bigfloat ay(a1y + avyt);
   const bigfloat avz(a2z - a1z);
   const bigfloat avzt(avz * at);
   const bigfloat az(a1z + avzt);
   const bigfloat bvx(b2x - b1x);
   const bigfloat bvxt(bvx * bt);
   const bigfloat bx(b1x + bvxt);
   const bigfloat bvy(b2y - b1y);
   const bigfloat bvyt(bvy * bt);
   const bigfloat by(b1y + bvyt);
   const bigfloat bvz(b2z - b1z);
   const bigfloat bvzt(bvz * bt);
   const bigfloat bz(b1z + bvzt);
   const bigfloat cvx(c2x - c1x);
   const bigfloat cvxt(cvx * ct);
   const bigfloat cx(c1x + cvxt);
   const bigfloat cvy(c2y - c1y);
   const bigfloat cvyt(cvy * ct);
   const bigfloat cy(c1y + cvyt);
   const bigfloat cvz(c2z - c1z);
   const bigfloat cvzt(cvz * ct);
   const bigfloat cz(c1z + cvzt);
   const bigfloat dvx(d2x - d1x);
   const bigfloat dvxt(dvx * dt);
   const bigfloat dx(d1x + dvxt);
   const bigfloat dvy(d2y - d1y);
   const bigfloat dvyt(dvy * dt);
   const bigfloat dy(d1y + dvyt);
   const bigfloat dvz(d2z - d1z);
   const bigfloat dvzt(dvz * dt);
   const bigfloat dz(d1z + dvzt);
   const bigfloat aex(ax - ex);
   const bigfloat aey(ay - ey);
   const bigfloat aez(az - ez);
   const bigfloat bex(bx - ex);
   const bigfloat bey(by - ey);
   const bigfloat bez(bz - ez);
   const bigfloat cex(cx - ex);
   const bigfloat cey(cy - ey);
   const bigfloat cez(cz - ez);
   const bigfloat dex(dx - ex);
   const bigfloat dey(dy - ey);
   const bigfloat dez(dz - ez);
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

inline int inSphere_LLLLE_exact(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double ex, double ey, double ez)
{
   double avx[2];
   expansionObject::two_Diff(a2x, a1x, avx);
   double avxt[4];
   expansionObject::Two_One_Prod(avx, at, avxt);
   double ax[5];
   int ax_len = expansionObject::Gen_Sum(1, &a1x, 4, avxt, ax);
   double avy[2];
   expansionObject::two_Diff(a2y, a1y, avy);
   double avyt[4];
   expansionObject::Two_One_Prod(avy, at, avyt);
   double ay[5];
   int ay_len = expansionObject::Gen_Sum(1, &a1y, 4, avyt, ay);
   double avz[2];
   expansionObject::two_Diff(a2z, a1z, avz);
   double avzt[4];
   expansionObject::Two_One_Prod(avz, at, avzt);
   double az[5];
   int az_len = expansionObject::Gen_Sum(1, &a1z, 4, avzt, az);
   double bvx[2];
   expansionObject::two_Diff(b2x, b1x, bvx);
   double bvxt[4];
   expansionObject::Two_One_Prod(bvx, bt, bvxt);
   double bx[5];
   int bx_len = expansionObject::Gen_Sum(1, &b1x, 4, bvxt, bx);
   double bvy[2];
   expansionObject::two_Diff(b2y, b1y, bvy);
   double bvyt[4];
   expansionObject::Two_One_Prod(bvy, bt, bvyt);
   double by[5];
   int by_len = expansionObject::Gen_Sum(1, &b1y, 4, bvyt, by);
   double bvz[2];
   expansionObject::two_Diff(b2z, b1z, bvz);
   double bvzt[4];
   expansionObject::Two_One_Prod(bvz, bt, bvzt);
   double bz[5];
   int bz_len = expansionObject::Gen_Sum(1, &b1z, 4, bvzt, bz);
   double cvx[2];
   expansionObject::two_Diff(c2x, c1x, cvx);
   double cvxt[4];
   expansionObject::Two_One_Prod(cvx, ct, cvxt);
   double cx[5];
   int cx_len = expansionObject::Gen_Sum(1, &c1x, 4, cvxt, cx);
   double cvy[2];
   expansionObject::two_Diff(c2y, c1y, cvy);
   double cvyt[4];
   expansionObject::Two_One_Prod(cvy, ct, cvyt);
   double cy[5];
   int cy_len = expansionObject::Gen_Sum(1, &c1y, 4, cvyt, cy);
   double cvz[2];
   expansionObject::two_Diff(c2z, c1z, cvz);
   double cvzt[4];
   expansionObject::Two_One_Prod(cvz, ct, cvzt);
   double cz[5];
   int cz_len = expansionObject::Gen_Sum(1, &c1z, 4, cvzt, cz);
   double dvx[2];
   expansionObject::two_Diff(d2x, d1x, dvx);
   double dvxt[4];
   expansionObject::Two_One_Prod(dvx, dt, dvxt);
   double dx[5];
   int dx_len = expansionObject::Gen_Sum(1, &d1x, 4, dvxt, dx);
   double dvy[2];
   expansionObject::two_Diff(d2y, d1y, dvy);
   double dvyt[4];
   expansionObject::Two_One_Prod(dvy, dt, dvyt);
   double dy[5];
   int dy_len = expansionObject::Gen_Sum(1, &d1y, 4, dvyt, dy);
   double dvz[2];
   expansionObject::two_Diff(d2z, d1z, dvz);
   double dvzt[4];
   expansionObject::Two_One_Prod(dvz, dt, dvzt);
   double dz[5];
   int dz_len = expansionObject::Gen_Sum(1, &d1z, 4, dvzt, dz);
   double aex[6];
   int aex_len = expansionObject::Gen_Diff(ax_len, ax, 1, &ex, aex);
   double aey[6];
   int aey_len = expansionObject::Gen_Diff(ay_len, ay, 1, &ey, aey);
   double aez[6];
   int aez_len = expansionObject::Gen_Diff(az_len, az, 1, &ez, aez);
   double bex[6];
   int bex_len = expansionObject::Gen_Diff(bx_len, bx, 1, &ex, bex);
   double bey[6];
   int bey_len = expansionObject::Gen_Diff(by_len, by, 1, &ey, bey);
   double bez[6];
   int bez_len = expansionObject::Gen_Diff(bz_len, bz, 1, &ez, bez);
   double cex[6];
   int cex_len = expansionObject::Gen_Diff(cx_len, cx, 1, &ex, cex);
   double cey[6];
   int cey_len = expansionObject::Gen_Diff(cy_len, cy, 1, &ey, cey);
   double cez[6];
   int cez_len = expansionObject::Gen_Diff(cz_len, cz, 1, &ez, cez);
   double dex[6];
   int dex_len = expansionObject::Gen_Diff(dx_len, dx, 1, &ex, dex);
   double dey[6];
   int dey_len = expansionObject::Gen_Diff(dy_len, dy, 1, &ey, dey);
   double dez[6];
   int dez_len = expansionObject::Gen_Diff(dz_len, dz, 1, &ez, dez);
   double aexbey_p[16], *aexbey = aexbey_p;
   int aexbey_len = expansionObject::Gen_Product_With_PreAlloc(aex_len, aex, bey_len, bey, &aexbey, 16);
   double bexaey_p[16], *bexaey = bexaey_p;
   int bexaey_len = expansionObject::Gen_Product_With_PreAlloc(bex_len, bex, aey_len, aey, &bexaey, 16);
   double ab_p[16], *ab = ab_p;
   int ab_len = expansionObject::Gen_Diff_With_PreAlloc(aexbey_len, aexbey, bexaey_len, bexaey, &ab, 16);
   double bexcey_p[16], *bexcey = bexcey_p;
   int bexcey_len = expansionObject::Gen_Product_With_PreAlloc(bex_len, bex, cey_len, cey, &bexcey, 16);
   double cexbey_p[16], *cexbey = cexbey_p;
   int cexbey_len = expansionObject::Gen_Product_With_PreAlloc(cex_len, cex, bey_len, bey, &cexbey, 16);
   double bc_p[16], *bc = bc_p;
   int bc_len = expansionObject::Gen_Diff_With_PreAlloc(bexcey_len, bexcey, cexbey_len, cexbey, &bc, 16);
   double cexdey_p[16], *cexdey = cexdey_p;
   int cexdey_len = expansionObject::Gen_Product_With_PreAlloc(cex_len, cex, dey_len, dey, &cexdey, 16);
   double dexcey_p[16], *dexcey = dexcey_p;
   int dexcey_len = expansionObject::Gen_Product_With_PreAlloc(dex_len, dex, cey_len, cey, &dexcey, 16);
   double cd_p[16], *cd = cd_p;
   int cd_len = expansionObject::Gen_Diff_With_PreAlloc(cexdey_len, cexdey, dexcey_len, dexcey, &cd, 16);
   double dexaey_p[16], *dexaey = dexaey_p;
   int dexaey_len = expansionObject::Gen_Product_With_PreAlloc(dex_len, dex, aey_len, aey, &dexaey, 16);
   double aexdey_p[16], *aexdey = aexdey_p;
   int aexdey_len = expansionObject::Gen_Product_With_PreAlloc(aex_len, aex, dey_len, dey, &aexdey, 16);
   double da_p[16], *da = da_p;
   int da_len = expansionObject::Gen_Diff_With_PreAlloc(dexaey_len, dexaey, aexdey_len, aexdey, &da, 16);
   double aexcey_p[16], *aexcey = aexcey_p;
   int aexcey_len = expansionObject::Gen_Product_With_PreAlloc(aex_len, aex, cey_len, cey, &aexcey, 16);
   double cexaey_p[16], *cexaey = cexaey_p;
   int cexaey_len = expansionObject::Gen_Product_With_PreAlloc(cex_len, cex, aey_len, aey, &cexaey, 16);
   double ac_p[16], *ac = ac_p;
   int ac_len = expansionObject::Gen_Diff_With_PreAlloc(aexcey_len, aexcey, cexaey_len, cexaey, &ac, 16);
   double bexdey_p[16], *bexdey = bexdey_p;
   int bexdey_len = expansionObject::Gen_Product_With_PreAlloc(bex_len, bex, dey_len, dey, &bexdey, 16);
   double dexbey_p[16], *dexbey = dexbey_p;
   int dexbey_len = expansionObject::Gen_Product_With_PreAlloc(dex_len, dex, bey_len, bey, &dexbey, 16);
   double bd_p[16], *bd = bd_p;
   int bd_len = expansionObject::Gen_Diff_With_PreAlloc(bexdey_len, bexdey, dexbey_len, dexbey, &bd, 16);
   double abc1_p[16], *abc1 = abc1_p;
   int abc1_len = expansionObject::Gen_Product_With_PreAlloc(aez_len, aez, bc_len, bc, &abc1, 16);
   double abc2_p[16], *abc2 = abc2_p;
   int abc2_len = expansionObject::Gen_Product_With_PreAlloc(bez_len, bez, ac_len, ac, &abc2, 16);
   double abc3_p[16], *abc3 = abc3_p;
   int abc3_len = expansionObject::Gen_Product_With_PreAlloc(cez_len, cez, ab_len, ab, &abc3, 16);
   double abc4_p[16], *abc4 = abc4_p;
   int abc4_len = expansionObject::Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 16);
   double abc_p[16], *abc = abc_p;
   int abc_len = expansionObject::Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 16);
   double bcd1_p[16], *bcd1 = bcd1_p;
   int bcd1_len = expansionObject::Gen_Product_With_PreAlloc(bez_len, bez, cd_len, cd, &bcd1, 16);
   double bcd2_p[16], *bcd2 = bcd2_p;
   int bcd2_len = expansionObject::Gen_Product_With_PreAlloc(cez_len, cez, bd_len, bd, &bcd2, 16);
   double bcd3_p[16], *bcd3 = bcd3_p;
   int bcd3_len = expansionObject::Gen_Product_With_PreAlloc(dez_len, dez, bc_len, bc, &bcd3, 16);
   double bcd4_p[16], *bcd4 = bcd4_p;
   int bcd4_len = expansionObject::Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 16);
   double bcd_p[16], *bcd = bcd_p;
   int bcd_len = expansionObject::Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 16);
   double cda1_p[16], *cda1 = cda1_p;
   int cda1_len = expansionObject::Gen_Product_With_PreAlloc(cez_len, cez, da_len, da, &cda1, 16);
   double cda2_p[16], *cda2 = cda2_p;
   int cda2_len = expansionObject::Gen_Product_With_PreAlloc(dez_len, dez, ac_len, ac, &cda2, 16);
   double cda3_p[16], *cda3 = cda3_p;
   int cda3_len = expansionObject::Gen_Product_With_PreAlloc(aez_len, aez, cd_len, cd, &cda3, 16);
   double cda4_p[16], *cda4 = cda4_p;
   int cda4_len = expansionObject::Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 16);
   double cda_p[16], *cda = cda_p;
   int cda_len = expansionObject::Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 16);
   double dab1_p[16], *dab1 = dab1_p;
   int dab1_len = expansionObject::Gen_Product_With_PreAlloc(dez_len, dez, ab_len, ab, &dab1, 16);
   double dab2_p[16], *dab2 = dab2_p;
   int dab2_len = expansionObject::Gen_Product_With_PreAlloc(aez_len, aez, bd_len, bd, &dab2, 16);
   double dab3_p[16], *dab3 = dab3_p;
   int dab3_len = expansionObject::Gen_Product_With_PreAlloc(bez_len, bez, da_len, da, &dab3, 16);
   double dab4_p[16], *dab4 = dab4_p;
   int dab4_len = expansionObject::Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 16);
   double dab_p[16], *dab = dab_p;
   int dab_len = expansionObject::Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 16);
   double al1_p[16], *al1 = al1_p;
   int al1_len = expansionObject::Gen_Product_With_PreAlloc(aex_len, aex, aex_len, aex, &al1, 16);
   double al2_p[16], *al2 = al2_p;
   int al2_len = expansionObject::Gen_Product_With_PreAlloc(aey_len, aey, aey_len, aey, &al2, 16);
   double al3_p[16], *al3 = al3_p;
   int al3_len = expansionObject::Gen_Product_With_PreAlloc(aez_len, aez, aez_len, aez, &al3, 16);
   double al4_p[16], *al4 = al4_p;
   int al4_len = expansionObject::Gen_Sum_With_PreAlloc(al1_len, al1, al2_len, al2, &al4, 16);
   double alift_p[16], *alift = alift_p;
   int alift_len = expansionObject::Gen_Sum_With_PreAlloc(al4_len, al4, al3_len, al3, &alift, 16);
   double bl1_p[16], *bl1 = bl1_p;
   int bl1_len = expansionObject::Gen_Product_With_PreAlloc(bex_len, bex, bex_len, bex, &bl1, 16);
   double bl2_p[16], *bl2 = bl2_p;
   int bl2_len = expansionObject::Gen_Product_With_PreAlloc(bey_len, bey, bey_len, bey, &bl2, 16);
   double bl3_p[16], *bl3 = bl3_p;
   int bl3_len = expansionObject::Gen_Product_With_PreAlloc(bez_len, bez, bez_len, bez, &bl3, 16);
   double bl4_p[16], *bl4 = bl4_p;
   int bl4_len = expansionObject::Gen_Sum_With_PreAlloc(bl1_len, bl1, bl2_len, bl2, &bl4, 16);
   double blift_p[16], *blift = blift_p;
   int blift_len = expansionObject::Gen_Sum_With_PreAlloc(bl4_len, bl4, bl3_len, bl3, &blift, 16);
   double cl1_p[16], *cl1 = cl1_p;
   int cl1_len = expansionObject::Gen_Product_With_PreAlloc(cex_len, cex, cex_len, cex, &cl1, 16);
   double cl2_p[16], *cl2 = cl2_p;
   int cl2_len = expansionObject::Gen_Product_With_PreAlloc(cey_len, cey, cey_len, cey, &cl2, 16);
   double cl3_p[16], *cl3 = cl3_p;
   int cl3_len = expansionObject::Gen_Product_With_PreAlloc(cez_len, cez, cez_len, cez, &cl3, 16);
   double cl4_p[16], *cl4 = cl4_p;
   int cl4_len = expansionObject::Gen_Sum_With_PreAlloc(cl1_len, cl1, cl2_len, cl2, &cl4, 16);
   double clift_p[16], *clift = clift_p;
   int clift_len = expansionObject::Gen_Sum_With_PreAlloc(cl4_len, cl4, cl3_len, cl3, &clift, 16);
   double dl1_p[16], *dl1 = dl1_p;
   int dl1_len = expansionObject::Gen_Product_With_PreAlloc(dex_len, dex, dex_len, dex, &dl1, 16);
   double dl2_p[16], *dl2 = dl2_p;
   int dl2_len = expansionObject::Gen_Product_With_PreAlloc(dey_len, dey, dey_len, dey, &dl2, 16);
   double dl3_p[16], *dl3 = dl3_p;
   int dl3_len = expansionObject::Gen_Product_With_PreAlloc(dez_len, dez, dez_len, dez, &dl3, 16);
   double dl4_p[16], *dl4 = dl4_p;
   int dl4_len = expansionObject::Gen_Sum_With_PreAlloc(dl1_len, dl1, dl2_len, dl2, &dl4, 16);
   double dlift_p[16], *dlift = dlift_p;
   int dlift_len = expansionObject::Gen_Sum_With_PreAlloc(dl4_len, dl4, dl3_len, dl3, &dlift, 16);
   double ds1_p[16], *ds1 = ds1_p;
   int ds1_len = expansionObject::Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 16);
   double ds2_p[16], *ds2 = ds2_p;
   int ds2_len = expansionObject::Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 16);
   double dl_p[16], *dl = dl_p;
   int dl_len = expansionObject::Gen_Diff_With_PreAlloc(ds2_len, ds2, ds1_len, ds1, &dl, 16);
   double dr1_p[16], *dr1 = dr1_p;
   int dr1_len = expansionObject::Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 16);
   double dr2_p[16], *dr2 = dr2_p;
   int dr2_len = expansionObject::Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 16);
   double dr_p[16], *dr = dr_p;
   int dr_len = expansionObject::Gen_Diff_With_PreAlloc(dr2_len, dr2, dr1_len, dr1, &dr, 16);
   double det_p[16], *det = det_p;
   int det_len = expansionObject::Gen_Sum_With_PreAlloc(dl_len, dl, dr_len, dr, &det, 16);

   double return_value = det[det_len - 1];
   if (det_p != det) FreeDoubles(det);
   if (dr_p != dr) FreeDoubles(dr);
   if (dr2_p != dr2) FreeDoubles(dr2);
   if (dr1_p != dr1) FreeDoubles(dr1);
   if (dl_p != dl) FreeDoubles(dl);
   if (ds2_p != ds2) FreeDoubles(ds2);
   if (ds1_p != ds1) FreeDoubles(ds1);
   if (dlift_p != dlift) FreeDoubles(dlift);
   if (dl4_p != dl4) FreeDoubles(dl4);
   if (dl3_p != dl3) FreeDoubles(dl3);
   if (dl2_p != dl2) FreeDoubles(dl2);
   if (dl1_p != dl1) FreeDoubles(dl1);
   if (clift_p != clift) FreeDoubles(clift);
   if (cl4_p != cl4) FreeDoubles(cl4);
   if (cl3_p != cl3) FreeDoubles(cl3);
   if (cl2_p != cl2) FreeDoubles(cl2);
   if (cl1_p != cl1) FreeDoubles(cl1);
   if (blift_p != blift) FreeDoubles(blift);
   if (bl4_p != bl4) FreeDoubles(bl4);
   if (bl3_p != bl3) FreeDoubles(bl3);
   if (bl2_p != bl2) FreeDoubles(bl2);
   if (bl1_p != bl1) FreeDoubles(bl1);
   if (alift_p != alift) FreeDoubles(alift);
   if (al4_p != al4) FreeDoubles(al4);
   if (al3_p != al3) FreeDoubles(al3);
   if (al2_p != al2) FreeDoubles(al2);
   if (al1_p != al1) FreeDoubles(al1);
   if (dab_p != dab) FreeDoubles(dab);
   if (dab4_p != dab4) FreeDoubles(dab4);
   if (dab3_p != dab3) FreeDoubles(dab3);
   if (dab2_p != dab2) FreeDoubles(dab2);
   if (dab1_p != dab1) FreeDoubles(dab1);
   if (cda_p != cda) FreeDoubles(cda);
   if (cda4_p != cda4) FreeDoubles(cda4);
   if (cda3_p != cda3) FreeDoubles(cda3);
   if (cda2_p != cda2) FreeDoubles(cda2);
   if (cda1_p != cda1) FreeDoubles(cda1);
   if (bcd_p != bcd) FreeDoubles(bcd);
   if (bcd4_p != bcd4) FreeDoubles(bcd4);
   if (bcd3_p != bcd3) FreeDoubles(bcd3);
   if (bcd2_p != bcd2) FreeDoubles(bcd2);
   if (bcd1_p != bcd1) FreeDoubles(bcd1);
   if (abc_p != abc) FreeDoubles(abc);
   if (abc4_p != abc4) FreeDoubles(abc4);
   if (abc3_p != abc3) FreeDoubles(abc3);
   if (abc2_p != abc2) FreeDoubles(abc2);
   if (abc1_p != abc1) FreeDoubles(abc1);
   if (bd_p != bd) FreeDoubles(bd);
   if (dexbey_p != dexbey) FreeDoubles(dexbey);
   if (bexdey_p != bexdey) FreeDoubles(bexdey);
   if (ac_p != ac) FreeDoubles(ac);
   if (cexaey_p != cexaey) FreeDoubles(cexaey);
   if (aexcey_p != aexcey) FreeDoubles(aexcey);
   if (da_p != da) FreeDoubles(da);
   if (aexdey_p != aexdey) FreeDoubles(aexdey);
   if (dexaey_p != dexaey) FreeDoubles(dexaey);
   if (cd_p != cd) FreeDoubles(cd);
   if (dexcey_p != dexcey) FreeDoubles(dexcey);
   if (cexdey_p != cexdey) FreeDoubles(cexdey);
   if (bc_p != bc) FreeDoubles(bc);
   if (cexbey_p != cexbey) FreeDoubles(cexbey);
   if (bexcey_p != bexcey) FreeDoubles(bexcey);
   if (ab_p != ab) FreeDoubles(ab);
   if (bexaey_p != bexaey) FreeDoubles(bexaey);
   if (aexbey_p != aexbey) FreeDoubles(aexbey);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int inSphere_LLLLE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double ex, double ey, double ez)
{
   int ret;
   ret = inSphere_LLLLE_filtered(a1x, a1y, a1z, a2x, a2y, a2z, at, b1x, b1y, b1z, b2x, b2y, b2z, bt, c1x, c1y, c1z, c2x, c2y, c2z, ct, d1x, d1y, d1z, d2x, d2y, d2z, dt, ex, ey, ez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = inSphere_LLLLE_interval(a1x, a1y, a1z, a2x, a2y, a2z, at, b1x, b1y, b1z, b2x, b2y, b2z, bt, c1x, c1y, c1z, c2x, c2y, c2z, ct, d1x, d1y, d1z, d2x, d2y, d2z, dt, ex, ey, ez);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return inSphere_LLLLE_exact(a1x, a1y, a1z, a2x, a2y, a2z, at, b1x, b1y, b1z, b2x, b2y, b2z, bt, c1x, c1y, c1z, c2x, c2y, c2z, ct, d1x, d1y, d1z, d2x, d2y, d2z, dt, ex, ey, ez);
}

