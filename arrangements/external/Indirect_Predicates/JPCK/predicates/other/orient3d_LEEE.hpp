#include "implicit_point.h"

#pragma intrinsic(fabs)

// Uncomment the following to activate overflow/underflow checks
#define CHECK_FOR_XYZERFLOWS

inline int orient3d_LEEE_filtered(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   const double pvx = p2x - p1x;
   const double pvxt = pvx * pt;
   const double px = p1x + pvxt;
   const double pvy = p2y - p1y;
   const double pvyt = pvy * pt;
   const double py = p1y + pvyt;
   const double pvz = p2z - p1z;
   const double pvzt = pvz * pt;
   const double pz = p1z + pvzt;
   const double qx_px = qx - px;
   const double qy_py = qy - py;
   const double rx_px = rx - px;
   const double ry_py = ry - py;
   const double rz_pz = rz - pz;
   const double qz_pz = qz - pz;
   const double sx_px = sx - px;
   const double sy_py = sy - py;
   const double sz_pz = sz - pz;
   const double tmp_a = qx_px * ry_py;
   const double tmp_b = qy_py * rx_px;
   const double m01 = tmp_a - tmp_b;
   const double tmq_a = qx_px * rz_pz;
   const double tmq_b = qz_pz * rx_px;
   const double m02 = tmq_a - tmq_b;
   const double tmr_a = qy_py * rz_pz;
   const double tmr_b = qz_pz * ry_py;
   const double m12 = tmr_a - tmr_b;
   const double mt1 = m01 * sz_pz;
   const double mt2 = m02 * sy_py;
   const double mt3 = m12 * sx_px;
   const double mtt = mt1 - mt2;
   const double m012 = mtt + mt3;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(p1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pt)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(rx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ry)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(rz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pvx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pvy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pvz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.718625242119744e-13;
   if (m012 > epsilon) return IP_Sign::POSITIVE;
   if (-m012 > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

inline int orient3d_LEEE_interval(interval_number p1x, interval_number p1y, interval_number p1z, interval_number p2x, interval_number p2y, interval_number p2z, interval_number pt, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number sx, interval_number sy, interval_number sz)
{
   setFPUModeToRoundUP();
   const interval_number pvx(p2x - p1x);
   const interval_number pvxt(pvx * pt);
   const interval_number px(p1x + pvxt);
   const interval_number pvy(p2y - p1y);
   const interval_number pvyt(pvy * pt);
   const interval_number py(p1y + pvyt);
   const interval_number pvz(p2z - p1z);
   const interval_number pvzt(pvz * pt);
   const interval_number pz(p1z + pvzt);
   const interval_number qx_px(qx - px);
   const interval_number qy_py(qy - py);
   const interval_number rx_px(rx - px);
   const interval_number ry_py(ry - py);
   const interval_number rz_pz(rz - pz);
   const interval_number qz_pz(qz - pz);
   const interval_number sx_px(sx - px);
   const interval_number sy_py(sy - py);
   const interval_number sz_pz(sz - pz);
   const interval_number tmp_a(qx_px * ry_py);
   const interval_number tmp_b(qy_py * rx_px);
   const interval_number m01(tmp_a - tmp_b);
   const interval_number tmq_a(qx_px * rz_pz);
   const interval_number tmq_b(qz_pz * rx_px);
   const interval_number m02(tmq_a - tmq_b);
   const interval_number tmr_a(qy_py * rz_pz);
   const interval_number tmr_b(qz_pz * ry_py);
   const interval_number m12(tmr_a - tmr_b);
   const interval_number mt1(m01 * sz_pz);
   const interval_number mt2(m02 * sy_py);
   const interval_number mt3(m12 * sx_px);
   const interval_number mtt(mt1 - mt2);
   const interval_number m012(mtt + mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return m012.sign();
}

inline int orient3d_LEEE_bigfloat(bigfloat p1x, bigfloat p1y, bigfloat p1z, bigfloat p2x, bigfloat p2y, bigfloat p2z, bigfloat pt, bigfloat qx, bigfloat qy, bigfloat qz, bigfloat rx, bigfloat ry, bigfloat rz, bigfloat sx, bigfloat sy, bigfloat sz)
{
   const bigfloat pvx(p2x - p1x);
   const bigfloat pvxt(pvx * pt);
   const bigfloat px(p1x + pvxt);
   const bigfloat pvy(p2y - p1y);
   const bigfloat pvyt(pvy * pt);
   const bigfloat py(p1y + pvyt);
   const bigfloat pvz(p2z - p1z);
   const bigfloat pvzt(pvz * pt);
   const bigfloat pz(p1z + pvzt);
   const bigfloat qx_px(qx - px);
   const bigfloat qy_py(qy - py);
   const bigfloat rx_px(rx - px);
   const bigfloat ry_py(ry - py);
   const bigfloat rz_pz(rz - pz);
   const bigfloat qz_pz(qz - pz);
   const bigfloat sx_px(sx - px);
   const bigfloat sy_py(sy - py);
   const bigfloat sz_pz(sz - pz);
   const bigfloat tmp_a(qx_px * ry_py);
   const bigfloat tmp_b(qy_py * rx_px);
   const bigfloat m01(tmp_a - tmp_b);
   const bigfloat tmq_a(qx_px * rz_pz);
   const bigfloat tmq_b(qz_pz * rx_px);
   const bigfloat m02(tmq_a - tmq_b);
   const bigfloat tmr_a(qy_py * rz_pz);
   const bigfloat tmr_b(qz_pz * ry_py);
   const bigfloat m12(tmr_a - tmr_b);
   const bigfloat mt1(m01 * sz_pz);
   const bigfloat mt2(m02 * sy_py);
   const bigfloat mt3(m12 * sx_px);
   const bigfloat mtt(mt1 - mt2);
   const bigfloat m012(mtt + mt3);
   return sgn(m012);
}

inline int orient3d_LEEE_exact(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   double pvx[2];
   expansionObject::two_Diff(p2x, p1x, pvx);
   double pvxt[4];
   expansionObject::Two_One_Prod(pvx, pt, pvxt);
   double px[5];
   int px_len = expansionObject::Gen_Sum(1, &p1x, 4, pvxt, px);
   double pvy[2];
   expansionObject::two_Diff(p2y, p1y, pvy);
   double pvyt[4];
   expansionObject::Two_One_Prod(pvy, pt, pvyt);
   double py[5];
   int py_len = expansionObject::Gen_Sum(1, &p1y, 4, pvyt, py);
   double pvz[2];
   expansionObject::two_Diff(p2z, p1z, pvz);
   double pvzt[4];
   expansionObject::Two_One_Prod(pvz, pt, pvzt);
   double pz[5];
   int pz_len = expansionObject::Gen_Sum(1, &p1z, 4, pvzt, pz);
   double qx_px[6];
   int qx_px_len = expansionObject::Gen_Diff(1, &qx, px_len, px, qx_px);
   double qy_py[6];
   int qy_py_len = expansionObject::Gen_Diff(1, &qy, py_len, py, qy_py);
   double rx_px[6];
   int rx_px_len = expansionObject::Gen_Diff(1, &rx, px_len, px, rx_px);
   double ry_py[6];
   int ry_py_len = expansionObject::Gen_Diff(1, &ry, py_len, py, ry_py);
   double rz_pz[6];
   int rz_pz_len = expansionObject::Gen_Diff(1, &rz, pz_len, pz, rz_pz);
   double qz_pz[6];
   int qz_pz_len = expansionObject::Gen_Diff(1, &qz, pz_len, pz, qz_pz);
   double sx_px[6];
   int sx_px_len = expansionObject::Gen_Diff(1, &sx, px_len, px, sx_px);
   double sy_py[6];
   int sy_py_len = expansionObject::Gen_Diff(1, &sy, py_len, py, sy_py);
   double sz_pz[6];
   int sz_pz_len = expansionObject::Gen_Diff(1, &sz, pz_len, pz, sz_pz);
   double tmp_a[72];
   int tmp_a_len = expansionObject::Gen_Product(qx_px_len, qx_px, ry_py_len, ry_py, tmp_a);
   double tmp_b[72];
   int tmp_b_len = expansionObject::Gen_Product(qy_py_len, qy_py, rx_px_len, rx_px, tmp_b);
   double m01_p[128], *m01 = m01_p;
   int m01_len = expansionObject::Gen_Diff_With_PreAlloc(tmp_a_len, tmp_a, tmp_b_len, tmp_b, &m01, 128);
   double tmq_a[72];
   int tmq_a_len = expansionObject::Gen_Product(qx_px_len, qx_px, rz_pz_len, rz_pz, tmq_a);
   double tmq_b[72];
   int tmq_b_len = expansionObject::Gen_Product(qz_pz_len, qz_pz, rx_px_len, rx_px, tmq_b);
   double m02_p[128], *m02 = m02_p;
   int m02_len = expansionObject::Gen_Diff_With_PreAlloc(tmq_a_len, tmq_a, tmq_b_len, tmq_b, &m02, 128);
   double tmr_a[72];
   int tmr_a_len = expansionObject::Gen_Product(qy_py_len, qy_py, rz_pz_len, rz_pz, tmr_a);
   double tmr_b[72];
   int tmr_b_len = expansionObject::Gen_Product(qz_pz_len, qz_pz, ry_py_len, ry_py, tmr_b);
   double m12_p[128], *m12 = m12_p;
   int m12_len = expansionObject::Gen_Diff_With_PreAlloc(tmr_a_len, tmr_a, tmr_b_len, tmr_b, &m12, 128);
   double mt1_p[128], *mt1 = mt1_p;
   int mt1_len = expansionObject::Gen_Product_With_PreAlloc(m01_len, m01, sz_pz_len, sz_pz, &mt1, 128);
   double mt2_p[128], *mt2 = mt2_p;
   int mt2_len = expansionObject::Gen_Product_With_PreAlloc(m02_len, m02, sy_py_len, sy_py, &mt2, 128);
   double mt3_p[128], *mt3 = mt3_p;
   int mt3_len = expansionObject::Gen_Product_With_PreAlloc(m12_len, m12, sx_px_len, sx_px, &mt3, 128);
   double mtt_p[128], *mtt = mtt_p;
   int mtt_len = expansionObject::Gen_Diff_With_PreAlloc(mt1_len, mt1, mt2_len, mt2, &mtt, 128);
   double m012_p[128], *m012 = m012_p;
   int m012_len = expansionObject::Gen_Sum_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 128);

   double return_value = m012[m012_len - 1];
   if (m012_p != m012) FreeDoubles(m012);
   if (mtt_p != mtt) FreeDoubles(mtt);
   if (mt3_p != mt3) FreeDoubles(mt3);
   if (mt2_p != mt2) FreeDoubles(mt2);
   if (mt1_p != mt1) FreeDoubles(mt1);
   if (m12_p != m12) FreeDoubles(m12);
   if (m02_p != m02) FreeDoubles(m02);
   if (m01_p != m01) FreeDoubles(m01);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

inline int orient3d_LEEE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   int ret;
   ret = orient3d_LEEE_filtered(p1x, p1y, p1z, p2x, p2y, p2z, pt, qx, qy, qz, rx, ry, rz, sx, sy, sz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient3d_LEEE_interval(p1x, p1y, p1z, p2x, p2y, p2z, pt, qx, qy, qz, rx, ry, rz, sx, sy, sz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_LEEE_exact(p1x, p1y, p1z, p2x, p2y, p2z, pt, qx, qy, qz, rx, ry, rz, sx, sy, sz);
}

