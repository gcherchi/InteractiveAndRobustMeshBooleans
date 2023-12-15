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

int orient3d_LEEE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz);
int orient3d_LLEE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double qt, double rx, double ry, double rz, double sx, double sy, double sz);
int orient3d_LLLE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double qt, double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double rt, double sx, double sy, double sz);
int orient3d_LLLL(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double qt, double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double rt, double s1x, double s1y, double s1z, double s2x, double s2y, double s2z, double st);
int inSphere_LEEEE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz, double ex, double ey, double ez);
int inSphere_LLEEE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double cx, double cy, double cz, double dx, double dy, double dz, double ex, double ey, double ez);
int inSphere_LLLEE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double dx, double dy, double dz, double ex, double ey, double ez);
int inSphere_LLLLE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double ex, double ey, double ez);
int inSphere_LLLLL(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double e1x, double e1y, double e1z, double e2x, double e2y, double e2z, double et);

inline int lnc_orient3d_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const explicitPoint3D &ae = a.toExplicit3D(), &be = b.toExplicit3D(), &ce = c.toExplicit3D(), &de = d.toExplicit3D();
	return orient3d(ae.X(), ae.Y(), ae.Z(), be.X(), be.Y(), be.Z(), ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC();
	const explicitPoint3D& be = b.toExplicit3D(), & ce = c.toExplicit3D(), & de = d.toExplicit3D();
	return orient3d_LEEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(), 
		be.X(), be.Y(), be.Z(), ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC(), &bl = b.toLNC();
	const explicitPoint3D& ce = c.toExplicit3D(), & de = d.toExplicit3D();
	return orient3d_LLEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(), 
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(), 
		ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC();
	const explicitPoint3D& de = d.toExplicit3D();
	return orient3d_LLLE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IIII(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC(), & dl = d.toLNC();
	return orient3d_LLLL(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		dl.P().X(), dl.P().Y(), dl.P().Z(), dl.Q().X(), dl.Q().Y(), dl.Q().Z(), dl.T());
}

// This version assumes that points are either explicit3D or LNC
inline int lnc_orient3D(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return lnc_orient3d_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit3D()) return lnc_orient3d_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return lnc_orient3d_IEEE(b, c, a, d);
		if (!c.isExplicit3D()) return lnc_orient3d_IEEE(c, d, a, b);
		return lnc_orient3d_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return lnc_orient3d_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return lnc_orient3d_IIEE(a, c, d, b);
		if (a.isExplicit3D() && d.isExplicit3D()) return lnc_orient3d_IIEE(b, c, a, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return lnc_orient3d_IIEE(d, a, c, b);
		if (a.isExplicit3D() && c.isExplicit3D()) return lnc_orient3d_IIEE(d, b, a, c);
		return lnc_orient3d_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return lnc_orient3d_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return lnc_orient3d_IIIE(d, b, a, c);
		if (b.isExplicit3D()) return lnc_orient3d_IIIE(a, c, d, b);
		return lnc_orient3d_IIIE(b, d, c, a);
	}

	return lnc_orient3d_IIII(a, b, c, d);
}

#include "lnc_predicates.hpp"
