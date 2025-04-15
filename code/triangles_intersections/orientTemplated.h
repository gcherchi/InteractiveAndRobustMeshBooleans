// orientTemplated.h
#ifndef ORIENT_TEMPLATED_H
#define ORIENT_TEMPLATED_H

#include <cinolib/predicates.h>
#include <CGAL/Gmpq.h>
#include <implicit_point.h>
#include <cinolib/rationals.h>
#include <cinolib/geometry/vec_mat.h>
#include <typeinfo>
#include <numerics.h>
template <typename T>
int orient3dT( const T * pa, const T * pb, const T * pc, const T * pd);

template <typename T>
int orient2dT( const T * pa, const T * pb, const T * pc);

#include "orientTemplated.cpp"

#endif // ORIENT_TEMPLATED_H
