The code-generation tool can parse files that strictly adhere to the following format.
Each line must be one of the following:
1) An explicit parameter declaration
2) An assignment
3) A comment
4) An empty line
5) A lambda declaration
6) A sign declaration
7) A nonzero declaration

1) Explicit parameter declaration (any number of different variable names)
v1
v1 v2 ... vn

2) Assignment
v1 = v2 + v3
v1 = v2 - v3
v1 = v2 * v3
v1 = 2 * v3

Any variable used (e.g. v2 and v3) must be declared beforehand. An assignment corresponds to an implicit declaration (e.g. of v1).
A previously declared variable may not be reassigned.
If a variable name starts with "lambda" (e.g. lambda_x) it will be considered as the output of a lambda function.
For the lambda determinant use the special variable name "lambda_d".
If a "lambda" variable is found, the whole code is assumed to be a lambda code.

3) Comment
// this is a comment

5) Lambda declaration
lambda_name(v1,v1,...,vn:lx;degree;size;err_bound;val_bound;ly;degree;size;err_bound;val_bound; ... d;degree;size;err_bound;val_bound)

where:
lambda_name is the name of the lambda function
v1..vn are the explicit parameters used by lambda (must be declared beforehand)
lx, ly, ..., d are the implicit parameters produced by lambda

degree is the polynomial degree for the parameter
size is the expansion size for the parameter
err_bound is the error bound for the parameter
val_bound is the value bound for the parameter
The aforementioned four values are produced when compiling the lambda code.

6) Sign declaration (deprecated)
SIGN d1 d2 ... dn
where di are normally lambda outputs that determine the final sign

7) Nonzero declaration
NONZERO v1 v2 ... vn
means that the predicate returns a boolean value if any of the v's is nonzero


**********************************

Source equations for indirect predicates
--

Lambdas
2d_SSI = 2D line-line intersection
3d_LPI = 3D line-plane intersection
3d_TPI = 3D three planes intersection
2d3d_LPI = XY coords of the 3D line-plane intersection
2d3d_TPI = XY coords of the 3D three planes intersection

