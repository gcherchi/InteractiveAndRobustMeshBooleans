----------------------------
NFG - Numbers For Geometry
----------------------------

by Marco Attene

Consiglio Nazionale delle Ricerche                                        
Istituto di Matematica Applicata e Tecnologie Informatiche                
Sezione di Genova                                                         
IMATI-GE / CNR                                                            

NFG is a standalone header-only C++ library providing useful number types for geometric computation.
Provides types and operations for arithmetic on:
- Intervals
- Floating point expansions
- Arbitrarily large natural numbers
- Arbitrarily precise floating point numbers
- Rational numbers

-------------------
Citation policy
--------------------
You are free to use this software according to the licensing terms specified at the end of this document.
If you use it for research purposes and produce publications, please cite the following paper 
that describes the underlying theory:

> M. Attene. Indirect Predicates for Geometric Constructions. In Elsevier Computer-Aided Design (2020).

-------------------
System Requirements
--------------------

The software has been tested on 64 bit PCs running:
 - Microsoft Windows OS with MSVC
 - Linux with standard gcc/g++ development environment
 - Mac OSX with both CLang and gcc

---------------------
Usage
---------------------

The repository provides a header-only C++ library.
To use in your code:
1) Add the "$(DIRECTORY_CONTAINING_THIS_README)/include" path to the list of 
   paths where your compiler searches header files
2) Include "numerics.h" in your code 
3) Tell your compiler to use the compilation directives reported at
   the beginning of include/numerics.h

As an example, check the CMakeLists.txt provided to compile the test.cpp code.

---------------------
Copyright and license
---------------------

NFG
Authors: Marco Attene                                                    

Copyright(C) 2019: IMATI-GE / CNR                                        

IMATI-GE / CNR is Consiglio Nazionale delle Ricerche                     
Istituto di Matematica Applicata e Tecnologie Informatiche               
Genova (Italy)                                                           

NFG is free software; you can redistribute it and/or modify     
it under the terms of the GNU Lesser General Public License as published 
by the Free Software Foundation; either version 3 of the License, or (at 
your option) any later version.                                          

NFG is distributed in the hope that it will be useful, but      
WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser 
General Public License for more details.                                 

You should have received a copy of the GNU Lesser General Public License 
along with the NFG. If not, see http://www.gnu.org/licenses/.
