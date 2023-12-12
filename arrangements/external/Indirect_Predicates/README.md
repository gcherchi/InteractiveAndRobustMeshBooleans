----------------------------
Indirect Predicates for Geometric Constructions
----------------------------

by Marco Attene

Consiglio Nazionale delle Ricerche                                        
Istituto di Matematica Applicata e Tecnologie Informatiche                
Sezione di Genova                                                         
IMATI-GE / CNR                                                            

This software implements fast and guaranteed floating point geometric predicates,
including novel indirect predicates as described in the following article:

M. Attene. Indirect Predicates for Geometric Constructions. In Elsevier Computer-Aided Design (2020, https://doi.org/10.1016/j.cad.2020.102856).

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
 - Mac OSX with CLANG.

---------------------
Usage
---------------------

The repository provides a header-only C++ library.
To use in your code:
1) Add the "Indirect_Predicates-master/include" path to the list of paths where your compiler searches header files
2) Include "implicit_point.h" in your code 
3) ALWAYS tell your compiler to use the following directives:
   MSVC: /fp:strict /Oi /STACK:8421376 /D _CRT_SECURE_NO_WARNINGS
   GCC/G++/CLANG: -frounding-math -O2 -Wl,-z,stacksize=8421376
4) Tell your compiler whether your CPU supports SSE2/AVX2 instructions
   MSVC: /arch:SSE2 or /arch:AVX2
   GCC/G++/CLANG: -msse2 or -mavx2
   
As an example, check the CMakeLists.txt provided to compile the test.cpp code.

---------------------
Copyright and license
---------------------

Indirect_Predicates
Authors: Marco Attene                                                    

Copyright(C) 2019: IMATI-GE / CNR                                        

IMATI-GE / CNR is Consiglio Nazionale delle Ricerche                     
Istituto di Matematica Applicata e Tecnologie Informatiche               
Genova (Italy)                                                           

Indirect_Predicates is free software; you can redistribute it and/or modify     
it under the terms of the GNU Lesser General Public License as published 
by the Free Software Foundation; either version 3 of the License, or (at 
your option) any later version.                                          

Indirect_Predicates is distributed in the hope that it will be useful, but      
WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser 
General Public License for more details.                                 

You should have received a copy of the GNU Lesser General Public License 
along with the Indirect_Predicates. If not, see http://www.gnu.org/licenses/.
