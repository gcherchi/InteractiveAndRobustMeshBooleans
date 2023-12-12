DEL indirect_predicates.h
DEL indirect_predicates.hpp
FOR %%f in (..\JPCK\predicates\direct\*.txt) DO ..\JPCK\x64\Release\converter.exe %%f -p
FOR %%f in (..\JPCK\predicates\indirect\*.txt) DO ..\JPCK\x64\Release\converter.exe %%f -p
ECHO.>> indirect_predicates.h
ECHO #include "indirect_predicates.hpp" >> indirect_predicates.h
COPY indirect_predicates.h ..\include
COPY indirect_predicates.hpp ..\include
