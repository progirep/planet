QMAKE_CXXFLAGS += -std=c++14 -g
DEFINES += USE_GLPK NDEBUG

TEMPLATE = app \
    console
CONFIG += release
CONFIG -= app_bundle
CONFIG -= qt
QT -= gui \
    core

HEADERS += minisat2/Alg.h minisat2/Heap.h minisat2/Map.h         minisat2/Queue.h   minisat2/SolverTypes.h  minisat2/Vec.h minisat2/Alloc.h   minisat2/IntMap.h    minisat2/Options.h     minisat2/Rnd.h     minisat2/Sort.h         minisat2/XAlloc.h minisat2/Dimacs.h  minisat2/IntTypes.h  minisat2/ParseUtils.h  minisat2/Solver.h  minisat2/System.h supersetdatabase.hpp

HEADERS += verifierContext.hpp lpSolvers.hpp config.h

SOURCES += minisat2/Options.cc minisat2/Solver.cc minisat2/System.cc main.cpp verifierContext.cpp supersetdatabase.cpp

TARGET = planet

INCLUDEPATH += /usr/include/lpsolve/

LIBS +=  -static -lglpk -lgmp -lumfpack -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lz -lltdl -ldl  #-static

PKGCONFIG += 

