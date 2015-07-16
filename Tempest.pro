TEMPLATE = app
CONFIG += console debug
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += $$PWD/include \
    lib/tioga/src

SOURCES += src/global.cpp \
    src/matrix.cpp \
    src/input.cpp \
    src/ele.cpp \
    src/operators.cpp \
    src/geo.cpp \
    src/output.cpp \
    src/face.cpp \
    src/flux.cpp \
    src/solver.cpp \
    lib/tioga/src/ADT.C \
    lib/tioga/src/MeshBlock.C \
    lib/tioga/src/parallelComm.C \
    lib/tioga/src/tioga.C \
    lib/tioga/src/tiogaInterface.C \
    lib/tioga/src/math.c \
    lib/tioga/src/utils.c \
    lib/tioga/src/cellVolume.f90 \
    lib/tioga/src/computeCellVolume.f90 \
    lib/tioga/src/kaiser.f \
    lib/tioga/src/median.F90 \
    src/solver_bounds.cpp \
    src/tempest.cpp
		   
HEADERS += include/global.hpp \
    include/matrix.hpp \
    include/input.hpp \
    include/geo.hpp \
    include/output.hpp \
    include/face.hpp \
    include/flux.hpp \
    include/solver.hpp \
    include/error.hpp \
    lib/tioga/src/ADT.h \
    lib/tioga/src/codetypes.h \
    lib/tioga/src/globals.h \
    lib/tioga/src/MeshBlock.h \
    lib/tioga/src/parallelComm.h \
    lib/tioga/src/tioga.h \
    lib/tioga/src/tiogaInterface.h \
    lib/tioga/src/utils.h \
    include/tempest.hpp

DISTFILES += \
    README.md \

OTHER_FILES += makefile \
    lib/tioga/src/makefile
