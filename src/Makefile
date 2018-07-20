# Makefile for Esh3D 

FFLAGS =
FPPFLAGS =
LDLIBS =

-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/conf/rules
-include ${PETSC_DIR}/lib/petsc/conf/variables
-include ${PETSC_DIR}/lib/petsc/conf/rules

ifeq ($(debug),true)
	FFLAGS += -O0 -g
else
	ifeq ($(CC),icc)
		FFLAGS += -ipo -diag-disable ipo -static_intel 
	else
		FFLAGS += -Wno-maybe-uninitialized -Wno-tabs
	endif
endif

OBJS = main.o m_esh3d.o m_global.o m_local.o m_elems.o m_utils.o

m_utils.o  : m_utils.F90
m_esh3d.o  : m_esh3d.F90 m_utils.o
m_elems.o  : m_elems.F90 m_utils.o
m_local.o  : m_local.F90 m_elems.o
m_global.o : m_global.F90 m_local.o
main.o     : main.F90 m_global.o m_esh3d.o

all: ${OBJS}
	-${FLINKER} ${OBJS} -o ../bin/esh3d ${PETSC_LIB} ${LDLIBS}
