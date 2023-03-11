#
# Paths to directories commands
#
include ../compiler.mk

DIR = .
INCDIR = ${DIR}/include
TSHDIR = ${DIR}/tsh

#
# Libraries
#
LIB = -mkl
ifdef LIB
LIB0 = ${LIB}
else
LIB = -L/home/softwares/intel/composer_xe_2013.1.117/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64
endif

#toolpath= /home/dhu13/Program/JADE_NAMD/src/dynamics/toolbox
#LT =  -L$(toolpath) -ltoolkit

#
# Compilers
#

F90C = ifort
F90FLAGS = -I${INCDIR}

CC = gcc
CFLAGS = -O2 -mcpu=i486 -I. -DLINUX -D__I8=0 -D__ENDIAN=0 

#
# Executable, source and objects
#
EXE  = main_sh.exe
MAIN = main_sh.f90

%.o: %.mod


#
# Dependencies
#

${EXE}:
	$(MAKE) -j2 -C ${TSHDIR} $(MAKECMDGOALS) 

	${F90C} -o $@ ${F90FLAGS} ${MAIN} ${TSHO} ${LIB}

TSHO = $(wildcard ./tsh/*o)

install: ${EXE} 
	mv ${EXE} ../../bin/jade.exe

clean:
	$(MAKE) -j2 -C ${TSHDIR} $(MAKECMDGOALS) 

#
# Suffix rules
#
.SUFFIXES: $(SUFFIXES) .f90
.f90.o:
	${F90C} -c ${F90FLAGS} $*.f90
.f.o:
	${FC} -c ${FFLAGS} $*.f
	
.c.o:
	${CC} -c ${CFLAGS} $*.c
