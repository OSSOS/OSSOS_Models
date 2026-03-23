U	= datadec ioutils elemutils rot
OBJU	= $(addsuffix .o,${U})
SRCU	= $(addsuffix .f95,${U})
FPPU	= $(addsuffix .fpp,${U})
MODU	= $(addsuffix .mod,${U})

M   = modelutils GiMeObj
OBJM	= $(addsuffix .o,${M})
SRCM	= $(addsuffix .f95,${M})
FPPM	= $(addsuffix .fpp,${M})
MODM	= $(addsuffix .mod,${M})


GOBJ	= GiMeObjF95

FC	= gfortran
FPP	= gfortran -E -x f95-cpp-input -fPIC -DPYTHON

FFLAGS	= -O3 -fPIC -x f95-cpp-input
FFLAGP	= -O3
CFLAGS  = 
CPPFLAGS = 

TOPTARGETS := clean test
GIMEOBJ ?= OSSOS-classical-model

all:	Driver

%.fpp :	%.f95
	$(FPP) $< -o $@

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

.PHONY: clean

clean: rmtmp
	rm -f _$(GOBJ)*.so Driver GiMeObj.f95

rmtmp:
	rm -f *~ *.o *.mod *.fpp f90wrap*f90
	rm -f $(GOBJ).py
	rm -f *.log LOG .f2py_f2cmap

.PHONY: link

GiMeObj.f95: link

link:
	\rm -f GiMeObj.f95
	\rm -f GiMeObj.o
	ln -s $(GIMEOBJ).f95 GiMeObj.f95

Driver: link Driver.f95 $(OBJU) ${SRCM} Makefile
	$(FC) $(FFLAGP) -o Driver $(OBJU) ${SRCM} Driver.f95

Modules: _$(GOBJ).so Makefile
	echo "Modules have been built"

_$(GOBJ).so: link $(FPPU) $(OBJU) $(GIMEOBJ).f95 Makefile
	\rm -f _$(GOBJ).so
	f90wrap -k kind_map -m $(GOBJ) $(FPPU) 1>f90wrap_GiMeObj.log 2>&1
	f2py-f90wrap --fcompiler=$(FC) -c -m _$(GOBJ) $(OBJU) f90wrap_*.f90 1>f2py_GiMeObj.log 2>&1
	\rm f90wrap_*.f90
	\rm *.fpp *.mod *.o
	sed -e 's/import _$(GOBJ)/from . import _$(GOBJ)/' -i.bck $(GOBJ).py
	rm $(GOBJ).py.bck
	rm .f2py_f2cmap
	mv _$(GOBJ)*.so _$(GOBJ).so

