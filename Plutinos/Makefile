# CFEPS Survey Simulator Makefile
#

### also works with ifort and f77
FC = gfortran
GIMEOBJ = OSSOS-plutinos-model

Plutinos: clean ../Commons/GiMeObj.f95 OSSOS-plutinos-model.f95 Makefile
	\rm -f Plutinos
	\rm -f ../Commons/$(GIMEOBJ).f95
	ln -s ../Plutinos/$(GIMEOBJ).f95 ../Commons/
	$(MAKE) -s -C ../Commons clean
	$(MAKE) -s -C ../Commons Driver GIMEOBJ=$(GIMEOBJ)
	\rm -f ../Commons/$(GIMEOBJ).f95 ../Commons/GiMeObj.f95
	cp ../Commons/Driver ./Plutinos

clean:
	\rm -f *.o core Plutinos OSSOS-model.out ModelUsed.dat
	$(MAKE) -s -C ../Commons clean

../Commons/GiMeObj.f95:
	\rm -f ../Commons/$(GIMEOBJ).f95 ../Commons/GiMeObj.f95
	ln -s ../Plutinos/$(GIMEOBJ).f95 ../Commons/
	ln -s $(GIMEOBJ).f95 ../Commons/GiMeObj.f95
