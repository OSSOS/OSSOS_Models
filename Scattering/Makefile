# CFEPS Survey Simulator Makefile
#

### also works with ifort and f77
FC = gfortran
GIMEOBJ = OSSOS-scattering-model

Scattering: clean ../Commons/GiMeObj.f95 OSSOS-scattering-model.f95 Makefile
	\rm -f Scattering
	\rm -f ../Commons/$(GIMEOBJ).f95
	ln -s ../Scattering/$(GIMEOBJ).f95 ../Commons/
	$(MAKE) -s -C ../Commons clean
	$(MAKE) -s -C ../Commons Driver GIMEOBJ=$(GIMEOBJ)
	\rm -f ../Commons/$(GIMEOBJ).f95 ../Commons/GiMeObj.f95
	cp ../Commons/Driver ./Scattering

clean:
	\rm -f *.o core Scattering OSSOS-model.out ModelUsed.dat
	$(MAKE) -s -C ../Commons clean

../Commons/GiMeObj.f95:
	\rm -f ../Commons/$(GIMEOBJ).f95 ../Commons/GiMeObj.f95
	ln -s ../Scattering/$(GIMEOBJ).f95 ../Commons/
	ln -s $(GIMEOBJ).f95 ../Commons/GiMeObj.f95
