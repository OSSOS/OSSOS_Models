# OSSOS model Makefile
#

### also works with ifort and f77
FC = gfortran
GIMEOBJ = OSSOS-classical-model

Classical: clean ../Commons/GiMeObj.f95 OSSOS-classical-model.f95 Makefile
	\rm -f Classical
	\rm -f ../Commons/$(GIMEOBJ).f95
	ln -s ../Classical/$(GIMEOBJ).f95 ../Commons/
	$(MAKE) -s -C ../Commons clean
	$(MAKE) -s -C ../Commons Driver GIMEOBJ=$(GIMEOBJ)
	\rm -f ../Commons/$(GIMEOBJ).f95 ../Commons/GiMeObj.f95
	cp ../Commons/Driver ./Classical

clean:
	\rm -f *.o core Classical OSSOS-model.out ModelUsed.dat
	$(MAKE) -s -C ../Commons clean

../Commons/GiMeObj.f95:
	\rm -f ../Commons/$(GIMEOBJ).f95 ../Commons/GiMeObj.f95
	ln -s ../Classical/$(GIMEOBJ).f95 ../Commons/
	ln -s $(GIMEOBJ).f95 ../Commons/GiMeObj.f95
