# Copyright (C) 2001-2016 Quantum ESPRESSO group
# Copyright (C) 2017 Filippo Spiga
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include make.inc

default :
	@echo 'to install GPU-accelerated Quantum ESPRESSO, type at the shell prompt:'
	@echo '  ./configure [--prefix=]'
	@echo '  make [-j] target'
	@echo ' '
	@echo 'where target identifies one or multiple CORE PACKAGES:'
	@echo '  pw           basic code for scf, structure optimization, MD'
#	@echo '  ph           phonon code, Gamma-only and third-order derivatives'
#	@echo '  pwcond       ballistic conductance'
#	@echo '  neb          code for Nudged Elastic Band method'
#	@echo '  pp           postprocessing programs'
#	@echo '  pwall        same as "make pw ph pp pwcond neb"'
#	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
#	@echo '  tddfpt       time dependent dft code'
#	@echo '  gwl          GW with Lanczos chains'
#	@echo '  ld1          utilities for pseudopotential generation'
#	@echo '  upf          utilities for pseudopotential conversion'
#	@echo '  xspectra     X-ray core-hole spectroscopy calculations'
#	@echo '  couple       Library interface for coupling to external codes'
#	@echo '  epw          Electron-Phonon Coupling with wannier functions'
#	@echo '  test-suite   run semi-automated test-suite for regression testing'
	@echo '  all          same as "make pw"'
	@echo ' '
	@echo 'where target is one of the following suite operation:'
#	@echo '  doc          build documentation'
	@echo '  links        create links to all executables in bin/'
	@echo '  tar          create a tarball of the source tree'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    remove files produced by "configure" as well'
	@echo '  distclean    revert distribution to the original status'

###########################################################
# Main targets
###########################################################

# The syntax "( cd PW ; $(MAKE) TLDEPS= all || exit 1)" below
# guarantees that error code 1 is returned in case of error and make stops
# If "|| exit 1" is not present, the error code from make in subdirectories
# is not returned and make goes on even if compilation has failed

pw : bindir libfft libla mods libs libiotk 
	if test -d PW ; then \
	( cd PW ; $(MAKE) TLDEPS= all || exit 1) ; fi

pw-lib : bindir libfft libla mods libs libiotk
	if test -d PW ; then \
	( cd PW ; $(MAKE) TLDEPS= pw-lib || exit 1) ; fi

#cp : bindir libfft libla mods liblapack libs libiotk
#	if test -d CPV ; then \
#	( cd CPV ; $(MAKE) TLDEPS= all || exit 1) ; fi

#ph : bindir libfft libla mods libs pw lrmods
#	if test -d PHonon; then \
#	(cd PHonon; $(MAKE) all || exit 1) ; fi

#neb : bindir libfft libla mods libs pw
#	if test -d NEB; then \
#  (cd NEB; $(MAKE) all || exit 1) ; fi

#tddfpt : bindir libfft libla mods libs pw
#	if test -d TDDFPT; then \
#	(cd TDDFPT; $(MAKE) all || exit 1) ; fi

#pp : bindir libfft libla mods libs pw
#	if test -d PP ; then \
#	( cd PP ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#pwcond : bindir libfft libla mods libs pw pp
#	if test -d PWCOND ; then \
#	( cd PWCOND ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#acfdt : bindir libfft libla mods libs pw ph
#	if test -d ACFDT ; then \
#	( cd ACFDT ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#gwl : ph
#	if test -d GWW ; then \
#	( cd GWW ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#ld1 : bindir liblapack libfft libla mods libs
#	if test -d atomic ; then \
#	( cd atomic ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#upf : libfft libla mods libs liblapack
#	if test -d upftools ; then \
#	( cd upftools ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#pw_export : libiotk bindir libfft mods libs pw
#	if test -d PP ; then \
#	( cd PP ; $(MAKE) TLDEPS= pw_export.x || exit 1 ) ; fi

#xspectra : bindir libfft mods libs pw
#	if test -d XSpectra ; then \
#	( cd XSpectra ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#couple : pw cp
#	if test -d COUPLE ; then \
#	( cd COUPLE ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

#all   : pwall cp ld1 upf tddfpt gwl xspectra
all	:	pw

###########################################################
# Auxiliary targets used by main targets:
# compile modules, libraries, directory for binaries, etc
###########################################################

libla : touch-dummy libeigensolver_gpu
	( cd LAXlib ; $(MAKE) TLDEPS= all || exit 1 )

libfft : touch-dummy
	( cd FFTXlib ; $(MAKE) TLDEPS= all || exit 1 )

mods : libiotk libla libfft
	( cd Modules ; $(MAKE) TLDEPS= all || exit 1 )

libs : mods
	( cd clib ; $(MAKE) TLDEPS= all || exit 1 )

lrmods : libs libla libfft 
	( cd LR_Modules ; $(MAKE) TLDEPS= all || exit 1 )

bindir :
	test -d bin || mkdir bin

#############################################################
# Targets for external libraries
############################################################

#libblas : liblapack

#liblapack: touch-dummy
#	cd install ; $(MAKE) -f extlibs_makefile $@

libeigensolver_gpu: touch-dummy
	cd install ; $(MAKE) -f extlibs_makefile $@

libiotk: libeigensolver_gpu touch-dummy
	cd install ; $(MAKE) -f extlibs_makefile $@

# In case of trouble with iotk and compilers, add
# FFLAGS="$(FFLAGS_NOOPT)" after $(MFLAGS)

touch-dummy :
	$(dummy-variable)

#########################################################
# "make links" produces links to all executables in bin/
#########################################################

links : bindir
	( cd bin/ ; \
	rm -f *.x ; \
	for exe in ../*/*/*.x ../*/bin/* ; do \
	    if test ! -L $$exe ; then ln -fs $$exe . ; fi \
	done ; 

#########################################################
# 'make install' works based on --with-prefix
# - If the final directory does not exists it creates it
#########################################################

install : touch-dummy
	@if test -d bin ; then mkdir -p $(PREFIX)/bin ; \
	for x in `find * ! -path "test-suite/*" -name *.x -type f` ; do \
		cp $$x $(PREFIX)/bin/ ; done ; \
	fi
	@echo 'GPU-accelerated Quantum ESPRESSO binaries installed in $(PREFIX)/bin'

#########################################################
# Run test-suite for numerical regression testing
# NB: it is assumed that reference outputs have been 
#     already computed once (usualy during release)
#########################################################

#test-suite: pw touch-dummy
#	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

#########################################################
# Other targets: clean up
#########################################################

# remove object files and executables
clean : 
	touch make.inc 
	for dir in \
		CPV LAXlib FFTXlib Modules PP PW EPW \
		NEB ACFDT COUPLE GWW XSpectra PWCOND \
		atomic clib LR_Modules pwtools upftools \
		dev-tools extlibs Environ TDDFPT PHonon GWW \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
	    fi \
	done
	- @(cd install ; $(MAKE) -f extlibs_makefile clean)
	- /bin/rm -rf bin/*.x tmp

# remove files produced by "configure" as well
veryclean : clean
	- @(cd install ; $(MAKE) -f extlibs_makefile veryclean)
	- rm -rf install/patch-plumed
	- cd install ; rm -f config.log configure.msg config.status \
		CPV/version.h ChangeLog* intel.pcl */intel.pcl
	- cd install ; rm -fr autom4te.cache
	- cd install; ./clean.sh ; cd -
	- rm -f espresso.tar.gz
	- rm -rf make.inc

distclean : veryclean

tar :
	@if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e /.svn/ -e'/\.' -e'\.o$$' -e'\.mod$$'\
		-e /.git/ -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'_tmp\.f90$$' -e'\.x$$' \
		-e'~$$' -e'\./GUI' -e '\./tempdir' | xargs tar rvf espresso.tar
	gzip espresso.tar

depend: libiotk 
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)
