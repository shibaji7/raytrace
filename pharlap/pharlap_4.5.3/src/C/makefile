ifeq (${OS}, Windows_NT)
  MAKEFLAG = -f makefile_cygwin
else
  MAKEFLAG = -f makefile_unix
endif
MEX_INSTALL_PATH = ../../mex

pharlap:
	$(MAKE) ${MAKEFLAG}

clean:
	rm -f *.o
	rm -f *.obj
	rm -f *.mex*
	rm -f .DS_Store
	rm -f *.exp
	rm -f *.map
	rm -f *.lib

distclean:
	rm -f *.o
	rm -f *.obj
	rm -f *~*
	rm -f *.mex*
	rm -f .DS_Store
	rm -f *.exp
	rm -f *.map
	rm -f *.lib

install: pharlap
	mkdir -p ${MEX_INSTALL_PATH}
	mv raytrace_2d.mex* ${MEX_INSTALL_PATH}
	mv raytrace_2d_sp.mex* ${MEX_INSTALL_PATH}
	mv raytrace_3d.mex* ${MEX_INSTALL_PATH}
	mv raytrace_3d_sp.mex* ${MEX_INSTALL_PATH}
	mv irreg_strength.mex* ${MEX_INSTALL_PATH}
	mv dop_spread_eq.mex* ${MEX_INSTALL_PATH}
	mv iri2007.mex* ${MEX_INSTALL_PATH}
	mv igrf2007.mex* ${MEX_INSTALL_PATH}
	mv iri2012.mex* ${MEX_INSTALL_PATH}
	mv igrf2011.mex* ${MEX_INSTALL_PATH}
	mv iri2016.mex* ${MEX_INSTALL_PATH}
	mv igrf2016.mex* ${MEX_INSTALL_PATH}
	mv iri2020.mex* ${MEX_INSTALL_PATH}
	mv igrf2020.mex* ${MEX_INSTALL_PATH}
	mv nrlmsise00.mex* ${MEX_INSTALL_PATH}
	mv ground_bs_loss.mex* ${MEX_INSTALL_PATH}	
	mv ground_fs_loss.mex* ${MEX_INSTALL_PATH}
	mv abso_bg.mex* ${MEX_INSTALL_PATH}

