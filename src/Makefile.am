lib_LTLIBRARIES = libmeep.la
include_HEADERS = meep.hpp
pkginclude_HEADERS = meep/mympi.hpp meep/vec.hpp meep/meep-config.h meepgeom.hpp material_data.hpp adjust_verbosity.hpp

AM_CPPFLAGS = -I$(top_srcdir)/src

BUILT_SOURCES = sphere-quad.h  meep/meep-config.h

HDRS = meep.hpp meep_internals.hpp meep/mympi.hpp meep/vec.hpp	\
bicgstab.hpp meepgeom.hpp material_data.hpp adjust_verbosity.hpp

libmeep_la_SOURCES = array_slice.cpp anisotropic_averaging.cpp 		\
bands.cpp boundaries.cpp bicgstab.cpp casimir.cpp 	\
cw_fields.cpp dft.cpp dft_ldos.cpp energy_and_flux.cpp 	\
fields.cpp fields_dump.cpp fix_boundary_sources.cpp loop_in_chunks.cpp h5fields.cpp h5file.cpp 	\
initialize.cpp integrate.cpp integrate2.cpp material_data.cpp monitor.cpp mympi.cpp 	\
multilevel-atom.cpp near2far.cpp output_directory.cpp random.cpp 	\
sources.cpp step.cpp step_db.cpp stress.cpp structure.cpp structure_dump.cpp		\
susceptibility.cpp time.cpp update_eh.cpp mpb.cpp update_pols.cpp 	\
vec.cpp step_generic.cpp step_generic_stride1.cpp meepgeom.cpp GDSIIgeom.cpp newton_raphson.cpp $(HDRS) $(BUILT_SOURCES)

SUBDIRS = support
libmeep_la_LIBADD = support/libsupport.la

libmeep_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@

PRELUDE = "/* This file was automatically generated --- DO NOT EDIT */"

noinst_PROGRAMS = sphere_quad
sphere_quad_SOURCES = sphere-quad.cpp
SPHERE_QUAD = $(top_builddir)/src/sphere_quad
sphere-quad.h:
	make sphere_quad
	(echo $(PRELUDE); echo; $(SPHERE_QUAD)) > $@

step_generic_stride1.cpp: step_generic.cpp
		(echo $(PRELUDE); echo; \
	 perl -0777 -pe 's/LOOP_OVER/S1LOOP_OVER/g; \
	                 s/\bstep_curl\b/step_curl_stride1/g; \
	                 s/\bstep_update_EDHB\b/step_update_EDHB_stride1/g; \
	                 s/\bstep_update_EDHB_NL\b/step_update_EDHB_NL_stride1/g; \
	                 s/\bstep_beta\b/step_beta_stride1/g; \
	                 s/\bstep_bfast\b/step_bfast_stride1/g' \
	     $(top_srcdir)/src/step_generic.cpp) > $@

MAINTAINERCLEANFILES = $(BUILT_SOURCES)
