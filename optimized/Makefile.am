AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov

bin_PROGRAMS = hd.opt helix.opt makeHexSquare.opt min.opt subdivide.opt scale.opt icosahedron.opt makeBud.opt makeBump.opt makeHexCylinder.opt  writeBiasPotential.opt writeCylRTraj.opt makePore.opt frc_sandbox.opt rd_jcp.opt
noinst_LIBRARIES = libinterp.a

AM_CPPFLAGS = -include ../config.h -include ../src/local_config.h $(AM_INCLUDE) 
AM_CFLAGS = -O2  
AM_CXXFLAGS = -O2

SRC_DIR = ../src
.cu.o:
	$(NVCC) -o $*.o -c $(SRC_DIR)/$*.cu

vpath %.C $(SRC_DIR)
vpath %.c $(SRC_DIR)
vpath %.cpp $(SRC_DIR)
vpath %.h  $(SRC_DIR)
vpath %.cu $(SRC_DIR)

complex_SRCS = NBAR.C dimer.C crowder.C syt.C ifitm3.C dynamin.C CNT.C dvl.C
util_SRCS    = mutil.C util.C
amoeba_SRCS  = amoeba_code/amotry.C amoeba_code/amoeba.C amoeba_code/nrutil.C 

rd_SRCS = fpr_subroutines/Faddeeva.cpp fpr_subroutines/passocF.cpp fpr_subroutines/pirr_pfree_ratio_psF.cpp fpr_subroutines/2D.cpp fpr_subroutines/DDmatrixcreate.cpp fpr_subroutines/DDpsur.cpp fpr_subroutines/TBLnorm.cpp fpr_subroutines/TBLpirr.cpp fpr_subroutines/TBLsur.cpp fpr_subroutines/fpir.cpp fpr_subroutines/fsur.cpp fpr_subroutines/integrator.cpp fpr_subroutines/fnorm.cpp rd.C fpr_subroutines/sizelookup.cpp fpr_subroutines/DDpirr_pfree_ratio_ps.cpp fpr_subroutines/pirr.cpp fpr_subroutines/pnorm.cpp

interp_SRCS = interpLib.C grad.C pboxing.C irr_kernel.C uv_map.C gauss.C modes.C write_structures.C min_distance.cpp GJK.cpp gjk_gudrun.cpp GJK_min2.cpp functions.cpp NearestSimplexR2.cpp NearestSimplexR.cpp  subdivide3.C fast_mm.C random_global.C height.C  gsl_random_globals.C sparse.C global_boxing.C lipid_density.C tachyon.C  lipid_composition.C vpass.C get_regions.C input.C l-bfgs.C srd.C compilation_time.C meshCollisionLib.C sans.C spline.C globals.C createAllAtom.C init.C fitRho.C 3d_interp.C M_matrix.C volume.C giftwrap.C aux_surfaces.C maxc.C face_mask.C aa_build_util.C protein_build_helper.C clathrin.C pdbFetch.C io_mol_read.C special_pore_plot_flips.C

sim_SRCS = simulation.C gather.C pcomplex.C p_p.C minimize.C dynamics.C npt.C parallel.C restart.C

io_SRCS = dcd.C pdb.C alignSet.C

global_SRCS = $(interp_SRCS) $(sim_SRCS) $(io_SRCS) $(amoeba_SRCS) $(rd_SRCS)

if CUDA
cuda_SRCS = cuda.cu
else
cuda_SRCS =
endif
hull_SRCS = dovoronoi.c clusterv.C

libinterp_a_SOURCES 	   = $(interp_SRCS) $(util_SRCS)
rd_jcp_opt_SOURCES        = hd.C $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS) $(hull_SRCS) 
hd_opt_SOURCES        = hd.C $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS) $(hull_SRCS) 
helix_opt_SOURCES        = helix.C $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS) $(hull_SRCS) 
frc_sandbox_opt_SOURCES        = frc_sandbox.C $(global_SRCS) $(complex_SRCS) $(hull_SRCS) $(util_SRCS) $(cuda_SRCS) 
min_opt_SOURCES       = min.C $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS) $(hull_SRCS)
subdivide_opt_SOURCES = subdivide.C $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS) $(hull_SRCS)
scale_opt_SOURCES	  = scale.C  $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS) $(hull_SRCS)
icosahedron_opt_SOURCES = icosahedron.C $(util_SRCS)
makeBud_opt_SOURCES   = makeBud.C $(hull_SRCS) $(util_SRCS)
makePore_opt_SOURCES   = makePore.C $(hull_SRCS) $(util_SRCS) getCylinderMesh.C
makeBump_opt_SOURCES  = makeBump.C $(hull_SRCS) $(util_SRCS)
makeHexCylinder_opt_SOURCES = makeHexCylinder.C $(hull_SRCS) $(util_SRCS)
makeHexSquare_opt_SOURCES = makeHexSquare.C $(hull_SRCS) $(util_SRCS)
writeBiasPotential_opt_SOURCES = writeBiasPotential.C $(util_SRCS)
writeCylRTraj_opt_SOURCES      = writeCylRTraj.C $(util_SRCS) 

rd_jcp_opt_LDFLAGS        = $(AM_LDFLAGS)
rd_jcp_opt_LDADD	  = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
hd_opt_LDFLAGS		  = $(AM_LDFLAGS) 
helix_opt_LDFLAGS		  = $(AM_LDFLAGS) 
hd_opt_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
helix_opt_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
frc_sandbox_opt_LDADD     = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
min_opt_LDADD             = $(AM_LDADD) $(AM_QHULL) -lstdc++
subdivide_opt_LDADD       = $(AM_LDADD) $(AM_QHULL) -lstdc++
scale_opt_LDADD           = $(AM_LDADD) $(AM_QHULL) -lstdc++
icosahedron_opt_LDADD     = $(AM_LDADD) -lstdc++
makeBud_opt_LDADD         = $(AM_LDADD) $(AM_QHULL) -lstdc++
makePore_opt_LDADD        = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeBump_opt_LDADD        = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeHexCylinder_opt_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeHexSquare_opt_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
writeBiasPotential_opt_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
writeCylRTraj_opt_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++

compilation_time.o: .FORCE
.PHONY: .FORCE


