AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov

bin_PROGRAMS = hd.dbg makeHexSquare.dbg min.dbg subdivide.dbg scale.dbg icosahedron.dbg makeBud.dbg makeBump.dbg makeHexCylinder.dbg writeBiasPotential.dbg writeCylRTraj.dbg makePore.dbg  frc_sandbox.dbg rd_jcp.dbg shiftMesh.dbg qii_simple.dbg qii_diamond.dbg qii_gyroid.dbg gyroid_ER.dbg join.dbg helix.dbg qii_x.dbg
noinst_LIBRARIES = libinterp.a

AM_CPPFLAGS = -include ../config.h -include ../src/local_config.h $(AM_INCLUDE) 
AM_CFLAGS = -g -pg
AM_CXXFLAGS = -g -pg

SRC_DIR = ../src

vpath %.C $(SRC_DIR)
vpath %.c $(SRC_DIR)
vpath %.cpp $(SRC_DIR)
vpath %.h  $(SRC_DIR)
vpath %.cu $(SRC_DIR)

.cu.o:
	$(NVCC) -o $*.o -c $(SRC_DIR)/$*.cu

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
hd_dbg_SOURCES        = hd.C $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
helix_dbg_SOURCES        = helix.C $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
rd_jcp_dbg_SOURCES        = rd_jcp.C $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
frc_sandbox_dbg_SOURCES        = frc_sandbox.C $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
shiftMesh_dbg_SOURCES        = shiftMesh.C $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) $(cuda_SRCS)
min_dbg_SOURCES       = min.C $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS) 
subdivide_dbg_SOURCES = subdivide.C $(global_SRCS)  $(hull_SRCS) $(complex_SRCS) $(util_SRCS)
scale_dbg_SOURCES	  = scale.C  $(global_SRCS) $(hull_SRCS) $(complex_SRCS) $(util_SRCS)
icosahedron_dbg_SOURCES = icosahedron.C $(util_SRCS)
join_dbg_SOURCES	= join.C $(global_SRCS) $(complex_SRCS) $(util_SRCS) $(hull_SRCS) getCylinderMesh.C
makeBud_dbg_SOURCES   = makeBud.C $(hull_SRCS) $(util_SRCS)
makePore_dbg_SOURCES   = makePore.C $(hull_SRCS) $(util_SRCS) getCylinderMesh.C
makeBump_dbg_SOURCES  = makeBump.C $(hull_SRCS) $(util_SRCS)
makeHexCylinder_dbg_SOURCES = makeHexCylinder.C $(hull_SRCS) $(util_SRCS)
makeHexSquare_dbg_SOURCES = makeHexSquare.C $(hull_SRCS) $(util_SRCS)
writeBiasPotential_dbg_SOURCES = writeBiasPotential.C $(util_SRCS)
writeCylRTraj_dbg_SOURCES      = writeCylRTraj.C $(util_SRCS) 
qii_diamond_dbg_SOURCES		= qii_diamond.C $(util_SRCS)
qii_gyroid_dbg_SOURCES		= qii_gyroid.C $(util_SRCS)
qii_simple_dbg_SOURCES		= qii_simple.C $(util_SRCS)
qii_x_dbg_SOURCES		= qii_x.C $(util_SRCS)
gyroid_ER_dbg_SOURCES		= gyroid_ER.C $(util_SRCS)

hd_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
helix_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
rd_jcp_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
frc_sandbox_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
shiftMesh_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) $(AM_CUDA) -lstdc++
join_dbg_LDADD              = $(AM_LDADD) $(AM_QHULL) -lstdc++
min_dbg_LDADD             = $(AM_LDADD) $(AM_QHULL) -lstdc++
subdivide_dbg_LDADD       = $(AM_LDADD) $(AM_QHULL) -lstdc++
scale_dbg_LDADD           = $(AM_LDADD) $(AM_QHULL) -lstdc++
icosahedron_dbg_LDADD     = $(AM_LDADD) -lstdc++
makeBud_dbg_LDADD         = $(AM_LDADD) $(AM_QHULL) -lstdc++
makePore_dbg_LDADD         = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeBump_dbg_LDADD        = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeHexCylinder_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
makeHexSquare_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
writeBiasPotential_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
writeCylRTraj_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
qii_simple_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
qii_diamond_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
qii_x_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
qii_gyroid_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++
gyroid_ER_dbg_LDADD = $(AM_LDADD) $(AM_QHULL) -lstdc++

compilation_time.o: .FORCE
.PHONY: .FORCE

