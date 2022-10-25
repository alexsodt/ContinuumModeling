#bilaye midplane: 
#makePore.opt 6 180 50 90 

set vexp = 1.25

# approximately matches boundary conditions.. sufficient density (6)
# Lx Lz Ly
makeHexSquare.opt 6 190 200 190

# mesh pt 72 is near the middle. I am going to change this so that it doesn't matter..
set pt1 = 72
set pt2 = 72
set rad = 50
set H = 100 

join.dbg planar.mesh $pt1 planar.mesh $pt2 $rad $H auto
subdivide.opt join.mesh
hd.opt min.inp mesh=join.mesh jobname=join1 nmin=10 movie=yes inner_volume_scale=$vexp fix_z_cut=0 outputMesh=yes > hd_opt.out 
mv join1.mesh pore_${vexp}.mesh
