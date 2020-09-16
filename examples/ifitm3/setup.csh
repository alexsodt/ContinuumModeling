makeHexSquare.opt 6 200 220 200
# join: mesh v1 mesh v2 radius length

set vexp = 1.25
set rad = 45
set H = 80 
join.dbg planar.mesh 50 planar.mesh 50 $rad $H auto
hd.opt min.inp mesh=join.mesh jobname=min nmin=3 movie=no inner_volume_scale=$vexp fix_z_cut=yes outputMesh=yes
#min.opt join.mesh
#min.opt min.mesh
#min.opt min.mesh
mv min.mesh pore.mesh
