makeHexSquare.opt 6 200 220 200
# join: mesh v1 mesh v2 radius length
join.dbg planar.mesh 50 planar.mesh 50 40 80 auto
min.opt join.mesh
min.opt min.mesh
min.opt min.mesh
mv min.mesh pore.mesh
