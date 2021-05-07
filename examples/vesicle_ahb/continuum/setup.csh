icosahedron.opt > icos.mesh
subdivide.opt icos.mesh
subdivide.opt subdiv.mesh
min.opt subdiv.mesh
scale.opt min.mesh 63
mv scale.mesh sphere.mesh
#mv sphere.inp_force sphere.inp

