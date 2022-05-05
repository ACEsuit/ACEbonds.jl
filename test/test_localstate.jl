
using ACE, JuLIP, ACEbonds, LinearAlgebra
using ACEbonds: bonds 

##

at = bulk(:Si, cubic=true) * 3 
set_pbc!(at, false )
rattle!(at, 0.3)


it = bonds(at, 3.0, 7.0, _ -> true) 

state = iterate(it)

# Js1, Rs1 = neigs(it.nlist_bond, 1)
# Js2, Rs2 = neigs(it.nlist_env, 1)