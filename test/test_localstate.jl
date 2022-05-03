
using ACE, JuLIP, ACEbonds, LinearAlgebra

##

at = bulk(:Si, cubic=true) * 3 
set_pbc!(at, false )

ACEbonds.get_bond(at, 1, 2)
