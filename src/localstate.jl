using JuLIP, StaticArrays, LinearAlgebra
using ACE: BondEnvelope, filter, State, CylindricalBondEnvelope

# function get_atoms(fd::HDF5.File; groupname=nothing)
#    groupname === nothing && (groupname = HDF5.name(first(fd)))
#    positions = HDF5.read(fd, string(groupname,"/positions"))
#    unitcell = HDF5.read(fd, string(groupname,"/unitcell"))
#    species = HDF5.read(fd, string(groupname,"/species"))
#    atoms = JuLIP.Atoms(; X = positions, Z = species,
#                         cell = unitcell,
#                         pbc = [true, true, true])
#    return [unitcell, species, positions], atoms
# end
# 
# function get_atoms(fname::String; groupname=nothing)
#    HDF5.h5open(fname, "r") do fd
#       get_atoms(fd; groupname=groupname)
#    end
# end

```
get_bond(at,i,j)
input:  an Atoms at, atomic index i and j
output: index of the nearest atom having index j, i.e., the one really forms a 
        bond with atom i
```
function get_bond(at::JuLIP.Atoms,i::Int64,j::Int64)
   nlist = JuLIP.neighbourlist(at,20.0)
   neigh_i = JuLIP.Potentials.neigsz(nlist,at,i)
   idx = findall(isequal(j),neigh_i[1])
   idx = idx[findmin(norm.(neigh_i[2][idx]))[2]]
   return idx
end

```
get_state(at,i,j,env)
input:  an Atoms at, atomic index i and j, and an local envelop env
output: local state of the bond formed by atom i and j


The variables i, j can be replaced by a vector of tuples, so that a set of states
can be found.
```
function get_state(at::JuLIP.Atoms,i::Int64,j::Int64,env::BondEnvelope; λ = 0)
   if i==j
      error("i,j should be distinct")
   end
   nlist = JuLIP.neighbourlist(at,20.0)
   neigh_i = JuLIP.Potentials.neigsz(nlist,at,i)
   idx = get_bond(at,i,j)
   # rr - vector form of the bond
   rr = neigh_i[2][idx]
   # rr_c - centred 
   rr_c = λ * rr
   st = State(rr = rr-rr_c, rr0 = rr, be=:bond)
   # TODO: include the atomic number
   for (jj, rj) in enumerate(neigh_i[2])
      st_temp = State(rr = rj-rr_c, rr0 = rr, be=:env)
      if rj≠rr && filter(env,st_temp)
         st = [st; st_temp]
      end
   end
   return st
end

get_state(at,index,env=CylindricalBondEnvelope(18.0,10.0,10.0)) = [ get_state(at,i,j,env) for (i,j) in index ]

# TODO: convert to cylinderical coordinate? Will need a default direction of ϕ...