using JuLIP, StaticArrays, LinearAlgebra
using ACE: State, filter 
# using ACE: BondEnvelope, filter, State, CylindricalBondEnvelope

struct BondsIterator 
   at
   env
end 

bonds(at::Atoms, env::BondEnvelope) = BondsIterator(at, env)

function iterate(iter::BondsIterator, state=(0,0))
   at = iter.at
   env = iter.env 
   i, q = state 
   Js, _ = neigs(nlist, i)
   Js = sort(Js) # sort in increasing order 

   # nothing left to do 
   if i >= length(at) && q >= length(Js)
      return nothing 
   end 

   # increment: 
   #   case 1: we haven't yet exhausted the current neighbours. 
   #           just increment the q index pointing into Js 
   if j < length(Js)
      q += 1
    #  case 2: if i < length(at) but q >= length(Js) then we need to 
    #          increment the i index and get the new neighbours 
   elseif 
      i += 1 
      Js, _ = neigs(nlist, i)
      Js = sort(Js) # sort in increasing order 
      q = 1 
      # skip all j = Js[q] for which we've already visited (i, j). Namely, if 
      # r_ij <= rcut and j < i then we've already found the bond (i, j) when we 
      # were sitting at j and looking at that neighbourhood. This should even 
      # work without minimum image convention. 
      # TODO: check this argument and then implement this skipping step!
   end 

   j = Js[q]
   
end


"""
get_bond(at,i,j)
input:  an Atoms at, atomic index i and j
output: index of the nearest atom having index j, i.e., the one really forms a 
        bond with atom i
"""
function get_bond(at::JuLIP.Atoms, i::Int64, j::Int64)
   nlist = JuLIP.neighbourlist(at,20.0)
   neigh_i = JuLIP.Potentials.neigsz(nlist,at,i)
   idx = findall(isequal(j),neigh_i[1])
   idx = idx[findmin(norm.(neigh_i[2][idx]))[2]]
   return idx
end

"""
get_state(at,i,j,env)
input:  an Atoms at, atomic index i and j, and an local envelop env
output: local state of the bond formed by atom i and j

The variables i, j can be replaced by a vector of tuples, so that a set of states
can be found.
"""
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

# get_state(at,index,env=CylindricalBondEnvelope(18.0,10.0,10.0)) = [ get_state(at,i,j,env) for (i,j) in index ]

# TODO: convert to cylinderical coordinate? Will need a default direction of ϕ...