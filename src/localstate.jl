using JuLIP, StaticArrays, LinearAlgebra
using ACE: State, filter 
# using ACE: BondEnvelope, filter, State, CylindricalBondEnvelope

struct BondsIterator 
   at
   nlist_bond
   rcutbond 
   nlist_env 
   rcutenv 
   filter 
end 

function bonds(at::Atoms, rcutbond, rcutenv, filter) 
   nlist_bond = neighbourlist(at, rcutbond) 
   nlist_env = neighbourlist(at, rcutenv)
   return BondsIterator(at, nlist_bond, rcutbond, nlist_env, rcutenv, filter)
end

function iterate(iter::BondsIterator, state=(0,0))
   i, q = state 
   Js, Rs = neigs(iter.nlist_bond, i)

   # nothing left to do 
   if i >= length(at) && q >= length(Js)
      return nothing 
   end 

   # increment: 
   #   case 1: we haven't yet exhausted the current neighbours. 
   #           just increment the q index pointing into Js 
   if q < length(Js)
      q += 1
      # here we could build in a rule to skip any pair for which we don't 
      # want to do the computation. 

    #  case 2: if i < length(at) but q >= length(Js) then we need to 
    #          increment the i index and get the new neighbours 
   elseif 
      i += 1 
      Js, Rs = neigs(iter.nlist_bond, i)
      Js = sort(Js) # sort in increasing order 
      q = 1 
   end 

   j = Js[q]   # index of neighbour (in central cell)
   rrj = Rs[q]  # position of neighbour (in shifted cell)
   ssj = Rs[q] - at.X[j]   # shift of atom j into shifted cell
   
   # now we construct the environment 

   
end

function _get_bond(nlist, i, rri, j, rrj, filter)
   Js_i, Rs_i = neigs(nlist, i)

   rrmid = 0.5 * (rri + rrj)
   rrbond = rrj - rri
   env = [] 
   for (q, rrq) in zip(Js_i, Rs_i) 
      rr = rrq + rri - rrmid 
      if rr ≈ rrj   # TODO: replace this with checking for j and shift!
         @assert Js_i[q] == j 
         X = State(rr = rrbond,  # or should it be zero? 
                   rr0 = rrbond, 
                   be = :bond)
         push!(env, X)
      elseif filter(rr)
         X = State(rr = rr, 
                   rr0 = rrbond, 
                   be = :env)
         push!(env, X)
      end
   end
   return identity.(env)
end
