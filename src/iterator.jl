using JuLIP, StaticArrays, LinearAlgebra
using ACE: State, filter 
using JuLIP.Potentials: neigsz
# using ACE: BondEnvelope, filter, State, CylindricalBondEnvelope

# TODO: make this type-stable
struct BondsIterator 
   at
   nlist_bond
   nlist_env 
   filter 
end 

"""
* rcutbond: include all bonds (i,j) such that rij <= rcutbond 
* `rcutenv`: include all bond environment atoms k such that `|rk - mid| <= rcutenv` 
* `filter` : `filter(X) == true` if particle `X` is to be included; `false` if to be discarded from the environment
"""
function bonds(at::Atoms, rcutbond, rcutenv, filter) 
   nlist_bond = neighbourlist(at, rcutbond; recompute=true, storelist=false) 
   nlist_env = neighbourlist(at, rcutenv; recompute=true, storelist=false)
   return BondsIterator(at, nlist_bond, nlist_env, filter)
end

function Base.iterate(iter::BondsIterator, state=(1,0))
   i, q = state 
   # store temporary arrays for those...
   Js, Rs = neigs(iter.nlist_bond, i)

   # nothing left to do 
   if i >= length(iter.at) && q >= length(Js)
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
   elseif i < length(iter.at) 
      i += 1 
      Js, Rs = neigs(iter.nlist_bond, i)
      q = 1 
   else 
      return nothing 
   end 

   j = Js[q]   # index of neighbour (in central cell)
   rr0 = rrij = Rs[q]  # position of neighbour (in shifted cell) relative to i
   # ssj = Rs[q] - iter.at.X[j]   # shift of atom j into shifted cell
   
   # now we construct the environment 
   Js_e, Rs_e, Zs_e = _get_bond_env(iter, i, j, rrij)

   return (i, j, rrij, Js_e, Rs_e, Zs_e), (i, q)
end


function _get_bond_env(iter::BondsIterator, i, j, rrij)
   # TODO: store temporary arrays 
   Js_i, Rs_i, Zs_i = neigsz(iter.nlist_env, iter.at, i)

   rri = iter.at.X[i]
   rrmid = rri + 0.5 * rrij
   Js = Int[]; sizehint!(Js,  length(Js_i) ÷ 4)
   Rs = typeof(rrij)[]; sizehint!(Rs,  length(Js_i) ÷ 4)
   Zs = AtomicNumber[]; sizehint!(Zs,  length(Js_i) ÷ 4)

   ŝ = rrij/norm(rrij) 
   
   # find the bond and remember it; 
   # TODO: this could now be integrated into the second loop 
   q_bond = 0 
   for (q, rrq) in enumerate(Rs_i)
      # rr = rrq + rri - rrmid 
      if rrq ≈ rrij   # TODO: replace this with checking for j and shift!
         @assert Js_i[q] == j
         q_bond = q 
         break 
      end
   end
   if q_bond == 0 
      error("the central bond neigbour atom j was not found")
   end

   # now add the environment 
   for (q, rrq) in enumerate(Rs_i)
      # skip the central bond 
      if q == q_bond; continue; end 
      # add the rest provided they fall within the provided filter 
      rr = rrq + rri - rrmid 
      z = dot(rr, ŝ)
      r = norm(rr - z * ŝ)
      if iter.filter(r, z)
         push!(Js, Js_i[q])
         push!(Rs, rr)
         push!(Zs, Zs_i[q])
      end
   end

   return Js, Rs, Zs 
end
