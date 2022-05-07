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
   transform 
end 

"""
* rcutbond: include all bonds (i,j) such that rij <= rcutbond 
* `rcutenv`: include all bond environment atoms k such that `|rk - mid| <= rcutenv` 
* `filter` : `filter(X) == true` if particle `X` is to be included; `false` if to be discarded from the environment
"""
function bonds(at::Atoms, rcutbond, rcutenv, filter, transform = identity) 
   nlist_bond = neighbourlist(at, rcutbond) 
   nlist_env = neighbourlist(at, rcutenv)
   return BondsIterator(at, nlist_bond, rcutbond, nlist_env, rcutenv, filter, transform)
end

function Base.iterate(iter::BondsIterator, state=(1,0))
   i, q = state 
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
   rrij = Rs[q]  # position of neighbour (in shifted cell) relative to i
   # ssj = Rs[q] - iter.at.X[j]   # shift of atom j into shifted cell
   
   # now we construct the environment 
   Js, env = _get_bond(iter.nlist_env, i, iter.at.X[i], j, rrij, 
                       iter.filter, iter.transform)

   return (i, Js, env), (i, q)
end


function _get_bond(nlist, i, rri, j, rrij, filter, transform)
   Js_i, Rs_i = neigs(nlist, i)
   @show rrij 
   @show Js_i[1], Rs_i[1] 

   rrmid = (rri + 0.5 * rrij)
   rrbond = rrij
   TX = typeof(transform( State(rr = rrbond, rr0 = rrbond, be = :bond) ))
   env = TX[]  ; sizehint!(env, length(Js_i) ÷ 4)
   Js = Int[]  ; sizehint!(Js,  length(Js_i) ÷ 4)
   

   # add the bond itself first.    
   q_bond = 0 
   for (q, rrq) in enumerate(Rs_i)
      # rr = rrq + rri - rrmid 
      if rrq ≈ rrij   # TODO: replace this with checking for j and shift!
         @assert Js_i[q] == j
         q_bond = q 
         X = State(rr = rrbond,  # or should it be zero? 
                   rr0 = rrbond, 
                   be = :bond) |> transform 
         push!(env, X)
         push!(Js, Js_i[q])
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
      X = State(rr = rr, 
                rr0 = rrbond, 
                be = :env) |> transform 
      if filter(X)
         push!(env, X)
         push!(Js, Js_i[q])
      end
   end

   return Js, env
end
