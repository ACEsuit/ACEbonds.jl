
"""
This implements a cylindrical cutoff for the bond environments: 
* A central bond (i,j) is within the cutoff if r_ij < rcutbond. 
* A neighbour atom at position rkij (relative position to midpoint) is within 
the environment if - after transformation to (r, θ, z) coordinates, it satisfies
`r <= rcutenv` and `abs(z) <= zcutenv`.

This struct implements the resulting filter under `env_filter`. 
"""
struct CylindricalCutoff{T}
   rcutbond::T 
   rcutenv::T
   zcutenv::T
end

env_filter(r, z, cutoff::CylindricalCutoff) = 
      (r <= cutoff.rcutenv) && (abs(z) <= cutoff.zcutenv)


struct ACEBondPotential{TM} <: AbstractCalculator
   models::Dict{Tuple{AtomicNumber, AtomicNumber}, TM}
   cutoff::CylindricalCutoff{Float64}
end


struct ACEBondPotentialBasis{TM} <: JuLIP.MLIPs.IPBasis
   models::Dict{Tuple{AtomicNumber, AtomicNumber}, TM}  # model = basis
   inds::Dict{Tuple{AtomicNumber, AtomicNumber}, UnitRange{Int}}
   cutoff::CylindricalCutoff{Float64}
end


ACEBondCalc = Union{ACEBondPotential, ACEBondPotentialBasis}


function _get_basisinds(V::ACEBondPotential)
   inds = Dict{Tuple{AtomicNumber, AtomicNumber}, UnitRange{Int}}()
   zz = sort(collect(keys(V.models)))
   i0 = 0
   for z in zz
      mo = V.models[z]
      len = length(mo.basis)
      inds[z] = (i0+1):(i0+len)   # to generalize for general models
      i0 += len
   end
   return inds 
end

_get_basisinds(V::ACEBondPotentialBasis) = V.inds

function basis(V::ACEBondPotential)
   models = Dict( [zz => model.basis for (zz, model) in V.models]... )
   inds = _get_basisinds(V)
   return ACEBondPotentialBasis{Base.valtype(models)}(models, inds)
end

# TODO: 
#   - nparams 
#   - get_params
#   - set_params! 

Base.length(basis::ACEBondPotentialBasis) = 
      sum(length(inds) for (_, inds) in basis.inds)

# --------------------------------------------------------

import JuLIP: energy, forces, virial 
import ACE: evaluate 

# overload the initiation of the bonds iterator to correctly extract the 
# right cutoffs. 
bonds(at::Atoms, calc::ACEBondCalc) = 
         bonds( at, calc.cutoff.rcutbond, 
                calc.cutoff.rcutbond/2 + calc.cutoff.zcutenv, 
                (r, z) -> env_filter(r, z, calc.cutoff) )


function energy(calc::ACEBondPotential, at::Atoms)
   E = 0.0 
   for (i, j, rrij, Js, Rs, Zs) in bonds(at, calc)
      # find the right ace model 
      is, js = min(i,j), max(i,j)
      ace = calc.models[(at.Z[is], at.Z[js])]
      # transform the euclidean to cylindrical coordinates
      env = eucl2cyl(rrij, at.Z[i], at.Z[j], Rs, Zs)
      # evaluate 
      E += evaluate(ace, env)
   end
   return E 
end


# function evaluate!(B, tmp, V::ACESitePotentialBasis, Rs, Zs, z0) 
#    # fill!(B, 0)
#    Bview = (@view B[V.inds[z0]])
#    evaluate!(Bview, V.models[z0], environment(V, Rs, Zs, z0))
#    return B 
# end

function energy(basis::ACEBondPotentialBasis, at::Atoms)
   E = zeros(Float64, length(basis))
   Et = zeros(Float64, length(basis))
   for (i, j, rrij, Js, Rs, Zs) in bonds(at, basis)
      # find the right ace model 
      is, js = min(i,j), max(i,j)
      ace = basis.models[(at.Z[is], at.Z[js])]
      # transform the euclidean to cylindrical coordinates
      env = eucl2cyl(rrij, at.Z[i], at.Z[j], Rs, Zs)
      # evaluate 
      ACE.evaluate!(Et, ace, ACE.ACEConfig(env))
      E += Et 
   end
   return E 
end
