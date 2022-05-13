
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
   return ACEBondPotentialBasis(models, inds, V.cutoff)
end

# TODO: 
#   - nparams 
#   - get_params
#   - set_params! 

Base.length(basis::ACEBondPotentialBasis) = 
      sum(length(inds) for (_, inds) in basis.inds)

# --------------------------------------------------------

import JuLIP: energy, forces, virial 
import ACE: evaluate, evaluate_d, grad_config

# overload the initiation of the bonds iterator to correctly extract the 
# right cutoffs. 
bonds(at::Atoms, calc::ACEBondCalc) = 
         bonds( at, calc.cutoff.rcutbond, 
                calc.cutoff.rcutbond/2 + calc.cutoff.zcutenv, 
                (r, z) -> env_filter(r, z, calc.cutoff) )

_get_model(calc::ACEBondCalc, zi, zj) = 
      calc.models[(min(zi, zj), max(zi,zj))]

function energy(calc::ACEBondPotential, at::Atoms)
   E = 0.0 
   for (i, j, rrij, Js, Rs, Zs) in bonds(at, calc)
      # find the right ace model 
      ace = _get_model(calc, at.Z[i], at.Z[j])
      # transform the euclidean to cylindrical coordinates
      env = eucl2cyl(rrij, at.Z[i], at.Z[j], Rs, Zs)
      # evaluate 
      Eij = evaluate(ace, env)
      E += Eij.val
   end
   return E
end


function energy(basis::ACEBondPotentialBasis, at::Atoms)
   E = zeros(Float64, length(basis))
   Et = zeros(Float64, length(basis))
   for (i, j, rrij, Js, Rs, Zs) in bonds(at, basis)
      # find the right ace model 
      ace = _get_model(basis, at.Z[i], at.Z[j])
      # transform the euclidean to cylindrical coordinates
      env = eucl2cyl(rrij, at.Z[i], at.Z[j], Rs, Zs)
      # evaluate 
      ACE.evaluate!(Et, ace, ACE.ACEConfig(env))
      E += Et 
   end
   return E 
end



function forces(calc::ACEBondPotential, at::Atoms)
   F = zeros(SVector{3, Float64}, length(at))
   for (i, j, rrij, Js, Rs, Zs) in bonds(at, calc)
      Zi, Zj = at.Z[i], at.Z[j]
      # find the right ace model 
      ace = _get_model(calc, Zi, Zj)
      # transform the euclidean to cylindrical coordinates
      env = eucl2cyl(rrij, Zi, Zj, Rs, Zs)
      # evaluate 
      dV_cyl = grad_config(ace, env)
      # transform back? 
      dV_drrij, dV_dRs = rrule_eucl2cyl(rrij::SVector, Zi, Zj, Rs, Zs, dV_cyl)
      # assemble the forces 
      F[i] -= dV_drrij 
      F[j] += dV_drrij 
      for (k, dv) in zip(Js, dV_dRs)
         F[k] -= dv
         F[i] += 0.5 * dv 
         F[j] += 0.5 * dv 
      end
   end
   return F 
end