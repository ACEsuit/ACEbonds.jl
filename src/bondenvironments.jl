module BondEnvironments

# export AbstractCutoff, EllipsoidCutoff, SphericalCutoff, DSphericalCutoff
# export env_filter, env_transform, env_cutoff

using StaticArrays
using JuLIP: AtomicNumber, chemical_symbol
using ACE
using LinearAlgebra: norm, I
using ACEbase: evaluate, evaluate!

abstract type AbstractBondEnvironment{T} end

include("cylindrical.jl")

"""
This implements a cylindrical cutoff for the bond environments: 
* A central bond (i,j) is within the cutoff if r_ij < rcutbond. 
* A neighbour atom at position rkij (relative position to midpoint) is within 
the environment if - after transformation to (r, θ, z) coordinates, it satisfies
`r <= rcutenv` and `abs(z) <= zcutenv`.

This struct implements the resulting filter under `env_filter`. 
"""
struct CylindricalCutoff{T} <: AbstractBondEnvironment{T}
   rcutbond::T 
   rcutenv::T
   zcutenv::T
end

env_filter(r, z, cutoff::CylindricalCutoff) = 
      (r <= cutoff.rcutenv) && (abs(z) <= cutoff.zcutenv)

env_transform(rrij::SVector, Zi, Zj, 
      Rs::AbstractVector{<: SVector}, 
      Zs::AbstractVector{<: AtomicNumber}, 
      ::CylindricalCutoff)  = eucl2cyl(rrij, Zi, Zj, Rs, Zs)

rrule_env_transform(rrij::SVector, Zi, Zj, 
                Rs::AbstractVector{<: SVector}, 
                Zs::AbstractVector{<: AtomicNumber}, 
                g_cyl::AbstractMatrix{<: DState},
                ::CylindricalCutoff )  = rrule_eucl2cyl(rrij::SVector, Zi, Zj, Rs, Zs, g_cyl)
rrule_env_transform(rrij::SVector, Zi, Zj, 
                Rs::AbstractVector{<: SVector}, 
                Zs::AbstractVector{<: AtomicNumber}, 
                dV_cyl::AbstractVector{<: DState}, 
                ::CylindricalCutoff) = rrule_eucl2cyl(rrij::SVector, Zi, Zj, Rs, Zs, dV_cyl)


struct EllipsoidCutoff{T} <: AbstractBondEnvironment{T}
    rcutbond::T 
    rcutenv::T
    zcutenv::T
    end

env_filter(r, z, cutoff::EllipsoidCutoff) = ((z/cutoff.zcutenv)^2 +(r/cutoff.rcutenv)^2 <= 1)
env_cutoff(ec::EllipsoidCutoff) = ec.zcutenv + ec.rcutenv 

# function env_transform(rrij::SVector, Zi, Zj, 
#     Rs::AbstractVector{<: SVector}, 
#     Zs::AbstractVector{<: AtomicNumber}, 
#     ec::EllipsoidCutoff)
#     rij = norm(rrij)

#     #Y0 = State( rr = rrij/ec.rcutbond, be = :bond,  mu = AtomicNumber(0)) # Atomic species of bond atoms does not matter at this stage.
#     Y0 = State( rr = rrij/ec.rcutbond, mube = :bond) # Atomic species of bond atoms does not matter at this stage.
#     cfg = Vector{typeof(Y0)}(undef, length(Rs)+1)
#     cfg[1] = Y0
#     trans = _ellipse_inv_transform(rrij,rij, ec)
#     for i = eachindex(Rs)
#         #cfg[i+1] = State(rr = trans(Rs[i]), be = :env,  mu = Zs[i])
#         cfg[i+1] = State(rr = trans(Rs[i]), mube = chemical_symbol(Zs[i]))
#     end
#     return cfg 
# end

# function _ellipsoid2sphere(rrij::SVector, rij::T, ec::EllipsoidCutoff) where {T<:Real}
#     rTr = rrij * transpose(rrij)/rij^2
#     G = SMatrix{3,3}(rTr/ec.zcutenv + (I - rTr)/ec.rcutenv)
#     return r -> G * r
# end

env_transform(rrij::SVector, Zi, Zj, 
    Rs::AbstractVector{<: SVector}, 
    Zs::AbstractVector{<: AtomicNumber}, 
    ec::EllipsoidCutoff) = ellipsoid2sphere(rrij, Zi, Zj, 
                            Rs,
                            Zs, 
                            ec.zcutenv, ec.rcutenv, ec.rcutbond)


rrule_env_transform(rrij::SVector, Zi, Zj, 
                Rs::AbstractVector{<: SVector}, 
                Zs::AbstractVector{<: AtomicNumber}, 
                g_ell::AbstractMatrix{<: DState},
                ec::EllipsoidCutoff )  = rrule_ellipsoid2sphere(rrij, Zi, Zj, Rs, Zs, g_ell, 
                                            ec.zcutenv, ec.rcutenv, ec.rcutbond)

rrule_env_transform(rrij::SVector, Zi, Zj, 
                Rs::AbstractVector{<: SVector}, 
                Zs::AbstractVector{<: AtomicNumber}, 
                dV_cyl::AbstractVector{<: DState}, 
                ec::EllipsoidCutoff) = rrule_ellipsoid2sphere(rrij, Zi, Zj, Rs, Zs, dV_cyl,
                                            ec.zcutenv, ec.rcutenv, ec.rcutbond)
            

function skewedhousholderreflection(rr0::SVector{3}, zc::T, rc::T) where {T<:Real}
    r02 = sum(rr0.^2)
    if r02 == 0
        return SMatrix{3,3}(1.0/rcutbond*I)
    end
    zc_inv, rc_inv = 1.0 ./zc, 1.0 ./ rc 
    return SMatrix{3,3}(rc_inv * I + (zc_inv - rc_inv)/r02 * rr0 * transpose(rr0) )
end

function pullback_skewedhousholderreflection(rr0::SVector{3}, zc::T, rc::T) where {T<:Real}
    H = skewedhousholderreflection(rr0,zc,rc)
    sf(rr) = skewedhousholderreflection(rr,zc,rc)
    dH = ForwardDiff.jacobian(rr -> sf(rr)[:], rr0)
    dHt = SMatrix{3,9}(dH)'
    return H, g -> SMatrix{3,3}(dHt * SVector{3}(g))
 end


function ellipsoid2sphere(rrij::SVector, Zi, Zj, 
    Rs::AbstractVector{<: SVector}, 
    Zs::AbstractVector{<: AtomicNumber}, zcutenv::T, rcutenv::T, rcutbond::T) where {T<:Real}
    @assert length(Rs) == length(Zs)
    G = skewedhousholderreflection(rrij,zcutenv, rcutenv)

    Y0 = State( rr = rrij/rcutbond, mube = :bond) # Atomic species of bond atoms does not matter at this stage.
    cfg = Vector{typeof(Y0)}(undef, length(Rs)+1)
    cfg[1] = Y0
    for i = eachindex(Rs)
        cfg[i+1] = State(rr = G * Rs[i], mube = chemical_symbol(Zs[i]))
    end
    return cfg 
end

function rrule_ellipsoid2sphere(rr0::SVector, Zi, Zj, 
    Rs::AbstractVector{<: SVector}, 
    Zs::AbstractVector{<: AtomicNumber}, 
    g_ell::AbstractMatrix{<: DState}, 
    zcutenv::T, rcutenv::T, rcutbond::T) where {T<:Real}
    lenB = size(g_ell, 1)
    lenR = length(Rs)
    @assert size(g_ell, 2) == lenR + 1 
    H = skewedhousholderreflection(rr0,zcutenv,rcutenv)
    H, pbH = pullback_skewedhousholderreflection(rr0, zcutenv, rcutenv)

    g_Rs = zeros(SVector{3, Float64}, lenB, lenR)
    # g_ell[:, 1] = derivative w.r.t. rr0 only
    g_rr0 = [ g_ell[n, 1].rr / rcutbond for n = 1:lenB ]

    for j = 1:lenR
        rrj = Rs[j] 
        for n = 1:lenB
            gj = g_ell[n, j+1].rr
            g_rr0[n] += pbH(gj)' * rrj
            g_Rs[n, j] = H' * gj
        end
    end

    return g_rr0, g_Rs 
end

function rrule_ellipsoid2sphere(rr0::SVector, Zi, Zj, 
                Rs::AbstractVector{<: SVector}, 
                Zs::AbstractVector{<: AtomicNumber}, 
                g_ell::AbstractVector{<: DState}, 
                zcutenv::T, rcutenv::T, rcutbond::T) where {T<:Real}
    lenR = length(Rs)
    @assert length(g_ell) == lenR + 1
    H = skewedhousholderreflection(rr0,zcutenv,rcutenv)
    H, pbH = pullback_skewedhousholderreflection(rr0, zcutenv, rcutenv)
    r̂0 = rr0 / norm(rr0) # ∇f(r0) = f'(r0) * r̂0

    g_Rs = zeros(SVector{3, Float64}, lenR)
    # deriv. of first element w.r.t. rr0 only 
    g_rr0 = g_ell[1].rr / rcutbond 
    for j = 1:lenR
        rrj = Rs[j] 
        gj = g_ell[j+1].rr
        #gj1 = J' * SVector(gj.r, gj.θ, gj.z)
        #@show gj
        g_rr0 += pbH(gj)' * rrj
        g_Rs[j] = H' * gj
    end

    return g_rr0, g_Rs 
end

end