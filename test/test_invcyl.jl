
using ACEbonds, ACE, ACEbase, Test, StaticArrays, LinearAlgebra, JuLIP
using ACEbase.Testing
using ACEbonds: CylindricalBondEnvelope
using ACE: discrete_jacobi, Rn1pBasis, scal1pbasis, Scal1pBasis,
           evaluate, evaluate_d, Trig1pBasis, @λ, 
           Product1pBasis, Categorical1pBasis, SimpleSparseBasis, 
           SymmetricBasis, Invariant, PIBasis

## basics 

@info("Test housholderreflection")
for _ = 1:15
   rr0 = randn(SVector{3, Float64})
   H = ACEbonds.housholderreflection(rr0)
   print_tf(@test( H * rr0 ≈ norm(rr0) * [0,0,1] ))
   print_tf(@test( H' * [0,0,1] ≈ rr0/norm(rr0) ))
   print_tf(@test H' * H ≈ I)
end


##

@info("test transformation from euclidean to cylindrical environments")

# TODO: HACK ACE States to achieve this without type piracy
#       introduce an `ace_isapprox` which defaults to `isapprox`
#       and then overload that one for Base types 
Base.isapprox(s1::Symbol, s2::Symbol) = (s1 == s2)

r0cut = 4.0 
rcut = 4.0 
zcut = 2.0 

for ntest = 1:30 
   rr0, Rs, Zs, Xs = ACEbonds.rand_env(r0cut, rcut, zcut)
   Xenv = ACEbonds.eucl2cyl(rr0, Rs, Zs)
   print_tf(@test( all(Xenv[2:end] .≈ Xs) ))
end

##


r0cut = 4.0
rcut = 4.0
zcut = 2.0

maxdeg = 8
maxL = 5 

Jz = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = r0cut+zcut, xin = -r0cut-zcut)
Jr = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = rcut, xin = -rcut)
Jr0 = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = r0cut, xin = -r0cut)

Cbe = Categorical1pBasis([:bond, :env], :be, :be, "Cbe")
Zm = Scal1pBasis(:z, nothing, :m, Jz, "Zm")
Rn = Scal1pBasis(:r, nothing, :n, Jr, "Rn")
El = Trig1pBasis(maxL; varsym = :θ, lsym = :l, label = "El")
Pk = Scal1pBasis(:r0, nothing, :k, Jr0, "Jk")

B1p = Product1pBasis((Cbe, Pk, Rn, El, Zm))


##
B1p = Product1pBasis((Cbe, Pk, Rn, El, Zm))
Bsel = ACEbonds.SparseBondBasis(; maxorder = 3, 
                                  default_maxdeg = maxdeg, 
                                  weight = Dict{Symbol, Float64}(:m => 1.0, :n => 1.0, :k => 1.0, :l => 1.0), 
                                 )
ACE.init1pspec!(B1p, Bsel)
length(B1p)

##

basis = PIBasis(B1p, Bsel; isreal=true)
@show length(basis)

##

rr0, Rs, Zs, Xs = ACEbonds.rand_env(r0cut, rcut, zcut)
Xenv = ACEbonds.eucl2cyl(rr0, Rs, Zs)

evaluate(B1p, Xenv)
evaluate(basis, ACEConfig(Xenv))

##

@info("Test invariance of the new basis")
for ntest = 1:30
   rr0, Rs, Zs, Xs = ACEbonds.rand_env(r0cut, rcut, zcut)
   Xenv = ACEbonds.eucl2cyl(rr0, Rs, Zs)
   B1 = evaluate(basis, ACEConfig(Xenv))
   Q = ACE.Random.rand_rot()
   rr0_Q = Q * rr0 
   Rs_Q = Ref(Q) .* Rs
   Xenv_Q = ACEbonds.eucl2cyl(rr0_Q, Rs_Q, Zs)
   B2 = evaluate(basis, ACEConfig(Xenv_Q))
   print_tf(@test( B1 ≈ B2 && !all(Xenv_Q .≈ Xenv) ))
end
