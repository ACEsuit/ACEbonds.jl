using ACEbonds, ACE, ACEbase, Test, StaticArrays, LinearAlgebra, JuLIP
using ACEbase.Testing
using ACEbonds: CylindricalCutoff
using ACE: discrete_jacobi, Rn1pBasis, scal1pbasis, Scal1pBasis,
           evaluate, evaluate_d, Trig1pBasis, @λ, 
           Product1pBasis, Categorical1pBasis, SimpleSparseBasis, 
           SymmetricBasis, Invariant, PIBasis

##

r0cut = 3.2
rcut = 4.0
zcut = 5.0

cutoff = CylindricalCutoff(r0cut, rcut, zcut)

maxdeg = 8
maxL = 5 

Jz = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = zcut, xin = -zcut)
Jr = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = rcut, xin = -rcut)
Jr0 = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = r0cut, xin = -r0cut)

Cbe = Categorical1pBasis([:bond, :env], :be, :be, "Cbe")
Zm = Scal1pBasis(:z, nothing, :m, Jz, "Zm")
Rn = Scal1pBasis(:r, nothing, :n, Jr, "Rn")
El = Trig1pBasis(maxL; varsym = :θ, lsym = :l, label = "El")

Pk = Scal1pBasis(:rij, nothing, :k, Jr0, "Jk")

##
B1p = Product1pBasis((Cbe, Pk, Rn, El, Zm))
Bsel = ACEbonds.SparseBondBasis(; maxorder = 3, 
                                  default_maxdeg = maxdeg, 
                                  weight = Dict{Symbol, Float64}(:m => 1.0, :n => 1.0, :k => 1.0, :l => 1.0), 
                                 )
ACE.init1pspec!(B1p, Bsel)
length(B1p)

basis = PIBasis(B1p, Bsel; isreal=true)
@show length(basis)

##

using ACEbonds: ACEBondPotentialBasis

zSi = AtomicNumber(:Si)
models = Dict((zSi, zSi) => basis)
inds = Dict((zSi, zSi) => 1:length(basis))
potbasis = ACEBondPotentialBasis(models, inds, cutoff)

at = rattle!(set_pbc!(bulk(:Si, cubic=true) * 3, false), 0.2)

energy(potbasis, at)