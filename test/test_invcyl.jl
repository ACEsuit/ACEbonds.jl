
using ACEbonds, ACE, ACEbase, Test, StaticArrays, LinearAlgebra
using ACEbase.Testing
using ACEbonds: CylindricalBondEnvelope
using ACE: discrete_jacobi, Rn1pBasis, scal1pbasis, Scal1pBasis,
           evaluate, evaluate_d, Trig1pBasis
##

r0cut = 4.0
rcut = 4.0
zcut = 2.0

maxdeg = 8
maxL = 6 

Jz = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = r0cut+zcut, xin = -r0cut-zcut)
Jr = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = rcut, xin = -rcut)
Jr0 = discrete_jacobi(maxdeg; pcut = 2, pin = 2, xcut = r0cut, xin = -r0cut)

Zm = Scal1pBasis(:rr, 3, :m, Jz, "Zm")
Rn = Scal1pBasis(:rr, nothing, :n, 
                 chain( (@λ rr -> 1/(rr[1]^2+rr[2]^2)), Jr ), "Rn")
El = Trig1pBasis(maxL; varsym = :rr, lsym = :l, label = "El")



