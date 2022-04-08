
using ACE: State 
using StaticArrays
# using NamedTupleTools 
using JuLIP: AtomicNumber 

function housholderreflection(rr0::SVector{3}) 
   r0 = norm(rr0)
   I3x3 = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
   if r0 == 0
      return I3x3
   end
   v = SVector(0, 0, 1) - rr0/r0 
   v̂ = v / norm(v)
   return I3x3 - 2 * v̂ * v̂' 
end

"""
transforms euclidean to cylindral bond environment 

### Input 
* `rr0` : central bond 
* `Rs` : list of environment atom positions, relative to rr0/2 
* `Zs` : associated list of atomic species 

### Output 
* An ACE configuration where each state represents either the bond in terms 
of just its length; or an atom from the environment in terms of its 
cylindrical coordinates (r, θ, z) and species mu, with origin at `rr0/2`.
"""
function eucl2cyl(rr0::SVector, Rs::AbstractVector{<: SVector}, 
                                Zs::AbstractVector{<: AtomicNumber})
   @assert length(Rs) == length(Zs)
   H = housholderreflection(rr0)
   r0 = norm(rr0)

   function _eucl2cyl(rr, mu)
      ss = H * rr 
      return State( mu = mu, 
                    r = sqrt(ss[1]^2 + ss[2]^2),
                    θ = atan(ss[2], ss[1]), 
                    z = ss[3], 
                    r0 = r0, 
                    be = :env )
   end

   Y0 = State( mu = AtomicNumber(0), r = 0.0, θ = 0.0, z = 0.0, 
               r0 = r0, be = :bond )
   cfg = Vector{typeof(Y0)}(undef, length(Rs)+1)
   cfg[1] = Y0
   for i = 1:length(Rs)
      cfg[i+1] = _eucl2cyl(Rs[i], Zs[i])
   end

   return cfg 
end

# monkey-path JuLIP since we aren't really updating it anymore 
Base.isapprox(a::AtomicNumber, b::AtomicNumber) = (a == b)

"""
For testing only. Generates a cylindrical bond environment, and returns it 
both in euclidean and cylindrical coordinates. The origin for both is the 
mid-point of the bond. 
"""
function rand_env(r0cut, rcut, zcut; Nenv = 10, species = [:Al, :Ti])
   species = AtomicNumber.(species)
   r0 = 1 + rand() 
   rr0 = randn(SVector{3, Float64})
   rr0 = r0 * (rr0/norm(rr0))
   H = ACEbonds.housholderreflection(rr0)
   Rs = SVector{3, Float64}[] 
   Zs = AtomicNumber[] 
   Xcyl = [] 
   for _ = 1:Nenv 
      z = (rand() - 0.5) * 2 * (zcut + r0cut)
      r = rand() * rcut 
      θ = rand() * 2 * π - π
      sθ, cθ = sincos(θ)
      rr = H' * SVector(cθ * r, sθ * r, z)
      mu = rand(species)
      push!(Zs, mu)
      push!(Rs, rr)
      push!(Xcyl, State(mu = mu, r = r, θ=θ, z = z, r0 = r0, be = :env))
   end 
   return rr0, Rs, Zs, Xcyl
end




