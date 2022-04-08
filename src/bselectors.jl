
# BondBasisSelector(Bsel::ACE.SparseBasis; 
#                   isym=:be, bond_weight = 1.0, env_weight = 1.0) = 
#    ACE.CategorySparseBasis(isym, [:bond, :env];
#             maxorder = ACE.maxorder(Bsel), 
#             p = Bsel.p, 
#             weight = Bsel.weight, 
#             maxlevels = Bsel.maxlevels,
#             minorder_dict = Dict( :bond => 1),
#             maxorder_dict = Dict( :bond => 1),
#             weight_cat = Dict(:bond => bond_weight, :env=> env_weight) 
#          )

# function SymmetricBond_basis(ϕ::ACE.AbstractProperty, env::ACE.BondEnvelope, Bsel::ACE.SparseBasis; RnYlm = nothing, bondsymmetry=nothing, kwargs...)
#    BondSelector =  BondBasisSelector(Bsel; kwargs...)
#    if RnYlm === nothing
#        RnYlm = RnYlm_1pbasis(;   r0 = ACE.cutoff_radialbasis(env), 
#                                            rin = 0.0,
#                                            trans = PolyTransform(2, ACE.cutoff_radialbasis(env)), 
#                                            pcut = 2,
#                                            pin = 0, 
#                                            kwargs...
#                                        )
#    end
#    filterfun = _->true
#    if bondsymmetry == "Invariant"
#       filterfun = ACE.EvenL(:be, [:bond])
#    end
#    if bondsymmetry == "Covariant"
#       filterfun = x -> !(ACE.EvenL(:be, [:bond])(x))
#    end
#    Bc = ACE.Categorical1pBasis([:bond, :env]; varsym = :be, idxsym = :be )
#    B1p =  Bc * RnYlm * env
#    return ACE.SymmetricBasis(ϕ, B1p, BondSelector; filterfun = filterfun)
# end


import ACE: AbstractSparseBasis, maxorder, Prodb, Onepb, OneParticleBasis, 
            degree 
import Base: filter 

struct SparseBondBasis <: AbstractSparseBasis
   maxorder::Int
   weight::Dict
   maxlevels::Dict
   p::Float64
   besym::Symbol
   bondsym::Symbol 
   envsym::Symbol 
   ksym::Symbol
   lsym::Symbol
   weight_cat::Dict 
end


@noinline function SparseBondBasis(; 
               besym = :be, bondsym = :bond, envsym = :env, ksym = :k,  
               lsym = :l, 
               maxorder = nothing, 
               p = 1, 
               weight = Dict{Symbol, Float64}(), 
               default_maxdeg = nothing,
               maxlevels = nothing,
               weight_cat = Dict(bondsym => 1.0, envsym => 1.0), 
               )
   @show maxorder, default_maxdeg
   if (default_maxdeg != nothing) && (maxlevels == nothing)
      return SparseBondBasis(maxorder, weight, 
                          Dict("default" => default_maxdeg), 
                          p, besym, bondsym, envsym, ksym, lsym, weight_cat)
   elseif (default_maxdeg == nothing) && (maxlevels != nothing)
      return SparseBondBasis(maxorder, weight, maxlevels, 
                          p, besym, bondsym, envsym, ksym, lsym, weight_cat)
   else
      @error """Either both or neither optional arguments `maxlevels` and 
                `default_maxdeg` were provided. To avoid ambiguity ensure that 
                exactly one of these arguments is provided."""
   end
end


function Base.filter(b::Onepb, Bsel::SparseBondBasis, basis::OneParticleBasis)
   d = degree(b, basis)
   k = b[Bsel.ksym]
   return (k == 1 && b[Bsel.besym] == Bsel.envsym) || 
            (d == k-1 && b[Bsel.besym] == Bsel.bondsym) 
end

function Base.filter(bb::Prodb, Bsel::SparseBondBasis, basis::OneParticleBasis)
   if length(bb) == 0; return true; end 
   has1bond = count((b[Bsel.besym] == Bsel.bondsym) for b in bb) == 1
   isinvariant = sum( b[Bsel.lsym] for b in bb ) == 0 
   return has1bond && isinvariant
end

# maxorder and maxlevel are inherited from the abstract interface 

level(b::Union{Prodb, Onepb}, Bsel::SparseBondBasis, basis::OneParticleBasis) =
      cat_weighted_degree(b, Bsel, basis)


# Category-weighted degree function
cat_weighted_degree(b::Onepb, Bsel::SparseBondBasis, basis::OneParticleBasis) =
      degree(b, basis, Bsel.weight) * Bsel.weight_cat[getproperty(b, Bsel.isym)]

cat_weighted_degree(bb::Prodb, Bsel::SparseBondBasis, basis::OneParticleBasis) = (
      length(bb) == 0 ? 0.0
                      : norm(cat_weighted_degree.(bb, Ref(Bsel), Ref(basis)), Bsel.p)
      )