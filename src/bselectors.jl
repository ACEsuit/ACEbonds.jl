
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
