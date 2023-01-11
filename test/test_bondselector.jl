using ACE
property = ACE.Invariant()
species_friction = [:H]
species_env = [:Cu]
maxorder_on=2
maxdeg_on=5 
rcut_on = 7.0
r0_on=.4*rcut_on
rin_on=.4
pcut_on=2
pin_on=2
p_sel_on = 2 
minorder_dict_on = Dict{Any, Float64}()
maxorder_dict_on = Dict{Any, Float64}()
weight_cat_on = Dict(c => 1.0 for c in hcat(species_friction,species_env))
maxorder_off=maxorder_on
maxdeg_off=maxdeg_on
maxorder_off=maxorder_on
maxdeg_off=maxdeg_on
rcut_off = rcut_on
r0_off=.4*rcut_off
rin_off=.4
pcut_off=2
pin_off=2
p_sel_off = 2
minorder_dict_off = Dict{Any, Float64}()
maxorder_dict_off = Dict{Any, Float64}()
weight_cat_off = Dict(c => 1.0 for c in hcat(species_friction,species_env,[:bond]))

species = vcat(species_friction,species_env)


import ACE: filter

function filter(b::ACE.Onepb, Bsel::ACE.CategorySparseBasis, basis::ACE.OneParticleBasis) 
    return true
end

function filter(bb, Bsel::ACE.CategorySparseBasis, basis::ACE.OneParticleBasis) 
    # auxiliary function to count the number of 1pbasis functions in bb 
    # for which b.isym == s.
    num_b_is_(s) = sum([(getproperty(b, Bsel.isym) == s) for b in bb])
 
    # Within category min correlation order constaint:
    cond_ord_cats_min = all( num_b_is_(s) >= ACE.minorder(Bsel, s)
                             for s in keys(Bsel.minorder_dict) )
    # Within category max correlation order constaint:   
    cond_ord_cats_max = all( num_b_is_(s) <= ACE.maxorder(Bsel, s)
                             for s in keys(Bsel.maxorder_dict) )
 
    return cond_ord_cats_min && cond_ord_cats_max
 end


@info "Generate offsite basis"

Bsel_off = ACE.SparseBasis(; maxorder=maxorder_off, p = p_sel_off, default_maxdeg = maxdeg_off ) 
RnYlm_off = ACE.Utils.RnYlm_1pbasis(;  r0 = r0_off, 
        rin = rin_off,
        trans = PolyTransform(2, r0_off), 
        pcut = pcut_off,
        pin = pin_off, 
        Bsel = Bsel_off, 
        rcut=rcut_off,
        maxdeg= maxdeg_off 
    );

@time offsite = SymmetricBondSpecies_basis(property, Bsel_off; 
    RnYlm=RnYlm_off, species=species,
    species_minorder_dict =  minorder_dict_off,
    species_maxorder_dict =  maxorder_dict_off,
    weight_cat = weight_cat_off
    );
length(offsite)
species = vcat(species_friction,species_env)

@info "Generate onsite basis"
env_on = SphericalCutoff(rcut_on)
Bsel_on = ACE.SparseBasis(; maxorder=maxorder_on, p = p_sel_on, default_maxdeg = maxdeg_on ) 
RnYlm_on = ACE.Utils.RnYlm_1pbasis(;  r0 = r0_on, 
        rin = rin_on,
        trans = PolyTransform(2, r0_on), 
        pcut = pcut_on,
        pin = pin_on, 
        Bsel = Bsel_on, 
        rcut=rcut_on,
        maxdeg= maxdeg_on * max(1,Int(ceil(1/minimum(values(weight_cat_on)))))
    );
Zk_on = ACE.Categorical1pBasis(species; varsym = :mu, idxsym = :mu) #label = "Zk"
Bselcat_on = ACE.CategorySparseBasis(:mu, species;
maxorder = ACE.maxorder(Bsel_on), 
    p = Bsel_on.p, 
    weight = Bsel_on.weight, 
    maxlevels = Bsel_on.maxlevels,
    minorder_dict = minorder_dict_on,
    maxorder_dict = maxorder_dict_on, 
    weight_cat = weight_cat_on
    )

    @time onsite = ACE.SymmetricBasis(property, RnYlm_on * Zk_on, Bselcat_on;);