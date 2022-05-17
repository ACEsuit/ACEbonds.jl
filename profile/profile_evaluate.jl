using BenchmarkTools

@btime energy($pot, $at)

##

@profview begin
   for n = 1:100
      energy(pot, at)
   end
end

##

cfg = ACE.ACEConfig([ State(be = :bond, rij = rand(), θ = 0.123, r = rand(), z = rand())
                     for _=1:10])

evaluator = pot.models[(zSi, zSi)].evaluator
evaluate(evaluator, cfg)

@btime evaluate($evaluator, $cfg)

@profview begin 
   for n = 1:10_000 
      evaluate(evaluator, cfg)
   end
end

## Basis 1p Components: 
X = State(be = :bond, rij = rand(), θ = 0.123, r = rand(), z = rand())

@btime for n = 1:1000; evaluate($El, $X); end 
@btime evaluate($Pk, $X)
@btime begin P = evaluate($Pk, $X); ACE.release!(P); end
@btime evaluate($Rn, $X)
@btime evaluate($Zm, $X)
@btime evaluate($Cbe, $X)

@code_warntype(evaluate(Pk, X))
@code_warntype(evaluate(Rn, X))
@code_warntype(evaluate(Zm, X))
@code_warntype(evaluate(El, X))
@code_warntype(evaluate(Cbe, X))

let Pk=Pk, Rn=Rn, Zm=Zm, Cbe=Cbe, X=X
   @profview begin
      for n = 1:1000_000
         P = evaluate(Pk, X)
         R = evaluate(Rn, X)
         Z = evaluate(Zm, X)
         C = evaluate(Cbe, X)
         ACE.release!(P)
         ACE.release!(R)
         ACE.release!(Z)
         ACE.release!(C)
      end
   end
end

## 

A = evaluate(B1p, X)
fill!(A, 0)

@profview begin
   for n = 1:100_000
      fill!(A, 0)
      ACE.add_into_A!(A, B1p, X)
   end
end