# plethysm


# adams(f::SymP,k) applies the k-th Adams operator to the p-symmetric function f.
function adams(f::SymP,k::Integer)::SymP
    n=f.degree
    g=SymP(n*k,zeros(eltype(f.coeffs),n_partitions(n*k)))
    for i in eachindex(f.coeffs)
        g.coeffs[partition2rank(k.*rank2partition(n,i))]=f.coeffs[i]
    end
    return g
end


# substitute_p(f::SymP,G::Vector{SymP}) substitute p_(k) in the p-polynomial f by G[k] for all k.
function substitute_p(f::SymP,G::Vector{SymP})::SymP
    nf=f.degree
    pf=length(f.coeffs)
    k=length(G)
    CoeffType=eltype(f.coeffs)
    if nf==0
        return f
    end
    ng=G[1].degree
    if k==0
        return SymP(nf*ng,zeros(CoeffType,n_partitions(nf*ng)))
    elseif k>nf
        return substitute_p(f::SymP,G[1:nf])
    elseif pf<=n_partitions(nf,k-1)
        return substitute_p(f,G[1:k-1])
    end
    r=n_partitions(nf,k-1)
    compute_h1_task=@spawn h1=G[k]*substitute_p(SymP(nf-k,f.coeffs[r+1:pf]),G)
    if k>1 
        h2=substitute_p(SymP(nf,f.coeffs[1:r]),G[1:k-1])
    else
        h2=SymP(nf*ng,zeros(CoeffType,n_partitions(nf*ng)))
    end
    h1=fetch(compute_h1_task)
    return h1+h2
end


# computation of plethysm in p-basis
"""
`plethysm(f,g)` composes the symmetric functions `f` and `g` (viewed as polynomial functors from the category\\
of vector spaces to itself).

# Examples
```julia-repl
julia> plethysm(s_(4),s_(3)) # computes the 4-th symmetric power of s_(3)
s_(4,4,4)+s_(5,4,2,1)+s_(6,2,2,2)+s_(6,4,2)+s_(6,6)+s_(7,3,2)+s_(7,4,1)+s_(8,2,2)+s_(8,4)+s_(9,3)+s_(10,2)+s_(12)

julia> plethysm(p_(4),5*p_(2,1)+7*p_(3)) # in the power sum basis, p_(k) acts as the k-th Adams operator
5*p_(8,4)+7*p_(12)
```
"""
function plethysm(f::SymP,g::SymP)::SymP
    nf=f.degree
    G=Vector{SymP}(undef,nf)
    for i in 1:nf
        G[i]=adams(g,i)
    end
    return substitute_p(f,G)
end


# computation of plethysm in arbitrary basis
function plethysm(f::Sym,g::Sym)::Sym
    gType=typeof(g)
    h=gType(plethysm(SymP(f),SymP(g)))
    CoeffType=common_coeff_type(f,g)
    return gType(h.degree,Vector{CoeffType}(h.coeffs))
end
