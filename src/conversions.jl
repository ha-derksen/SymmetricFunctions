# Conversions between different types of symmetric functions

# SymC(f::SymS) converts a Symmetric function of degree n in s-basis to a symmetric function in c-basis.
# SymC(A::SymS,k::Integer=A.degree) does the same, but only computes the n_partitions(n,k) first coefficients
# in the c-basis, which are the partitions of n with largest part <=k. (This is used for the recursion.)
"""
`SymC(f)` converts a symmetric function `f` to type `SymC` and writes it as a linear combination in the c-basis.

#Examples

```julia-repl    
julia> SymC(s_(3,2,1)) # computes the irreducible S_6 character s_(3,2,1)
16*c_(1,1,1,1,1,1)-2*c_(3,1,1,1)-2*c_(3,3)+c_(5,1)
```
"""
function SymC(f::SymS,k::Integer=f.degree)::SymC
    n=f.degree
    CoeffType=eltype(f.coeffs)
    if k>n
        return SymC(SymS(n,f.coeffs))
    elseif n==0
        return SymC(n,f.coeffs[1:1])
    elseif k==0   
        return SymC(n,zeros(CoeffType,0))
    end
    C=zeros(CoeffType,n_partitions(n-k))
    for i in eachindex(f.coeffs)
        if f.coeffs[i]!=0
            MN2=murnaghan_nakayama_rule2(n,i,k)
            for j in eachindex(MN2)
                C[abs(MN2[j])]+=f.coeffs[i]*sign(MN2[j])
            end
        end
    end            
    computeCC_task=@spawn CC=SymC(SymS(n-k,C),min(n-k,k))
    return SymC(n,[SymC(SymS(n,f.coeffs),k-1).coeffs;fetch(computeCC_task).coeffs])
end 





# computation of the produc of an s-symmetric function with p_(k)
function multiply_by_p(f::SymS,k::Integer)::SymS
    n=f.degree
    CoeffType=eltype(f.coeffs)
    g=SymS(n+k,zeros(CoeffType,n_partitions(n+k)))
    for i in eachindex(f.coeffs)
        if f.coeffs[i]!=0
            MN1=murnaghan_nakayama_rule1(n,i,k)
            for j in eachindex(MN1)
                g.coeffs[abs(MN1[j])]+=f.coeffs[i]*sign(MN1[j])
            end
        end
    end
    return g
end


    

# Computation of the product of an s-symmetric function with a p-symmetric function.
# The result is an s-symmetric function.
function multiply_by_SymP(f::SymS,g::SymP,k::Integer=g.degree)::SymS
    nf=f.degree
    ng=g.degree
    pg=length(g.coeffs)
    CoeffType=eltype(f.coeffs)
    if ng==0
        return g.coeffs[1]*f
    elseif k==0
        return SymS(nf+ng,zeros(CoeffType,n_partitions(ng)))
    elseif k>ng
        return multiply_by_SymP(f,g)
    elseif pg<=n_partitions(ng,k-1)
        return multiply_by_SymP(f,g,k-1)
    end
    r=n_partitions(ng,k-1)
    compute_h1_task=@spawn h1=multiply_by_SymP(f,SymP(ng,g.coeffs[1:r]),k-1)
    h=multiply_by_p(multiply_by_SymP(f,SymP(ng-k,g.coeffs[r+1:pg]),k),k)
    h1=fetch(compute_h1_task)
    h+=h1
    return h    
end




# conversion from SymP to SymS
"""
`SymS(f)` converts a symmetric function to type `SymS` and writes it as a linear combination of Schur functions.

# Examples
```julia-repl
julia> SymS(p_(3,2,1))
-s_(1,1,1,1,1,1)-s_(2,2,2)+s_(3,1,1,1)+s_(3,3)-s_(4,1,1)+s_(6)

julia> SymS(c_(2,2))
1//8*s_(1,1,1,1)-1//8*s_(2,1,1)+1//4*s_(2,2)-1//8*s_(3,1)+1//8*s_(4)
```
"""
function SymS(f::SymP)::SymS
    return multiply_by_SymP(s_(),f)
end





# z_(lambda) is a scalar that counts the number of elements in the stabilizer of an element of cycle type
# lambda in S_n,  under the conjugation action of S_n on itself.
"""
`z_(lambda)` computes the size of the centralizer group of an element in S_n with cycle type `lambda`.

# Examples
```julia-repl
julia> z_(7,7,7,5,5,5,5) # centralizer has 7^3*5^4*3!*4! elements
30870000
```
"""
function z_(args...)::Integer
      if isa(args[1],Vector)
        lambda=args[1]
    else
        lambda=collect(args)
    end
    k=prod(Vector{BigInt}(lambda))
       mu=conjugate_partition(lambda)
    r=length(mu)
    for i=1:r-1
        k*=factorial(BigInt(mu[i]-mu[i+1]))
    end
    k*=factorial(BigInt(mu[r]))
    return k
end


# Find the type of f.coeffs after division with an integer
function quotient_coeff_type(f::Sym)
    CoeffType=eltype(f.coeffs)
    return typeof(CoeffType(1)//CoeffType(2))
end

# Find the type after elements of f.coeffs are multiplied by elements of g.coeffs.
function common_coeff_type(f::Sym,g::Sym)
      return typeof(eltype(f.coeffs)(1)*eltype(g.coeffs)(1))
end

# conversion from SymC to SymP
"""
`SymP(f)` converts a symmetric function to type `SymP` and writes it as a linear combination in the power sum basis.

# Examples
```julia-repl
julia> SymP(s_(3,2,1))
1//45*p_(1,1,1,1,1,1)-1//9*p_(3,1,1,1)-1//9*p_(3,3)+1//5*p_(5,1)

julia> SymP(c_(5,5,2,2,2))
1//2400*p_(5,5,2,2,2)
```
""" 
function SymP(f::SymC)::SymP
    n=f.degree
    p=length(f.coeffs)
    QuotientType=quotient_coeff_type(f)
    a=zeros(QuotientType,p)
    for i in 1:p
       lambda=rank2partition(n,i)
        a[i]=f.coeffs[i]//z_(lambda)
    end
    return SymP(n,a)
end


 
# conversion from SymP to SymC
function SymC(f::SymP)::SymC
    n=f.degree
    p=length(f.coeffs)
    a=zeros(BigInt,p)
    for i in 1:p
       lambda=rank2partition(n,i)
        a[i]=f.coeffs[i]*z_(lambda)
    end
    return SymC(n,a)
end

#conversion from SymS to SymP
function SymP(f::SymS)::SymP
    return SymP(SymC(f))
end

# conversion from SymC to SymS
function SymS(f::SymC)::SymS
    return SymS(SymP(f))
end




# identity conversions

function SymP(f::SymP)::SymP
    return f
end


function SymC(f::SymC)::SymC
    return f
end


function SymS(f::SymS)::SymS
    return f
end
