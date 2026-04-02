
# symmetric function structures: Sym, SymS, SymP SymC


# Create structure for symmetric functions in s-basis (Schur basis).
"""
`SymS` is the type of a homogeneous symmetric function expressed in Schur basis.

# Fields: 
`degree::Integer`: the degree of the symmetric function\\
`coeffs::Vector`: stores the coefficients when written in the s-basis
"""
mutable struct SymS
    degree::Integer
    coeffs::Vector
end

# Create structure for symmetric functions in p-basis (power sum basis).
"""
`SymP` is a homogeneous symmetric function expressed in the power sum basis.

# Fields: 
`degree::Integer`: the degree of the symmetric function\\
`coeffs::Vector`: stores the coefficients when written in the p-basis
"""
mutable struct SymP
    degree::Integer
    coeffs::Vector
end

# Create structure for symmetric functions in c-basis (conjugacy class basis).
# The c-basis differs from the p-basis by scalars.
"""
`SymC` is the type of a homogeneous symmetric function expressed in the conjugacy class basis.

# Fields: 
`degree::Integer`: the degree of the symmetric function\\
`coeffs::Vector`: stores the coefficients when written in the c-basis
"""
mutable struct SymC
    degree::Integer
    coeffs::Vector
end


# Sym is the abstract type of symmetric functions in any basis
"""
`Sym` is the abstract type for symmetric functions of type `SymS`, `SymP` or `SymC`.
"""
Sym=Union{SymS,SymP,SymC}



# s_(lambda) creates the corresponding Schur symmetric function (with BigInt coefficients).
"""
`s_(lambda)` creates a Schur symmetric function of type `SymS` for the partition `lambda`.

# Examples
```julia-repl
julia> s_(3,3,2,1)
s_(3,3,2,1)

julia> s_([5,5,4])
s_(5,5,4)
```
"""
function s_(args...)
    if isa(args[1],Vector)
        lambda=args[1]
    else
        lambda=collect(args)
    end
    n=sum(lambda)
    f=zeros(BigInt,n_partitions(n))
    f[partition2rank(lambda)]=1
    return SymS(n,f)
end

# exception when lambda is empty
function s_()
    return SymS(0,ones(BigInt,1))
end






# p_(lambda) creates the corresponding power sum monomial (with BigInt coefficients).
# p_(5,5,3,1,1,1)=p_(5)^2*p_(3)*p_(1)^3
"""
`p_(lambda)` creates a monomial in the power sums of type `SymP`. \\
Note that `p_(3,3,2,1)` is the same as `p_(3)*p_(3)*p_(2)*p_(1)`.

# Examples
```julia-repl
julia> p_(3,3,2,1)
p_(3,3,2,1)

julia> p_([5,5,4])
p_(5,5,4)
```
"""
function p_(args...)
     if isa(args[1],Vector)
        lambda=args[1]
    else
        lambda=collect(args)
    end
    n=sum(lambda)
    f=zeros(BigInt,n_partitions(n))
    f[partition2rank(lambda)]=1
    return SymP(n,f)
end

# exception when lambda is empty
function p_()
    return SymP(0,ones(BigInt,1))
end





# c_(lambda) creates the symmetric function corresponding to the conjugacy class of c
# (with BigInt coefficients).
# c_(lambda)=z_(lambda)*p_(lambda), where z_(lambda) is the size of the conjugacy class 
# corresponding to lambda.
"""
`c_(lambda)` creates a symmetric function of type `SymC` corresponding to a conjugacy class.\\
A schur function expressed in c-basis gives the S_n character.

# Examples
```julia-repl
julia> p_(3,3,2,1)
p_(3,3,2,1)

julia> p_([5,5,4])
p_(5,5,4)

julia> SymC(s_(3,1)) # gives values of the irreducible S_4 character corresponding to s_(3,1)
3*c_(1,1,1,1)+c_(2,1,1)-c_(2,2)-c_(4)
```
"""
function c_(args...)
     if isa(args[1],Vector)
        lambda=args[1]
    else
        lambda=collect(args)
    end
    n=sum(lambda)
    f=zeros(BigInt,n_partitions(n))
    f[partition2rank(lambda)]=1
    return SymC(n,f)
end

# exception when lambda is empty
function c_()
    return SymC(0,ones(BigInt,1))
end



# This function defines how symmetric functions are displayed.
function Base.show(io::IO, f::Sym)
    if f.degree==0
        print(io,f.coeffs[1])
    else
        j=0
        for i in eachindex(f.coeffs)
            if f.coeffs[i]!=0
                if f.coeffs[i]<0
                    print(io,"-")
                elseif j>0
                    print(io,"+")
                end
                if abs(f.coeffs[i])!=1
                    print(io,"$(abs(f.coeffs[i]))*")
                end
                sym_type=typeof(f)
                if sym_type==SymC
                    print(io,"c_(")
                elseif sym_type==SymP
                    print(io,"p_(")
                elseif sym_type==SymS
                    print(io,"s_(")
                end
                lambda=rank2partition(f.degree,i)
                print(io,lambda[1])
                k=2
                while k<=length(lambda)
                    print(io,",$(lambda[k])")
                    k+=1
                end
                print(io,")")
                j+=1
            end
        end
        if j==0
            print(io,"0")
        end
    end
end

# coeff(f,lambda) gives the coefficient of the symmetric function f with respect to
# the partition lambda.
"""
`coeff(f,lambda)` gives the coefficient of the symmetric function `f` with respect to the partition `lambda`.

#Examples
```
julia> coeff(5*s_(3,2,1)-4*s_(4,2),[4,2])
-4

julia> coeff(p_(2,2)+6*p_(3,1),[2,2])
1
```
"""
function coeff(f::Sym,lambda::Vector{<:Integer})
    return f.coeffs[partition2rank(lambda)]
end

# redefine equality

import Base: ==

function ==(f::SymS,g::SymS)
    return f.degree==g.degree && f.coeffs==g.coeffs
end

function ==(f::SymP,g::SymP)
    return f.degree==g.degree && f.coeffs==g.coeffs
end

function ==(f::SymC,g::SymC)
    return f.degree==g.degree && f.coeffs==g.coeffs
end
