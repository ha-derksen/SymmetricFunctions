# partitions



# make_n_partitions_table(n) creates a global (n+1)x(n+1) table  n_partitions_table
# n_partitions_table[i+1,j+1] is the number of integer partitions of i for whose parts are at most j.
function make_n_partitions_table(n::Integer)
    if n>121
        error("maximum partition size exceeded")
    end
    global n_partitions_table=Matrix{Int64}(undef,n+1,n+1)
    n_partitions_table[1,1:n+1].=1
    for i in 2:n+1
        n_partitions_table[i,1]=0
        for j in 2:n+1
            if i<j
                n_partitions_table[i,j]=n_partitions_table[i,j-1]
            else
                n_partitions_table[i,j]=n_partitions_table[i,j-1]+n_partitions_table[i-j+1,j]
            end
        end
    end
end


make_n_partitions_table(121)


# n_partitions(n) counts the number of integer partitions of n.
# n_partitions(n,k) counts the number of integer partitions of n for which all parts are at most k.
# n_partitions uses the table 'n_partitions_table'.
"""
`n_partitions(n,k)` gives the number of partitions of `n` with parts <=`k`.

# Examples
```julia-repl

julia> n_partitions(4) # there are 5 partitions of 4, namely [1,1,1,1],[2,1,1],[2,2],[3,1],[4]
5

julia> n_partitions(5,2) # there are 3 partitions of 5 with parts <=2, namely [1,1,1,1,1],[2,1,1,1],[2,2,1]
3
```
"""
function n_partitions(n::Integer,k::Integer=n)::Integer
    if n>=size(n_partitions_table,1)
        println("increasing the size of n_partitions_table to $n")
        make_n_partitions_table(n)
    end
    return n_partitions_table[n+1,k+1]
end





# rank2partition(n,rk) gives the rk-th partition of n.
# Partitions of n are ordered from small to large according to the pure lexicographic order.
# For example [1,1,1,1,1] < [2,1,1,1] < [2,2,1] < [3,1,1] < [3,2] < [4,1] < [5] if n=5.
function rank2partition(n::Integer,rk::Integer)::Vector{<:Integer}
    if rk>n_partitions(n)
        error("rank of partition exceeds the number of partitions of $n")
    elseif rk==1
        return ones(Int,n)
    end
    lambda=zeros(Int,n)
    nn=copy(n)
    i=0
    r=copy(rk)
    while nn>0
        kk=0
        while n_partitions(nn,kk)<r
            kk+=1
        end
        i+=1
        lambda[i]=kk
        r-=n_partitions(nn,kk-1)
        nn-=kk
    end
    return lambda[1:i]
end





# partition2rank(lambda::Vector{<:Integer}) gives the rank of lambda among all integer partitions of n=|lambda|.
# The partitions of n are ordered from small to large according to the pure lexicographic order.
# For example [1,1,1,1,1] < [2,1,1,1] < [2,2,1] < [3,1,1] < [3,2] < [4,1] < [5] if n=5.
function partition2rank(lambda::Vector{<:Integer})::Integer
    n=sum(lambda)
    k=length(lambda)
    j=1
    rk=1
    while j<=k && lambda[j]>1
        rk+=n_partitions(n,lambda[j]-1)
        n=n-lambda[j]

        j+=1
    end
    return rk
end


# merge_partition(lambda,mu) merges the partitions lambda and mu.
function merge_partitions(lambda::Vector{<:Integer},mu::Vector{<:Integer})::Vector{<:Integer}
    a=length(lambda)
    b=length(mu)
    if a==0
        return mu
    elseif b==0
        return lambda
    end
    nu=zeros(eltype(lambda),a+b)
    i=1
    j=1
    while i+j<=a+b+1
        if j>b || (i<=a && lambda[i]>mu[j])
            nu[i+j-1]=lambda[i]
            i+=1
        else
            nu[i+j-1]=mu[j]
            j+=1
        end
    end
    return nu
end


# conjugate_partition(lambda) computes the conjugate partition of lambda
"""
`conjugate_partition(lambda)` computes the conjugate partition of `lambda`.

# Examples
```julia-repl
julia> conjugate_partition([5,5,2,1])
5-element Vector{Int64}:
 4
 3
 2
 2
 2
```
"""
function conjugate_partition(lambda::Vector{<:Integer})::Vector{<:Integer}
    i=0
    k=length(lambda)
    if k==0
        return(lambda)
    end
    mu=zeros(eltype(lambda),lambda[1])
     mu[1:lambda[k]].=k
    for j in 1:k-1
        mu[lambda[j+1]+1:lambda[j]].=j
    end
    return mu
end