
"""
A module for computations with symmetric functions.
 
# Public names: 
`SymS`, `SymP`, `SymC`, `Sym`,`s_` ,`p_` ,`c_` ,`z_` ,`coeff`, `plethysm`, `n_partitions`, `conjugate partition`
"""
module SymmetricFunctions


using Base.Threads
# using LinearAlgebra

include("partitions.jl") 
include("murnaghan-nakayama.jl")
include("types.jl")
include("linear_arithmetic.jl")
include("conversions.jl")
include("multiplications.jl")
include("plethysm.jl")

export SymS,SymP,SymC,Sym
export z_,s_,p_,c_,coeff,plethysm,n_partitions, conjugate_partition

end
