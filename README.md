# The `SymmetricFunctions` Module for `julia`

`SymmetricFunctions` is module for doing computations with symmetric functions. It can be used to computations such as:

1. irreducible characters of the symmetric group
2. tensor product decompositions of representations of the symmetric group
3. tensor product decompositions of representations of the general linear group
4. plethysms of Schur functors

## Symmetric Functions

Symmetric functions can be expressed as linear combinations in the power sum basis, or the Schur basis. The core of the module consists of fast algorithms that use the Murnaghan-Nakayama rule, recursion and (if there is more than one core) multi-threading to go from the power sum basis to the Schur basis and back. 

## Jupyter notebook

Basic functions are explained on the `Jupyter` notebook `src/SFnotbook.`.
