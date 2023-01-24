# IntegerSmithNormalForm

[![CI](https://github.com/dmerkert/IntegerSmithNormalForm.jl/workflows/CI/badge.svg)](https://github.com/dmerkert/IntegerSmithNormalForm.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![codecov](https://codecov.io/gh/dmerkert/IntegerSmithNormalForm.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dmerkert/IntegerSmithNormalForm.jl)

This package computes the Smith Normal Form on integer matrices and provides the respective transformation matrices.

## Example

A simple example is the matrix
```julia-REPL
julia> A = [4 2; 0 8]
2×2 Matrix{Int64}:
 4  2
 0  8
```

where

```julia-REPL
julia> using IntegerSmithNormalForm

julia> S, B, T = snf(A)
([1 0; 4 -1], [2 0; 0 16], [0 1; 1 -2])
```

and we have

```julia-REPL
julia> S*A*T
2×2 Matrix{Int64}:
 2   0
 0  16

julia> elementary_divisors(A)
2-element Vector{Int64}:
  2
 16
```
