"""
    polyeval(z, p, N)

Evaluate the polynomial ``\\sum_k c[k] z^{k-1}`` for the coefficients `c[1]`, `c[2]`, ..., `c[N]`;
that is, the coefficients are given in ascending order by power of `z`.  This macro expands
to efficient inline code that uses either Horner's method.

```jldoctest

julia> p = randn(5)
5-element Array{Float64,1}:
  0.106455 
  0.0351716
 -0.204764 
 -0.87642  
  0.260407 

julia> polyeval(0.5, p, 5)
-0.020427265391245605

```
"""
function polyeval(x, p::Vector, N)
    return evalpoly(x,@view p[1:N])
end

function polyeval(x, p::Tuple, N)
    return evalpoly(x,p[1:N])
end



"""
    calcz(vm, b0, c0, [EPS, [MAXITER, [relax]]])

Calculates the compressibility factor using a virial equation:

`z = 1 + b₀/z + c₀/vₘ²`

The algorithm transforms the values to a cubic polynomial and solves the roots using 
cardan's method.

 * `b0` b₀ parameter of the virial equation
 * `c0` c₀ parameter of the virial equation
"""
function calcz(b0, c0)
    #transformation:
    #z3 = z2 + b₀z + c0
    #z3-z2- b₀z- c0
    _c0,_b0,_1 = promote(c0,b0,1.0)
    poly =(-_c0,-_b0,-_1,_1)
    f0(x) = evalpoly(x,poly)
    xs = _1 #supposing ideal gas
    return Roots.secant_method(f0,xs)
    
end

