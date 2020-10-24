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
    poly =(-c0,-b0,-1.0,1.0)
    res =  cardan(poly)
    real_res = cardan_reals(res)
    return maximum(real_res)
end


function cardan(poly::NTuple{4,T}) where {T<:AbstractFloat}
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    third = T(1//3)
    _1 = T(1.0)
    a1  =  complex(_1/poly[4])
    E1  = -complex(poly[3])*a1
    E2  =  complex(poly[2])*a1
    E3  = -complex(poly[1])*a1
    s0  =  E1
    E12 =  E1*E1
    A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
    B   =  E12 - 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = sqrt(A*A - 4*B*B*B)
    if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s1 = exp(log(0.5 * (A + Δ)) * third)
    else
        s1 = exp(log(0.5 * (A - Δ)) * third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2)
end



function cardan_reals(res::NTuple{3,Complex{T}}) where T
    nan = T(NaN)
    a1, a2, a3 = res
    isreal1 = abs(imag(a1) < 16*eps(imag(a1)))
    isreal2 = abs(imag(a2) < 16*eps(imag(a2)))
    isreal3 = abs(imag(a3) < 16*eps(imag(a3)))
    r1 = real(a1)
    r2 = real(a2)
    r3 = real(a3)
    if isreal1 & isreal2 & isreal3
        return (r1,r2,r3)
    elseif !isreal1 & !isreal2
        return (r3,r3,r3)
    elseif  !isreal3 & !isreal2
        return (r1,r1,r1)
    elseif  !isreal3 & !isreal1
        return (r2,r2,r2)
    elseif !isreal1
        return (r3,r2,r2)
    elseif !isreal2
        return (r1,r3,r3)
    elseif !isreal3
        return (r1,r3,r3)
    else
        return (nan,nan,nan)
    end
end