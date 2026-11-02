# Path to dll

dqawse = normpath(joinpath(@__DIR__, "../../deps/src/libdqawse.dll"))

# Try calling the test triple_() function.  Fortran is pass-by-reference, so triple_ expects a Ref{Cdouble}.

x = Float64(π)

(@ccall dqawse.triple_(x::Ref{Cdouble})::Cdouble), 3x


# # Integrating an endpoint singularity with DQAWSE
# Compute:
# ```math
#   \int_0^1 \frac{x}{\sqrt{x(1-x)}}dx=\frac{\pi}{2}
# ```
# Here, `integr = 1` selects (x-a)^alfa * (b-x)^beta, so the callback supplies only $f(x) = x$.

# Fortran passes x by reference, including when calling back into Julia.
smooth_part(x::Ptr{Cdouble}) = unsafe_load(x)
f_callback = @cfunction(smooth_part, Cdouble, (Ptr{Cdouble},))

a, b = 0.0, 1.0 # Limits of integration
alfa, beta = -0.5, -0.5 # α > -1 and β > -1
integr = Cint(1) # see Fortran file, indicates the weight function
epsabs, epsrel = 1.0e-14, 1.0e-14 # absolute accuracy requested, relative accuracy requested
limit = Cint(100) # upper bound on number of subintervals

result = Ref{Cdouble}(0.0) # Fortran puts the result in this variable
abserr = Ref{Cdouble}(0.0) # estimate of the absolute error magnitude
neval = Ref{Cint}(0) # number of integrand evals
ier = Ref{Cint}(0) # 0: normal termination, otherwise error code (see Fortran)
last = Ref{Cint}(0) # Number of subintervals that was actually produced
# Below: other details (see Fortran)
alist = zeros(Cdouble, limit)
blist = zeros(Cdouble, limit)
rlist = zeros(Cdouble, limit)
elist = zeros(Cdouble, limit)
iord = zeros(Cint, limit)

@ccall dqawse.dqawse_(
    f_callback::Ptr{Cvoid},
    a::Ref{Cdouble}, b::Ref{Cdouble},
    alfa::Ref{Cdouble}, beta::Ref{Cdouble}, integr::Ref{Cint},
    epsabs::Ref{Cdouble}, epsrel::Ref{Cdouble}, limit::Ref{Cint},
    result::Ref{Cdouble}, abserr::Ref{Cdouble},
    neval::Ref{Cint}, ier::Ref{Cint},
    alist::Ptr{Cdouble}, blist::Ptr{Cdouble},
    rlist::Ptr{Cdouble}, elist::Ptr{Cdouble}, iord::Ptr{Cint},
    last::Ref{Cint},
)::Cvoid

# Compare the numerical result and its estimated error with the exact value.
exact = π / 2
@show result[] exact abserr[] abs(result[] - exact) neval[] ier[] last[]
@assert abs(result[] - exact) <= max(epsabs, epsrel * abs(exact))

# # Integrating a logarithmic endpoint singularity
# Compute
# ```math
#   \int_0^1 x\cdot\frac{\log(1-x)}{\sqrt{1-x}}dx=-\frac{32}{9}
# ```
# `integr = 3` selects (x-a)^alfa * (b-x)^beta * log(b-x).
# The Julia callback function supplies $f(x) = x$.

alfa, beta = 0.0, -0.5
integr = Cint(3)

@ccall dqawse.dqawse_(
    f_callback::Ptr{Cvoid},
    a::Ref{Cdouble}, b::Ref{Cdouble},
    alfa::Ref{Cdouble}, beta::Ref{Cdouble}, integr::Ref{Cint},
    epsabs::Ref{Cdouble}, epsrel::Ref{Cdouble}, limit::Ref{Cint},
    result::Ref{Cdouble}, abserr::Ref{Cdouble},
    neval::Ref{Cint}, ier::Ref{Cint},
    alist::Ptr{Cdouble}, blist::Ptr{Cdouble},
    rlist::Ptr{Cdouble}, elist::Ptr{Cdouble}, iord::Ptr{Cint},
    last::Ref{Cint},
)::Cvoid

exact = -32 / 9
@show result[] exact abserr[] abs(result[] - exact) neval[] ier[] last[]
@assert abs(result[] - exact) <= max(epsabs, epsrel * abs(exact))

# # Benchmark it

using BenchmarkTools

bm = @benchmark begin
    @ccall dqawse.dqawse_(
        $f_callback::Ptr{Cvoid},
        $a::Ref{Cdouble}, $b::Ref{Cdouble},
        $alfa::Ref{Cdouble}, $beta::Ref{Cdouble}, $integr::Ref{Cint},
        $epsabs::Ref{Cdouble}, $epsrel::Ref{Cdouble}, $limit::Ref{Cint},
        $result::Ref{Cdouble}, $abserr::Ref{Cdouble},
        $neval::Ref{Cint}, $ier::Ref{Cint},
        $alist::Ptr{Cdouble}, $blist::Ptr{Cdouble},
        $rlist::Ptr{Cdouble}, $elist::Ptr{Cdouble}, $iord::Ptr{Cint},
        $last::Ref{Cint},
    )::Cvoid
end

@show bm

# ### How `dqawse` works
# The performance is a bit disappointing, about 4.2us/evaluation on my machine. (1e-12 precision is achieved in 1.7us, so lower precision is already a big performance win.)
# A lot of the overhead comes from computation of quadrature rules, which can be cached.
# Implementation details can be found in [Piessens1983](@cite), specifically subsection 2.2.3.5 for the algebraic-logarithmic singularity case,
# and 2.2.3.3 for the overview of Clenshaw-Curtis quadrature involving a weight function.
# Clenshaw-Curtis quadrature uses a Chebysev polynomial basis.
# The difficulty of computing the quadrature with a weight function is avoided by calculating "modified Chebyshev moments"
# that can be relatively cheaply obtained using recurrence relations.
# Modified Chebysev moments of the specified weight function are calculated using `subroutine dqmomo`.
# `dqmomo` fills arrays `ri`, `rj`, `rg`, `rh` depending on the specified weight function.
# Once these are filled in, `subroutine dqc25s` is used for the actual integration, with the modified Chebyshev moments passed.
# `dqc25s` is a 25-point Clenshaw-Curtis quadrature routine, with a fallback to 15-point Gauss-Kronrod
# when the endpoint singularity is not included.
# The main dqawse function then bisects the integration region to h-refine until the desired accuracy is achieved.
# Precomputation of the quadrature weights is an easy way to write a speedy version of this later.
