# Path to dll

dqawse = normpath(joinpath(@__DIR__, "../../deps/src/libdqawse.dll"))

# Try calling the test triple_() function.  Fortran is pass-by-reference, so triple_ expects a Ref{Cdouble}.

x = Float64(π)

(@ccall dqawse.triple_(x::Ref{Cdouble})::Cdouble), 3x
