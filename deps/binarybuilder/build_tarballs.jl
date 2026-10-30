# 2026 10 30, Daniel Knapp, Amsterdam: Initial version (based on https://github.com/JuliaPackaging/Yggdrasil/blob/master/O/OpenSpecFun/build_tarballs.jl)

using BinaryBuilder

name = "dqawse"
version = v"0.1.0"

sources = [
    DirectorySource(
        joinpath(@__DIR__, "..", "src");
        target = "dqawse",
    ),
]

# Bash recipe for building across all platforms
script = raw"""
cd "${WORKSPACE}/srcdir/dqawse"

# It needs to be told it's on Windows
if [[ ${target} == *mingw* ]]; then
    OS=WINNT
elif [[ ${target} == *darwin* ]]; then
    OS=Darwin
fi

# Build it
make OS=${OS} CC="$CC" CXX="$CXX" FC="$FC" -j${nproc}

# Install it
make install OS=${OS} prefix=$prefix
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line. Since openspecfun uses
# Fortran for AMOS, we need the combinatorial explosion of platforms
# and GCC versions.
platforms = expand_gfortran_versions(supported_platforms(; experimental = true))

# The products that we will ensure are always built
products = [
    LibraryProduct("libdqawse", :libdqawse),
]

# Dependencies that must be installed before this package can be built
dependencies = [
]

# Build the tarballs.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; julia_compat = "1.6")

# Build trigger: 1
