# type in Julia: Pkg.test("emABCDQR")

using emABCDQR
using Base.Test

include("ABCDQR_kfilter_test.jl")
include("ABCDQR_find_test.jl")
