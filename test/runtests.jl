# type in Julia: Pkg.test("emABCDQR")

using emABCDQR
using Base.Test

include("ABCDQR_kfilter_test.jl")
include("ABCDQR_em_test.jl")
include("ABCDQR_em1_test.jl")
