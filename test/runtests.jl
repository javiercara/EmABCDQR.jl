# type in Julia: Pkg.test("EmABCDQR")

using EmABCDQR
using Base.Test

include("../src/ABCDQR_simula_test.jl")
include("../src/ABCDQR_kfilter_test.jl")
include("../src/ABCDQR_em_test.jl")
include("../src/ABCDQR_em1_test.jl")
