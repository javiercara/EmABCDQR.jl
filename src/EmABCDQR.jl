__precompile__()

module EmABCDQR

type ABCDQR
	A
	B
	C
	D
	Q
	R
	m1
	P1
	loglik::Float64
	aic::Float64
end

using EmACQR

include("ABCDQR_simula.jl")
include("ABCDQR_kfilter.jl")
include("ABCDQR_kfilter_s.jl")
include("ABCDQR_em.jl")
include("ABCDQR_em_s.jl")
include("ABCDQR_em1.jl")
include("ABCDQR_em_s1.jl")

export
	ABCDQR_simula,
	ABCDQR_kfilter, ABCDQR_kfilter_s,
	ABCDQR_em, ABCDQR_em_s,
	ABCDQR, ABCDQR_em1, ABCDQR_em_s1

end # module
