function ABCDQR_em1(y,u,nx::Int;max_iter::Int=100,tol::Float64=1e-6,txo::Bool=true)
	#=
	estimate A, C, Q, R, m1, P1 using the EM algorithm
	and the subspace algorithm for the starting point
	for model
	
	x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	
	cov(w_{t},v_{t}) = [Q 0;0 R]
	x1 -> N(m1,P1)
	
	javier.cara@upm.es, 2017-01
	
	ss=true : estationary algorithm
	
	=#	
		
	# data as a matrix and by rows
	y,ny,nt = byrow(y)
	u,nu,aux = byrow(u)
	
	# starting values from the ssi algorithm
	i = nx+1
	Ai,Ci,Qi,Ri,Si = ACQRS_sub(y,nx,i)
	Bi = zeros(nx,nu)
	Di = zeros(ny,nu)	
	
	# x11 = C^{-1}(y_1 - v_1)
	m1i = Ci\y[:,1:1]
	P1i = zeros(nx,nx)

	A,B,C,D,Q,R,m1,P1,loglikv,aic = ABCDQR_em(y,u,Ai,Bi,Ci,Di,Qi,Ri,m1i,P1i,max_iter,tol,txo)
	
	# output
	return ABCDQR(A,B,C,D,Q,R,m1,P1,loglikv[end],aic)

end

