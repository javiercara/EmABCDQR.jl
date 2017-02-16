function ABCDQR_em(y::Array{Float64,2},u::Array{Float64,2},Ai::Array{Float64,2},Bi::Array{Float64,2},Ci::Array{Float64,2},Di::Array{Float64,2},
                   Qi::Array{Float64,2},Ri::Array{Float64,2},m1i::Array{Float64,1},P1i::Array{Float64,2},max_iter::Int,tol::Float64,txo::Bool)
	#
	# estimate A, B, C, D, Q, R, m1, P1 using the EM algorithm for model
	# 
	# x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	# y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	#
	# cov(w_{t},v_{t}) = [Q 0;0 R]
	# x1 -> N(m1,P1)
	#
	# javier.cara@upm.es, 2016-02
	#

	# number of states
	nx = size(Ai,1)

	# data by rows
	y,ny,nt = byrow(y)
	u,nu,nt = byrow(u)

	# initial values
	A = Ai
	B = Bi
	C = Ci
	D = Di
	Q = Qi
	R = Ri
	m1 = m1i
	P1 = P1i

	# log-likelihood values
	loglikv = zeros(max_iter)

	# Syy, Suu and Syu do not depend on the iterations
	Syy = zeros(ny,ny)
	Suu = zeros(nu,nu)
	Syu = zeros(ny,nu)
	for t in 1:nt
		Syy = Syy + y[:,t]*y[:,t]'
		Suu = Suu + u[:,t]*u[:,t]'
		Syu = Syu + y[:,t]*u[:,t]'
	end

	tol1 = 1.0
	iter = 1
	while (iter <= max_iter) && (tol1 > tol)
		tic()

		# E-step
		# ---------------------------------------------------------------------------------
		# Kalmanfilter
		(xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik) = ABCDQR_kfilter(y,u,A,B,C,D,Q,R,m1,P1)
		(xtN,PtN,Pt1tN) = ACQR_ksmoother(A,xtt,Ptt,xtt1,Ptt1)		

		loglikv[iter] = loglik
		if iter > 1
			tol1 = abs( (loglikv[iter] - loglikv[iter-1])/loglikv[iter-1] )
		end	

		# initial values
		Sxx = zeros(nx,nx)
		Sx1x = zeros(nx,nx)
		Syx = zeros(ny,nx)
		Sx1u = zeros(nx,nu)
		Sxu = zeros(nx,nu)

		# matrices Sxx, Sx1x, Syx, Sx1u, Sxu, Sx1x1
		for t in 1:nt
			Sxx = Sxx + PtN[:,:,t] + xtN[:,t]*xtN[:,t]'
			Sx1x = Sx1x + Pt1tN[:,:,t] + xtN[:,t+1]*xtN[:,t]'
			Syx = Syx + y[:,t]*xtN[:,t]'
		  Sx1u= Sx1u + xtN[:,t+1]*u[:,t]'
		  Sxu = Sxu + xtN[:,t]*u[:,t]'
		end
		Sx1x1 = Sxx - (PtN[:,:,1] + xtN[:,1]*xtN[:,1]') + (PtN[:,:,nt+1] + xtN[:,nt+1]*xtN[:,nt+1]')

		# M-step
		# -------------------------------------------------------------------------------------
		# Matrices x0e y P0e
		m1 = xtN[:,1]
		P1 = PtN[:,:,1]

		# AB
		AB = [Sx1x Sx1u]/[Sxx Sxu;Sxu' Suu]

		# Matrix A
		A = AB[:,1:nx]

		# Matrix B
		B = AB[:,nx+1:nx+nu]

		# Matrix Q
		M1 = Sx1x*A'
		M2 = Sx1u*B'
		M3 = A*Sxu*B'
		Q = Sx1x1 - M1 - M1' - M2 - M2' + M3 + M3' + A*Sxx*A' + B*Suu*B'
		Q = 1/nt*Q
		Q = (Q + Q')/2 # to make sure it's a symmetric matrix
		
		# CD
		CD = [Syx Syu]/[Sxx Sxu;Sxu' Suu]

		# Matrix Ae
		C = CD[:,1:nx]

		# Matrix Be
		D = CD[:,nx+1:nx+nu]

		# Matrix Re
		M1 = Syx*C'
		M2 = Syu*D'
		M3 = C*Sxu*D'
		R = Syy - M1' - M1 - M2 - M2' + M3 + M3' + C*Sxx*C' + D*Suu*D'
		R = 1/nt*R
		R = (R + R')/2 # to make sure it's a symmetric matrix
		
		etime = toq()
		if txo
			println( "Iter " * @sprintf("%3d",iter) * ",   @time = " * @sprintf("%5.2f",etime) * ",   logLik = " * @sprintf("%.6E",loglik) * 
			         ",   tol = ", @sprintf("%.2E",tol1) )
		end
		
		iter += 1

	end

	loglikv = loglikv[1:(iter-1)]	
	
	# Akaike Information Criterion
	P = ny*nx + nx*nx + nx*(nx+1)/2 + ny*(ny+1)/2
	aic = -2*loglikv[end] + 2*P

	# output
	return A,B,C,D,Q,R,m1,P1,loglikv,aic

end


