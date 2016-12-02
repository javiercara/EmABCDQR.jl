function ABCDQR_em_s(y,u,Ai,Bi,Ci,Di,Qi,Ri,kmax,txo=false)
	#
	# estimate A, B, C, D, Q, R, m1, P1 using the STATIONARY EM algorithm for model
	# 
	# x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	# y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	#
	# cov(w_{t},v_{t}) = [Q 0;0 R]
	# x1 -> N(m1,P1)
	#
	# javier.cara@upm.es, 2106-02
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
	# x11 = C^{-1}y_1 - C^{-1}Du_1 - C^{-1}v_1
	m1 = C\(y[:,1] - D*u[:,1])
	P1 = zeros(nx,nx)

	# log-likelihood values
	loglikv = zeros(kmax)

	# Syy, Suu and Syu do not depend on the iterations
	Syy = zeros(ny,ny)
	Suu = zeros(nu,nu)
	Syu = zeros(ny,nu)
	for t in 1:nt
		Syy = Syy + y[:,t]*y[:,t]'
		Suu = Suu + u[:,t]*u[:,t]'
		Syu = Syu + y[:,t]*u[:,t]'
	end

	# total time
	tcalculo=0.0

	for k in 1:kmax
		tic()

		# E-step
		# ---------------------------------------------------------------------------------
		# Kalmanfilter
		(xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik) = ABCDQR_kfilter_s(y,u,A,B,C,D,Q,R,m1)
		(xtN,PtN,Pt1tN) = ACQR_ksmoother_s(A,xtt,Ptt,xtt1,Ptt1)		

		loglikv[k] = loglik

		# inutial values
		Sxx = zeros(nx,nx)
		Sx1x = zeros(nx,nx)
		Syx = zeros(ny,nx)
		Sx1u = zeros(nx,nu)
		Sxu = zeros(nx,nu)

		# matrices Sxx, Sx1x, Syx, Sx1u, Sxu, Sx1x1
		for t in 1:nt
			Sxx = Sxx + xtN[:,t]*xtN[:,t]'
			Sx1x = Sx1x + xtN[:,t+1]*xtN[:,t]'
			Syx = Syx + y[:,t]*xtN[:,t]'
			Sx1u= Sx1u + xtN[:,t+1]*u[:,t]'
			Sxu = Sxu + xtN[:,t]*u[:,t]'
		end
		Sxx = Sxx + nt*PtN
		Sx1x = Sx1x + nt*Pt1tN
		Sx1x1 = Sxx - xtN[:,1]*xtN[:,1]' + xtN[:,nt+1]*xtN[:,nt+1]'

		# M-step
		# -------------------------------------------------------------------------------------
		# Matrices x0e y P0e
		m1 = xtN[:,1]
		P1 = PtN

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
		tcalculo = tcalculo + etime
		if txo
			println(string("Iteration ",k,".   Elapsed time = ", etime))	
		end

	end

	# output
	return A,B,C,D,Q,R,m1,P1,loglikv

end
