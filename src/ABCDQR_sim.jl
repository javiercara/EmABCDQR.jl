function ABCDQR_sim(nt,no,ni)
	#
	# simulated model to test the functions
	#
	# nt: number of time instants
	# no: number of outputs, no in [1,2]
	# ni: number of inputs, ni in [1,2]
		
	# random seed for the simulations
	srand(99)
	
	# simulated mechanical system
	#
	# |    k1    ------     k2    ------
	# |---###---|  m1  |---###---|  m2  |
	# |    c1    ------     c2    ------
	# |           @  @             @  @
	# |========================================
  
  # vibration parameters
  # ----------------------------
	# nyquist frequency (Hz)
	fnq = 8
	# sampling frequency (Hz)
	fs = 2*fnq
	# time interval (s)
	dt = 1/fs
  
  # matrices for the mechanical system 
  # ----------------------------------------
  # mass matrix 
  m1 = 35
  m2 = 17.5
	M = [m1 0;0 m2]
	
	# stiffness matrix
	k1 = 8750
	k2 = 3500
	K = [k1+k2 -k2;-k2 k2]
	# theoretical natural frecuencies: 1.71 Hz, 3.31 Hz
	
	# eigenvalues and eigenectors
	D,V = eig(K,M)

	# modal matrices
	Mm = V'*M*V
	Km = V'*K*V
	
	# damping ratios
	zm = [0.02 0.02]
	
	# modal damping matrix
	Gm = zeros(2,2)
	for i in 1:2
		Gm[i,i] = 2*sqrt(Mm[i,i]*Km[i,i])*zm[i]
	end
	
	# damping matrix
	G = (V')\Gm/V
	
	# matrices for the discrete state space model
	# --------------------------------------------------------------
	# Matrix A
	Ac = [zeros(2,2) eye(2);-M\K -M\G] # Continuous Matrix A 
	A = expm(Ac*dt) # discrete Matrix A

	# Matrix C (continuous = discrete)
	C = [-M\K -M\G]

	# Matrix B
	# Js is the input location matrix
	if ni == 1
		J = [1;0]
	else
		J = eye(2)
	end
	MinvJ = M\J
	Bc = [zeros(2,ni);MinvJ] # continuous matrix B
	B = (A-eye(4,4))*inv(Ac)*Bc # discrete matrix B

	# matrix D (continuous = discrete)
	D = MinvJ
	
	# input definition
	# ---------------------------------------
	# system input u_{t} = a_{t} + a_{t-1} (MA1)
	a = 100*randn(ni,nt+1) # if 100 is not used, Q and R are too little
	u = a[:,2:end] + a[:,1:(end-1)]
	
	# state equation, no noise
	# ---------------------------------
	x = zeros(4,nt)
	for t=1:nt-1
    x[:,t+1] = A*x[:,t] + B*u[:,t]
	end
	P = cov(x')
	Q = 0.01*P # std approx 0.1
	U = chol(Q)
	w = U'*randn(4,nt)
	
	# state equation with noise
	# ------------------------------------------
	x = zeros(4,nt)
	for t = 1:nt-1
    x[:,t+1] = A*x[:,t] + B*u[:,t] + w[:,t]
	end
	
	# observation equation
	# -------------------------------------
	y = C*x + D*u
	L = cov(y')
	R = 0.01*L # std approx 0.1
	U = chol(R)
	v = U'*randn(2,nt)
	y = y + v
	y = y[1:no,:]
	
	#
	S = cov(w',v')
    
	return Dict("y"=>y,"u"=>u,"A"=>A,"B"=>B,"C"=>C,"D"=>D,"Q"=>Q,"R"=>R,"S"=>S)
    
end
