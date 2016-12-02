function ABCDQR_kfilter(y,u,A,B,C,D,Q,R,x11,P11)
	# 
	# Kalman filter for model
	# 
	# x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
	# y_{t}   = C*x_{t} + D*u_{t} + v_{t}
	# 
	# cov(w_{t},v_{t}) = [Q 0;0 R]
	#
	# javier.cara@upm.es, 2015-10 
	# 

	(ny,nt) = size(y)
	nx = size(A,1)	

	# allocation
	xtt = zeros(nx,nt)
	Ptt = zeros(nx,nx,nt)
	xtt1 = zeros(nx,nt+1)
	Ptt1 = zeros(nx,nx,nt+1)
	et = zeros(ny,nt)
	St = zeros(ny,ny,nt)
	Kt = zeros(nx,ny,nt)
	loglik = 0.0

	# Filter 
	xtt[:,1] = x11 
	Ptt[:,:,1] = P11
	for t in 2:nt
		
		# one-step ahead prediction
		xtt1[:,t] = A*xtt[:,t-1] + B*u[:,t-1]
		Ptt1[:,:,t] = A*Ptt[:,:,t-1]*A' + Q
		
		#  innovations
		et[:,t] = y[:,t] - C*xtt1[:,t] - D*u[:,t]
		St[:,:,t] = C*Ptt1[:,:,t]*C' + R # et variance
		
		# Kalman gain
		Stinv = eye(ny,ny)/St[:,:,t] # numerically preferible to Stinv=inv(St)
		Kt[:,:,t] = Ptt1[:,:,t]*C'*Stinv # kalman gain		
			
		# filtered values
		xtt[:,t] = xtt1[:,t] + Kt[:,:,t]*et[:,t]
		Ptt[:,:,t] = (eye(nx) - Kt[:,:,t]*C)*Ptt1[:,:,t]
				
		# likelihood
		l0 = et[:,t]'*Stinv*et[:,t] # typeof(l0) = Array{Float64,1}
		loglik = loglik + log(det(St[:,:,t])) + l0[1]	
	end
	# E[x_{N+1}|y_N], Var[x_{N+1}|y_N]
	xtt1[:,nt+1] = A*xtt[:,nt] + B*u[:,nt]
	Ptt1[:,:,nt+1] = A*Ptt[:,:,nt]*A' + Q
	
	loglik =  - ny*nt/2*log(2*pi) - 0.5*loglik

	return xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik

end

