function modalparami(A,B,C,dt,dofu)
	# 
	# Compute natural frequencies, damping ratios, mode shapes
	# and modal masses from A, C, D and dt
	#
	# dofu: degrees of freedon where inputs are applied
	# 	example : consider we measure at points (a,b,c,d,e) and 
	# 	the inputs are located at points b,e => dofu = [2,5]
	#
	# javier.cara@upm.es, 2016-02 
	#  
	
	# size
	no,ns = size(C)

	# invariants
	# ---------------------------------
	# eigenvalues and eigenvectors
	D,V = eig(A)
	# L
	Lambda = log(D)/dt

	# invariante V^{-1}*B
	VinvB = V\B

	# invariante C*V
	CV = C*V

	# modal parameters
	# ---------------------------------
	# allocating
	wm1 = zeros(ns)
	zm1 = zeros(ns)
	vm1 = zeros(no,ns)
	mm1 = zeros(ns)

	s = 1
	ss = 1
	while s < ns
		 if abs(imag(D[s])) > 1e-10 # complex number
				wm1[ss] = abs(Lambda[s])
				zm1[ss] = -real(Lambda[s])/wm1[ss]
				     
				vm0 = CV[:,s]/Lambda[s]^2
				# maximum element = 1    			
				maxval = vm0[1]
				for j in 2:no
					if abs(vm0[j]) > abs(maxval)
						maxval = vm0[j]
					end
				end
				vm1[:,ss] = real(vm0/maxval) # in theory, it should be a real vector
				
				c = (D[s] - 1.0)/( Lambda[s]*(Lambda[s] - conj(Lambda[s])) ) * (vm1[dofu,ss:ss]')/(VinvB[s:s,:]*maxval) # c is a vector with 1 component
				mm1[ss] = abs(real( c[1] ) )
		     
				# next value
				s = s+2
				ss = ss+1

		 else
		     # real eigenvalue
		     s = s+1
		 end
	end
	
	# deleting zero values
	wm = wm1[1:ss-1]
	zm = zm1[1:ss-1]    	
	vm = vm1[:,1:ss-1]        
	mm = mm1[1:ss-1] 

	# sorting frequencies			
	pos = sortperm(wm)
	wm = wm[pos]
	zm = zm[pos]	
	vm = vm[:,pos]
	mm = mm[pos]

	return wm,zm,vm,mm
	
end

# ======================================================================
function modalparami_test()
	#
	
	M = [35.0 0.0;0.0 17.5] # mass matrix
	K = [12250.0 -3500.0;-3500.0 3500.0] # stiffness matrix
	g = [0.02;0.02] # damping ratios
	
	# teorethical eigenvalues and eigenvectors
	W2,V = eig(K,M)
	W = sqrt(W2)
	
	# damping matrix
	Mm = V'*M*V
	Gm = 2*Mm*diagm(W)*diagm(g)
	G = inv(V')*Gm*inv(V) # damping matrix
	
	# state-space matrices
	Minv = [ 1/M[1,1] 0.0;0.0 1/M[2,2] ]
	Ac = [zeros(2,2) eye(2,2);-Minv*K -Minv*G] # continuous A matrix
	dt = 0.02 # time step
	A = expm(Ac*dt) # discrete A matrix
	C = [-Minv*K -Minv*G] # discrete C matrix
	Bc = [zeros(2,2);Minv]
	B = (A - eye(4))*inv(Ac)*Bc
	
	# computing the modal parameters through the state-space matrices
	dofu = [1,2]
	wm,zm,vm,mm = modalparami(A,B,C,dt,dofu)
	
	# eigenvectors with max. component = 1
	V1 = zeros(2,2)
	for j in 1:2
		maxval = V[1,j]
		if abs(V[2,j]) > abs(maxval)
			maxval = V[2,j]
		end
		V1[:,j] = V[:,j]/maxval
	end	
	# updated mass masses
	Mm1 = V1'*M*V1
	mm1 = diag(Mm1)
	
	return W,g,V1,mm1,wm,zm,vm,mm
	
end

