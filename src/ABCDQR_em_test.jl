nt = 1000
m = ABCDQR_simula(nt)

y = m["y"]
u = m["u"]

############################### starting point
nx = 4
i = nx+1
Ai,Ci,Qi,Ri,Si = EmACQR.ACQRS_sub(y,nx,i)
Bi = zeros(nx,2)
Di = zeros(2,2)
m1i = Ci\y[:,1:1]
P1i = 0*Ai

#######################
max_iter=10
tol=1e-6
txo=false
print("Testing ABCDQR_em: ")
A,B,C,D,Q,R,m1,P1,loglik,aic = ABCDQR_em(y,u,Ai,Bi,Ci,Di,Qi,Ri,m1i,P1i,max_iter,tol,txo)
@test sum(A.^2) > 1e-6
println("OK")

#######################
max_iter=10
tol=1e-10
txo=false
print("Testing ABCDQR_em_s: ")
A,B,C,D,Q,R,m1,P1,loglik,aic = ABCDQR_em_s(y,u,Ai,Bi,Ci,Di,Qi,Ri,m1i,max_iter,tol,txo)
@test sum(A.^2) > 1e-6
println("OK")
