m = ABCDQR_sim(1000,2,2)

y = m["y"]
u = m["u"]
A = m["A"]
B = m["B"]
C = m["C"]
D = m["D"]
Q = m["Q"]
R = m["R"]
x10 = 0*A[:,1]
P10 = 0*A

###########################
print("Testing ABCDQR_kfilter: ")
xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik = ABCDQR_kfilter(y,u,A,B,C,D,Q,R,x10,P10)
@test sum(xtt.^2) > 1e-6
println("OK")

#######################
print("Testing ABCDQR_kfilter_s: ")
xtt,Ptt,xtt1,Ptt1,et,St,Kt,loglik = ABCDQR_kfilter_s(y,u,A,B,C,D,Q,R,x10)
@test sum(xtt.^2) > 1e-6
println("OK")

