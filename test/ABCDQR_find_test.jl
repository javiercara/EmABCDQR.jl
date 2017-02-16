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
print("Testing ABCDQR_find: ")
nx = 4
mest = ABCDQR_find(y,u,nx,max_iter=10,ss=false,txo=false)
mest1 = ABCDQR_find(y,u,nx,max_iter=10,ss=true,txo=false)
@test sum(mest.A.^2) > 1e-6
@test sum(mest1.A.^2) > 1e-6
println("OK")


