
def root_newton_raphson(x0, f, dfdx):
    eps_s = 1e-9
    eps_a = 1
    k=0
    xr=x0
    while eps_a > eps_s:
        dx = -f(xr)/dfdx(xr)
        xr+=dx
        eps_a = abs(dx/xr)
        k+=1
    return xr,k,eps_a
    
    

