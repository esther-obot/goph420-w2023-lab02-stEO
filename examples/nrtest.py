import matplotlib.pyplot as plt
import numpy as np

from lab02.root import root_newton_raphson

def main():
    def f(x):
        return x**2-1
    def dfdx(x):
        return 2*x
    
    xr0,k0,eps_a0 = root_newton_raphson(1.5, f, dfdx)
    print(f'xr0:{xr0}, k:{k0}, eps_a0:{eps_a0}')

    xr1,k1,eps_a1 = root_newton_raphson(-1.5, f, dfdx)
    print(f'xr1:{xr1}, k:{k1}, eps_a1:{eps_a1}')
    

if __name__ == '__main__':
    main()