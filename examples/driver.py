#Name: Obot Esther
#UCID: 30124569

import matplotlib.pyplot as plt
import numpy as np

from lab02.root import root_newton_raphson

def main():
    P1 =  1800
    P2 = 2500
    B1 = 1900
    B2 = 3200
    H  = 4000
    zeta_max= np.sqrt(H**2*(1/B1**2-1/B2**2))
    M = P2/P1
    zeta_list =[]
    f_list =[]
    nmax=10
    for n in range(nmax):
        f=0.25*(2*n+1)/zeta_max
        f_list.append(f)
        def g(zeta):
            return M*np.sqrt((zeta_max/zeta)**2-1)-np.tan(2*np.pi*f*zeta)
        
        
        def dgdzeta(zeta):
            return (-M*zeta_max**2/ zeta**3/np.sqrt((zeta_max/zeta)**2-1)
                    -2*np.pi*f/np.cos(2*np.pi*f*zeta)**2)
        
        zeta_r_list =[]
        for k in range(n+1):
            zeta_a = 0.25*(2*k+1)/f
            zeta_0 = zeta_a - 1e-6
            zeta_r = root_newton_raphson( zeta_0,g,dgdzeta)[0]
            zeta_r_list.append(zeta_r)
        zeta_list.append(zeta_r_list)

    f_modes =[]
    zeta_modes = []
    cl_modes= []
    lamda_modes =[]
    for m in range(4):
        f_mode =[]
        zeta_mode = []
        for k in range(m,nmax):
            f_mode.append(f_list[k])
            zeta_mode.append(zeta_list[k][m])
        f_mode = np.array(f_mode)
        zeta_mode = np.array(zeta_mode)
        cl_mode= 1/np.sqrt(1/B1**2-zeta_mode**2/H**2)
        lamda_mode = cl_mode/f_mode
        f_modes.append(f_mode)
        zeta_modes.append(zeta_mode)
        cl_modes.append(cl_mode)
        lamda_modes.append(lamda_mode)

    plt.figure(figsize=(6,8))
    plt.subplot(2,1,1)
    for f,cl in zip(f_modes,cl_modes):
        plt.plot(f,cl)
    plt.xlabel('f[Hz]')
    plt.ylabel('c_L[m/s]')
    plt.legend([f'mode {k}' for k in range(4)])

    plt.subplot(2,1,2)
    for f,lamda in zip(f_modes,lamda_modes):
        plt.plot(f,lamda)
    plt.xlabel('f[Hz]')
    plt.ylabel('lamda_L[m]')
    plt.legend([f'mode {k}' for k in range(4)])

    plt.savefig("modes.png")



    
if __name__ == '__main__':
    main()

