import numpy as np
import numba as nb

@nb.njit
def two_stream(tau, w0, u0, Rsfc):
    nz = len(tau)
    
    gt = np.zeros((nz),np.float64)
    gam1 = np.empty((nz),np.float64)
    gam2 = np.empty((nz),np.float64)
    gam3 = np.empty((nz),np.float64)
    gam4 = np.empty((nz),np.float64)
    lamda = np.empty((nz),np.float64)
    cap_gam = np.empty((nz),np.float64)
    e1 = np.empty((nz),np.float64)
    e2 = np.empty((nz),np.float64)
    e3 = np.empty((nz),np.float64)
    e4 = np.empty((nz),np.float64)
    tauc = np.empty((nz+1),np.float64)
    direct = np.empty((nz+1),np.float64)
    cp0 = np.empty((nz),np.float64)
    cpb = np.empty((nz),np.float64)
    cm0 = np.empty((nz),np.float64)
    cmb = np.empty((nz),np.float64)
    A = np.empty((nz*2),np.float64)
    B = np.empty((nz*2),np.float64)
    D = np.empty((nz*2),np.float64)
    E = np.empty((nz*2),np.float64)
    y1 = np.empty((nz),np.float64)
    y2 = np.empty((nz),np.float64)
    
    gam1 = np.sqrt(3.0)*(2.0-w0*(1+gt))/2.0
    gam2 = np.sqrt(3.0)*w0*(1.0-gt)/2.0
    gam3 = (1.0-np.sqrt(3.0)*gt*u0)/2.0
    gam4 = 1.0 - gam3
    u1 = 1.0/np.sqrt(3.0)

    lamda = (gam1**2.0 - gam2**2.0)**(0.50)
    
    cap_gam = gam2 / (gam1 + lamda)
    
    e1 = 1.0 + cap_gam*np.exp(-lamda*tau)
    e2 = 1.0 - cap_gam*np.exp(-lamda*tau)
    e3 = cap_gam + np.exp(-lamda*tau)
    e4 = cap_gam - np.exp(-lamda*tau)
    
    tauc[0] = 0.0
    for i in range(1,nz+1):
        tauc[i] = tauc[i-1] + tau[i-1]
        
    Fs_pi = 1.0
    direct[0] = u0*Fs_pi
    
    for i in range(nz):
        facp = w0[i]*Fs_pi*((gam1[i]-1.0/u0)*gam3[i]+gam4[i]*gam2[i])
        facm = w0[i]*Fs_pi*((gam1[i]+1.0/u0)*gam4[i]+gam2[i]*gam3[i])
        et0 = np.exp(-tauc[i]/u0)
        etb = et0*np.exp(-tau[i]/u0)
        denom = lamda[i]**2.0 - 1.0/u0**2.0

        direct[i+1] = u0*Fs_pi*etb
        cp0[i] = et0*facp/denom
        cpb[i] = etb*facp/denom
        cm0[i] = et0*facm/denom
        cmb[i] = etb*facm/denom
        
    Ssfc = Rsfc*direct[nz]
        
    A[0] = 0.0
    B[0] = e1[0]
    D[0] = -e2[0]
    E[0] = 0.0 - cm0[0]
    
    # Odd coeficients (Equation 41)
    for i in range(nz-1):
        l = 2*i + 2
        A[l] = e2[i]*e3[i] - e4[i]*e1[i]
        B[l] = e1[i]*e1[i+1] - e3[i]*e3[i+1]
        D[l] = e3[i]*e4[i+1] - e1[i]*e2[i+1]
        E[l] = e3[i]*(cp0[i+1] - cpb[i]) + e1[i]*(cmb[i] - cm0[i+1])
    
    # Even coefficients (Equation 42)
    for i in range(nz-1):
        l = 2*i + 1
        A[l] = e2[i+1]*e1[i] - e3[i]*e4[i+1]
        B[l] = e2[i]*e2[i+1] - e4[i]*e4[i+1]
        D[l] = e1[i+1]*e4[i+1] - e2[i+1]*e3[i+1]
        E[l] = e2[i+1]*(cp0[i+1] - cpb[i]) - e4[i+1]*(cm0[i+1] - cmb[i])
    
    l = 2*nz - 1
    A[l] = e1[nz-1] - Rsfc*e3[nz-1]
    B[l] = e2[nz-1] - Rsfc*e4[nz-1]
    D[l] = 0.0
    E[l] = Ssfc - cpb[nz-1] + Rsfc*cmb[nz-1]
    
    E = solve_tridiag(A, B, D, E)
    
    for i in range(nz):
        l = 2*i + 1
        y1[i] = E[l-1]
        y2[i] = E[l]
    
    amean = np.empty((nz+1),np.float64)
    amean[0] = (1.0/u1)*(y1[0]*e3[0]-y2[0]*e4[0]+cp0[0]) + direct[0]/u0;
    for i in range(nz):
        # J_n*4*pi = amean (Equation 49)
        amean[i+1] = (1.0/u1)*(y1[i]*(e1[i]+e3[i])+y2[i]*(e2[i]+e4[i]) + cpb[i]+cmb[i]) + direct[i+1]/u0

    return amean
    
@nb.njit
def solve_tridiag(a, b, c, d):
    n = len(a)
    # d = np.empty((n),np.float64)
    n = n - 1
    
    c[0] /= b[0]
    d[0] /= b[0]

    for i in range(1,n):
        c[i] /= b[i] - a[i]*c[i-1]
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1])

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1])

    for i in range(n-1, -1, -1):
        d[i] -= c[i]*d[i+1]
    return d
    
    
    
    

    