#include <vector>
#include <cmath>
#include "twostream.h"

void two_stream(std::vector<double>& tau, std::vector<double>& w0, 
               double& u0, double& Rsfc, double& surface_radiance, 
               std::vector<double>& amean)
{
  
  int nz = tau.size();
  double u1;
  
  std::vector<double> gam1(nz);
  std::vector<double> gam2(nz);
  std::vector<double> gam3(nz);
  std::vector<double> gam4(nz);
  std::vector<double> gt(nz,0.0);
  std::vector<double> lambda(nz);
  std::vector<double> cap_gam(nz);
  std::vector<double> e1(nz), e2(nz), e3(nz), e4(nz);
  std::vector<double> tauc(nz+1);
  std::vector<double> direct(nz+1);
  std::vector<double> cp0(nz), cpb(nz), cm0(nz), cmb(nz);
  std::vector<double> A(nz*2), B(nz*2), D(nz*2), E(nz*2);
  std::vector<double> y1(nz), y2(nz);
  
  // Eddington Two-Stream coefficients (Table 1)
  for (int i = 0; i < nz; i++)
  {
    gam1[i] = std::sqrt(3.0)*(2.0-w0[i]*(1.0+gt[i]))/2.0;
    gam2[i] = std::sqrt(3.0)*w0[i]*(1.0-gt[i])/2.0;
    gam3[i] = (1.0-std::sqrt(3.0)*gt[i]*u0)/2.0;
    gam4[i] = 1.0 - gam3[i];
  }
  u1 = 1.0/std::sqrt(3.0);
  
  // lambda, and capital Gamma (Equations 21, 22)
  for (int i = 0; i < nz; i++)
  {
    lambda[i] = std::pow(std::pow(gam1[i],2.0) 
                - std::pow(gam2[i],2.0),0.5);
    cap_gam[i] = gam2[i] / (gam1[i] + lambda[i]);
  }
  
  // e's (Equation 44)
  for (int i = 0; i < nz; i++)
  {
    double wrk_real = std::exp(-lambda[i]*tau[i]);
    e1[i] = 1.0 + cap_gam[i]*wrk_real;
    e2[i] = 1.0 - cap_gam[i]*wrk_real;
    e3[i] = cap_gam[i] + wrk_real;
    e4[i] = cap_gam[i] - wrk_real;
  }
  
  // tauc - cumulative optical depth at the top of each layer
  tauc[0] = 0.0;
  for (int i = 1; i < nz+1; i++)
    tauc[i] = tauc[i-1] + tau[i-1];
    
  // C+ and C- (Equation 23, 24)
  double Fs_pi = 1.0; // We take the solar flux at the top of the atmosphere to be = 1
  // Note: Toon et al. 1989 defines Fs * pi = solar flux. So Fs = solar flux / pi
  direct[0] = u0*Fs_pi;
  for (int i = 0; i < nz; i++)
  {
    double facp = w0[i]*Fs_pi*((gam1[i]-1.0/u0)*gam3[i]+gam4[i]*gam2[i]);
    double facm = w0[i]*Fs_pi*((gam1[i]+1.0/u0)*gam4[i]+gam2[i]*gam3[i]);
    double et0 = std::exp(-tauc[i]/u0);
    double etb = et0*std::exp(-tau[i]/u0);
    double denom = std::pow(lambda[i],2.0) - 1.0/std::pow(u0,2.0);
    
    direct[i+1] = u0*Fs_pi*etb;
    cp0[i] = et0*facp/denom;
    cpb[i] = etb*facp/denom;
    cm0[i] = et0*facm/denom;
    cmb[i] = etb*facm/denom;
  }
  double Ssfc = Rsfc*direct[nz];
  
  // Coefficients of tridiagonal linear system (Equations 39 - 43)
  // Odd coeficients (Equation 41)
  A[0] = 0.0;
  B[0] = e1[0];
  D[0] = -e2[0];
  E[0] = 0.0 - cm0[0]; // assumes no downward diffuse flux at top-of-atmosphere
  for (int i = 0; i < nz-1; i++)
  {
    int l = 2*i + 2;
    A[l] = e2[i]*e3[i] - e4[i]*e1[i];
    B[l] = e1[i]*e1[i+1] - e3[i]*e3[i+1];
    D[l] = e3[i]*e4[i+1] - e1[i]*e2[i+1];
    E[l] = e3[i]*(cp0[i+1] - cpb[i]) + e1[i]*(cmb[i] - cm0[i+1]);
  }
  
  // Even coefficients (Equation 42)
  for (int i = 0; i < nz-1; i++)
  {
    int l = 2*i + 1;
    A[l] = e2[i+1]*e1[i] - e3[i]*e4[i+1];
    B[l] = e2[i]*e2[i+1] - e4[i]*e4[i+1];
    D[l] = e1[i+1]*e4[i+1] - e2[i+1]*e3[i+1];
    E[l] = e2[i+1]*(cp0[i+1] - cpb[i]) - e4[i+1]*(cm0[i+1] - cmb[i]);
  }
  int l = 2*nz - 1;
  A[l] = e1[nz-1] - Rsfc*e3[nz-1];
  B[l] = e2[nz-1] - Rsfc*e4[nz-1];
  D[l] = 0.0;
  E[l] = Ssfc - cpb[nz-1] + Rsfc*cmb[nz-1];
  
  solve_tridiag(A, B, D, E, 2*nz);
  
  for (int i = 0; i < nz; i++)
  {
    int l = 2*i + 1;
    y1[i] = E[l-1];
    y2[i] = E[l];
  }
  
  // amean = integral(J_n d_Omega) = J_n*4*pi
  // Above is an integration of J_n over a complete sphere, 
  // and has units of flux density, or irradiance.
  // amean(i) is for at top of layer i. amean(nz + 1) is the ground.
  
  // very top edge of atmosphere (Not in paper. Derive from Equation 17 and 31. bit confusing)
  amean[0] = (1.0/u1)*(y1[0]*e3[0]-y2[0]*e4[0]+cp0[0]) + direct[0]/u0;
  for (int i = 0; i < nz; i++)
  {
    // J_n*4*pi = amean (Equation 49)
    amean[i+1] = (1.0/u1)*(y1[i]*(e1[i]+e3[i])+y2[i]*(e2[i]+e4[i]) 
                 + cpb[i]+cmb[i]) + direct[i+1]/u0;
  }
  
  int i = nz-1;
  surface_radiance = (y1[i]*e3[i]+y2[i]*e4[i]+cmb[i])/u1 + std::exp(-tauc[i+1]/u0);
}

void solve_tridiag(std::vector<double>& a, std::vector<double>& b, 
                   std::vector<double>& c, std::vector<double>& d, int n) {
  n--;
  c[0] /= b[0];
  d[0] /= b[0];

  for (int i = 1; i < n; i++) {
    c[i] /= b[i] - a[i]*c[i-1];
    d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
  }

  d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

  for (int i = n; i-- > 0;) {
    d[i] -= c[i]*d[i+1];
  }
}

