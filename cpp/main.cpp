#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include "twostream.h"

int main()
{
  int nz = 200;
  std::vector<double> tau(nz);
  std::vector<double> w0(nz);
  double u0 = 0.6427876096865394;
  double Rsfc = 0.25; 
  std::vector<double> amean(nz);
  double surface_radiance;
  
  // read in tau and w0 //
  std::ifstream myfile;
  myfile.open("tau_and_w0.txt");
  for (int i = 0; i < nz; i++)
  {
    myfile >> tau[i];
    myfile >> w0[i];
  }
  myfile.close();
  
  two_stream(tau, w0, u0, Rsfc, surface_radiance, amean);
  
  // write the solution //
  std::ofstream outfile;
  outfile.open("results/amean_cpp.txt");
  for (int i = 0; i < nz+1; i++)
  {
    outfile << amean[i];
    outfile << "\n";
  }
  outfile.close();
  
  // timing //
  
  int nt = 10000;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i<nt; i++)
    two_stream(tau, w0, u0, Rsfc, surface_radiance, amean);
  auto t2 = std::chrono::high_resolution_clock::now();
  
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "One run in C++ takes ";
  std::cout << ms_double.count()/1000.0/nt << " seconds\n";

}