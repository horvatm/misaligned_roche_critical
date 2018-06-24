/*
  Searching for critical points of the Kopal potential Omega of the
  misaligned binary star:

    Omega(x,y,z,params) = 
      1/r1 + q(1/r2 - x/delta^2) + 
      1/2 (1 + q) F^2 [(x cos theta' - z sin theta')^2 + y^2]
      
    r1 = sqrt(x^2 + y^2 + z^2)
    r2 = sqrt((x-delta)^2 + y^2 + z^2)

  The critical point  r is defined as
    
    Nabla Omega(r) = 0;
  
  Note:
  *  Sometimes we referee to theta as beta.
  
  Author: Martin Horvat, Jan 2018
*/ 

#include <iostream>
#include <cmath>
#include <fstream>

#include "main.h"

int main(){
  
  std::cout.precision(16);
  std::cout << std::scientific;
   
  std::ofstream fl("res1.dat"); 
 
  fl.precision(16);
  fl << std::scientific;
   
  double
    X0 = -10000,        // limits of the bounding box [X0,Y0] , [X1,Y1]
    Y0 = -10000,
    X1 = +10000,
    Y1 = +10000,
    
    Nx = 1000000,       // maximal number of steps in X direction
    Ny = 1000000,       // maximal number of steps in Y direction
    
    Nl = 50,            // maximal number of recursions
      
    Nth = 400,          // steps in theta, theta in [th_min, th_max]
    th,
    th_min = -M_PI/2,
    th_max = +M_PI/2,
    dth = (th_max - th_min)/(Nth-1),
    
    Nq = 50,            // steps in q, q in [q_min, q_max]
    q,
    q_min = 0.1,
    q_max = 10,
    dq = (q_max - q_min)/(Nq-1),
    
    NF = 50,            // steps in F, F in [F_min, F_max]
    F,
    F_min = 0.1,
    F_max = 10,
    dF = (F_max - F_min)/(NF-1),
    
    p[4], H[2][2];
    
  for (int i = 0; i < Nq; ++i) {
     p[0] = q = q_min + dq*i;
    
    for (int j = 0; j < NF; ++j) {
      F = F_min + dF*j;
      p[1] = (1 + q)*F*F;
      
      for (int k = 0; k < Nth; ++k) { 
        th = th_min + k*dth;
        
        sincos(th, p+2, p+3);
        
        Troots <double> roots;
        
        //scan_with_lines<double>(X0, Y0, X1, Y1, p, roots, Nx, Ny, Nl);
        scan_with_image<double>(X0, Y0, X1, Y1, p, roots, Nx, Ny, Nl);
        
        for (auto && x: roots)
          if (!(x[0] == 0 && x[1] == 0) && !(x[0] == 1 && x[1] == 0)) {
            Dfun(x.ptr(), H, p);
            fl << q << ' ' << F << ' ' << th << ' ' << x[0] << ' ' << x[1] << ' ' << extreme_type(H) << '\n';
          }
        
      }
     fl.flush(); 
    }
  }
  return 0;
}
