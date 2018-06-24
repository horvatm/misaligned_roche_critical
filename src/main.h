#pragma once

#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

/* 
 Return sign of x
*/ 
template <class T>
int sgn(const T & x) {
  if (x > 0) return +1; 
  if (x < 0) return -1;
  return 0;
}

/*
  Class for storing a root [x,y] representing x + I y in complex plane.  
*/

template <typename T>
struct Troot {
  
  T x[2];
  
  Troot(T y[2]) {
    x[0] = y[0];
    x[1] = y[1];
  }
  
  Troot(const T & x0, const T & x1) {
    x[0] = x0;
    x[1] = x1;  
  }
  
  T * ptr() {return x;}
  
  
  T& operator[](const std::size_t & idx)       { return x[idx]; }
  const T& operator[](const std::size_t & idx) const { return x[idx]; }
};

/*
  Class for storing roots and checking theirs uniqueness.
*/ 
template <class T>
struct Troots: public std::vector<Troot<T>> {
  
  void add(T x[2], const T & eps = 0, const T & min = 0){
    
    bool ok = false;
    
    for (auto && t: *this) 
      if ( 
          std::abs(t[0] - x[0]) <= eps*std::max(std::abs(x[0]), std::abs(t[0])) + min &&
          std::abs(t[1] - x[1]) <= eps*std::max(std::abs(x[1]), std::abs(t[1])) + min
         ){ 
        ok = true;
        break;
      }
      
    if (!ok) this->emplace_back(x);
  }
  
  void add(const T &x0, const T &x1, const T & eps = 0, const T & min = 0){
    
    bool ok = false;
        
    for (auto && t: *this) 
      if ( 
          std::abs(t[0] - x0) <= eps*std::max(std::abs(x0), std::abs(t[0])) + min &&
          std::abs(t[1] - x1) <= eps*std::max(std::abs(x1), std::abs(t[1])) + min
         ){ 
        ok = true;
        break;
      }
    
    if (!ok) this->emplace_back(x0,x1);
  }
  
  void print(){
    for (auto  && x: *this) std::cout << x[0] << ' ' << x[1] << '\n';  
  }
};

/* 
  System of two functions f - gradient of the Kopal potential Omega
    f := (fx,fy) = grad(Omega)
*/ 

template  <class T>
void fun(T *x, T *f, T *p, unsigned choice = 3){
  T 
    q = p[0], 
    b = p[1],
    s = p[2], // sin(beta) 
    c = p[3], // cos(beta)
    
    u = x[0] - 1,
    
    x2 = x[0]*x[0],
    u2 = u*u,
    z2 = x[1]*x[1],
    
    v = b*(x[0]*c - x[1]*s),
    
    f1 = 1/(x2 + z2),
    f2 = 1/(u2 + z2);
  
  f1 *= std::sqrt(f1);
  f2 *= std::sqrt(f2);
     
  if ((choice  & 1u) == 1u) f[0]  = -f1*x[0] - q*(f2*u + 1) + c*v;
  if ((choice  & 2u) == 2u) f[1]  = -(f1 + q*f2)*x[1] - s*v;  
}

/*
  Gradient of the two functions -- Hessian of the Kopal potential Omega
    grad(f) = Hessian(Omega)
*/
template  <class T>
void Dfun(T x[2], T H[2][2], T *p, T *f = 0){
  
  T 
    q = p[0], 
    b = p[1],
    s = p[2], 
    c = p[3],
    
    u = x[0] - 1,
    
    x2 = x[0]*x[0],
    u2 = u*u,
    z2 = x[1]*x[1],
    
    v = b*(x[0]*c - x[1]*s),
    
    f12 = 1/(x2 + z2),
    f22 = 1/(u2 + z2),
    
    f11 = std::sqrt(f12),
    f21 = std::sqrt(f22),
    
    f13 = f12*f11,
    f23 = f22*f21,
    
    f15 = f13*f12,
    f25 = f23*f22;
    
  if (f) {
    f[0]  = -f13*x[0] - q*(f23*u + 1) + c*v;
    f[1]  = -(f13 + q*f23)*x[1] - s*v;
  }
  
  H[0][0] = 3*x2*f15 - f13 + q*(2*u2 - z2)*f25 + b*c*c;

  H[0][1] = H[1][0] = 3*x[1]*(x[0]*f15 + q*u*f25) - b*c*s;

  H[1][1] = -f13 -q*f23 + 3*z2*(f15 + q*f25) + b*s*s;
}

template <class T>
int extreme_type(T H[2][2]){
  
  T h = H[0][0] + H[1][1],
    det = H[0][0]*H[1][1] - H[0][1]*H[0][1];

  int s[2];
 
  if (h != 0) { 
    s[0] = sgn(h);
    s[1] = s[0]*sgn(det);
  } else {
    s[0] = +1;
    s[1] = -1;  
  }
  
  if (s[1] == 0) std::cerr << "Singular hessian\n";
  
  return s[0] + s[1];
}

/* 
  Regularized system of two functions fR -- without singularities
    
    fR = ((x0-1)^2 + x1^2)^(3/2) (x0^2 + x1^2)^(3/2) f
    
  If could be simpler pre-factor, but it designed to simplify expressions.
*/  
template  <class T>
void funR(T *x, T *f, T *p, unsigned choice = 3){
  T 
    q = p[0], 
    b = p[1],
    s = p[2], 
    c = p[3],
    
    u = x[0] - 1,
    u2 = u*u,
    x2 = x[0]*x[0],
    z2 = x[1]*x[1],
    
    v = b*(x[0]*c - x[1]*s),
    
    f1 = x2 + z2,
    f2 = u2 + z2;
  
  f1 *= std::sqrt(f1);
  f2 *= std::sqrt(f2);
  
  T  f12 = f1*f2;
  
  if ((choice  & 1u) == 1u) f[0]  = -f2*x[0] - q*(f1*u + f12) + c*v*f12;
  if ((choice  & 2u) == 2u) f[1]  = -(f2 + q*f1)*x[1] - s*v*f12;  
}



/*
  Solving Eq. A x = b
*/

template <class T> bool solve2D(T A[2][2], T b[2], T x[2]){
  T det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
      
  if (det == 0) return false;
  
  x[0] = (A[1][1]*b[0] - A[0][1]*b[1])/det;
  x[1] = (A[0][0]*b[1] - A[1][0]*b[0])/det;
  return true;
}


/* 
  Newton-Raphson method for finding roots of the system of equations
  around the point x:
  
    x _{n+1} = x_n - Hess(x_n)^-1 f(x_n)
*/
template <class T>
bool newt(T *x, T *p, 
  int max_iter = 100, 
  const T &eps = 10*std::numeric_limits<T>::epsilon(),
  const T & min = 10*std::numeric_limits<T>::min()) {
  

  long double 
    y[2] = {x[0], x[1]},
    q[4] = {p[0],p[1],p[2],p[3]},
    t[2], dx[2], f[2], H[2][2];
    
  
  do {
    
    Dfun(y, H, q, f);
    
    solve2D(H, f, dx);
    
    for (int i = 0; i < 2; ++i) y[i] -= dx[i];
    
    t[0] = std::max(std::abs(dx[0]), std::abs(dx[1]));
    t[1] = std::max(std::abs(y[0]), std::abs(y[1]));
    
  } while (t[0] > t[1]*eps + min && --max_iter > 0);

  if (max_iter == 0) {
    std::cerr 
      << "Too many iterations\n"
      << "p=" << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << p[3] << '\n'; 
    return false;
  }
  
  x[0] = y[0], x[1] = y[1];
  
  return true;
}

/* Recursive 2D bisection of the triangle for nr levels on an rectangle:
 
  (x0,y1) ----- (x1,y1)      (s01, q01) ----- (s11, q11)
    |              |              |                |
    |              |              |                |
  (x0,y0) ----- (x1,y0)      (s00, q00) ----- (s01, q01)
  
  At parameters p.
*/
template <class T>
void bis(
  const T & x0, const T & y0, 
  const T & x1, const T & y1, 
  const int &s00, const int &q00,
  const int &s01, const int &q01,
  const int &s10, const int &q10,
  const int &s11, const int &q11,
  T *p, const int & nr,
  Troots <T> & roots) {
  
  const T eps = 10*std::numeric_limits<T>::epsilon();
  const T min = 10*std::numeric_limits<T>::min();
  
  // End of the road
  //
  bool 
    d0 = std::abs(x0- x1) <= eps*std::max(std::abs(x0), std::abs(x1)) + min,
    d1 = std::abs(y0- y1) <= eps*std::max(std::abs(y0), std::abs(y1)) + min;
  
  if ( d0 || d1 || nr == 0) {
    
    #if 0
    std::cout 
        << "B:" 
        << x0 << ' ' << x1 << ' ' << y0  << ' ' << y1 << '|'
        << s00 << ' ' << q00 << '|' 
        << s01 << ' ' << q01 << '|'
        << s10 << ' ' << q10 << '|'
        << s11 << ' ' << q11 << '\n';
    #endif 
    
    // Discussing corners
    if (s00 == 0 && q00 == 0) {
      roots.add(x0, y0, eps, min);
    } else if (s11 == 0 && q11 == 0) {
      roots.add(x1, y1, eps, min);
    } else if (s01 == 0 && q01 == 0) {
      roots.add(x0, y1, eps, min);
    } else if (s10 == 0 && q10 == 0) {
      roots.add(x1, y0, eps, min);
    } else { 
      // discussing special known roots
      if ( x0 <= 0 && 0 <= x1 && y0 <= 0 && 0 <= y1) {
        roots.add(0., 0., eps, min);
      } else if ( x0 <= 1 && 1 <= x1 && y0 <= 0 && 0 <= y1) {
        roots.add(1., 0., eps, min);
      } else {   
        T x[2] = {(x0+x1)/2, (y0+y1)/2};
        bool st = newt(x, p, eps, min);
        
       // std::cout << "N:" << x[0] << ' ' << x[1] << "|" << st << '\n';
        if (st) roots.add(x, eps, min);
      }
    }
    return;
  }

  
  //
  // Divide rectangle into 4 x Rectangles
  //
  
  
  int i, j,
      s[3][3],
      q[3][3];
  
  T h[3] = {x0, (x0 + x1)/2, x1},
    v[3] = {y0, (y0 + y1)/2, y1},
    x[2], f[2];
  
  s[0][0] = s00; q[0][0] = q00;
  s[0][2] = s01; q[0][2] = q01;
  s[2][0] = s10; q[2][0] = q10;
  s[2][2] = s11; q[2][2] = q11;
  
  i = 0; j = 1;
  x[0] = h[i]; x[1] = v[j]; funR(x, f, p);
  s[i][j] = sgn(f[0]); q[i][j] = sgn(f[1]);
  
  i = 1; j = 0;  
  x[0] = h[i]; x[1] = v[j]; funR(x, f, p);
  s[i][j] = sgn(f[0]); q[i][j] = sgn(f[1]);

  i = 1; j = 1;       
  x[0] = h[i]; x[1] = v[j]; funR(x, f, p);
  s[i][j] = sgn(f[0]); q[i][j] = sgn(f[1]);

  i = 1; j = 2;  
  x[0] = h[i]; x[1] = v[j]; funR(x, f, p);
  s[i][j] = sgn(f[0]); q[i][j] = sgn(f[1]);
  
  i = 2; j = 1;
  x[0] = h[i]; x[1] = v[j]; funR(x, f, p);
 
  s[i][j] = sgn(f[0]); q[i][j] = sgn(f[1]);
  
  //    
  // Bisection on 4 rectangles
  //
  
  for (i = 0; i < 2; ++i)
    for (j = 0; j < 2; ++j)
        
        // do bisection of there are some signs in both functions
        if (
            std::abs(s[i][j]+s[i+1][j]+s[i][j+1]+s[i+1][j+1]) != 4 && 
            std::abs(q[i][j]+q[i+1][j]+q[i][j+1]+q[i+1][j+1]) != 4
           )
          bis(
            h[i], v[j], h[i+1], v[j+1],
            s[i][j], q[i][j], s[i][j+1], q[i][j+1],
            s[i+1][j], q[i+1][j], s[i+1][j+1], q[i+1][j+1],
            p, nr-1, roots);
}


#if defined(GPU_ENABLED)
__global__ void signs_vline_gpu(
  const double  q, 
  const double  b, 
  const double  s,
  const double  c, 
  const double  x0, 
  const double  y0, 
  const double  dy, 
  const int  Ny, 
  int *sn, int *qn){

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
 
  if (idx < Ny) {
 
    double 
      x[2] = {x0, y0 + idx*dy},
      
      u = x[0] - 1,
      u2 = u*u,
      x2 = x[0]*x[0],
      z2 = x[1]*x[1],
    
      v = b*(x[0]*c - x[1]*s),
    
      f1 = x2 + z2,
      f2 = u2 + z2;
  
    f1 *= std::sqrt(f1);
    f2 *= std::sqrt(f2);
  
    double f12 = f1*f2, f[2];
  
    f[0]  = -f2*x[0] - q*(f1*u + f12) + c*v*f12;
    f[1]  = -(f2 + q*f1)*x[1] - s*v*f12;  
    
    int *P[2] = {sn + idx, qn + idx};
    
    for (int j = 0; j < 2; ++j) 
      if (f[j] > 0) 
        *(P[j]) = 1;  
      else if (f[j] < 0)
        *(P[j]) = -1;
      else 
        *(P[j]) = 0;
  }
}


__global__ void signs_point_gpu(
  const double  q, 
  const double  b, 
  const double  s,
  const double  c, 
  const double  x0, 
  const double  y0, 
  const double  dx,
  const double  dy,
  const int  Nx,
  const int  Ny,
  unsigned char *m){

  int 
    idx = blockIdx.x * blockDim.x + threadIdx.x,
    idy = blockIdx.y * blockDim.y + threadIdx.y;
 
  if (idx < Nx && idy < Ny) {
 
    double 
      x[2] = {x0 + idx*dx, y0 + idy*dy},
      
      u = x[0] - 1,
      u2 = u*u,
      x2 = x[0]*x[0],
      z2 = x[1]*x[1],
    
      v = b*(x[0]*c - x[1]*s),
    
      f1 = x2 + z2,
      f2 = u2 + z2;
  
    f1 *= std::sqrt(f1);
    f2 *= std::sqrt(f2);
  
    double f12 = f1*f2, f[2];
  
    f[0]  = -f2*x[0] - q*(f1*u + f12) + c*v*f12;
    f[1]  = -(f2 + q*f1)*x[1] - s*v*f12;  
    
    *(m + idx*Ny + idy) = 
      (f[0] < 0 ? 0 : (f[0] > 0 ? 2 : 0)) +
      ((f[1] < 0 ? 0 : (f[1] > 0 ? 2 : 0)) << 4);
  }
}

template <class T>
void calc_vline_gpu(
  T p[4], 
  const T & x0, 
  const T & y0, 
  const T & dy, 
  const int & Ny, 
  int * sn, 
  int * qn){

  int size = Ny*sizeof(int);
  
  int *sn_, *qn_;
  
  cudaMalloc(&sn_, size);
  cudaMalloc(&qn_, size);

  int blockSize = 1000;
  int numBlocks = (Ny + blockSize - 1) / blockSize;

  signs_vline_gpu<<<numBlocks, blockSize>>>(p[0], p[1], p[2], p[3], x0, y0, dy, Ny, sn_, qn_);
  //cudaDeviceSynchronize();  
  
  cudaMemcpy(sn, sn_, size, cudaMemcpyDeviceToHost);
  cudaMemcpy(qn, qn_, size, cudaMemcpyDeviceToHost);

  cudaFree(sn_);
  cudaFree(qn_);
}

#else

template <class T>
void calc_vline(
  T p[4], 
  const T & x0, 
  const T & y0, 
  const T & dy, 
  const int & Ny, 
  int * sn, 
  int * qn){

  #pragma omp parallel for
  for (int i = 0; i < Ny; ++i) {    
    T y[2] = {x0, y0 + i*dy}, g[2];
    
    funR(y, g, p);
    
    sn[i] = sgn(g[0]);
    qn[i] = sgn(g[1]);
  }
}

#endif

  
/*
  Scanning the rectangle from left to right:

  (x0,y1) ----- (x1,y1)
    |              |
    |              |
  (x0,y0) ----- (x1,y0)
  
  
  Nx - number of bins in x direction
  Ny - number of bins in y direction
  p - parameters of the map
*/ 

template <class T>
void scan_with_lines(
  const T & x0, const T & y0, 
  const T & x1, const T & y1,
  T p[4],  Troots <T> & roots, 
  int Nx = 1000, int Ny = 1000, 
  int max_levels = 40) {

  T dx = (x1 - x0)/Nx,
    dy = (y1 - y0)/Ny;
      
  int *s = new int [2*Ny], *sn = s + Ny, 
      *q = new int [2*Ny], *qn = q + Ny;
 
  // calculating signs for the first two columns (two x's)  
  #if defined(GPU_ENABLED)
  calc_vline_gpu(p, x0, y0, dy, Ny, sn, qn);
  #else
  calc_vline(p, x0, y0, dy, Ny, sn, qn);    
  #endif

  for (int i = 0; i < Nx - 1; ++i) {
   
    // calculating sign for the next columns
    for (int j = 0; j < Ny; ++j) {
      s[j] = sn[j];
      q[j] = qn[j];
    }
    
    #if defined(GPU_ENABLED)
    calc_vline_gpu(p, x0 + (i+1)*dx, y0, dy, Ny, sn, qn);
    #else
    calc_vline(p, x0 + (i+1)*dx, y0, dy, Ny, sn, qn);    
    #endif
    
    for (int j = 0; j < Ny - 1; ++j)
    
      // try bisection
      if (
          std::abs(s[j]+s[j+1]+sn[j]+sn[j+1]) != 4 && 
          std::abs(q[j]+q[j+1]+qn[j]+qn[j+1]) != 4
         ) 
        bis(
          x0 + i*dx, y0 + j*dy, 
          x0 + (i+1)*dx, y0 + (j+1)*dy, 
          s[j], q[j], s[j+1], q[j+1],
          sn[j], qn[j], sn[j+1], qn[j+1],
          p, max_levels, roots);
  }
  
  delete [] s;
  delete [] q;
}

/*
  Scanning the rectangle from left to right:

  (x0,y1) ----- (x1,y1)
    |              |
    |              |
  (x0,y0) ----- (x1,y0)
  
  
  Nx - number of bins in x direction
  Ny - number of bins in y direction
  p - parameters of the map
*/ 

template <class T>
void scan_with_image(
  const T & x0, const T & y0, 
  const T & x1, const T & y1,
  T p[4],  Troots <T> & roots, 
  int Nx = 1000, int Ny = 1000, 
  int max_levels = 40) {

  T dx = (x1 - x0)/Nx,
    dy = (y1 - y0)/Ny;
  
  // generate image  
  #if defined(GPU_ENABLED)
  int 
    BLOCKSIZE_x = 16,
    BLOCKSIZE_y = 16;
  
  unsigned char *r;
  
  cudaMallocManaged(&r, Nx*Ny);
  
  dim3 gridSize((Nx + BLOCKSIZE_x - 1)/BLOCKSIZE_x, (Ny + BLOCKSIZE_y - 1)/BLOCKSIZE_y);
  dim3 blockSize(BLOCKSIZE_x, BLOCKSIZE_y);

  signs_point_gpu<<<gridSize, blockSize>>>(p[0], p[1], p[2], p[3], x0, y0, dx, dy, Nx, Ny, r);
  cudaDeviceSynchronize();
  
  #else 
  
  unsigned char *r = new unsigned char [unsigned(Nx)*Ny];

  #pragma omp parallel for
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j){
      
      T y[2] = {x0 + i*dx, y0 + j*dy}, g[2];
    
      funR(y, g, p);
    
      r[i*Ny + j] = 
        (g[0] < 0 ? 0 : (g[0] > 0 ? 2 : 1)) + 
        ((g[1] < 0 ? 0 : (g[1] > 0 ? 2 : 1)) << 4);
    }
  }
  #endif
  
  // analyse image
  unsigned char *o = r, *o_ = r + Ny, *v;
     
  int s[2][2], q[2][2];

  for (int i = 0; i < Nx - 1; ++i, o = o_, o_ += Ny) {
       
    for (int j = 0; j < Ny - 1; ++j){
      v = o + j;
      
      s[0][0] = int(*v & 3) - 1;
      q[0][0] = int(*v >> 4) - 1;
    
      ++v;
      s[0][1] = int(*v & 3) - 1;
      q[0][1] = int(*v >> 4) - 1;
      
      v = o_ + j;
      s[1][0] = int(*v & 3) - 1;
      q[1][0] = int(*v >> 4) - 1;
      
      ++v;
      s[1][1] = int(*v & 3) - 1;
      q[1][1] = int(*v >> 4) - 1;
      
      
      // try bisection
      if (
        std::abs(s[0][0]+s[0][1]+s[1][0]+s[1][1]) != 4 && 
        std::abs(q[0][0]+q[0][1]+q[1][0]+q[1][1]) != 4
        ) 
        bis(
          x0 + i*dx, y0 + j*dy, 
          x0 + (i+1)*dx, y0 + (j+1)*dy, 
          s[0][0], q[0][0], s[0][1], q[0][1],
          s[1][0], q[1][0], s[1][1], q[1][1],
          p, max_levels, roots);
      
    }
  }
  
  #if defined(GPU_ENABLED)
  cudaFree(r); 
  #else 
  delete [] r;
  #endif
}

