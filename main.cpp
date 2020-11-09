#include <iostream>
#include <cmath>
#include <fstream>

#include <omp.h>
#include <time.h>

#include <fstream>
#include <armadillo>

#include "MC_solver.hpp"

using namespace std;
using namespace arma;

int main()
{
  srand(clock());

  int L  = 5;
  int J = 1;
  int mc_cycles = 1E6;

  float T = 2.3;
  float k = 1;


  MC_solver *test;
  test = new MC_solver(L,J,mc_cycles,T,k);
  test->initialize_matrix_ordered();
  test->MC_solve();

}
