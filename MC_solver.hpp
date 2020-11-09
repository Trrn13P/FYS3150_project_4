#include <armadillo>

using namespace arma;

#ifndef MC_SOLVER_HPP
#define MC_SOLVER_HPP

class MC_solver{
public:
  mat spin_matrix;
  int L, N, J, mc_cycles;

  float T, k, beta;

  float exp_precalc[17];

  int PBC_index(int i);
  float atom_energy(int i, int j);
  float local_energy(int i, int j);
  float system_energy();
  float random_between_zero_and_one();
  int minus_one_or_one();
  int random_index(int start, int finish);
  void initialize_matrix_random();
  void initialize_matrix_ordered();
  float MC_step();
  void MC_solve();

  MC_solver(int L_,int J_,int mc_cycles_, float T_, float k_){
    spin_matrix = zeros(L,L);
    N = int(L*L);
    L = L_;
    J = J_;
    mc_cycles = mc_cycles_;

    T = T_; k = k_;
    beta = k*T;


    float deltaE_precalc[17] = {-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16};
    for (int i=0;i<17;i++){
      exp_precalc[i]=exp(-beta*deltaE_precalc[i]);
    }
  }
};
#endif
