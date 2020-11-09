#include "MC_solver.hpp"

#include <iostream>
#include <cmath>
#include <fstream>

#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

/* This function returns N if i=-1, returns 0 if i=N
and else returns i, for PBC.
*/
int MC_solver::PBC_index(int i){
  return (L-(L-i%L)%L)%L+(L-i-(L-i)%(L+1))%(L-i%L);
}

/*
This function returns the energy of one atom i,j.
*/
//not returning -J, so remember to scale back
float MC_solver::atom_energy(int i, int j){
  return spin_matrix(i,j)*(spin_matrix(PBC_index(i),PBC_index(j-1))
  +spin_matrix(PBC_index(i),PBC_index(j+1))
  +spin_matrix(PBC_index(i-1),PBC_index(j))
  +spin_matrix(PBC_index(i+1),PBC_index(j)));
}


/*
This function returns the energy from a point i,j to it's neighbours and again
from the neighbours to itself it the neigbour is not a PBC point. This is needed
to calculate the DeltaE if we flip a spin.
*/

float MC_solver::local_energy(int i, int j){
  float energy = atom_energy(i,j)

  +spin_matrix(PBC_index(i),PBC_index(j))*(
    spin_matrix(PBC_index(i+1),PBC_index(j))* ((i+1-i%(L-1))%(L))
    +spin_matrix(PBC_index(i),PBC_index(j+1))* ((j+1-j%(L-1))%(L))

    +spin_matrix(PBC_index(i-1),PBC_index(j)) * ((1+i-i%L)%(i+1))
    +spin_matrix(PBC_index(i),PBC_index(j-1)) * ((1+j-j%L)%(j+1)));

  return energy;
}


/*
returning the total energy of the system
*/
float MC_solver::system_energy(){
  float Energy = 0;
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      Energy+=-J*atom_energy(i,j);
    }
  }
  return Energy;
}

/*
returns a random number between 0 and 1.
*/
float MC_solver::random_between_zero_and_one(){
  return rand()*1./RAND_MAX;
}

/*
This function returns a -1 or 1 randomly
*/
int MC_solver::minus_one_or_one(){
  float number = random_between_zero_and_one();
  return 2*round(number)-1;
}

/*
this function returns a random index betwen start and finish
*/
int MC_solver::random_index(int start, int finish){
  return rand() % (finish+1-start) + start;

}


void MC_solver::initialize_matrix_random(){
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
        spin_matrix(i,j)= minus_one_or_one();;
      }
    }
}

void MC_solver::initialize_matrix_ordered(){
    spin_matrix = ones(L,L);

}


float MC_solver::MC_step(){

  float r = random_between_zero_and_one();
  int i = random_index(0,L-1);
  int j = random_index(0,L-1);

  float E_j = local_energy(i,j);
  spin_matrix(i,j) = minus_one_or_one();

  float deltaE = local_energy(i,j)-E_j;

  int index = round((deltaE+16)/2);
  float exponential_func = exp_precalc[index];

  /* Metropolis test */
  float a = r-exponential_func;
  float b = a/abs(a);
  float c = (1-b)/2;
  /* Flips the matrix element back if spin is rejected*/
  spin_matrix(i,j) = spin_matrix(i,j)*(-1)*b;
  return c*deltaE;
}


void MC_solver::MC_solve(){

  float deltaE;
  vec energies = zeros(mc_cycles);
  energies(0) = system_energy();

  for(int i=1;i<mc_cycles;i++){
    cout << i << endl;
    deltaE = 0;
    for(int j=0;j<int(L*L);j++){
      deltaE += MC_step();
    }
    energies(i)=energies(i-1)+ deltaE;
  }
  ofstream outfile("test2.txt");
  outfile << energies << endl;
  outfile.close();
}
