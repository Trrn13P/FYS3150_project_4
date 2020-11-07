#include <iostream>
#include <cmath>

#include <omp.h>
#include <time.h>

#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

/*
IKKE VITS MED METROPOLIS OM MAN IKKE FLIPPER
kjør ca 1 mill MC cycles
*/


/* This function returns N if i=-1, returns 0 if i=N
and else returns i, for PBC.
*/
int PBC_index(int i, int L){
  return (L-(L-i%L)%L)%L+(L-i-(L-i)%(L+1))%(L-i%L);
}

/*
This function returns the energy of one atom i,j.
*/
//not returning -J, so remember to scale back
float atom_energy(mat spin_matrix, int i, int j, int L){
  return spin_matrix(i,j)*(spin_matrix(PBC_index(i,L),PBC_index(j-1,L))
  +spin_matrix(PBC_index(i,L),PBC_index(j+1,L))
  +spin_matrix(PBC_index(i-1,L),PBC_index(j,L))
  +spin_matrix(PBC_index(i+1,L),PBC_index(j,L)));
}


/*
This function returns the energy from a point i,j to it's neighbours and again
from the neighbours to itself it the neigbour is not a PBC point. This is needed
to calculate the DeltaE if we flip a spin.
*/

float local_energy(mat spin_matrix,int i, int j, int L){
  float energy = atom_energy(spin_matrix,i,j,L)

  +spin_matrix(PBC_index(i,L),PBC_index(j,L))*(
    spin_matrix(PBC_index(i+1,L),PBC_index(j,L))* ((i+1-i%(L-1))%(L))
    +spin_matrix(PBC_index(i,L),PBC_index(j+1,L))* ((j+1-j%(L-1))%(L))

    +spin_matrix(PBC_index(i-1,L),PBC_index(j,L)) * ((1+i-i%L)%(i+1))
    +spin_matrix(PBC_index(i,L),PBC_index(j-1,L)) * ((1+j-j%L)%(j+1)));

  return energy;
}


/*
returning the total energy of the system
*/
float system_energy(mat spin_matrix,int L,float J){
  float Energy = 0;
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      Energy+=-J*atom_energy(spin_matrix,i,j,L);
    }
  }
  return Energy;
}

/*
returns a random number between 0 and 1.
*/
float random_between_zero_and_one(){
  return rand()*1./RAND_MAX;
}

/*
This function returns a -1 or 1 randomly
*/
int minus_one_or_one(){
  float number = random_between_zero_and_one();
  return 2*round(number)-1;
}

/*
this function returns a random index betwen start and finish
*/
int random_index(int start, int finish){
  return rand() % (finish+1-start) + start;

}


mat initialize_matrix_random(int L){
  mat spin_matrix = zeros(L,L);

  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
        spin_matrix(i,j)= minus_one_or_one();;
      }
    }
  return spin_matrix;
}

mat initialize_matrix_ordered(int L){
  mat spin_matrix = ones(L,L);
  return spin_matrix;
}


float MC_step(mat spin_matrix,int L, float J){
  float T = 1.2; float k = 1;
  float beta = k*T;

  float r = random_between_zero_and_one();
  int i = random_index(0,L-1);
  int j = random_index(0,L-1);

  float E_j = local_energy(spin_matrix,i,j,L);
  spin_matrix(i,j) = minus_one_or_one();

  float deltaE = local_energy(spin_matrix,i,j,L)-E_j;
  cout << deltaE << endl;


  //if(r<=exp(-beta*deltaE)){
  //  E_i = E_i+deltaE;
  //}
  return 1.0;
  //return deltaE;
}

/*
void MC_solve(mat spin_matrix,int L,float J){
  float E = 0;
  float M = 0;
  float deltaE;
  for(int i=0;i<mc_cycles;i++){
    deltaE = MC_step(spin_matrix,L,J);
    E+=deltaE;
  }
}
*/



int main()
{
  srand(clock());
  //srand(1);

  int L  = 100;
  int J = 1;
  mat spin_matrix = initialize_matrix_random(L);

  for(int i=0;i<1000;i++){
    MC_step(spin_matrix,L,J);
    //cout << minus_one_or_one() << endl;
  //

  //cout << atom_energy(spin_matrix,-1,2,L) << endl;
  //spin_matrix.print();
  //cout << system_energy(spin_matrix,L,J) << endl;
}


  //cout << system_energy(spin_matrix,L,J) << endl;
}
