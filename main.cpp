#include <iostream>
#include <cmath>
#include <fstream>

#include <omp.h>
#include <time.h>

#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;


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


float MC_step(mat spin_matrix,int L, float J, float exp_precalc[17]){


  float r = random_between_zero_and_one();
  int i = random_index(0,L-1);
  int j = random_index(0,L-1);

  float E_j = local_energy(spin_matrix,i,j,L);
  spin_matrix(i,j) = minus_one_or_one();

  float deltaE = local_energy(spin_matrix,i,j,L)-E_j;

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


void MC_solve(mat spin_matrix,int L,float J,int mc_cycles){
  float T = 1.2; float k = 1;
  float beta = k*T;

  float exp_precalc[17];
  float deltaE_precalc[17]={-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16};
  for (int i=0;i<17;i++){
    exp_precalc[i]=exp(-beta*deltaE_precalc[i]);
  }


  float deltaE;
  vec energies = zeros(mc_cycles);
  energies(0) = system_energy(spin_matrix,L,J);

  for(int i=1;i<mc_cycles;i++){
    deltaE = MC_step(spin_matrix,L,J,exp_precalc);
    energies(i)=energies(i-1)+ deltaE;
    cout << system_energy(spin_matrix,L,J) << endl;
  }
  ofstream outfile("test2.txt");
  outfile << energies << endl;
  outfile.close();
}




int main()
{
  srand(clock());

  int L  = 5;
  int J = 1;
  mat spin_matrix = initialize_matrix_ordered(L);


  int mc_cycles = 1E4;
  MC_solve(spin_matrix,L,J,mc_cycles);
  for(int i=0;i<1;i++){
    //cout << i << endl;

    //MC_step(spin_matrix,L,J);
    //cout << MC_step(spin_matrix,L,J) << endl;

    //cout << minus_one_or_one() << endl;
  //

  //cout << atom_energy(spin_matrix,-1,2,L) << endl;
  //spin_matrix.print();
  //cout << system_energy(spin_matrix,L,J) << endl;
}



  //cout << system_energy(spin_matrix,L,J) << endl;
}
