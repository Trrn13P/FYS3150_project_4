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


float atom_energy(mat spin_matrix, int i, int j){
  //not returning -J, so remember to scale back
  return spin_matrix(i,j)*(spin_matrix(i,j-1)+spin_matrix(i,j+1)+spin_matrix(i-1,j)+spin_matrix(i+1,j));
}

float local_energy_change(mat spin_matrix,int i, int j){
  
}

float system_energy(mat spin_matrix,int L,float J){
  /*
  Denne maa endres to only look at close neighbors, så funksjonen må ta inn
  argument i,j. Må ha gammel og ny energi i disse punktene
  */
  float Energy = 0;
  for(int i=1;i<L+1;i++){
    for(int j=1;j<L+1;j++){
      Energy+=-J*spin_matrix(i,j)*(spin_matrix(i,j-1)+spin_matrix(i,j+1)+spin_matrix(i-1,j)+spin_matrix(i+1,j));
    }
  }
  return Energy;
}

float random_between_zero_and_one(){
  return rand()*1./RAND_MAX;
}

int minus_one_or_one(){
  float number = random_between_zero_and_one();
  if(number<0.5){
    return -1;
  }
  else{
    return 1;
  }
}

int random_index(int start, int finish){
  return rand() % (finish+1-start) + start;

}


mat initialize_matrix(int L){
  mat spin_matrix = zeros(L+2,L+2);
  float number;

  for(int i=1;i<L+1;i++){
    for(int j=1;j<L+1;j++){
        spin_matrix(i,j)= minus_one_or_one();;
      }
    }

    spin_matrix.row(0)=spin_matrix.row(L);
    spin_matrix.row(L+1)=spin_matrix.row(1);

    spin_matrix.col(0)=spin_matrix.col(L);
    spin_matrix.col(L+1)=spin_matrix.col(1);

  return spin_matrix;
}

mat flip_random_spin(mat spin_matrix,int L){
  int i = random_index(1,L+1);
  int j = random_index(1,L+1);
  spin_matrix(i,j) = minus_one_or_one();
  return spin_matrix;
}

float MC_step(mat spin_matrix,int L, float J){
  float T = 1.2; float k = 1;
  float beta = k*T;
  float r = random_between_zero_and_one();

  float E_j = system_energy(spin_matrix,L,J);

  spin_matrix = flip_random_spin(spin_matrix,L);
  E_i = system_energy(spin_matrix,L,J);

  float deltaE = E_i-E_j;

  if(r<=exp(-beta*deltaE)){
    E_i = E_i+deltaE;
  }
  return deltaE;
}

void MC_solve(mat spin_matrix,int L,float J){
  float E = 0;
  float M = 0;
  float deltaE;
  for(int i=0;i<mc_cycles;i++){
    deltaE = MC_step(spin_matrix,L,J);
    E+=deltaE;
  }
}




int main()
{
  srand(clock());

  int L  = 2;
  mat spin_matrix = initialize_matrix(L);

  int magnetization = 0;
  int J = 1;

  for(int i=1;i<L+1;i++){
    for(int j=1;j<L+1;j++){
      magnetization += spin_matrix(i,j);
    }
  }
  cout << system_energy(spin_matrix,L,J) << endl;
}
