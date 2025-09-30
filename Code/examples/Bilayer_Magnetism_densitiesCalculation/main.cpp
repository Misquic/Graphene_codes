#include <iostream>
#include <vector>
#include "funkcje.h" // functions just for this project
#include "Utilities.h" // my functions for std::vector and other functionalities
#include "Constants.h" // freqquently and widely used variables used as parameters
        
int main(int argc, char* argv[]){
    
  dmsg("debug\n");

  return 0;

  size_t M = 10;
  if(argc >= 2)
  {
    M = std::atoi(argv[1]);
  }
  Array2D<double> n(3,M*M);

  double n0_min = -8e11*Const::inv_cmsq2au;
  double n0_max =  8e11*Const::inv_cmsq2au;
  double delta_n0 = (n0_max - n0_min)/M;
    

  double B_min = 0.1*Const::T2au;
  double B_max = 2*Const::T2au;
  double delta_B = (B_max - B_min)/M;

  std::vector<double> E_nl_vec(std::abs(Const::L_max));
  for(size_t i = 0; i < M; i++)
  {
    for(size_t j = 0; j < M; j++)
    {
      double n0 = n0_min + i*delta_n0; 
      double B = B_min + j*delta_B; 

      if(n0 > 0)
      {
        E_nl_vec = prepare_E_nl(B, std::abs(Const::L_max));
      }
      else
      {
        E_nl_vec = prepare_E_nl(B, -std::abs(Const::L_max));
      }

      n(0,i*M+j) = n0/Const::inv_cmsq2au;
      n(1,i*M+j) = B/Const::T2au;
      n(2,i*M+j) = count_E_0(E_nl_vec, n0, B, prepare_part_sum(E_nl_vec))/Const::eV2au;
    }
  }
  
  save(n, "results/n_map.csv");

  return 0;
}
