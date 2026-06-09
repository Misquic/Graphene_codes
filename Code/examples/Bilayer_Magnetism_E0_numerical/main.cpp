#include <iostream>
#include <vector>
#include "Bilayer.h" // functions just for this project
#include "Utilities.h" // my functions for std::vector and other functionalities
#include "Constants.h" // frequently and widely used variables used as parameters
#include <string.h>

int main(int argc, char* argv[]){

  std::cout << "Program calculates E0 numerically in Graphene in Magnetic field (eq.14)\n";

  size_t N = 100, M = 100;
  if (argc >= 2)
  {
    N = std::atoi(argv[1]);
  }
  char resultsFile[512] = "results/n_map.csv";
  if (argc >= 3)
  {
    strncpy(resultsFile, argv[2], 511);
    resultsFile[511] = '\0';
  }
  double Barg = 4;
  if (argc >= 4)
  {
    Barg = std::atof(argv[3]);
    M = 1;
  }

  std::cout << resultsFile << '\n';
  std::cout << Barg << '\n';

  double n0_min = -8e11;
  double n0_max =  8e11;
  double delta_n0 = (n0_max - n0_min)/N;

  double B_min = 0.1;
  double B_max = 2;
  if (argc >= 4)
  {
    B_min = Barg;
    B_max = Barg;
  }
  double delta_B = (B_max - B_min)/M;

  // Array with parameters and results [density, Magnetic field, base energy]
  Array2D<double> n(6,M*N);

  for (size_t i = 0; i < N; i++)
  {
    double n0 = n0_min + i*delta_n0;

    progressBar(n0, n0_min, n0_max);

    for(size_t j = 0; j < M; j++)
    {
      double B = B_min + j*delta_B;

      Bilayer bilayer(n0, n0);
      Bilayer::resultsB results = bilayer.countDensitiesAndPotential(-20., 0., B);

      n(0,i*M+j) = results.nt;
      n(1,i*M+j) = B;
      n(2,i*M+j) = results.E0t;
      n(3,i*M+j) = n0;
      n(4,i*M+j) = results.nt - n0;
      n(5,i*M+j) = results.Vgt;
    }
  }

  progressBar(n0_max, n0_min, n0_max);

  save(n, resultsFile);

  std::cout << "end\n";
  return 0;
}
