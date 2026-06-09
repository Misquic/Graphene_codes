#include <iostream>
#include <vector>
#include "Bilayer.h" // functions just for this project
#include "Utilities.h" // my functions for std::vector and other functionalities
#include "Constants.h" // frequently and widely used variables used as parameters

int main(int argc, char* argv[]){

  std::cout << "Program calculates Voltages numerically in Graphene in Magnetic field (eq.16)\n";

  double Vt = 0;
  if(argc >= 2)
  {
    Vt = std::atof(argv[1]);
  }

  Bilayer bilayer;

  const double BMin = 5.f;
  const double BMax = 8.f;
  const std::vector<double> BTab = linspace<double>(BMin, BMax, .05f);

  const double VbMin = -20.f;
  const double VbMax = 20.f;
  const std::vector<double> VbTab = linspace<double>(VbMin, VbMax, .125f);

  save(BTab, "./results/B.csv");
  save(VbTab, "./results/Vb.csv");

  Array2D<double> resultsVgt(BTab.size(), VbTab.size());
  Array2D<double> resultsVgb(BTab.size(), VbTab.size());

  std::cout << '\n';

  size_t VbIdx = 0;
  for (const double Vb: VbTab)
  {
    progressBar(Vb, VbMin, VbMax);
    size_t BIdx = 0;
    for (const double B: BTab)
    {
      Bilayer::resultsB res = {};

      res = bilayer.countDensitiesAndPotential(Vt, Vb, B);

      resultsVgt(BIdx, VbIdx) = res.E0t;
      resultsVgb(BIdx, VbIdx) = res.E0b;

      BIdx++;
    }
    VbIdx++;
  }

  save(resultsVgt, "./results/E0t.csv");
  save(resultsVgb, "./results/E0b.csv");

  // print_hist(resultsVgb);
  // print_hist(resultsVgt);
}
