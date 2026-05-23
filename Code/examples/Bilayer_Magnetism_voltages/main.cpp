#include <iostream>
#include <vector>
#include "Bilayer.h" // functions just for this project
#include "Utilities.h" // my functions for std::vector and other functionalities
#include "Constants.h" // frequently and widely used variables used as parameters

int main(int argc, char* argv[]){

  std::cout << "Program calculates Voltages numerically in Graphene in Magnetic field (eq.16)\n";

  double Vt = 3;
  if(argc >= 2)
  {
    Vt = std::atof(argv[1]);
  }

  Bilayer bilayer;

  const double BMin = 1.f;
  const double BMax = 8.f;
  const std::vector<double> BTab = linspace<double>(BMin, BMax, .05f);

  const double VbMin = -60.f;
  const double VbMax = 30.f;
  const std::vector<double> VbTab = linspace<double>(VbMin, VbMax, .5f);

  save(BTab, "./Vt_" + str(Vt) + "/B.csv");
  save(VbTab, "./Vt_" + str(Vt) + "/Vb.csv");

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

      try
      {
        res = bilayer.countDensitiesAndPotential(Vt, Vb, B);
      }
      catch(const std::runtime_error& e)
      {
        std::cerr << e.what() << '\n';
        res.Vgb = -0.2;
        res.Vgt = -0.2;
      }

      resultsVgt(BIdx, VbIdx) = res.Vgt;
      resultsVgb(BIdx, VbIdx) = res.Vgb;

      BIdx++;
    }
    VbIdx++;
  }

  save(resultsVgt, "./Vt_" + str(Vt) + "/Vgt.csv");
  save(resultsVgb, "./Vt_" + str(Vt) + "/Vgb.csv");

  print_hist(resultsVgb);
  print_hist(resultsVgt);
}
