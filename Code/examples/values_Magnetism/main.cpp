#include <iostream>
#include <vector>
#include "Bilayer.h" // functions just for this project
#include "Utilities.h" // my functions for std::vector and other functionalities
#include "Constants.h" // frequently and widely used variables used as parameters

int main()
{
  Bilayer bilayer(3 * 1e11, 3 * 1e11);

  const double dVb = 0.5;
  const double dB = 4.59999993 - 4.49999994;

  //std::vector<T> linspace(const T start, const T end, const T step)
  double minB = 0.49999999;
  double maxB = 7.99999989 + 0.5* dB;
  const std::vector<double> BTab = linspace(minB, maxB , dB);
  const std::vector<double> VbTab = linspace(-30., 20., dVb);

  FILE* file = fopen("./out.csv", "w");

  assertVerbose(file, "file");

  for (double B: BTab)
  {
    progressBar(B, minB, maxB);
    for (double Vb: VbTab)
    {
      Bilayer::resultsB results = bilayer.countDensitiesAndPotential(0., Vb, B);
      //           Vb,    B,    nt,   nb,   o_t,  o_b,  E0t, E0b
      std::fprintf(file,
                   "%04.6f,\t %04.6f,\t %04.6f,\t %04.6f,\t %04.6f,\t %04.6f,\t %04.6f,\t %04.6f\n",
                   Vb,
                   B,
                   results.nt / 1e11,
                   results.nb / 1e11,
                   -results.Vgt - results.E0t,
                   -results.Vgb - results.E0b,
                   results.E0t,
                   results.E0b);
    }
    std::fprintf(file, "\n");
  }
  fclose(file);
}
