#include <iostream>
#include <vector>
#include "Bilayer.h" // functions just for this project
#include "Utilities.h" // my functions for std::vector and other functionalities
#include "Constants.h" // freqquently and widely used variables used as parameters
        
int main(int argc, char* argv[]){
  std::cout << "Program calculates densities for bilayer Graphene (eq.3 & eq.5)\n"
               "use: ./main dG show N\n"
               "\t     dG: thickness of graphene layer, possible values:\n"
               "\t\t      0: counting from 0.1 to 0.9\n"
               "\t\t<value>: counting for only <value>\n"
               "\t   show: printing plot in terminal. For huge N use 'ctrl'+'-' to view whole image\n"
               "\t      N: number of pixels/points for which to count densities\n\n";

  double dt = 21.;
  double dGmin = 0.10;
  double dGmax = 0.90;
  double db = 53.;

  double et = 3.3;
  double eG = 1.;
  double eb = 3.3;

  size_t N = 150;
  if (argc >= 4)
  {
    N = std::atoi(argv[3]);
  }
  printf("N = %lu\n", N);
  Array2D<double> nt(N, N);
  Array2D<double> nb(N, N);

  double dG = dGmin;
  if (argc >= 2)
  {
    dG = std::atof(argv[1]);
    if(dG == 0)
    {
      dG = dGmin;
    }
    else{
      dGmax = dG;
    }
    printf("dG = %f\n", dG);
  }
  else{
    printf("dG = %f\n", 0.);
  }

  bool show = false;
  if (argc >= 3)
  {
    show = bool(std::atoi(argv[2]));
  }
  printf("show = %s\n", show?"true":"false");

  for(; dG <= dGmax; dG += (dGmax-dGmin)/50)
  {
    Bilayer system(0, 0,
                   dt, dG, db,
                   et, eG, eb);
    
    double Vt_min = -6, Vt_max = 6;
    double delta_Vt = (Vt_max - Vt_min)/(N+1);
    
    double Vb_min = -10, Vb_max = 10;
    double delta_Vb = (Vb_max - Vb_min)/(N+1);
    
    for (size_t i = 0; i < N; i++)
    {
      double Vb = Vb_max - delta_Vb*i;
      for (size_t j = 0; j < N; j++)
      {
        double Vt = Vt_min + delta_Vt*j;
        auto[nt_temp, nb_temp] = system.countDensitiesAndPotential(Vt, Vb);
        nt(i, j) = nt_temp/Const::inv_cmsq2au; 
        nb(i, j) = nb_temp/Const::inv_cmsq2au;
      }
    }
    
    std::cout << str(int(dG*100)) << "\n";
    save(nt, "results/bilayer_nt/dG_0." + str(int(dG*100)) + ".csv");
    save(nb, "results/bilayer_nb/dG_0." + str(int(dG*100)) + ".csv");
    
    if(show)
    {
      // printf("\e[0;0H");
      nb.apply([](double a){
        if(std::abs(a) < 1e10){
          return 1.;
        }
        return 0.;
      });
      
      nt.apply([](double a){
        if(std::abs(a) < 1e10){
          return 1.;
        }
        return 0.;
      });
      print_hist(nt+nb);
    }
  }
  printf("\nend\n");
}
