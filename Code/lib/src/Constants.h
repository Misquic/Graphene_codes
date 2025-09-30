#ifndef CONSTANTS_H
#define CONSTANTS_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <limits>

// inline constexpr double m_eff2 = 0.067;

double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
{
  return curr == prev
    ? curr
    : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);    
}

double constexpr const_sqrt(double x)
{
  return x >= 0 && x < std::numeric_limits<double>::infinity()
    ? sqrtNewtonRaphson(x,x,0)
    : std::numeric_limits<double>::quiet_NaN();
};

class Const
{
public:
  inline static constexpr double nm2au = 1/0.0529; // nm*nm2au = au
  inline static constexpr double m_eff = 0.067; 
  inline static constexpr double E_h = 27.211; //eV
  inline static constexpr double eV2au = 1./E_h; // eV*eV2au = au
  inline static constexpr double eps_0 = 1./(4.*M_PI); // au (e^2 * (eV->au)^-1 * (m->au)) 
  inline static constexpr double cm2au = 1e-2*1e9*nm2au; 
  inline static constexpr double inv_cmsq2au = 1./cm2au/cm2au;
  inline static constexpr double t = 3.*eV2au;
  inline static constexpr double a = 1./(4.*const_sqrt(3.))*nm2au;
  inline static constexpr double hv = 3./2.*t*a;
  inline static constexpr double T2au = 4.254382E-6;
  inline static constexpr double h = 1*2*M_PI; // normal planc constant

  // parameters constants only for now it will change
  inline static constexpr double Lorentzian_par = 0.4e-3*eV2au;  // witdh in approx for delta function
  inline static constexpr double inv_Lorentzian_par = 1./Lorentzian_par;  

  static int L_max;

  Const(const Const& other) = delete;
  Const(Const&& other) = delete;
  Const operator=(const Const& other) = delete;
  Const operator=(Const&& other) = delete;

};

#endif