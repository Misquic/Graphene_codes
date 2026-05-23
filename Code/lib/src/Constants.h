#ifndef CONSTANTS_H
#define CONSTANTS_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <limits>
#include <cstdint>

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
  inline static constexpr double e     = 1;                // electron charge
  inline static constexpr double eps_0 = 1./(4.*M_PI);     // au (e^2 * (eV->au)^-1 * (m->au))
  inline static constexpr double hbar  = 1;                // reduced planc constant
  inline static constexpr double h     = hbar * 2 * M_PI;  // normal planc constant

  inline static constexpr double nm2au       = 1 / 0.0529;          // nm*nm2au = au
  inline static constexpr double m_eff       = 0.067;
  inline static constexpr double E_h         = 27.211;              //eV
  inline static constexpr double eV2au       = 1. / E_h;            // eV*eV2au = au
  inline static constexpr double V2au        = eV2au / Const::e;    // V*V2au = au
  inline static constexpr double cm2au       = 1e-2 * 1e9 * nm2au;
  inline static constexpr double inv_cmsq2au = 1. / cm2au / cm2au;
  inline static constexpr double T2au        = 4.254382E-6;

  inline static constexpr double t           = 3. * eV2au;
  inline static double a                     = 1. / (4. * std::sqrt(3.)) * nm2au;
  inline static double hbar_v                = 3. / 2. * t * a;
  inline static double v                     = hbar_v / hbar;

  // parameters constants only for now it will change
  inline static constexpr double Lorentzian_par     = 1e-3 * eV2au;         // witdh in approx for delta function
  inline static constexpr double inv_Lorentzian_par = 1. / Lorentzian_par;
  inline static constexpr double VgAbsTol           = 1e-6 * V2au;          // Tolerance of absolute error for Vg
  inline static constexpr double VgRelTol           = 1e-6;                 // Tolerance of relative error for Vg
  inline static constexpr double smallNumber        = 1e-14 * V2au;         // Don't devide by zero
  inline static constexpr double VgAlpha            = .3;                   // relaxation, 1 -> only new Value 0 -> only previous value;
  inline static constexpr uint8_t maxIterationsB    = 30 / VgAlpha;         // max iterations of Bilayer::iterationB

  //default values fo Bilayer
  inline static constexpr double nit_default = 8.12/2 * 1e11 * inv_cmsq2au;
  inline static constexpr double nib_default = 8.12/2 * 1e11 * inv_cmsq2au;
  inline static constexpr double dt_default  = 330.0 * nm2au; // SiO2
  inline static constexpr double dg_default  = 0.5 * nm2au;
  inline static constexpr double db_default  = 40.0 * nm2au; // hBN
  inline static constexpr double et_default  = 3.9; // SiO2
  inline static constexpr double eg_default  = 1.0;
  inline static constexpr double eb_default  = 3.7; // hBN
  inline static constexpr double Ct_default  = eps_0 * et_default / dt_default;
  inline static constexpr double Cg_default  = eps_0 * eg_default / dg_default;
  inline static constexpr double Cb_default  = eps_0 * eb_default / db_default;
  inline static constexpr int L_max = 500;
  inline static constexpr size_t Nl = 2 * static_cast<size_t>(L_max) + 1;

  inline static double E0Min = -30 * eV2au;
  inline static double E0Max = 30 * eV2au;
  inline static double VgMin = -20 * V2au;
  inline static double VgMax = 20 * V2au;


  Const(const Const& other) = delete;
  Const(Const&& other) = delete;
  Const operator=(const Const& other) = delete;
  Const operator=(Const&& other) = delete;

};

#endif