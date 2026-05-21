#include "Bilayer.h"
#include "Constants.h" // frequently and widely used variables used as parameters
#include <array>
#include "Utilities.h"

extern Rnd rnd;
//////////////////////////// local function declarations ///////////////////////

static inline std::array<double, 2> prepareTop(
  const double topValue, const double grapheneValue);

static inline std::array<double, 2> prepareBot(
  const double grapheneValue, const double bottomValue);

static inline double prepareSumAtan_EnlTimesInvLorenzianPar(
  const std::array<double, Const::Nl>& landauEnergies);

static inline size_t getLevelIndex(const int level);

static inline std::array<double, Const::Nl> prepareLandauLevelEnergies(
  const double Bz);

static inline std::array<double, Const::Nl>
  preparelandauLevelsEnergiesTimesInvLorenzianPar(
    std::array<double, Const::Nl> LandauLevelEnergies);

static inline double countN(
  const double n0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double Vg);

static inline double landauLevelEnergy(const int level, const double Bz);

static inline double count4eBz_hPi(const double Bz);

static double countE0(
  const double Bz,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergies,
  const double sumAtan_EnlTimesInvLorenzianPar);

static double countVg(
  const double Bz,
  const double E0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar,
  const double sumAtan_EnlTimesInvLorenzianPar);

///////////////////////////// method declarations //////////////////////////////

/* returns {nt, nb, Vgt, Vgb, E0t, E0b} -> densities calculated with external magnetic field,
    provide Vt and Vb in eV and B in T */
// TODO compare equations with original paper not tdmac5536supp1.pdf
Bilayer::resultsB Bilayer::countDensitiesAndPotential(
  double Vt, double Vb, double B) const
{
  const double Vt_au = Vt * Const::V2au;
  const double Vb_au = Vb * Const::V2au;
  const double B_au = B * Const::T2au;

  // initial values and conversion to au
  double Vgt = clamp((Vb_au + 2./3. * (Vt_au - Vb_au)) * 0.01,
                     Const::VgMin,
                     Const::VgMax);
  double Vgb = clamp((Vb_au + 1./3. * (Vt_au - Vb_au)) * 0.01,
                     Const::VgMin,
                     Const::VgMax);

  // precalculate some values
  const std::array<double, Const::Nl> landauLevelEnergies = prepareLandauLevelEnergies(B_au);
  const std::array<double, Const::Nl> landauLevelsEnergiesTimesInvLorenzianPar =
    preparelandauLevelsEnergiesTimesInvLorenzianPar(landauLevelEnergies);
  const double sumAtan = prepareSumAtan_EnlTimesInvLorenzianPar(
    landauLevelsEnergiesTimesInvLorenzianPar);

  const double E0t = countE0(B_au, _nit, landauLevelEnergies, sumAtan);
  const double E0b = countE0(B_au, _nib, landauLevelEnergies, sumAtan);

  // ASK jeśli obliczamy tutaj E0 to czy ma sens później obliczać dla innego E0 w fortranie?
  // Jeśli nie to co zrobić?

  bool converged = false;
  uint8_t it = 0;

  do
  {
    it++;
    auto[newVgt, newVgb, newConverged] =
      iterationB(
        Vt_au, Vb_au, Vgt, Vgb, B_au, E0t, E0b, landauLevelsEnergiesTimesInvLorenzianPar, sumAtan);

    Vgt = Const::VgAlpha * newVgt + (1. - Const::VgAlpha) * Vgt;
    Vgb = Const::VgAlpha * newVgb + (1. - Const::VgAlpha) * Vgb;
    converged = newConverged;

  } while (!converged && it < Const::maxIterationsB);

  assertVerbose(converged, "IterationB didn't converge, Vgt=" << Vgt << " Vgb=" << Vgb << '\n');

  // count densities for final voltages
  const std::array<double, 2> capacitiesT = prepareTop(_Ct, _Cg);
  const std::array<double, 2> capacitiesB = prepareBot(_Cg, _Cb);
  const std::array<double, 2> voltagesT = prepareTop(Vt_au, Vgb);
  const std::array<double, 2> voltagesB = prepareBot(Vgt, Vb_au);

  const double nt = countN(_nit, capacitiesT, voltagesT, Vgt);
  const double nb = countN(_nib, capacitiesB, voltagesB, Vgb);

  // return results and log them

  resultsB results = {
    .nt = nt / Const::inv_cmsq2au,
    .nb = nb / Const::inv_cmsq2au,
    .Vgt = Vgt / Const::V2au,
    .Vgb = Vgb / Const::V2au,
    .E0t = E0t / Const::eV2au,
    .E0b = E0b / Const::eV2au
  };

  return results;
};

// returns Vgt, Vgb, converged
std::tuple<double, double, bool> Bilayer::iterationB(
  const double Vt,
  const double Vb,
  const double prevVgt,
  const double prevVgb,
  const double Bz,
  const double E0t,
  const double E0b,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar,
  const double sumAtan_EnlTimesInvLorenzianPar) const
{
  const std::array<double, 2> capacitiesT = prepareTop(_Ct, _Cg);
  const std::array<double, 2> capacitiesB = prepareBot(_Cg, _Cb);

  std::array<double, 2> voltagesT = prepareTop(Vt, prevVgb);
  std::array<double, 2> voltagesB = prepareBot(prevVgt, Vb);

  double Vgt = 0;
  double Vgb = 0;

  // if (r < 0.333) // both from prev // ASK czy to ma sens to mieszanie? chyba nie wpływa na wyniki bardzo
  // {

  //   Vgt =
  //   countVg(
  //     Bz,
  //     E0t,
  //     capacitiesT,
  //     voltagesT,
  //     _nit,
  //     landauLevelEnergies,
  //     sumAtan_EnlTimesInvLorenzianPar);
  //   Vgb =
  //     countVg(
  //       Bz,
  //       E0b,
  //       capacitiesB,
  //       voltagesB,
  //       _nib,
  //       landauLevelEnergies,
  //       sumAtan_EnlTimesInvLorenzianPar);
  // }
  // else if (r < 0.666) // count Vgt first and take it to calculation of Vgb
  // {
  Vgt =
    countVg(
      Bz,
      E0t,
      capacitiesT,
      voltagesT,
      _nit,
      landauLevelsEnergiesTimesInvLorenzianPar,
      sumAtan_EnlTimesInvLorenzianPar);

  voltagesB[0] = Vgt;
  Vgb =
    countVg(
      Bz,
      E0b,
      capacitiesB,
      voltagesB,
      _nib,
      landauLevelsEnergiesTimesInvLorenzianPar,
      sumAtan_EnlTimesInvLorenzianPar);
  // }
  // else // count Vgb first and take it to calculation of Vgt
  // {

  //   Vgb =
  //     countVg(
  //       Bz,
  //       E0b,
  //       capacitiesB,
  //       voltagesB,
  //       _nib,
  //       landauLevelEnergies,
  //       sumAtan_EnlTimesInvLorenzianPar);

  //   voltagesT[1] = Vgb;
  //   Vgt =
  //     countVg(
  //       Bz,
  //       E0t,
  //       capacitiesT,
  //       voltagesT,
  //       _nit,
  //       landauLevelEnergies,
  //       sumAtan_EnlTimesInvLorenzianPar);
  // }

  const double absoluteErrorT = std::abs(Vgt - prevVgt);
  const double absoluteErrorB = std::abs(Vgb - prevVgb);

  const bool convergedT =
    absoluteErrorT < Const::VgAbsTol + Const::VgRelTol * std::abs(Vgt);
  const bool convergedB =
    absoluteErrorB < Const::VgAbsTol + Const::VgRelTol * std::abs(Vgb);

  const bool converged = convergedT && convergedB;

  return {Vgt, Vgb, converged};
}

//////////////////////////// local function definitions ////////////////////////

// above top layer, below top layer
static inline std::array<double, 2> prepareTop(
  const double topValue, const double grapheneValue)
{
  return {topValue, grapheneValue};
}

// above bottom layer, below bottom layer
static inline std::array<double, 2> prepareBot(
  const double grapheneValue, const double bottomValue)
{
  return {grapheneValue, bottomValue};
}

static inline double prepareSumAtan_EnlTimesInvLorenzianPar(
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar)
{
  double sum = 0;

  for (double EnlTimesInvLorenzianPar: landauLevelsEnergiesTimesInvLorenzianPar)
  {
    sum += std::atan(EnlTimesInvLorenzianPar);
  }

  return sum;
}

static inline size_t getLevelIndex(const int level)
{
  assert(level >= -Const::L_max && level <= Const::L_max);
  return Const::L_max + level;
}

static inline std::array<double, Const::Nl> prepareLandauLevelEnergies(
  const double Bz)
{
  std::array<double, Const::Nl> landauLevelEnergies{0};

  for (int level = -Const::L_max; level <= Const::L_max; level++)
  {
    landauLevelEnergies[getLevelIndex(level)] = landauLevelEnergy(level, Bz);
  }

  return landauLevelEnergies;
}

static inline std::array<double, Const::Nl>
  preparelandauLevelsEnergiesTimesInvLorenzianPar(
    std::array<double, Const::Nl> landauLevelEnergies)
{
  std::array<double, Const::Nl> landauLevelsEnergiesTimesInvLorenzianPar{0};

  for (size_t level = 0; level < landauLevelEnergies.size(); level++)
  {
    landauLevelsEnergiesTimesInvLorenzianPar[level] =
      landauLevelEnergies[level] * Const::inv_Lorentzian_par;
  }

  return landauLevelsEnergiesTimesInvLorenzianPar;
}

static inline double countN(
  const double n0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double Vg)
{
  const double nG = (capacities[0] * (voltages[0] - Vg) +
                     capacities[1] * (voltages[1] - Vg)) / Const::e;

  return n0 + nG;
}

// level = nl in paper
static inline double landauLevelEnergy(const int level, const double Bz)
{
  return signum(level) *
    std::sqrt(
      2 * Const::e * Bz * Const::hbar * Const::v * Const::v * std::abs(level));
}

static inline double count4eBz_hPi(const double Bz)
{
  return 4. * Const::e * Bz / (Const::h * M_PI);
}

static double countE0(
  const double Bz,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergies,
  const double sumAtan_EnlTimesInvLorenzianPar)
{
  const double expr = count4eBz_hPi(Bz);

  auto func =
  [n0, Bz, expr, &landauLevelsEnergies, sumAtan_EnlTimesInvLorenzianPar]
  (const double E0)
  {
    double sum = 0;
    for (double Enl: landauLevelsEnergies)
    {
      sum += std::atan((E0 - Enl) * Const::inv_Lorentzian_par);
    }

    sum += sumAtan_EnlTimesInvLorenzianPar;

    return expr * sum - n0;
  };

  return bisection(Const::E0Min, Const::E0Max, func);
}

static double countVg(
  const double Bz,
  const double E0,
  const std::array<double, 2>& capacities,
  const std::array<double, 2>& voltages,
  const double n0,
  const std::array<double, Const::Nl>& landauLevelsEnergiesTimesInvLorenzianPar,
  const double sumAtan_EnlTimesInvLorenzianPar)
{
  const double expr = count4eBz_hPi(Bz);

  auto func =
  [n0, &capacities, &voltages, expr, Bz, E0,
    &landauLevelsEnergiesTimesInvLorenzianPar,
    sumAtan_EnlTimesInvLorenzianPar]
  (double Vg)
  {
    const double left = countN(n0, capacities, voltages, Vg);
    const double shiftedE = (E0 + Const::e * Vg) * Const::inv_Lorentzian_par;

    double right = 0;
    for (double EnlTimesInvLorenzianPar: landauLevelsEnergiesTimesInvLorenzianPar)
    {
      right += std::atan(shiftedE - EnlTimesInvLorenzianPar);
    }

    right += sumAtan_EnlTimesInvLorenzianPar;

    return expr*right - left;
  };

  return bisection(Const::VgMin, Const::VgMax, func);
}



