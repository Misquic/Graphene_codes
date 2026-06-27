#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <string.h>
#include <charconv>

#include "resistances.h"

////////////////////////////////// Static //////////////////////////////////////

std::ifstream openIFile(const std::string& path, bool& valid)
{
  dmsg("Reading file " << path << '\n');

  std::ifstream file(path);
  if (!file.is_open())
  {
    std::cerr << "Couldn't open file " << path << " terminating\n";
    valid = false;
  };

  return file;
}

std::ofstream openOFile(const std::string& path, bool& valid)
{
  dmsg("Opening file " << path << '\n');

  std::ofstream file(path);
  if (!file.is_open())
  {
    std::cerr << "Couldn't open file " << path << " terminating\n";
    valid = false;
  };

  return file;
}

// read while **current_pp is not c, stop before c
static inline bool readUntil(char c, char** const current_pp, size_t& remainingSize)
{
  while ((**current_pp != c) && (remainingSize > 0))
  {
    (*current_pp)++;
    remainingSize--;
  }
  return (**current_pp) == c;
}

// read while **current_pp is c, stop at first after not c
static inline bool readWhile(char c, char** const current_pp, size_t& remainingSize)
{
  while ((**current_pp == c) && (remainingSize > 0))
  {
    (*current_pp)++;
    remainingSize--;
  }
  return true;
}

static bool readNumber(int& number, char** const current_pp, size_t& remainingSize, const bool last)
{
  readWhile(' ', current_pp, remainingSize); // trim spaces
  char* start = *current_pp;
  char* stop = NULL;
  if (!last)
  {
    readUntil(' ', current_pp, remainingSize);
    stop = (*current_pp) + 1; // set stop to just after number
    // readWhile(' ', current_pp, remainingSize); // go to next
  }
  else
  {
    readUntil('\n', current_pp, remainingSize);
    stop = (*current_pp) + 1; // set stop to just after number
  }

  const auto res = std::from_chars(start, stop, number, 10);
  if ((int)res.ec != 0) return false;
  return true;
}

static bool readNumber(float& number, char** const current_pp, size_t& remainingSize, const bool last)
{
  readWhile(' ', current_pp, remainingSize); // trim spaces
  char* start = *current_pp;

  char* stop = NULL;
  if (!last)
  {
    readUntil(' ', current_pp, remainingSize);
    stop = (*current_pp) + 1; // set stop to just after number
    // readWhile(' ', current_pp, remainingSize); // go to next
  }
  else
  {
    readUntil('\n', current_pp, remainingSize);
    stop = (*current_pp) + 1; // set stop to just after number
  }

  const auto res = std::from_chars(start, stop, number, std::chars_format::scientific);
  if ((int)res.ec != 0)
  {
    dmsg((int)res.ec);
    return false;
  }
  return true;
}

// c -> what char is past desired stop
static inline bool readPast(char c, char** const current_pp, size_t& remainingSize)
{
  if (!readUntil(c, current_pp, remainingSize)) return false;

  if (remainingSize > 0)
  {
    remainingSize--;
    (*current_pp)++;
    return true;
  }
  return false;
}

static inline bool convertAndCheckNumber(double& number, const char* const start_p, const char* const end_p)
{
  // std::printf("start_p: %p, *start_p: %c, end_p: %p, *end_p: %c, diff: %d\n",
              // start_p, *start_p, end_p, *end_p, (int)(end_p - start_p));

  std::from_chars_result res = std::from_chars(start_p, end_p, number, std::chars_format::general);

  if ((int)res.ec != 0)
  {
    std::cerr << "Error in reading number, ec = " << (int)res.ec << " read: " << number << '\n';
    return false;
  }
  return true;

}

static inline bool convertAndCheckNumber(int& number, const char* const start_p, const char* const end_p)
{
  // std::printf("start_p: %p, *start_p: %c, end_p: %p, *end_p: %c, diff: %d\n",
              // start_p, *start_p, end_p, *end_p, (int)(end_p - start_p));

  std::from_chars_result res = std::from_chars(start_p, end_p, number);
  if ((int)res.ec != 0)
  {
    std::cerr << "Error in reading number, ec = " << (int)res.ec << " read: " << number << '\n';
    return false;
  }
  return true;

}

static void remapLeads(int& to, int& from, const leadInfoS& leadInfo, const leadInfoS& leadDesired)
{
  if (to == leadInfo.currentLeadFrom) to = leadDesired.currentLeadFrom;
  else if (to == leadInfo.currentLeadTo) to = leadDesired.currentLeadTo;
  else if (to == leadInfo.voltageLeadHigh) to = leadDesired.voltageLeadHigh;
  else to = leadDesired.voltageLeadLow;

  if (from == leadInfo.currentLeadFrom) from = leadDesired.currentLeadFrom;
  else if (from == leadInfo.currentLeadTo) from = leadDesired.currentLeadTo;
  else if (from == leadInfo.voltageLeadHigh) from = leadDesired.voltageLeadHigh;
  else from = leadDesired.voltageLeadLow;
};
////////////////////////////////// Global //////////////////////////////////////

size_t getNumOfSubDirs(const std::filesystem::path& dirs)
{
  size_t getNumOfSubDirs = 0;

  for (auto& entry: std::filesystem::directory_iterator(dirs))
  {
    if (entry.is_directory()) getNumOfSubDirs++;
  }

  return getNumOfSubDirs;
}

bool check_(int result, int line)
{
  if (result != true)
  {
    std::cerr << "Error " << result << " in line " << line << '\n';
  }
  return result;
};

double calculateR(double* const T)
{
  // T_to_from
  // G_12 =/= G_21 due to magnetic field
  double G[9] =
  {
    T[index4(0, 1)] + T[index4(0, 2)] + T[index4(0, 3)],
      -T[index4(0, 1)],
        -T[index4(0, 2)],
    -T[index4(1, 0)],
      T[index4(1, 0)] + T[index4(1, 2)] + T[index4(1, 3)],
        -T[index4(1, 2)],
    -T[index4(2, 0)],
      -T[index4(2, 1)],
        T[index4(2, 0)] + T[index4(2, 1)] + T[index4(2, 3)]
  };

  double R[9];
  inverseMatrix(G, R);
  // printMatD(R);

  // R_kl,mn = (V_m - V_n) / I_k = R_mk - R_nk
  // V_i = sum_{j=1, j!=i}^N R_ij * Ij
  //

  double R4t = R[index3(1, 0)] - R[index3(2, 0)];
  // std::printf("current leads: 1 & 4, voltage leads 2 & 3, 3 higher, R: %f\n", R4t);

  return R4t;
}

void printMat(const double* const mat)
{
  std::printf("[ % 5.3f, % 5.3f, % 5.3f ]\n"\
              "[ % 5.3f, % 5.3f, % 5.3f ]\n"\
              "[ % 5.3f, % 5.3f, % 5.3f ]\n\n",
               mat[0 * 3 + 0], mat[0 * 3 + 1], mat[0 * 3 + 2],
               mat[1 * 3 + 0], mat[1 * 3 + 1], mat[1 * 3 + 2],
               mat[2 * 3 + 0], mat[2 * 3 + 1], mat[2 * 3 + 2]);
}

void inverseMatrix(const double mat[9], double invMat[9])
{
  // mat = {
  //  a b c
  //  d e f
  //  g h i
  // }
  // mat = [a d g b e h c f i]
  const double& a = mat[0 * 3 + 0];
  const double& b = mat[0 * 3 + 1];
  const double& c = mat[0 * 3 + 2];
  const double& d = mat[1 * 3 + 0];
  const double& e = mat[1 * 3 + 1];
  const double& f = mat[1 * 3 + 2];
  const double& g = mat[2 * 3 + 0];
  const double& h = mat[2 * 3 + 1];
  const double& i = mat[2 * 3 + 2];

  double det = a*e*i +
               d*h*c +
               g*b*f -
               g*e*c -
               d*b*i -
               a*h*f;


  // double invMat[9] = {
  //    (e*i - h*f), -(b*i - h*c),  (b*f - e*c),
  //   -(d*i - g*f),  (a*i - g*c), -(a*f - d*c),
  //    (d*h - g*e), -(a*h - g*b),  (a*e - d*b)
  // };
  double invDet = 1. / det;

  invMat[0] =  (e*i - h*f) * invDet;
    invMat[1] = -(b*i - h*c) * invDet;
      invMat[2] =  (b*f - e*c) * invDet;

  invMat[3] = -(d*i - g*f) * invDet;
    invMat[4] =  (a*i - g*c) * invDet;
      invMat[5] = -(a*f - d*c) * invDet;

  invMat[6] =  (d*h - g*e) * invDet;
    invMat[7] = -(a*h - g*b) * invDet;
      invMat[8] =  (a*e - d*b) * invDet;

  if (std::abs(det) < 1e-16)
  {
    std::cerr << "det = " << det << '\n';
    printMatD(mat);
    printMatD(invMat);
    return;
  }

};

void TRlead(double* const TMat, uint8_t lead)
{
  double T = 0;
  double R = 0;
  for (uint8_t i = 0; i < 4; i++)
  {
    if (i != lead - 1)
    {
      // T += TMat[lead - 1][i];
      T += TMat[index4(i, lead - 1)];
    }
    else
    {
      // R += TMat[lead - 1][i];
      R += TMat[index4(i, lead - 1)];
    }
    /* code */
  }
  std::printf("T: % 5.5f, R: % 5.5f\n", T, R);
}

int readNumLeads(const std::filesystem::path& path)
{
  bool valid = true;
  std::ifstream file = openIFile(path, valid);
  if (!valid) return false;

  constexpr size_t buffSize = 128;
  char buff[buffSize] = {};

  file.read(buff, buffSize);

  char* current_p = buff;
  size_t remainingSize = buffSize;
  valid &= check(readPast('=', &current_p, remainingSize));
  valid &= check(readWhile(' ', &current_p, remainingSize));

  // Save start of number
  char* start_p = current_p;

  valid &= check(readPast('\n', &current_p, remainingSize));
  char* end_p = current_p;

  int numLeads = -1;

  valid &= check(convertAndCheckNumber(numLeads, start_p, end_p));
  // dmsg("start: " << (size_t)start << " end: " << (size_t)end << " diff " << (size_t)(end - start));

  return numLeads;
}

inline bool operator!=(const leadInfoS& a, const leadInfoS& b)
{
  return !(a.currentLeadFrom == b.currentLeadFrom &&
           a.currentLeadTo == b.currentLeadTo &&
           a.voltageLeadHigh == b.voltageLeadHigh &&
           a.voltageLeadLow == b.voltageLeadLow);
}

bool readMatrix(const std::filesystem::path& path,
                double* const mat,
                const leadInfoS& leadInfo)
{
  constexpr uint numLeads = 4;
  constexpr uint numElements = numLeads * numLeads;

  // remapping is done to keep way of calculating of resistanced based on Datta (eradicate row and column 4)
  constexpr leadInfoS leadDesired =
  {
    .currentLeadFrom = 2,
    .currentLeadTo = 3,
    .voltageLeadHigh = 1,
    .voltageLeadLow = 4
  };

  const bool needsRemap = leadDesired != leadInfo;

  bool valid = true;
  std::ifstream file = openIFile(path, valid);

  if (!valid) return false;

  constexpr size_t buffSize = 1 << 12;
  char buff[buffSize] = {};
  size_t remainingSize = buffSize - 1; //1 less to keep last byte 0
  file.read(buff, remainingSize);

  // read header
  char* current_p = buff;
  valid &= check(readPast('=', &current_p, remainingSize));
  valid &= check(readWhile(' ', &current_p, remainingSize));

  char* start_p = current_p;
  valid &= check(readPast('\n', &current_p, remainingSize));
  char* end_p = current_p;

  int numLeadsCheck = -1;
  valid &= check(convertAndCheckNumber(numLeadsCheck, start_p, end_p));
  if (numLeadsCheck != numLeads)
  {
    std::cout << "Wrong number of leads in a file: " << numLeadsCheck;
    return false;
  }

  // read elements
  for (int matLine = 0; matLine < numElements; matLine++)
  {
    // read 'to'
    valid &= check(readWhile(' ', &current_p, remainingSize));
    char* start_p = current_p;
    valid &= check(readPast(',', &current_p, remainingSize));
    char* end_p = current_p;
    int to = -1;
    valid &= check(convertAndCheckNumber(to, start_p, end_p));

    // read 'from'
    valid &= check(readWhile(' ', &current_p, remainingSize));
    start_p = current_p;
    valid &= check(readPast(',', &current_p, remainingSize));
    end_p = current_p;
    int from = -1;
    valid &= check(convertAndCheckNumber(from, start_p, end_p));

    // read transmission
    // read 1 after ','
    // current_p++;
    // remainingSize--;
    start_p = current_p;
    valid &= check(readPast('\n', &current_p, remainingSize));
    end_p = current_p;
    double t = -1;
    valid &= check(convertAndCheckNumber(t, start_p, end_p));

    current_p++;
    // current_p++;
    remainingSize--;
    // remainingSize--;

    // std::printf("to: %d, from %d ", to, from);
    valid &= (to - 1 < numLeads);
    valid &= (from - 1 < numLeads);
    valid &= (to - 1 >= 0);
    valid &= (from - 1 >= 0);
    valid &= (numLeads == 4);

    if (needsRemap)
    {
      // remapping is done to keep way of calculating of resistanced based on Datta (eradicate row and column 4)
      remapLeads(to, from, leadInfo, leadDesired);
    }
    // std::printf("to: %d, from %d, t: %f\n", to, from, t);

    mat[index4(to - 1, from - 1)] = t;
  }

  return valid;
}
