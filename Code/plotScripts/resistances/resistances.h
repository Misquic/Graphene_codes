#pragma once

// #include <string_view>
#include <inttypes.h>
#include <filesystem>

bool check_(int result, int line);
#ifdef DEBUG
  #define printMatD(mat) printMat(mat)
  #define printMat4D(mat) printMat4(mat)
  #define dmsg(x) std::cerr << x << '\n'
  #define check(result) check_((result), __LINE__);
#else
  #define printMatD(mat)
  #define printMat4D(mat)
  #define dmsg(x)
  #define check(result) (result)
#endif

#define index4(i, j) ((i) * 4 + (j))
#define index3(i, j) ((i) * 3 + (j))


/*
      [2] --V-> [3]
       |_________|
  [1]--[         ]--[4]
   |   [_________]   |
   |                 |
   |  <-----I-----   |

   4 - currentLeadFrom   // e- flows to
   1 - currentLeadTo     // e- flows from
   3 - voltageLeadHigh
   2 - voltageLeadLow
*/

typedef struct LeadInfoS
{
  uint currentLeadFrom; //
  uint currentLeadTo;   //
  uint voltageLeadHigh; //
  uint voltageLeadLow;  //
} LeadInfoS;

constexpr LeadInfoS leadDesired =
{
  .currentLeadFrom = 4,
  .currentLeadTo = 1,
  .voltageLeadHigh = 3,
  .voltageLeadLow = 2
};

std::ifstream openIFile(const std::string& path, bool& valid);
std::ofstream openOFile(const std::string& path, bool& valid);

size_t getNumOfSubDirs(const std::filesystem::path& dirs);

double calculateR(double* const T);

void printMat(const double* const mat);
void printMat4(const double* const mat);

void inverseMatrix(const double mat[9], double invMat[9]);

int readNumLeads(const std::filesystem::path& path);

bool readMatrix(const std::filesystem::path& path,
                double* const mat,
                const LeadInfoS& leadInfo,
                bool needsRemap);

void TRlead(double* TMat, uint8_t lead);

inline bool operator!=(const LeadInfoS& a, const LeadInfoS& b)
{
  return !(a.currentLeadFrom == b.currentLeadFrom &&
           a.currentLeadTo == b.currentLeadTo &&
           a.voltageLeadHigh == b.voltageLeadHigh &&
           a.voltageLeadLow == b.voltageLeadLow);
}
