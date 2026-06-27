#pragma once

// #include <string_view>
#include <inttypes.h>
#include <filesystem>

bool check_(int result, int line);

#ifdef NDEBUG
  #define printMatD(mat)
  #define dmsg(x)
  #define check(result) (result)
#else
  #define printMatD(mat) printMat(mat)
  #define dmsg(x) std::cerr << x << '\n'
  #define check(result) check_((result), __LINE__);
#endif

#define index4(i, j) ((i) * 4 + (j))
#define index3(i, j) ((i) * 3 + (j))

typedef struct leadInfoS
{
  uint currentLeadFrom;
  uint currentLeadTo;
  uint voltageLeadHigh;
  uint voltageLeadLow;
} leadInfoS;

std::ifstream openIFile(const std::string& path, bool& valid);
std::ofstream openOFile(const std::string& path, bool& valid);

size_t getNumOfSubDirs(const std::filesystem::path& dirs);

double calculateR(double* const T);

void printMat(const double* const mat);

void inverseMatrix(const double mat[9], double invMat[9]);

int readNumLeads(const std::filesystem::path& path);

bool readMatrix(const std::filesystem::path& path,
                double* const mat,
                const leadInfoS& leadInfo);

void TRlead(double* TMat, uint8_t lead);
