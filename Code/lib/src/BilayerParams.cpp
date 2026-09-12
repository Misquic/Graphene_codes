#include "BilayerPositionGeneration.h"

//////////////// local function declarations ////////////////

static void displayUsage(const char* const* const argv);

static void readArgsAndSetParams(const int argc, const char* const* const argv, ParamsS& params);

static void printParams(const ParamsS& params);

//////////////// global function definitions ////////////////

void readArgs(const int argc, const char* const* const argv, ParamsS& params)
{
  // setDefaultParams(params);

  if (argc == 1)
  {
    printParams(params);
    return;
  }

  std::string firstArg = argv[1];

  // display help and
  if (firstArg == "help")
  {
    displayUsage(argv);
    std::exit(0);
  }
  else
  {
    readArgsAndSetParams(argc, argv, params);
  }

  printParams(params);
}

void printPositions(const std::vector<Vec3S>& positions)
{
  std::printf("  {\n");
  for (const Vec3S& pos: positions)
  {
    std::printf("    {%f, %f, %f},\n", pos.x, pos.y, pos.z);
  }
  std::printf("  };\n");
}

void savePositions(const std::string_view name, const std::vector<Vec3S>& positions)
{
  std::ofstream outFile(name.data(), std::ios::out | std::ios::binary);
  if (!outFile.is_open())
  {
    std::cerr << "Couldn't open a file " << name << "\n Saving to ./temp.csv";
    outFile.clear();
    outFile.open("./temp.csv", std::ios::out | std::ios::binary);
  }

  for (const Vec3S& pos: positions)
  {
    outFile << pos.x << ',' << pos.y << ',' << pos.z << '\n';
  }
}

//////////////// local function definitions ////////////////


static void displayUsage(const char* const* const argv)
{
  std::string appName = argv[0];
  size_t slashPos = appName.find_last_of('\\');
  size_t dotPos = appName.find_last_of('.');
  appName = appName.substr(slashPos + 1, dotPos - slashPos - 1);

  std::printf("Usage: %s Nx Ny R ",
              appName.c_str());
}

static void readArgsAndSetParams(const int argc, const char* const* const argv, ParamsS& params)
{
  char* end = nullptr;
  int currentArg = 2;
  if (argc >= currentArg)
  {
    params.nX = std::strtoull(argv[currentArg - 1], &end, 10);
  }
  currentArg++;

  if (argc >= currentArg)
  {
    params.nY = std::strtoull(argv[currentArg - 1], &end, 10);
  }
  currentArg++;

  if (argc >= currentArg)
  {
    params.foldRadius = std::strtof(argv[currentArg - 1], &end);
  }
  currentArg++;

  if (argc >= currentArg)
  {
    params.leadWidth = std::strtoull(argv[currentArg - 1], &end, 10);
  }
  currentArg++;


  makeParamsRight(params);
}

void makeParamsRight(ParamsS& params)
{
  // make it even
  if (params.nY % 2 == 1)
  {
    params.nY++;
  }

  // make sure leadWidth is not too big
  if (params.leadWidth > params.nY - params.cutLead)
  {
    params.leadWidth = params.nY - params.cutLead - 1;
  }

  // make it odd
  if (params.leadWidth % 2 == 0)
  {
    params.leadWidth++;
  }

}


static void printParams(const ParamsS& params)
{
  std::printf("Params has been set: \n"
              "Nx = %lu "
              "Ny = %lu "
              "foldRadius = %f "
              "leadWidth = %lu "
              "saveFileName = %s"
              "\n",
              params.nX,
              params.nY,
              params.foldRadius,
              params.leadWidth,
              params.saveFileName.data());
}

