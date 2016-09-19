//
// File:    Input.h
//
#if !defined(INPUT_H)
#define INPUT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

// vector
#ifdef HAVE_VECTOR
#include <vector>
#else
#ifdef HAVE_VECTOR_H
#include <vector.h>
#else
#error No vector or vector.h available
#endif // HAVE_VECTOR_H
#endif // HAVE_VECTOR

#ifdef STDC_HEADERS
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#include "DataTypes.h"

//
//
//
namespace quasicontinuum {

/**
 * @brief Complex input class
 */
class Input {

  //
  // public data types
  //
public:
  //
  // public methods
  //

  //  material info
  char d_initFile[256];
  char d_dataFile[256];
  char d_material[20];
  char d_potentialName[100];
  char d_potentialType[20];
  char d_materialsFile[256];

  //  restart file info
  int d_restartFlag;
  char d_restartFile[256];

  //  minimization method related data
  int d_minMethod;
  int d_minMethodMaxIter;
  double d_minMethodTol;

  //  test flag
  int d_testFlag;

  //  QC option
  struct qc_options_t d_qcOption;

  //  lattice info
  int d_numQuasi;
  char d_latticeType[10];
  std::vector<std::vector<double>> d_latticeConstants;
  std::pair<std::vector<int>, std::vector<int>> d_latticeMinMax;
  std::pair<std::vector<double>, std::vector<double>> d_configurationMinMax;
  std::vector<std::vector<double>> d_shift;

  // initial config info
  int d_initialConfigurationFlag;

  // temperature info
  int d_tempFlag;
  double d_temperature;
  double d_tau;
  double d_boltzmanConstant;
  std::vector<double> d_initialFrequency;
  int d_quadratureMethod;

  // electrostatics info
  int d_electrostaticFlag;
  int d_externalChargeFlag;
  std::vector<std::pair<double, std::vector<double>>> d_externalCharges;
  std::vector<std::pair<int, double>> d_electroBCs;
  double d_electrostaticCutoffRadius;
  int d_electrostaticIntegration;

  // external loading info
  int d_externalLoadFlag;
  int d_numLoads;
  std::vector<std::vector<double>> d_deformationMatrix;

  // residual force info
  std::vector<int> d_residualForceRemoveFlag;

  // void info
  int d_voidFlag;
  char d_voidType[20];
  std::vector<double> d_voidCenter;
  int d_voidNumParams;
  std::vector<double> d_voidParams;

  //  indent info
  int d_indentFlag;
  std::vector<double> d_indentInitialLocation;
  std::vector<double> d_indentIncrement;
  double d_indentRadius;
  double d_indentConstant;

  // cg info
  double d_cgTolerance;
  int d_cgMaxIterations;
  int d_debugLevel;
  double d_cgLineSearchTolerance;
  int d_cgLineSearchMaxIterations;
  double d_remeshTolerance;

  //
  //  threading info
  //
  int d_numberThreads;

  //
  //  Output directory
  //
  char *d_outputDirectory;

  //
  //  lattice info
  //
  std::pair<int, std::vector<std::pair<int, std::pair<int, double>>>>
      d_latticeInfo;

  //
  //  uyniversal constants
  //
  std::pair<int, std::vector<double>> d_universalConst;

public:
  /**
   * @brief getInstance.
   */
  static Input *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  //
  //  parameterInit()
  //
  void parameterInit();

  /**
   * @brief Get main input
   *
   * @param numQuasi Number of Quasicontinuum in problem.
   * @param numLoads Number of load steps in problem.
   * @param boundaryFlag Handle body boundary conditions.
   * @param loadDeformation Incremental load to apply.
   * @param tolerance CG convergence tolerance.
   * @param maxIterations Max number of CG iterations.
   * @param debugLevel CG debug output level.
   * @param lineTolerance Line search tolerance.
   * @param lineIterations Max number of line search iterations.
   * @param getDefects Process file for defects option.
   * @param addNodeCutoff Add node cutoff total.
   * @param addNodeFlag Add node flag options.
   * @param addNodeRestartNum Add node restart number.
   */
  void mainInput(int &numQuasi, int &numLoads, int &shearFlag,
                 std::vector<std::vector<double>> &loadDeformation,
                 double &tolerance, int &maxIterations, int &debugLevel,
                 double &lineTolerance, int &lineIterations, int &getDefects,
                 int &addNodeCutoff, int &addNodeFlag, int &addNodeRestartNum,
                 std::vector<double> &defectCenter,
                 std::vector<double> &defectBox, double &defectLat,
                 std::vector<std::pair<double, std::vector<double>>> &charges,
                 double &temperature, int &statistics, int &quadMethod,
                 int &outputFlag, char *&outputDirectory, int &testFlag,
                 int &minMethod, int &minMethodMaxIter, double &minMethodTol,
                 int &initialConfigFlag, int ac, char **av);

  /**
   * @brief Get quasi input
   * @param data types of quasicontinuum
   */

  struct qc_options_t quasiInput(const int iQuasi, struct lattice_t *P_lattice,
                                 struct all_node_list_t *P_node_list,
                                 struct element_list_t *P_element_list,
                                 struct indentor_t *P_indentor,
                                 char *materialName, double &mass, int ac,
                                 char **av);

  //
  //  systemInfoOutput()
  //
  void systemInfoOutput();

private:
  /**
   * @brief Constructor.
   */
  Input();

  /**
   * @brief Copy constructor.
   */
  Input(Input const &);

  /**
   * @brief Assignment operator.
   */
  const Input &operator=(const Input &);

  /**
   * @brief Destructor.
   */
  ~Input();

  //
  // private methods
  //

private:
  /**
   * @brief Load potential.
   *
   * @param potentialType Potential to load.
   */
  void LoadPotential(const int potentialType);

  /**
   * @brief read quasi data.
   *
   */
  void ReadQuasiData(const int iQuasi, struct lattice_t *P_lattice,
                     struct all_node_list_t *P_node_list, char *materialName,
                     char *input_filename);

  //
  //  ReadRestartFile()
  //
  void ReadRestartFile(int iQuasi, const char *filename,
                       struct all_node_list_t *P_all_node_list,
                       struct element_list_t *P_element_list,
                       struct indentor_t *P_indentor,
                       struct lattice_t *P_lattice, char *materialName);

  //
  //  ReadAtomicMass()
  //
  void ReadAtomicMass(int iQuasi, double &mass);

  //
  //  check if Quasi is core
  //
  int IsQuasiCore(int iQuasi);

  /**
   * @brief ProcessOptions()
   *
   */
  void ProcessOptions(char **init_file, char **data_file, int *number_threads,
                      enum switch_t *restart, char **restart_file,
                      double *temperature, double *tau, int argc, char **argv);

  /**
   * @brief create defaul init filename
   *
   */
  int CreateDefaultInitFilename(char *init_file, size_t n, char *name,
                                char const *ext);

  //
  // display command line input error
  //
  void ErrUsage(char *pname);

  //
  // private data types
  //
private:
  static Input *_instance;
};
}

#endif // INPUT_H
