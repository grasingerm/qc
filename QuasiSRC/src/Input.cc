//
//  Input.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
/* no _REENTRANT */
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#ifdef STDC_HEADERS
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_CTYPE_H
#include <ctype.h>
#else
#error ctype.h not found.
#endif /* HAVE_CTYPE_H */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found.
#endif /* HAVE_UNISTD_H */

#ifdef HAVE_LIBGEN_H
#include <libgen.h>
#else
#error libgen.h not found.
#endif /* HAVE_LIBGEN_H */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#error limits.h not found.
#endif /* HAVE_LIMITS_H */

#ifdef HAVE_RPC_RPC_H
#include <rpc/rpc.h>
#else
#error rpc/rpc.h not found.
#endif /* HAVE_RPC_RPC_H */

#ifdef HAVE_RPC_XDR_H
#include <rpc/xdr.h>
#else
#error rpc/xdr.h not found.
#endif /* HAVE_RPC_XDR_H */

#ifdef HAVE_VECTOR
#include <vector>
#else
#ifdef HAVE_VECTOR_H
#include <vector.h>
#else
#error No vector or vector.h available
#endif // HAVE_VECTOR_H
#endif // HAVE_VECTOR

#if defined(__QC_SGI)
#include <invent.h>
//#include <limits.h>  Already included above
#include <sys/types.h>
#endif /* sgi */

#include <math.h>

#include "C_Interface.h"
#include "CreateMesh.h"
#include "DataTypes.h"
#include "Electrostatics.h"
#include "Element.h"
#include "Error.h"
#include "ForceEnergyCalculation.h"
#include "Indent.h"
#include "Input.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "Node.h"
#include "Output.h"
#include "PairPotentials.h"
#include "QuadraturePoints.h"
#include "Quasicontinua.h"
#include "Void.h"

#include "monitor.h"
#include "read_pipe.h"
#include "threads.h"

#if !defined(lint)
static char rcsid[] = "$Id: input.c,v 1.14 2002/12/06 01:49:45 fago Exp $";
#endif /* !lint */

static const char *keywords[] = {
#define NUMBER_THREADS 0
    "number_threads",
#define PERIODIC 1
    "periodic",
#define DATA_FILE 2
    "data_file",
#define MATERIALS_FILE 3
    "materials_file",
#define INDENT_RADIUS 4
    "radius",
#define INDENT_CONSTANT 5
    "constant",
#define INDENT_DISPLACEMENT 6
    "displacement",
#define INDENT_POSITION 7
    "position",
#define INDENT_ENABLE 8
    "indentEnable",
#define VOID_ENABLE 9
    "voidEnable",
#define VOID_CENTER 10
    "voidCenter",
#define VOID_NUM_PARAMS 11
    "voidNumParams",
#define VOID_PARAMS 12
    "voidParams",
#define VOID_TYPE 13
    "voidType",
#define NUM_LOADS 14
    "numLoads",
#define BOUNDARY_FLAG 15
    "boundaryFlag",
#define DEF1 16
    "F1",
#define DEF2 17
    "F2",
#define DEF3 18
    "F3",
#define CHARGE 19
    "charge",
#define NUM_LATTICE 20
    "numLattice",
#define SHIFT 21
    "shift",
#define POTENTIAL 22
    "potential",
#define TOLERANCE 23
    "tolerance",
#define MAXITERATIONS 24
    "maxIterations",
#define DEBUGLEVEL 25
    "debugLevel",
#define LINETOLERANCE 26
    "lineTolerance",
#define LINEITERATIONS 27
    "lineIterations",
#define REMESHTOLERANCE 28
    "remeshTolerance",
#define CROSSNEIGHLISTFLAG 29
    "crossNeighListFlag",
#define GETDEFECTSFLAG 30
    "getDefectsFlag",
#define ADDNODEFLAG 31
    "addNodeFlag",
#define ADDNODERESTARTNUM 32
    "addNodeRestartNum",
#define ADDNODECUTOFF 33
    "addNodeCutoff",
#define DEFECT_CENTER 34
    "defectCenter",
#define DEFECT_BOX 35
    "defectBox",
#define DEFECT_LAT 36
    "defectLat",
#define X_BEGIN_E_BCS 37
    "xBeginElectroBCS",
#define X_END_E_BCS 38
    "xEndElectroBCS",
#define Y_BEGIN_E_BCS 39
    "yBeginElectroBCS",
#define Y_END_E_BCS 40
    "yEndElectroBCS",
#define Z_BEGIN_E_BCS 41
    "zBeginElectroBCS",
#define Z_END_E_BCS 42
    "zEndElectroBCS",
#define ALL_RESIDUAL 43
    "allResidual",
#define X_BEGIN_RESIDUAL 44
    "xBeginResidual",
#define X_END_RESIDUAL 45
    "xEndResidual",
#define Y_BEGIN_RESIDUAL 46
    "yBeginResidual",
#define Y_END_RESIDUAL 47
    "yEndResidual",
#define Z_BEGIN_RESIDUAL 48
    "zBeginResidual",
#define Z_END_RESIDUAL 49
    "zEndResidual",
#define EAM_METHOD 50
    "eam_method",
#define QUAD_METHOD 51
    "quadMethod",
#define STATISTICS 52
    "statistics",
#define TEMPERATURE 53
    "temperature",
#define BOLTZMAN_CONSTANT 54
    "boltzman_constant",
#define TAU 55
    "tau",
#define ELECTROSTATICS 56
    "electrostatics",
#define INITIAL_CONFIGURATION_FLAG 57
    "initial_config_flag",
#define OUTPUT_FLAG 58
    "output_flag",
#define OUTPUT_DIRECTORY 59
    "output_directory",
#define LATTICE_CONSTANT_1 60
    "a1",
#define LATTICE_CONSTANT_2 61
    "a2",
#define LATTICE_CONSTANT_3 62
    "a3",
#define LATTICE_MIN 63
    "l_min",
#define LATTICE_MAX 64
    "l_max",
#define PHYS_CONST_FLAG 65
    "phys_const_flag",
#define MAXPLANCK_CONST 66
    "maxplanck_constant",
#define ELECTRIC_CONST 67
    "electric_constant",
#define ATOMIC_MASS_FLAG 68
    "atomic_mass_flag",
#define RESTART_INFO 69
    "restart_flag",
#define TEST_FLAG 70
    "test_flag",
#define MIN_METHOD 71
    "min_method_flag",
#define ITER_MIN_METHOD 72
    "max_iter_min_method",
#define TOL_MIN_METHOD 73
    "tol_min_method",
// required in quasiInput()
#define EMPTY_KEY_MASK 0x000
#define NUMBER_THREADS_MASK 0x001
#define DATA_FILE_MASK 0x004
#define RADIUS_MASK 0x010
#define CONSTANT_MASK 0x020
#define DISPLACEMENT_MASK 0x040
#define POSITION_MASK 0x080
#define TEMPERATURE_MASK 0x100
};

// required in quasiInput()
#define DEFAULT_INIT_FILE_EXT ".ini"
#define NODE_RESTART_FILE "node_restart.plt"
#define ATOM_RESTART_FILE "atom_restart.plt"
#define MAX_NODE_NUMBER 128

//
//  boundary_fixity_flag
//  -1 - default
//   0 - fix x,y,z in z=0 plane
//   1 - fix x in x=0 plane, y in y=0 plane, z in z=0 plane
//   2 - all free
//
static int boundary_fixity_flag = 2;

// required in mainInput()

#define XDR_ERR(msg)                                                           \
  do {                                                                         \
    fprintf(stderr, "[%s:%d] %s\n", __FILE__, __LINE__, msg);                  \
    exit(EXIT_FAILURE);                                                        \
  } while (0)

//
//
//
namespace quasicontinuum {
// struct qc_options_t qc = QC_OPTIONS_INITIALIZER;
// static struct qc_options_t global_qc_options = qc;

//
//  initializing static _instance
//
Input *Input::_instance = NULL;

//
// constructor
//

Input::Input() {
  //
  //  call funtion to initialize the parameters
  //
  parameterInit();

  //
  //
  //
  return;
}

//
// destructor
//

Input::~Input() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Input *Input::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Input();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Input::destroyInstance() {

  //
  // delete instance
  //
  delete _instance;

  //
  //
  //
  return;
}

//
//  parameterInit()
//
void Input::parameterInit() {
  //
  //  material info
  //
  //  d_material = "NiAl" "GaN"
  //  d_potentialName = "MishinNiAl" "Ceria"
  //  d_potentialType = "pairiwise" "eam"
  //
  strcpy(d_material, "");
  strcpy(d_potentialName, "");
  strcpy(d_potentialName, "");

  //
  //  restart file handle
  //
  d_restartFlag = 0;

  //
  //  minimization related data
  //
  d_minMethod = 0; // default
  d_minMethodMaxIter = 1;
  d_minMethodTol = 0.00001;

  //
  //  test flag
  //
  d_testFlag = 0;

  //
  //  QC option
  //
  d_qcOption = QC_OPTIONS_INITIALIZER;

  //
  //  lattice info
  //
  //  d_numQuasi
  //  d_latticeType : "sc" "bcc" "wuzerite"
  //  d_latticeConstants : a1, a2, a3
  //  d_latticeMinMax : lattice coordinates beginning and ending
  //  d_configurationMinMax : initial config min and max
  //  d_shift
  //
  d_numQuasi = 0;
  strcpy(d_latticeType, "");

  d_latticeConstants.clear();
  d_latticeConstants.resize(3);
  for (int dof = 0; dof < 3; dof++) {
    d_latticeConstants[dof].push_back(0.0);
    d_latticeConstants[dof].push_back(0.0);
    d_latticeConstants[dof].push_back(0.0);
  }

  d_latticeMinMax.first.push_back(0);
  d_latticeMinMax.first.push_back(0);
  d_latticeMinMax.first.push_back(0);

  d_latticeMinMax.second.push_back(0);
  d_latticeMinMax.second.push_back(0);
  d_latticeMinMax.second.push_back(0);

  d_configurationMinMax.first.push_back(0.0);
  d_configurationMinMax.first.push_back(0.0);
  d_configurationMinMax.first.push_back(0.0);

  d_configurationMinMax.second.push_back(0.0);
  d_configurationMinMax.second.push_back(0.0);
  d_configurationMinMax.second.push_back(0.0);

  d_shift.clear();

  //
  //  initial config info
  //  d_initialConfigurationFlag = 0 - crystal lattice configuration
  //                               1 - equilibrated from crystal lattice
  //                               configuration
  //                               2 - specified by input file
  //
  d_initialConfigurationFlag = 0;

  //
  //  temperature info
  //
  //  d_tempFlag = 0 - zero temp
  //               1 - finite temp
  //               2 - non-equilibrium
  //  d_temperature
  //  d_tau
  //  d_boltzmanConstant
  //  d_initialFrequency
  //  d_quadratureMethod
  //
  d_tempFlag = 1;
  d_temperature = 300.0;
  d_tau = 0.1;
  d_boltzmanConstant = 1.0;
  d_initialFrequency.clear();
  d_quadratureMethod = 3;

  //
  //  electrostatics info
  //
  //  d_electrostaticFlag : 0 - no
  //                        1 - yes
  //  d_externalChargeFlag : 0 - no
  //                         1 - yes
  //  d_externalCharges : first - incremet
  //                      second - location
  //  d_electroBCs : first - 0 - zero charge density
  //                         1 - free surface
  //                         2 - specified charge density
  //  d_electrostaticCutoffRadius
  //  d_electrostaticIntegration
  //
  d_electrostaticFlag = 0;
  d_externalChargeFlag = 0;
  d_externalCharges.clear();

  d_electroBCs.clear();
  std::pair<int, double> dummyPair(0, 0.0);
  for (int i = 0; i < 6; ++i)
    d_electroBCs.push_back(dummyPair);

  d_electrostaticCutoffRadius = 0.0;
  d_electrostaticIntegration = 1;

  //
  //  external loading info
  //
  //  d_externalLoadFlag : -1 - no
  //                       0 - bulk
  //                       1 - surface
  //  d_numLoads
  //  d_deformationMatrix
  //
  d_externalLoadFlag = -1;
  d_numLoads = 0;
  d_deformationMatrix.clear();
  d_deformationMatrix.resize(3);
  for (int dof = 0; dof < 3; dof++) {
    d_deformationMatrix[dof].push_back(0.0);
    d_deformationMatrix[dof].push_back(0.0);
    d_deformationMatrix[dof].push_back(0.0);
  }

  //
  //  residual force info
  //
  //  d_residualForceRemoveFlag
  //
  d_residualForceRemoveFlag.clear();
  d_residualForceRemoveFlag.push_back(0);
  d_residualForceRemoveFlag.push_back(0);
  d_residualForceRemoveFlag.push_back(0);
  d_residualForceRemoveFlag.push_back(0);
  d_residualForceRemoveFlag.push_back(0);
  d_residualForceRemoveFlag.push_back(0);
  d_residualForceRemoveFlag.push_back(0);

  //
  //  void info
  //
  //  d_voidFlag : 0 - no
  //               1 - yes
  //  d_voidType : "circle" "vpitt"
  //  d_voidParams
  //
  d_voidFlag = 0;
  strcpy(d_voidType, "");

  d_voidCenter.clear();
  d_voidCenter.push_back(0.0);
  d_voidCenter.push_back(0.0);
  d_voidCenter.push_back(0.0);

  d_voidNumParams = 0;
  d_voidParams.clear();

  //
  //  indent info
  //
  //  d_indentFlag : 0 - no
  //                 1 - yes
  //  d_indentInitialLocation
  //
  d_indentFlag = 0;
  d_indentInitialLocation.clear();
  d_indentInitialLocation.push_back(0.0);
  d_indentInitialLocation.push_back(0.0);
  d_indentInitialLocation.push_back(0.0);

  d_indentIncrement.clear();
  d_indentIncrement.push_back(0.0);
  d_indentIncrement.push_back(0.0);
  d_indentIncrement.push_back(0.0);

  d_indentRadius = 0.0;
  d_indentConstant = 0.0;

  //
  //  cg info
  //
  d_cgTolerance = 0.0;
  d_cgMaxIterations = 0;
  d_debugLevel = 1;
  d_cgLineSearchTolerance = 0.0;
  d_cgLineSearchMaxIterations = 0;
  d_remeshTolerance = 0.0;

  //
  //  threading info
  //
  d_numberThreads = 0;

  //
  //  atomic mass
  //
  d_latticeInfo.first = 0;

  //
  //  universal const
  //
  d_universalConst.first = 0;
  d_universalConst.second.clear();
  d_universalConst.second.push_back(1.0);
  d_universalConst.second.push_back(1.0);
  d_universalConst.second.push_back(1.0);

  //
  //  Output Directory
  //
  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);

  d_outputDirectory = (char *)malloc(path_max * sizeof(char));

  //
  //  done
  //
  return;
}

//
// read data for main.cc
//
void Input::mainInput(
    int &numQuasi, int &numLoads, int &boundaryFlag,
    std::vector<std::vector<double>> &loadDeformation, double &tolerance,
    int &maxIterations, int &debugLevel, double &lineTolerance,
    int &lineIterations, int &getDefects, int &addNodeCutoff, int &addNodeFlag,
    int &addNodeRestartNum, std::vector<double> &defectCenter,
    std::vector<double> &defectBox, double &defectLat,
    std::vector<std::pair<double, std::vector<double>>> &charges,
    double &temperature, int &statistics, int &quadMethod, int &outputFlag,
    char *&outputDirectory, int &testFlag, int &minMethod,
    int &minMethodMaxIter, double &minMethodTol, int &initialConfigFlag, int ac,
    char **av) {
  //
  //  local parameters
  //
  FILE *init;
  long line_max;
  char *data;
  char *line = NULL;

  char dummy[20];
  int numShifts = 0;

  //
  //  flags that indicates whether electroBCS and d_electrostaticFlag
  //  data has been read before setting up the interatomic
  //  potential
  //
  std::vector<int> flags_for_LoadPotential;
  flags_for_LoadPotential.resize(2);
  flags_for_LoadPotential[0] = 0;
  flags_for_LoadPotential[1] = 0;

  //
  //  universal constant read flag
  //  value of this flag will show whether we have read all
  //  universal constants
  //
  int universalConstFlag = 0;

  long path_max = pathconf("/", _PC_PATH_MAX);

  if (path_max == -1L) {
    fprintf(stderr, "%s %d: Path max failure.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  //
  //  ProcessOptions()  : to read from command line input
  //
  char *init_file = NULL;
  char *data_file = NULL;
  char *restart_file = NULL;
  double temperature_cmd;
  double tau_cmd;
  int number_threads = 0;

  //
  temperature_cmd = -0.1;
  tau_cmd = -0.1;

  ProcessOptions(&init_file, &data_file, &number_threads, &d_qcOption.restart,
                 &restart_file, &temperature_cmd, &tau_cmd, ac, av);

  //
  //  set the Input data according to run time input
  //

  // initial node input data file
  if (data_file != NULL) {
    strcpy(d_dataFile, data_file);
  }

  // init data file
  if (init_file != NULL) {
    // copy the command line filename
    strcpy(d_initFile, init_file);
  } else {
    if ((init_file = (char *)malloc(path_max * sizeof(char))) == NULL) {
      fprintf(stderr, "%s %d: Malloc error.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    // CreateDefaultInitFilename( init_file, path_max, av[0],
    // DEFAULT_INIT_FILE_EXT );

    // Why did someone decide this was a good idea?
    // If you're reading this, it's too late -- Drake
    strcpy(init_file, "quasi.ini");

    // copy the default init filename
    strcpy(d_initFile, init_file);
  }

  //  number of threads
  if (number_threads > 0) {
    d_numberThreads = number_threads;

    // initialize in this case
    qc_thread_init(d_numberThreads);
  }

  //
  //  if d_restartFlag is 1, but d_qcOption.restart is OFF
  //  change it to ON and copy the d_dataFile to d_restartFile
  //
  //  if d_qcOption is ON put 1 in d_restartFlag
  //
  if (d_qcOption.restart == ON) {
    // restart filename is provided in command line, therefore
    // set restart flag to 1, and also copy the restart filename
    // to d_restartFile
    d_restartFlag = 1;

    strcpy(d_restartFile, restart_file);
  }

  //
  // open init file
  //
  // printf("d_initFile = %s\n", d_initFile);
  if ((init = fopen(d_initFile, "r")) == NULL) {
    fprintf(stderr, "%s %d: Unable to open init file: %s\n", __FILE__,
            __LINE__, d_initFile);
    exit(EXIT_FAILURE);
  }

  // allocate space for line buffer
  if ((line_max = sysconf(_SC_LINE_MAX)) == -1L) {
    fprintf(stderr, "%s %d: Line buffer space allocation failed.\n", 
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  if ((line = (char *)malloc(line_max * sizeof(char))) == NULL) {
    fprintf(stderr, "%s %d: Malloc failure.\n", 
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  while (fgets(line, (int)(line_max - 1L), init) != NULL) {
    // skip comments or empty lines
    if (line[0] == '#' || line[0] == '\n')
      continue;

    // match keywords

    //
    // read number of thread and call qc_thread_init()
    //
    if (strstr(line, keywords[NUMBER_THREADS]) != NULL) {
      // proceed only if command line input argv did not
      // give number of threads data
      if (d_numberThreads == 0) {
        // set if according to the quasi.ini file
        sscanf(line, "%*s %d", &d_numberThreads);

        // also initialize threads
        qc_thread_init(d_numberThreads);
      }
      // else
      // {
      //   // d_numberThreads was fixed by command line argument
      //   // nothing to do
      // }
    }

    //
    //  read data file
    //
    if (strstr(line, keywords[DATA_FILE]) != NULL) {
      sscanf(line, "%*s %s", d_dataFile);
    }

    //
    //  read data file
    //
    if (strstr(line, keywords[MATERIALS_FILE]) != NULL) {
      sscanf(line, "%*s %s", d_materialsFile);
    }

    //
    // restart flag
    //
    if (strstr(line, keywords[RESTART_INFO]) != NULL) {
      // when restart flag is 1, we expect the data file to be
      // restart file

      // only read the restart flag value if it is not been set
      // from command line input
      // As command line input must override the
      // quasi.ini file content.

      if (d_restartFlag == 0) {
        // command line did not set restart flag, and hence,
        // read from quasi.ini
        sscanf(line, "%*s %d", &d_restartFlag);

        if (d_restartFlag == 1) {
          // value is 1, meaning we are reading from restart
          // data file
          d_qcOption.restart = ON;
          strcpy(d_restartFile, d_dataFile);
        }
      }
    }

    //
    //  minimization method related input parameters
    //
    if (strstr(line, keywords[MIN_METHOD]) != NULL) {
      sscanf(line, "%*s %d", &d_minMethod);
    }

    if (strstr(line, keywords[ITER_MIN_METHOD]) != NULL) {
      sscanf(line, "%*s %d", &d_minMethodMaxIter);
    }

    if (strstr(line, keywords[TOL_MIN_METHOD]) != NULL) {
      sscanf(line, "%*s %lf", &d_minMethodTol);
    }

    //
    // test flag
    //
    if (strstr(line, keywords[TEST_FLAG]) != NULL) {
      sscanf(line, "%*s %d", &d_testFlag);
    }

    //
    // read temperature from quasi.ini
    //
    if (strstr(line, keywords[TEMPERATURE]) != NULL) {
      //
      //  since temperature may also be provided through command
      //  line, and, command line data will override the data in
      //  quasi init file.
      //
      //  check temperature_cmd
      if (temperature_cmd < 0.0) {
        // temperature is not provided
        // read temperature from quasi init file
        sscanf(line, "%*s %lf", &d_temperature);
      } else {
        // use command line input for temp
        d_temperature = temperature_cmd;
      }

      // set temperature of quasicontinua
      Quasicontinua::getInstance()->setTemperature(d_temperature);
    }

    //
    // read whether it is zero temp problem or finite temp
    // zero temp is not implemented at present
    //
    if (strstr(line, keywords[STATISTICS]) != NULL) {
      sscanf(line, "%*s %d", &statistics);

      if (statistics == 1)
        d_tempFlag = 1;
    }

    //
    // read tau from quasi.ini
    //
    if (strstr(line, keywords[TAU]) != NULL) {
      //
      //  since tau may also be provided through command
      //  line, and, command line data will override the data in
      //  quasi init file.
      //
      //  check tau_cmd
      if (tau_cmd < 0.0) {
        // tau is not provided
        // read tau from quasi init file
        sscanf(line, "%*s %lf", &d_tau);
        // printf("here tau = %f\n",d_tau);
      } else {
        // use command line input for tau
        d_tau = tau_cmd;
      }
    }

    //
    //
    //
    if (strstr(line, keywords[PHYS_CONST_FLAG]) != NULL) {
      sscanf(line, "%*s %d", &d_universalConst.first);
    }

    //
    //
    //
    if (d_universalConst.first == 1) {
      if (strstr(line, keywords[BOLTZMAN_CONSTANT]) != NULL) {
        sscanf(line, "%*s %lf", &d_universalConst.second[0]);
        universalConstFlag = universalConstFlag + 1;
      }

      if (strstr(line, keywords[MAXPLANCK_CONST]) != NULL) {
        sscanf(line, "%*s %lf", &d_universalConst.second[1]);
        universalConstFlag = universalConstFlag + 1;
      }

      if (strstr(line, keywords[ELECTRIC_CONST]) != NULL) {
        sscanf(line, "%*s %lf", &d_universalConst.second[2]);
        universalConstFlag = universalConstFlag + 1;
      }

      // if we have read all the data, modify it in PairPotentials class
      if (universalConstFlag == 3) {
        PairPotentials::getInstance()->setUniversalConstants(
            d_universalConst.second);
      }
    }

    //
    // read electrostatics flag. If 0 it's not enabled, 1 enabled
    //
    if (strstr(line, keywords[ELECTROSTATICS]) != NULL) {
      sscanf(line, "%*s %d", &d_electrostaticFlag);
      flags_for_LoadPotential[0] = 1;
    }

    //
    //  read outputFlag
    //
    if (strstr(line, keywords[OUTPUT_FLAG]) != NULL) {
      sscanf(line, "%*s %d", &outputFlag);
    }

    //
    //  based on output flag either read the directory from quasi.ini
    //  or create new output directory
    //
    if (strstr(line, keywords[OUTPUT_DIRECTORY]) != NULL) {
      // check the value of flag
      if (outputFlag == 1) {
        // read directory from quasi.ini
        sscanf(line, "%*s %s", outputDirectory);

        // call Output class function to create directory
        Output::getInstance()->createOutputDirectory(outputFlag,
                                                     outputDirectory);
      } else {
        // create output directory where quasi executable is stored (default)
        Output::getInstance()->createOutputDirectory(outputFlag,
                                                     outputDirectory);
      }

      strcpy(d_outputDirectory, outputDirectory);
    }

    //
    // read initial_config_flag from quasi.ini
    //
    if (strstr(line, keywords[INITIAL_CONFIGURATION_FLAG]) != NULL) {
      sscanf(line, "%*s %d", &d_initialConfigurationFlag);
    }

    //
    // read nano indentation flag
    //
    if (strstr(line, keywords[INDENT_ENABLE]) != NULL) {
      sscanf(line, "%*s %d", &d_indentFlag);
      // set the indent flag in Indent class
      Indent::getInstance()->setNanoIndentation(d_indentFlag);
    }

    if (strstr(line, keywords[INDENT_RADIUS]) != NULL) {
      sscanf(line, "%*s %lf", &d_indentRadius);
    }

    if (strstr(line, keywords[INDENT_CONSTANT]) != NULL) {
      sscanf(line, "%*s %lf", &d_indentConstant);
    }

    if (strstr(line, keywords[INDENT_DISPLACEMENT]) != NULL) {
      sscanf(line, "%*s %lf %lf %lf", &d_indentIncrement[0],
             &d_indentIncrement[1], &d_indentIncrement[2]);
    }

    if (strstr(line, keywords[INDENT_POSITION]) != NULL) {
      sscanf(line, "%*s %lf %lf %lf", &d_indentInitialLocation[0],
             &d_indentInitialLocation[1], &d_indentInitialLocation[2]);
    }

    //
    // read void flag and void parameters
    //
    if (strstr(line, keywords[VOID_ENABLE]) != NULL) {
      sscanf(line, "%*s %d", &d_voidFlag);

      // set the flag in Void class
      Void::getInstance()->setVoid(d_voidFlag);
    }

    if (strstr(line, keywords[VOID_CENTER]) != NULL) {
      sscanf(line, "%*s %lf %lf %lf", &d_voidCenter[0], &d_voidCenter[1],
             &d_voidCenter[2]);
    }

    if (strstr(line, keywords[VOID_NUM_PARAMS]) != NULL) {
      sscanf(line, "%*s %d", &d_voidNumParams);

      // resize voidParams
      d_voidParams.resize(d_voidNumParams);
    }

    if (strstr(line, keywords[VOID_PARAMS]) != NULL) {
      if (d_voidNumParams == 1) {
        sscanf(line, "%*s %lf", &d_voidParams[0]);
      } else if (d_voidNumParams == 2) {
        sscanf(line, "%*s %lf %lf", &d_voidParams[0], &d_voidParams[1]);
      } else {
        d_print("Error mainInput()\n");
        fprintf(stderr, "%s %d: I have no idea what is suppose to happen here\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
    }

    if (strstr(line, keywords[VOID_TYPE]) != NULL) {
      sscanf(line, "%s %s", dummy, d_voidType);

      // set void parameters in Void class
      Void *voidC = Void::getInstance();

      if (voidC->isVoidEnable() == 1) {
        double center[3];
        for (int dof = 0; dof < 3; dof++)
          center[dof] = d_voidCenter[dof];

        if (d_voidType[0] == 'C') {
          enum void_t type = CIRCLE;
          voidC->setVoidParameters(center, type, d_voidNumParams, d_voidParams);
        } else if (d_voidType[0] == 'V') {
          enum void_t type = VPIT;
          voidC->setVoidParameters(center, type, d_voidNumParams, d_voidParams);
        } else {
          d_print("Error in mainInput\n");
          fprintf(stderr, "%s %d: I have no idea what is suppose to happen here\n",
                  __FILE__, __LINE__);
          exit(EXIT_FAILURE);
        }
      }
    }

    //
    //  read number of loads
    //
    if (strstr(line, keywords[NUM_LOADS]) != NULL) {
      sscanf(line, "%*s %d", &d_numLoads);
    }

    //
    //  read boundary flag
    //
    if (strstr(line, keywords[BOUNDARY_FLAG]) != NULL) {
      sscanf(line, "%*s %d", &d_externalLoadFlag);
    }

    //
    //  read load deformation
    //
    if (strstr(line, keywords[DEF1]) != NULL) {
      sscanf(line, "%*s %lf %lf %lf", &d_deformationMatrix[0][0],
             &d_deformationMatrix[0][1], &d_deformationMatrix[0][2]);
    }

    if (strstr(line, keywords[DEF2]) != NULL) {
      sscanf(line, "%*s %lf %lf %lf", &d_deformationMatrix[1][0],
             &d_deformationMatrix[1][1], &d_deformationMatrix[1][2]);
    }

    if (strstr(line, keywords[DEF3]) != NULL) {
      sscanf(line, "%*s %lf %lf %lf", &d_deformationMatrix[2][0],
             &d_deformationMatrix[2][1], &d_deformationMatrix[2][2]);
    }

    //
    //  read charges
    //
    if (strstr(line, keywords[CHARGE]) != NULL) {
      std::vector<double> charge_location(3, 0.0);
      double charge_increment;
      int chargeNum;

      sscanf(line, "%*s %d %lf %lf %lf %lf", &chargeNum, &charge_increment,
             &charge_location[0], &charge_location[1], &charge_location[2]);

      if (chargeNum == d_externalCharges.size())
        d_externalCharges.push_back(
            std::make_pair(charge_increment, charge_location));
    }

    //
    //  read atomic mass flag
    //
    if (strstr(line, keywords[ATOMIC_MASS_FLAG]) != NULL) {
      sscanf(line, "%*s %d", &d_latticeInfo.first);
    }

    //
    //  read number of lattices in multilattice system
    //
    if (strstr(line, keywords[NUM_LATTICE]) != NULL) {
      sscanf(line, "%*s %d", &d_numQuasi);
    }

    //
    //  read shift vector
    //
    if (strstr(line, keywords[SHIFT]) != NULL) {
      std::vector<double> shifts(3, 0.0);
      int shiftNum;
      int shiftFlag;
      double atomic_mass;

      sscanf(line, "%*s %d %lf %lf %lf %d %lf", &shiftNum, &shifts[0],
             &shifts[1], &shifts[2], &shiftFlag, &atomic_mass);

      std::pair<int, std::pair<int, double>> data;
      data.first = shiftNum;
      data.second.first = shiftFlag;
      data.second.second = atomic_mass;

      d_latticeInfo.second.push_back(data);

      if (shiftNum == numShifts) {
        Quasicontinua::getInstance()->insertShift(shiftNum, shifts, shiftFlag);
        numShifts++;
      }

      d_shift.push_back(shifts);
    }

    //
    //  read boundary conditions
    //
    if (strstr(line, keywords[X_BEGIN_E_BCS]) != NULL) {
      sscanf(line, "%*s %d %lf", &d_electroBCs[0].first,
             &d_electroBCs[0].second);
      flags_for_LoadPotential[1] += 1;
    }

    if (strstr(line, keywords[X_END_E_BCS]) != NULL) {
      sscanf(line, "%*s %d %lf", &d_electroBCs[1].first,
             &d_electroBCs[1].second);
      flags_for_LoadPotential[1] += 1;
    }

    if (strstr(line, keywords[Y_BEGIN_E_BCS]) != NULL) {
      sscanf(line, "%*s %d %lf", &d_electroBCs[2].first,
             &d_electroBCs[2].second);
      flags_for_LoadPotential[1] += 1;
    }

    if (strstr(line, keywords[Y_END_E_BCS]) != NULL) {
      sscanf(line, "%*s %d %lf", &d_electroBCs[3].first,
             &d_electroBCs[3].second);
      flags_for_LoadPotential[1] += 1;
    }

    if (strstr(line, keywords[Z_BEGIN_E_BCS]) != NULL) {
      sscanf(line, "%*s %d %lf", &d_electroBCs[4].first,
             &d_electroBCs[4].second);
      flags_for_LoadPotential[1] += 1;
    }

    if (strstr(line, keywords[Z_END_E_BCS]) != NULL) {
      sscanf(line, "%*s %d %lf", &d_electroBCs[5].first,
             &d_electroBCs[5].second);
      flags_for_LoadPotential[1] += 1;
    }

    //
    //  read potential data
    //
    if (strstr(line, keywords[POTENTIAL]) != NULL) {
      sscanf(line, "%s %s %lf %d", dummy, d_potentialName,
             &d_electrostaticCutoffRadius, &d_electrostaticIntegration);

      if (flags_for_LoadPotential[0] != 1) {
        fprintf("%s %d: Need to read electrostaics flag before calling "
                "LoadPotential()\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }

      if (flags_for_LoadPotential[1] != 6) {
        fprintf("%s %d: Need to read electro boundary conditions before "
                "calling LoadPotential()\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }

      if (strcmp(d_potentialName, "MishinNiAl") == 0) {
        int a = 0;
        LoadPotential(a);
      }

      if (strcmp(d_potentialName, "Ceria") == 0) {
        int a = 1;
        LoadPotential(a);
      }

      if (strcmp(d_potentialName, "NiMn_LJ") == 0) {
        int a = 2;
        LoadPotential(a);
      }

      if (strcmp(d_potentialName, "GaNCoreShell") == 0) {
        int a = 3;
        LoadPotential(a);
      }

      if (strcmp(d_potentialName, "GaNCoreShellCorrect") == 0) {
        int a = 4;
        LoadPotential(a);
      }

      if (strcmp(d_potentialName, "PTO_Shimada") == 0) {
        int a = 5;
        LoadPotential(a);
      }

      if (strcmp(d_potentialName, "ArLJ") == 0) {
        int a = 6;
        LoadPotential(a);
      }
    }

    //
    // read residual flag
    //
    if (strstr(line, keywords[ALL_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[0]);
    }

    if (strstr(line, keywords[X_BEGIN_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[1]);
    }

    if (strstr(line, keywords[X_END_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[2]);
    }

    if (strstr(line, keywords[Y_BEGIN_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[3]);
    }

    if (strstr(line, keywords[Y_END_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[4]);
    }

    if (strstr(line, keywords[Z_BEGIN_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[5]);
    }

    if (strstr(line, keywords[Z_END_RESIDUAL]) != NULL) {
      sscanf(line, "%*s %d", &d_residualForceRemoveFlag[6]);
    }

    if (strstr(line, keywords[QUAD_METHOD]) != NULL) {
      sscanf(line, "%*s %d", &d_quadratureMethod);
    }

    //
    //  read tolerance and other data for CGSolver
    //
    if (strstr(line, keywords[TOLERANCE]) != NULL) {
      sscanf(line, "%*s %lf", &d_cgTolerance);
    }

    if (strstr(line, keywords[MAXITERATIONS]) != NULL) {
      sscanf(line, "%*s %d", &d_cgMaxIterations);
    }

    if (strstr(line, keywords[DEBUGLEVEL]) != NULL) {
      sscanf(line, "%*s %d", &d_debugLevel);
    }

    if (strstr(line, keywords[LINETOLERANCE]) != NULL) {
      sscanf(line, "%*s %lf", &d_cgLineSearchTolerance);
    }

    if (strstr(line, keywords[LINEITERATIONS]) != NULL) {
      sscanf(line, "%*s %d", &d_cgLineSearchMaxIterations);
    }

    //
    // read remesh tolerance
    //
    if (strstr(line, keywords[REMESHTOLERANCE]) != NULL) {
      double tolerance_remesh;
      sscanf(line, "%*s %lf", &d_remeshTolerance);

      // set remesh tolerance in CreateMesh
      CreateMesh::getInstance()->setRemeshTolerance(d_remeshTolerance);
    }

    //
    //  read location flag, but don't use it.
    //  New CrossNeighborList.cc doesn't have option of save_location_flag
    //
    if (strstr(line, keywords[CROSSNEIGHLISTFLAG]) != NULL) {
      int saveLocationsFlag;
      sscanf(line, "%*s %d", &saveLocationsFlag);
    }

    //
    // read get defects flag
    //
    if (strstr(line, keywords[GETDEFECTSFLAG]) != NULL) {
      sscanf(line, "%*s %d", &getDefects);
    }

    if (strstr(line, keywords[DEFECT_CENTER]) != NULL) {
      defectCenter.resize(3, 0.0);
      sscanf(line, "%*s %lf %lf %lf", &defectCenter[0], &defectCenter[1],
             &defectCenter[2]);
    }

    if (strstr(line, keywords[DEFECT_BOX]) != NULL) {
      defectBox.resize(3, 0.0);
      sscanf(line, "%*s %lf %lf %lf", &defectBox[0], &defectBox[1],
             &defectBox[2]);
    }

    if (strstr(line, keywords[DEFECT_LAT]) != NULL) {
      sscanf(line, "%*s %lf", &defectLat);
    }

    //
    //  read add node flag. This is used during remeshing.
    //  If sufficient nodes are not added after remeshing, we stop
    //  further remeshing.
    //
    if (strstr(line, keywords[ADDNODEFLAG]) != NULL) {
      sscanf(line, "%*s %d", &addNodeFlag);
    }

    //
    //  read add node restart num, i.e. max limit of remeshing
    //
    if (strstr(line, keywords[ADDNODERESTARTNUM]) != NULL) {
      sscanf(line, "%*s %d", &addNodeRestartNum);
    }

    //
    //  read add node cut off
    //
    if (strstr(line, keywords[ADDNODECUTOFF]) != NULL) {
      sscanf(line, "%*s %d", &addNodeCutoff);
    }
  } // end of loop over lines of quasi.ini

  //
  //
  //
  if (d_universalConst.first == 0) {
    d_universalConst.second[0] = 0.831445986;
    d_universalConst.second[1] = 39.903127109;
    d_universalConst.second[2] = 0.0055263488697;
  }

  //
  // set the values
  //
  numQuasi = d_numQuasi;
  numLoads = d_numLoads;
  boundaryFlag = d_externalLoadFlag;

  std::vector<double> zeros(3, 0.0);
  loadDeformation.clear();
  loadDeformation.push_back(zeros);
  loadDeformation.push_back(zeros);
  loadDeformation.push_back(zeros);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      loadDeformation[i][j] = d_deformationMatrix[i][j];

  tolerance = d_cgTolerance;
  maxIterations = d_cgMaxIterations;
  debugLevel = d_debugLevel;
  lineTolerance = d_cgLineSearchTolerance;
  lineIterations = d_cgLineSearchMaxIterations;

  charges.clear();
  if (d_externalCharges.size() != 0) {
    d_externalChargeFlag = 1;

    charges.resize(d_externalCharges.size());
    for (int i = 0; i < d_externalCharges.size(); i++) {
      charges[i].first = d_externalCharges[i].first;
      charges[i].second.resize(3);
      for (int dof = 0; dof < 3; dof++)
        charges[i].second[dof] = d_externalCharges[i].second[dof];
    }
  }

  temperature = d_temperature;
  quadMethod = d_quadratureMethod;

  testFlag = d_testFlag;
  minMethod = d_minMethod;
  minMethodMaxIter = d_minMethodMaxIter;
  minMethodTol = d_minMethodTol;

  initialConfigFlag = d_initialConfigurationFlag;

  //
  //  Also modify the PairPotentials data
  //
  PairPotentials::getInstance()->d_minMethod = d_minMethod;
  PairPotentials::getInstance()->d_quadMethod = d_quadratureMethod;
  PairPotentials::getInstance()->d_statistics = statistics;

  //
  //  put the residualForceFlags into ForceEnergyCalculation
  //
  if (d_residualForceRemoveFlag.size() != 0)
    ForceEnergyCalculation::getInstance()->setRemoveResidualFlags(
        d_residualForceRemoveFlag);

  // close init file
  if (fclose(init)) {
    fprintf(stderr, "%s %d: File closure failed.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // cleanup
  free(line);

  return;
} // end of mainInput()

//
// quasiInput()
//
struct qc_options_t Input::quasiInput(const int iQuasi,
                                      struct lattice_t *P_lattice,
                                      struct all_node_list_t *P_node_list,
                                      struct element_list_t *P_element_list,
                                      struct indentor_t *P_indentor,
                                      char *materialName, double &mass, int ac,
                                      char **av) {
  struct qc_options_t qc_options = QC_OPTIONS_INITIALIZER;

  int quasicontinuum_id = iQuasi + 1;

  FILE *init;

  char data_file[256];
  char restart_file[256];

  double initial_freq = 1.0;

  // temp int for number_threads. Note : We have read number_threads in
  //  mainInput() and there we have already called qc_thread_init()

  /* counter to adjust data file name*/
  int data_file_count = 0;
  int restart_file_count = 0;

  /* variable to hold size of data file*/
  int data_file_size;
  int restart_file_size;

  /* variable to hold number of digits, +1 spot for _*/
  int quasi_id_digit = 2;

  unsigned key_flag = 0;

  // system dependent variables
  long path_max = pathconf("/", _PC_PATH_MAX);

  if (path_max == -1L) {
    fprintf("%s %d: pathconf error.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  /**
    * in order to ensure proper placement of data number of threads is
    * necessary so reading should be done when the number of threads is known
    */

  // compute initial frequency from tau, Boltzman_constant, and temperature
  double atomicMass = 1.0;
  ReadAtomicMass(iQuasi, atomicMass);
  mass = atomicMass;

  initial_freq =
      sqrt(atomicMass * d_universalConst.second[0] * d_temperature) / d_tau;
  d_initialFrequency.push_back(initial_freq);

  int coreFlag = -1; // -1 : shell,    0 : core
  coreFlag = IsQuasiCore(iQuasi);

  //
  // read input file
  //
  if (d_qcOption.restart == OFF) {
    qc_options.restart = OFF;

    //
    //  create data_file for current quasicontinuum
    //
    if (quasicontinuum_id > 1) {
      /**
        * dummy strings for making new data_file string
        */
      char *str1 = NULL;
      char *str2 = NULL;

      /**
        * initialize data_file_count and data_file_size
        */
      data_file_count = 0;
      data_file_size = strlen(d_dataFile);

      /**
        * find where first '.' is reached
        */
      while (d_dataFile[data_file_count] != '.') {
        data_file_count++;

        if (data_file_count == data_file_size) {
          fprintf("%s %d: something about data files and too many instances\n.",
                  __FILE__, __LINE__);
          exit(EXIT_FAILURE);
        }
      }

      /**
        * allocate sizes for data
        */
      // if( (str1 = malloc(data_file_count * sizeof(char))) == NULL ){
      if ((str1 = (char *)malloc(path_max * sizeof(char))) == NULL) {
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      // if( (str2 = malloc((strlen(data_file) - data_file_count)*sizeof(har)))
      // == NULL ){
      if ((str2 = (char *)malloc(path_max * sizeof(char))) == NULL) {
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      /**
        * copy tail end of data_file string
        */
      strcpy(str2, &d_dataFile[data_file_count]);

      // printf("str2 : %s\n", str2);
      /**
        * copy beginning of data_file string
        */
      strncpy(str1, d_dataFile, data_file_count);
      str1[data_file_count] = '\0';

      // printf("str1 : %s\n", str1);
      /**
        * determine number of digits in quasicontinuum_id
        */
      if (quasicontinuum_id > 9) {
        quasi_id_digit++;
        if (quasicontinuum_id > 99) {
          quasi_id_digit++;
          if (quasicontinuum_id > 999) {
            ERROR("too many quasicontinua for data file input");
            exit(EXIT_FAILURE);
          }
        }
      }

      /**
        * write new data file name
        */
      sprintf(data_file, "%s_%i%s", str1, quasicontinuum_id, str2);

      // printf("data_file : %s\n", data_file);
      /**
        * free dummy strings
        */
      free(str1);
      free(str2);
    } else {
      strcpy(data_file, d_dataFile);
    }

    ReadQuasiData(iQuasi, P_lattice, P_node_list, materialName, data_file);

    /**
      * after reading the node data from input and restart file, put fix_w_mask,
      * initial frequency, temperature.
      */

    int fixity;
    fixity = FREE_MASK;
    fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

    int i_node;
    for (i_node = 0; i_node < P_node_list->node_list.number_nodes; i_node++) {
      struct node_t *P_node = P_node_list->node_list.nodes[i_node];

      if (coreFlag == -1)
        P_node->fix_w_mask = FREE_MASK; // for shell keep w free
      else
        P_node->fix_w_mask = FIX_W_MASK; // for core w is fixed

      P_node->initial_temperature = d_temperature;

      P_node->temperature = P_node->initial_temperature;

      P_node->initial_frequency = initial_freq;

      P_node->frequency = P_node->initial_frequency;

      P_node->initial_tau = d_tau;

      P_node->tau = d_tau;

      //
      // depending on boundary_fixity_flag modify the fixity
      //
      if (boundary_fixity_flag == 0) {
        // be default set if free and modify it only if
        // node satisfies any condition below
        P_node->fix_mask = FREE_MASK;

        // fix XYZ at z=0 plane
        if (P_node->l[2] == 0)
          P_node->fix_mask = FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
      } else if (boundary_fixity_flag == 1) {
        // be default set if free and modify it only if
        // node satisfies any condition below
        P_node->fix_mask = FREE_MASK;

        // fix z at z=0 plane, y at y=0 plane, x at x=0 plane
        if (P_node->l[0] == 0 && P_node->l[1] != 0 && P_node->l[2] != 0) {
          // x=0 plane
          P_node->fix_mask |= FIX_X_MASK;
        }

        if (P_node->l[1] == 0 && P_node->l[0] != 0 && P_node->l[2] != 0) {
          // y=0 plane
          P_node->fix_mask |= FIX_Y_MASK;
        }

        if (P_node->l[2] == 0 && P_node->l[0] != 0 && P_node->l[1] != 0) {
          // z=0 plane
          P_node->fix_mask |= FIX_Z_MASK;
        }

        if (P_node->l[0] == 0 && P_node->l[1] == 0 && P_node->l[2] != 0) {
          // intersection of x=0 and y=0 plane
          P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK;
        }

        if (P_node->l[0] == 0 && P_node->l[2] == 0 && P_node->l[1] != 0) {
          // intersection of x=0 and z=0 plane
          P_node->fix_mask |= FIX_X_MASK | FIX_Z_MASK;
        }

        if (P_node->l[1] == 0 && P_node->l[2] == 0 && P_node->l[0] != 0) {
          // intersection of z=0 and y=0 plane
          P_node->fix_mask |= FIX_Z_MASK | FIX_Y_MASK;
        }

        if (P_node->l[0] == 0 && P_node->l[1] == 0 && P_node->l[2] == 0) {
          // intersection of x=0, y=0, z=0 plane
          P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
        }
      } else if (boundary_fixity_flag == 2) {
        // all free
        P_node->fix_mask = FREE_MASK;
      }
    }

    /**
      * to process indentor position nodal coordinates are necessary
      */
    P_indentor->displacement_incr[0] = d_indentIncrement[0];
    P_indentor->displacement_incr[1] = d_indentIncrement[1];
    P_indentor->displacement_incr[2] = d_indentIncrement[2];

    P_indentor->initial_position[0] = d_indentInitialLocation[0];
    P_indentor->initial_position[1] = d_indentInitialLocation[1];
    P_indentor->initial_position[2] = d_indentInitialLocation[2];

    P_indentor->position[0] = P_indentor->initial_position[0];
    P_indentor->position[1] = P_indentor->initial_position[1];
    P_indentor->position[2] = P_indentor->initial_position[2];

    P_indentor->radius = d_indentRadius;
    P_indentor->constant = d_indentConstant;
  } else {
    qc_options.restart = ON;

    if (quasicontinuum_id > 1) {
      /**
       * dummy strings for making new data_file string
       */
      char *str1 = NULL;
      char *str2 = NULL;

      /**
       * initialize restart_file_count
       */
      restart_file_count = 0;
      restart_file_size = strlen(d_restartFile);

      /**
       * find where first '.' is reached
       */
      while (d_restartFile[restart_file_count] != '.') {
        restart_file_count++;

        if (restart_file_count == restart_file_size) {
          D_ERROR("Restart file input for multiple instances");
          exit(EXIT_FAILURE);
        }
      }

      /**
       * allocate sizes for data
       */
      if ((str1 = (char *)malloc(restart_file_count * sizeof(char))) == NULL) {
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      if ((str2 = (char *)malloc((strlen(d_restartFile) - restart_file_count) *
                                 sizeof(char))) == NULL) {
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      /**
       * copy tail end of restart_file string
       */
      strcpy(str2, &d_restartFile[restart_file_count]);

      /**
       * copy beginning of restart_file string
       */
      strncpy(str1, d_restartFile, restart_file_count);
      str1[restart_file_count] = '\0';

      /**
       * determine number of digits in quasicontinuum_id
       */
      if (quasicontinuum_id > 9) {
        quasi_id_digit++;
        if (quasicontinuum_id > 99) {
          quasi_id_digit++;
          if (quasicontinuum_id > 999) {
            D_ERROR("too many quasicontinua for restart_file input");
            exit(EXIT_FAILURE);
          }
        }
      }

      /**
       * write new restart_file name
       */
      sprintf(restart_file, "%s_%i%s", str1, quasicontinuum_id, str2);

      /**
       * free dummy strings
       */
      free(str1);
      free(str2);
    } // if quasicontinua_ id > 1
    else {
      strcpy(restart_file, d_restartFile);
    }

    //
    //  read restart file
    //
    ReadRestartFile(iQuasi, restart_file, P_node_list, P_element_list,
                    P_indentor, P_lattice, materialName);
  }

  /**
    * initialize indenter  position
    */

  if (qc_options.restart == OFF) {
    memcpy(&P_indentor->position[0], &P_indentor->initial_position[0],
           sizeof(double[3]));
  }
  /**
    * close init file
    */

  /**
    * cleanup
    */

  return (qc_options);
} // end of quasiInput()

//
// load potential
//
void Input::LoadPotential(const int potentialType) {
  //
  //  potential type
  //  0 : LJ
  //  1 : EAM Johnson
  //  2 : EAM data
  //  3 : LJ NiMn
  //  4 : Buckingham
  //  5 : Harmonic
  //  6 : Anharmonic
  //  7 : Rydberg (not yet added)
  //
  long double dielectricConstant =
      4.0 * M_PI * PairPotentials::getInstance()->getElectricConstant();

  //
  // switch between cases
  //
  switch (potentialType) {
  //
  // Mishin02NiAl
  //
  case 0: {
    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      for (int i = 0; i < d_latticeInfo.second.size(); i++) {
        if (d_latticeInfo.second[i].first == 0)
          d_latticeInfo.second[i].second.second = 1.0; // put Ni mass here
        else
          d_latticeInfo.second[i].second.second = 1.0; // put Al mass here
      }
    }

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    int nGridPoints;
    bool EAMFlag;

    //
    // populate first potential between Quasicontinuum 0 and itself
    //
    interactionData.first = 0;
    interactionData.second = 0;

    potentialType = 2;

    double cutoff_radius = 5.95414;

    parameterData.clear();
    parameterData.push_back(cutoff_radius); // cutoff radius

    nGridPoints = 301;

    EAMFlag = true;

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertEAMPotential(interactionData, potentialType,
                                       parameterData, "Mishin02_NiNi.eam",
                                       nGridPoints, EAMFlag);

    //
    // populate first potential between Quasicontinuum 0 and 1
    //
    interactionData.first = 0;
    interactionData.second = 1;

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertEAMPotential(interactionData, potentialType,
                                       parameterData, "Mishin02_NiAl.eam",
                                       nGridPoints, EAMFlag);

    //
    // populate first potential between Quasicontinuum 1 and 0
    //
    interactionData.first = 1;
    interactionData.second = 0;

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertEAMPotential(interactionData, potentialType,
                                       parameterData, "Mishin02_AlNi.eam",
                                       nGridPoints, EAMFlag);

    //
    // populate first potential between Quasicontinuum 1 and 1
    //
    interactionData.first = 1;
    interactionData.second = 1;

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertEAMPotential(interactionData, potentialType,
                                       parameterData, "Mishin02_AlAl.eam",
                                       nGridPoints, EAMFlag);

    //
    //  process electrostatics data if it's enabled
    //
    if (d_electrostaticFlag == 1) {
      //
      // insert charge data
      //
      Electrostatics *electrostaticsC = Electrostatics::getInstance();

      //
      // set charges
      //
      std::vector<double> charges;
      charges.push_back(-0.08);
      charges.push_back(0.08);

      //
      // insert data for electrostatic forces
      //
      const double elementSize = 0.35;
      const double integration_error = 0.001;
      const int num_shells = 8;

      // if using electric constant other than from quasi.ini or
      //  actual electric const, modify the d_universalConst data
      const double electricConstant = 1.0 / dielectricConstant;
      electrostaticsC->insertInputData(charges, d_electrostaticCutoffRadius,
                                       elementSize, electricConstant,
                                       d_electrostaticIntegration, d_electroBCs,
                                       integration_error, num_shells);
    }

    //
    //
    //
    return;

    break;
  }

  //
  // Ceria (Gotte 2007)
  //
  case 1: {
    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      for (int i = 0; i < d_latticeInfo.second.size(); i++) {
        if (d_latticeInfo.second[i].first % 2 == 0) {
          // core lattice
          if (d_latticeInfo.second[i].first < 12) {
            // Ce atom
            d_latticeInfo.second[i].second.second = 1.0; // put Ce mass here
          } else {
            // O atom
            d_latticeInfo.second[i].second.second = 1.0; // put O mass here
          }
        } else
          d_latticeInfo.second[i].second.second =
              1.0; // core-mass is irrelevant in present code
      }
    }

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    bool EAMFlag = false;

    //
    // populate potentials between cores and shells
    //
    for (int n = 0; n < 12; ++n) {
      interactionData.first = 2 * n;
      interactionData.second = 2 * n + 1;
      potentialType = 5;
      parameterData.clear();
      parameterData.push_back(1.0); // cutoff radius
      if (n < 4)
        parameterData.push_back(43.451); // k
      else
        parameterData.push_back(1759.8); // k
      parameterData.push_back(0.0);      // C

      //
      // insert 1st potential information into singleton
      //
      pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                            parameterData, EAMFlag);
    }

    //
    // populate potentials between shells of Ce4+ and O
    //
    for (int i = 0; i < 4; ++i) {
      interactionData.first = 2 * i;
      for (int j = 0; j < 8; ++j) {
        interactionData.second = 2 * j + 8;
        potentialType = 4;
        parameterData.clear();
        parameterData.push_back(14.0);         // cutoff radius
        parameterData.push_back(755.1311);     // A
        parameterData.push_back(0.429);        // rho
        parameterData.push_back(0.0);          // C
        parameterData.push_back(5.07273e-12);  // vC
        parameterData.push_back(-1.18246e-10); // vF

        //
        // insert 1st potential information into singleton
        //
        pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                              parameterData, EAMFlag);
      }
    }

    //
    // populate potentials between shells of O and O
    //
    for (int i = 8; i <= 22; i = i + 2) {
      interactionData.first = i;
      for (int j = i; j <= 22; j = j + 2) {
        interactionData.second = j;
        potentialType = 4;
        parameterData.clear();
        parameterData.push_back(14.0);          // cutoff radius
        parameterData.push_back(9533.421);      // A
        parameterData.push_back(0.234);         // rho
        parameterData.push_back(224.88);        // C
        parameterData.push_back(-0.0000298664); // vC
        parameterData.push_back(0.0000127999);  // vF

        //
        // insert 1st potential information into singleton
        //
        pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                              parameterData, EAMFlag);
      }
    }

    //
    //  process electrostatics information
    //
    if (d_electrostaticFlag == 1) {
      //
      // insert charge data
      //
      Electrostatics *electrostaticsC = Electrostatics::getInstance();

      //
      // set charges
      //
      std::vector<double> charges;
      charges.push_back(4.6475);
      charges.push_back(-0.6475);
      charges.push_back(4.6475);
      charges.push_back(-0.6475);
      charges.push_back(4.6475);
      charges.push_back(-0.6475);
      charges.push_back(4.6475);
      charges.push_back(-0.6475);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      charges.push_back(-6.5667);
      charges.push_back(4.5667);
      const double elementSize = 27.0;
      //
      // insert data for electrostatic forces
      //
      const double integration_error = 0.001;
      const int num_shells = 8;

      // if using electric constant other than from quasi.ini or
      //  actual electric const, modify the d_universalConst data
      const double electricConstant = 1.0 / dielectricConstant;
      electrostaticsC->insertInputData(charges, d_electrostaticCutoffRadius,
                                       elementSize, electricConstant,
                                       d_electrostaticIntegration, d_electroBCs,
                                       integration_error, num_shells);
    }

    //
    //
    //
    return;

    /* NOTREACHED */
    break;
  }

  //
  // NiMn_LJ
  //
  case 2: {
    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      for (int i = 0; i < d_latticeInfo.second.size(); i++) {
        if (d_latticeInfo.second[i].first == 0)
          d_latticeInfo.second[i].second.second = 1.0; // put Ni mass here
        else
          d_latticeInfo.second[i].second.second = 1.0; // put Mn mass here
      }
    }

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    bool EAMFlag = false;

    //
    // populate first potential between Quasicontinuum 0 and itself
    //
    interactionData.first = 0;
    interactionData.second = 0;

    double sigma = 2.2808;
    double epsilon = 0.518883157; // in unit eV

    double sigma_2 = sigma * sigma;
    double sigma_6 = sigma_2 * sigma_2 * sigma_2;
    double sigma_12 = sigma_6 * sigma_6;

    potentialType = 3;
    parameterData.clear();
    parameterData.push_back(8.8377200); // cutoff radius
    parameterData.push_back(sigma);     // sigma NiNi
    parameterData.push_back(sigma_6);   // sigma^6 NiNi
    parameterData.push_back(sigma_12);  // sigma^12 NiNi
    parameterData.push_back(epsilon);   // epsilon NiNi

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate second potential between Quasicontinuum 0 and 1
    //
    interactionData.first = 0;
    interactionData.second = 1;

    sigma = 2.30415;
    epsilon = 0.513931525; // in eV unit

    sigma_2 = sigma * sigma;
    sigma_6 = sigma_2 * sigma_2 * sigma_2;
    sigma_12 = sigma_6 * sigma_6;

    potentialType = 3;
    parameterData.clear();
    parameterData.push_back(8.8377200); // cutoff radius
    parameterData.push_back(sigma);     // sigma NiNi
    parameterData.push_back(sigma_6);   // sigma^6 NiNi
    parameterData.push_back(sigma_12);  // sigma^12 NiNi
    parameterData.push_back(epsilon);   // epsilon NiNi

    //
    // insert 2nd potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate third potential between Quasicontinuum 1 and 1
    //
    interactionData.first = 1;
    interactionData.second = 1;

    sigma = 2.3275;
    epsilon = 0.509026419; // in eV unit

    sigma_2 = sigma * sigma;
    sigma_6 = sigma_2 * sigma_2 * sigma_2;
    sigma_12 = sigma_6 * sigma_6;

    potentialType = 3;
    parameterData.clear();
    parameterData.push_back(8.8377200); // cutoff radius
    parameterData.push_back(sigma);     // sigma NiNi
    parameterData.push_back(sigma_6);   // sigma^6 NiNi
    parameterData.push_back(sigma_12);  // sigma^12 NiNi
    parameterData.push_back(epsilon);   // epsilon NiNi

    //
    // insert 3rd potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  process electrostatics data if it's enabled
    //
    if (d_electrostaticFlag == 1) {
      //
      // insert charge data
      //
      Electrostatics *electrostaticsC = Electrostatics::getInstance();

      //
      // set charges
      //
      std::vector<double> charges;
      charges.push_back(-1.0);
      charges.push_back(1.0);

      //
      // insert data for electrostatic forces
      //
      const double elementSize = 0.35;
      const double integration_error = 0.001;
      const int num_shells = 8;

      // if using electric constant other than from quasi.ini or
      //  actual electric const, modify the d_universalConst data
      const double electricConstant = 1.0; // or dielectricConstant
      d_universalConst.second[2] = 1.0 / electricConstant;

      electrostaticsC->insertInputData(charges, d_electrostaticCutoffRadius,
                                       elementSize, electricConstant,
                                       d_electrostaticIntegration, d_electroBCs,
                                       integration_error, num_shells);
    }

    return;

    /* NOTREACHED */
    break;
  }

  //
  // GaNCoreShell
  //
  case 3: {
    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      for (int i = 0; i < d_latticeInfo.second.size(); i++) {
        if (d_latticeInfo.second[i].first == 0 ||
            d_latticeInfo.second[i].first == 1)
          d_latticeInfo.second[i].second.second = 1.0; // put Ga mass here
        else
          d_latticeInfo.second[i].second.second = 1.0; // put N mass here
      }
    }

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    bool EAMFlag = false;

//
// start off with just shell-shell interactions
//
#if 0
        //
        // populate potentials between N cores and shells
        //
        for(n = 2; n <= 4; n = n + 2){
          interactionData.first = n;
          interactionData.second = n + 1;
          potentialType = 5;
          parameterData.clear();
          parameterData.push_back(1.0); // cutoff radius
          parameterData.push_back(11.7092); // k
          parameterData.push_back(0.0); // C

          //
          // insert 1st potential information into singleton
          //
          pairPotentials->insert(interactionData,
               potentialType,
               parameterData,
               EAMFlag);
        }
#endif
    //
    // populate potential between Ga shells
    //
    interactionData.first = 0;
    interactionData.second = 1;
    potentialType = 4;
    parameterData.clear();
    parameterData.push_back(10.0);     // cutoff radius
    parameterData.push_back(6068.14);  // A
    parameterData.push_back(0.31846);  // rho
    parameterData.push_back(250.0);    // C
    parameterData.push_back(-0.00025); // vC
    parameterData.push_back(0.00015);  // vF

    //
    // insert potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate potential between N shells
    //
    interactionData.first = 2;
    interactionData.second = 3;
    potentialType = 4;
    parameterData.clear();
    parameterData.push_back(10.0);     // cutoff radius
    parameterData.push_back(4115.42);  // A
    parameterData.push_back(0.31949);  // rho
    parameterData.push_back(280.0);    // C
    parameterData.push_back(-0.00028); // vC
    parameterData.push_back(0.000168); // vF

    //
    // insert potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate potentials between shells of Ga and N
    //
    potentialType = 4;
    parameterData.clear();
    parameterData.push_back(10.0);         // cutoff radius
    parameterData.push_back(872.42);       // A
    parameterData.push_back(0.31318);      // rho
    parameterData.push_back(0.0);          // C
    parameterData.push_back(1.184e-11);    // vC
    parameterData.push_back(-3.78167e-11); // vF

    //
    // insert potential information into singleton
    //
    interactionData.first = 0;
    interactionData.second = 2;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // insert potential information into singleton
    //
    interactionData.first = 0;
    interactionData.second = 3;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // insert potential information into singleton
    //
    interactionData.first = 1;
    interactionData.second = 2;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // insert potential information into singleton
    //
    interactionData.first = 1;
    interactionData.second = 3;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  process electrostatics data if it's enabled
    //
    if (d_electrostaticFlag == 1) {
      //
      // insert charge data
      //
      Electrostatics *electrostaticsC = Electrostatics::getInstance();

      //
      // set charges
      //
      std::vector<double> charges;
      charges.push_back(2.0);
      charges.push_back(2.0);
      charges.push_back(-2.0);
      charges.push_back(-2.0);

      //
      // insert data for electrostatic forces
      //
      const double elementSize = 8.0;
      const double integration_error = 0.001;
      const int num_shells = 8;

      // if using electric constant other than from quasi.ini or
      //  actual electric const, modify the d_universalConst data
      const double electricConstant = 14.39964; // or dielectricConstant

      electrostaticsC->insertInputData(charges, d_electrostaticCutoffRadius,
                                       elementSize, electricConstant,
                                       d_electrostaticIntegration, d_electroBCs,
                                       integration_error, num_shells);
    }

    //
    //
    //
    return;

    /* NOTREACHED */
    break;
  }

  //
  // GaNCoreShell with actual core shell
  //
  case 4: {
    // 0 - Ga core
    // 1 - Ga core
    // 2 - N core
    // 3 - N core
    // 4 - N shell
    // 5 - N shell

    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      for (int i = 0; i < d_latticeInfo.second.size(); i++) {
        if (d_latticeInfo.second[i].first == 0 ||
            d_latticeInfo.second[i].first == 1)
          d_latticeInfo.second[i].second.second = 1.0; // put Ga mass here
        else
          d_latticeInfo.second[i].second.second = 1.0; // put N mass here
      }
    }

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    bool EAMFlag = false;

    //
    // populate potentials between N cores and shells
    //
    for (int n = 2; n < 4; ++n) {
      interactionData.first = n;
      interactionData.second = n + 2;
      potentialType = 5;
      parameterData.clear();
      parameterData.push_back(0.5);     // cutoff radius
      parameterData.push_back(11.7092); // k
      parameterData.push_back(0.0);     // C

      //
      // insert 1st potential information into singleton
      //
      pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                            parameterData, EAMFlag);
    }

    //
    // populate potential between Ga cores
    //
    interactionData.first = 0;
    interactionData.second = 1;
    potentialType = 4;
    parameterData.clear();
    parameterData.push_back(10.0);     // cutoff radius
    parameterData.push_back(6068.14);  // A
    parameterData.push_back(0.31846);  // rho
    parameterData.push_back(250.0);    // C
    parameterData.push_back(-0.00025); // vC
    parameterData.push_back(0.00015);  // vF

    //
    // insert potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate potential between N shells
    //
    interactionData.first = 2;
    interactionData.second = 3;
    potentialType = 4;
    parameterData.clear();
    parameterData.push_back(10.0);     // cutoff radius
    parameterData.push_back(4115.42);  // A
    parameterData.push_back(0.31949);  // rho
    parameterData.push_back(280.0);    // C
    parameterData.push_back(-0.00028); // vC
    parameterData.push_back(0.000168); // vF

    //
    // insert potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate potentials between cores of Ga and N shell
    //
    potentialType = 4;
    parameterData.clear();
    parameterData.push_back(10.0);         // cutoff radius
    parameterData.push_back(872.42);       // A
    parameterData.push_back(0.31318);      // rho
    parameterData.push_back(0.0);          // C
    parameterData.push_back(1.184e-11);    // vC
    parameterData.push_back(-3.78167e-11); // vF

    //
    // insert potential information into singleton
    //
    interactionData.first = 0;
    interactionData.second = 2;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // insert potential information into singleton
    //
    interactionData.first = 0;
    interactionData.second = 3;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // insert potential information into singleton
    //
    interactionData.first = 1;
    interactionData.second = 2;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // insert potential information into singleton
    //
    interactionData.first = 1;
    interactionData.second = 3;
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    const double element_size = 8.0;
    //
    //  process electrostatics data if it's enabled
    //
    if (d_electrostaticFlag == 1) {
      //
      // insert charge data
      //
      Electrostatics *electrostaticsC = Electrostatics::getInstance();

      //
      // set charges
      //
      std::vector<double> charges;
      charges.push_back(2.0);
      charges.push_back(2.0);
      charges.push_back(-2.5);
      charges.push_back(-2.5);
      charges.push_back(0.5);
      charges.push_back(0.5);
#if 0
          charges.push_back(0.0);
          charges.push_back(0.0);
          charges.push_back(0.0);
          charges.push_back(0.0);
          charges.push_back(0.0);
          charges.push_back(0.0);
#endif
      //
      // insert data for electrostatic forces
      //
      const double integration_error = 0.001;
      const int num_shells = 8;

      // if using electric constant other than from quasi.ini or
      //  actual electric const, modify the d_universalConst data
      const double electricConstant = 14.39964; // or dielectricConstant

      electrostaticsC->insertInputData(charges, d_electrostaticCutoffRadius,
                                       element_size, electricConstant,
                                       d_electrostaticIntegration, d_electroBCs,
                                       integration_error, num_shells);
    }

    // set atomistic dead load forces
    // num_cells = 8
    const int num_quasi = 6;
    std::vector<std::vector<double>> forces;

    forces.resize(num_quasi);

    forces[0].push_back(-2.2321e-02);
    forces[0].push_back(-1.2887e-02);
    forces[0].push_back(-1.2492e-02);

    forces[1].push_back(8.5015e-03);
    forces[1].push_back(4.9083e-03);
    forces[1].push_back(-1.4966e-03);

    forces[2].push_back(1.2308e-03);
    forces[2].push_back(7.1058e-04);
    forces[2].push_back(2.7238e-02);

    forces[3].push_back(-3.7470e-02);
    forces[3].push_back(-2.1633e-02);
    forces[3].push_back(1.3656e-02);

    forces[4].push_back(-2.5439e-04);
    forces[4].push_back(-1.4687e-04);
    forces[4].push_back(-5.4459e-03);

    forces[5].push_back(7.4856e-03);
    forces[5].push_back(4.3218e-03);
    forces[5].push_back(-2.7295e-03);

    ForceEnergyCalculation *forceC = ForceEnergyCalculation::getInstance();

    forceC->setAtomisticDeadLoads(forces, element_size);

    //
    return;

    // NOTREACHED
    break;
  }

  //
  // PTO_Shimada (Shimada, Wakahara, Umeno and Kitamura 2008)
  //
  case 5: {
    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      for (int i = 0; i < d_latticeInfo.second.size(); i++) {
        if (d_latticeInfo.second[i].first % 2 == 0) {
          // shell lattice
          if (d_latticeInfo.second[i].first == 0)
            d_latticeInfo.second[i].second.second = 1.0; // put Pb mass here
          else if (d_latticeInfo.second[i].first == 2)
            d_latticeInfo.second[i].second.second = 1.0; // put Ti mass here
          else
            d_latticeInfo.second[i].second.second = 1.0; // put O mass here
        } else {
          // core lattice; atomic mass is immaterial for core lattice
          d_latticeInfo.second[i].second.second = 1.0;
        }
      }
    }

    //
    //  Species details
    //
    //  (0,1) : (shell, core) of Pb
    //  (2,3) : (shell, core) of Ti
    //  (4,5) : (shell, core) of O1
    //  (6,7) : (shell, core) of O2
    //  (8,9) : (shell, core) of O3
    //

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    bool EAMFlag = false;

    //
    // populate potentials between cores and shells
    //

    //
    //  Pb --> (shell, core)
    //
    interactionData.first = 0;
    interactionData.second = 1;
    potentialType = 6; // anharmonic
    parameterData.clear();
    parameterData.push_back(0.1);      // r_cut
    parameterData.push_back(154.8713); // k2
    parameterData.push_back(22416.67); // k4

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  Ti --> (shell, core)
    //
    interactionData.first = 2;
    interactionData.second = 3;
    potentialType = 6; // anharmonic
    parameterData.clear();
    parameterData.push_back(0.1);        // r_cut
    parameterData.push_back(8829.4096);  // k2
    parameterData.push_back(1928581.70); // k4

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  O1 --> (shell, core)
    //
    interactionData.first = 4;
    interactionData.second = 5;
    potentialType = 6; // anharmonic
    parameterData.clear();
    parameterData.push_back(0.1);      // r_cut
    parameterData.push_back(180.9134); // k2
    parameterData.push_back(6945.78);  // k4

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  O2 --> (shell, core)
    //
    interactionData.first = 6;
    interactionData.second = 7;
    potentialType = 6; // anharmonic
    parameterData.clear();
    parameterData.push_back(0.1);      // r_cut
    parameterData.push_back(180.9134); // k2
    parameterData.push_back(6945.78);  // k4

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  O3 --> (shell, core)
    //
    interactionData.first = 8;
    interactionData.second = 9;
    potentialType = 6; // anharmonic
    parameterData.clear();
    parameterData.push_back(0.1);      // r_cut
    parameterData.push_back(180.9134); // k2
    parameterData.push_back(6945.78);  // k4

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    // populate potentials between shells of Pb, Ti, O
    //

    //
    //  shells of Pb - Ti
    //
    interactionData.first = 0;
    interactionData.second = 2;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(387.7316);  // A
    parameterData.push_back(0.394957);  // rho
    parameterData.push_back(223.24409); // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of Pb - O1
    //
    interactionData.first = 0;
    interactionData.second = 4;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(2538.4110); // A
    parameterData.push_back(0.300698);  // rho
    parameterData.push_back(2.61676);   // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of Pb - O2
    //
    interactionData.first = 0;
    interactionData.second = 6;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(2538.4110); // A
    parameterData.push_back(0.300698);  // rho
    parameterData.push_back(2.61676);   // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of Pb - O3
    //
    interactionData.first = 0;
    interactionData.second = 8;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(2538.4110); // A
    parameterData.push_back(0.300698);  // rho
    parameterData.push_back(2.61676);   // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of Ti - O1
    //
    interactionData.first = 2;
    interactionData.second = 4;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(2555.2075); // A
    parameterData.push_back(0.278391);  // rho
    parameterData.push_back(2.25557);   // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of Ti - O2
    //
    interactionData.first = 2;
    interactionData.second = 6;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(2555.2075); // A
    parameterData.push_back(0.278391);  // rho
    parameterData.push_back(2.25557);   // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of Ti - O3
    //
    interactionData.first = 2;
    interactionData.second = 8;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(2555.2075); // A
    parameterData.push_back(0.278391);  // rho
    parameterData.push_back(2.25557);   // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of O1 - O2
    //
    interactionData.first = 4;
    interactionData.second = 6;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(1698.6653); // A
    parameterData.push_back(0.271756);  // rho
    parameterData.push_back(61.84354);  // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of O1 - O3
    //
    interactionData.first = 4;
    interactionData.second = 8;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(1698.6653); // A
    parameterData.push_back(0.271756);  // rho
    parameterData.push_back(61.84354);  // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  shells of O2 - O3
    //
    interactionData.first = 6;
    interactionData.second = 8;
    potentialType = 4; // buckingham
    parameterData.clear();
    parameterData.push_back(8.0);       // cutoff radius
    parameterData.push_back(1698.6653); // A
    parameterData.push_back(0.271756);  // rho
    parameterData.push_back(61.84354);  // C
    parameterData.push_back(0.0);       // vC -- for now keeping it zero
    parameterData.push_back(0.0);       // vF -- for now keeping it zero

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    //
    //  process electrostatics information
    //
    if (d_electrostaticFlag == 1) {
      //
      // insert charge data
      //
      Electrostatics *electrostaticsC = Electrostatics::getInstance();

      //
      // set charges
      //
      std::vector<double> charges;
      //  Pb Atom - (shell, core)
      charges.push_back(-3.63322);
      charges.push_back(5.49586);
      //  Ti Atom - (shell, core)
      charges.push_back(-16.27952);
      charges.push_back(19.36901);
      // O1 Atom - (shell, core)
      charges.push_back(-4.19917);
      charges.push_back(2.54843);
      // O2 Atom - (shell, core)
      charges.push_back(-4.19917);
      charges.push_back(2.54843);
      // O3 Atom - (shell, core)
      charges.push_back(-4.19917);
      charges.push_back(2.54843);

      const double elementSize = 27.0;
      //
      // insert data for electrostatic forces
      //
      const double integration_error = 0.001;
      const int num_shells = 8;
      const double electricConstant = 1.0 / dielectricConstant;
      electrostaticsC->insertInputData(charges, d_electrostaticCutoffRadius,
                                       elementSize, electricConstant,
                                       d_electrostaticIntegration, d_electroBCs,
                                       integration_error, num_shells);
    }

    //
    //
    //
    return;

    /* NOTREACHED */
    break;
  }

  //
  //  ArLJ : Argon Lennard-Jonnes
  //
  case 6: {
    //
    //  modify the atomic mass data if the flag is 0
    //
    if (d_latticeInfo.first == 0) {
      d_latticeInfo.second[0].second.second = 1.0; // put Ar mass here
    }

    //
    // get pair potentials singleton
    //
    PairPotentials *pairPotentials = PairPotentials::getInstance();

    //
    // create potential variables
    //
    std::pair<int, int> interactionData;
    int potentialType;
    std::vector<double> parameterData;
    bool EAMFlag = false;

    //
    // populate first potential between Quasicontinuum 0 and itself
    //
    interactionData.first = 0;
    interactionData.second = 0;

    double sigma = 3.4000;
    double epsilon = 0.0104;

    double sigma_2 = sigma * sigma;
    double sigma_6 = sigma_2 * sigma_2 * sigma_2;
    double sigma_12 = sigma_6 * sigma_6;

    potentialType = 3;
    parameterData.clear();
    parameterData.push_back(8.50000);  // cutoff radius
    parameterData.push_back(sigma);    // sigma NiNi
    parameterData.push_back(sigma_6);  // sigma^6 NiNi
    parameterData.push_back(sigma_12); // sigma^12 NiNi
    parameterData.push_back(epsilon);  // epsilon NiNi

    //
    // insert 1st potential information into singleton
    //
    pairPotentials->insertNonEAMPotential(interactionData, potentialType,
                                          parameterData, EAMFlag);

    return;

    /* NOTREACHED */
    break;
  }

  default:

    //
    // error
    //
    d_print("Exit: mainInput failure\n");
    exit(EXIT_FAILURE);
  }

  //
  //
  //
  return;
} // end of LoadPotential()

//
//  ProcessOptions()
//
void Input::ProcessOptions(char **init_file, char **data_file,
                           int *number_threads, enum switch_t *restart,
                           char **restart_file, double *temperature,
                           double *tau, int argc, char **argv) {
  // extern char *optarg;
  int c;

  // i: quasi init file name
  // d: input data file name
  // n: number of threads
  // r: restart flag and restart file name
  // t: temperature value
  // w: std deviation of position, i.e. tau, value

  while ((c = getopt(argc, argv, "i:d:n:r:t:w:")) != -1)
    switch (c) {

    case 'i':

      *init_file = optarg;
      printf("*************** init_file = %s\n", *init_file);
      break;

    case 'd':

      *data_file = optarg;
      break;

    case 'n':

      *number_threads = atoi(optarg);
      break;

    case 'r':
      *restart = ON;
      *restart_file = optarg;

      d_qcOption.restart = ON;
      break;

    case 't':

      *temperature = atof(optarg);
      break;

    case 'w':

      *tau = atof(optarg);
      break;

    case '?':

      ErrUsage(argv[0]);
      break;
    }
  return;
} // end of ProcessOptions()

//
// create default init filename
//
int Input::CreateDefaultInitFilename(char *init_file, size_t n, char *name,
                                     char const *ext) {

  char *base_name = basename(name);
  char local_base_name[_POSIX_PATH_MAX];
  char *ptr;

  /**
   * make local copy of base_name
   */

  strcpy(local_base_name, base_name);

  /**
   * trim local_base_name to the first "."
   */

  ptr = strchr(local_base_name, '.');

  if (ptr != NULL)
    *ptr = '\0';

  /**
   * check if there is enough space in the buffer
   */

  if (strlen(local_base_name) + strlen(ext) >= n)
    return (1);

  /**
   * fill in init_file
   */

  sprintf(init_file, "%s%s", local_base_name, ext);

  return (0);

} // end of CreateDefautlInitFilename()

//
//  read quasi node data, lattice data
//
void Input::ReadQuasiData(const int quasicontinuum_id,
                          struct lattice_t *P_lattice,
                          struct all_node_list_t *P_node_list,
                          char *materialName, char *input_filename) {

  FILE *input_file;

  XDR xdrs;

  struct node_t *P_node;

  int i_node = -1;

  char *name;

  /**
    * open data file
    */

  input_file = open_data_file(input_filename, "r");
  if (input_file == NULL) {
    ERROR(input_filename);
    exit(EXIT_FAILURE);
  }

  /**
    * create input XDR stream
    */
  xdrstdio_create(&xdrs, input_file, XDR_DECODE);

  printf("data file = %s\n", input_filename);

  /**
    * read material name
    */
  char name_name[256];

  //
  // if(xdr_string(&xdrs, &(name_name[0]), _POSIX_NAME_MAX) == FALSE)
  if (xdr_string(&xdrs, &materialName, _POSIX_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  /**
    * read lattice
    */

  name = P_lattice->name;
  if (xdr_string(&xdrs, &name, LATTICE_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  if (xdr_enum(&xdrs, (enum_t *)&(P_lattice->type)) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_lattice->a1[0]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a1[1]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a1[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_lattice->a2[0]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a2[1]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a2[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_lattice->a3[0]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a3[1]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a3[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_int(&xdrs, &P_lattice->l_start[0]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_start[1]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_start[2]) == FALSE)
    XDR_ERR("xdr_int()");

  if (xdr_int(&xdrs, &P_lattice->l_end[0]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_end[1]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_end[2]) == FALSE)
    XDR_ERR("xdr_int()");

  //
  //  update the data of Input
  //
  d_latticeConstants[0][0] = P_lattice->a1[0];
  d_latticeConstants[0][1] = P_lattice->a1[1];
  d_latticeConstants[0][2] = P_lattice->a1[2];

  d_latticeConstants[1][0] = P_lattice->a2[0];
  d_latticeConstants[1][1] = P_lattice->a2[1];
  d_latticeConstants[1][2] = P_lattice->a2[2];

  d_latticeConstants[2][0] = P_lattice->a3[0];
  d_latticeConstants[2][1] = P_lattice->a3[1];
  d_latticeConstants[2][2] = P_lattice->a3[2];

  d_latticeMinMax.first[0] = P_lattice->l_start[0];
  d_latticeMinMax.first[1] = P_lattice->l_start[1];
  d_latticeMinMax.first[2] = P_lattice->l_start[2];

  d_latticeMinMax.second[0] = P_lattice->l_end[0];
  d_latticeMinMax.second[1] = P_lattice->l_end[1];
  d_latticeMinMax.second[2] = P_lattice->l_end[2];

  /**
    * read number of nodes
    */

  if (xdr_int(&xdrs, &P_node_list->node_list.number_nodes) == FALSE)
    XDR_ERR("xdr_int()");

  /**
    * initialize P_node_list
    */
  memset((void *)&P_node_list->new_node_list, 0, sizeof(struct node_list_t));

  pthread_mutex_init(&P_node_list->node_list.lock, (pthread_mutexattr_t *)NULL);
  pthread_mutex_init(&P_node_list->new_node_list.lock,
                     (pthread_mutexattr_t *)NULL);

  /**
    * read in nodes
    */
  Node *nodeClass = Node::getInstance();

  for (i_node = 0; i_node < P_node_list->node_list.number_nodes; i_node++) {

    /**
      * create a node
      */
    nodeClass->makeActiveNode(&P_node_list->node_list, i_node);

    P_node = P_node_list->node_list.nodes[i_node];

    /**
      * read in lattice coordinates of the node
      */

    if (xdr_int(&xdrs, &(P_node->l[0])) == FALSE)
      XDR_ERR("xdr_int()");
    if (xdr_int(&xdrs, &(P_node->l[1])) == FALSE)
      XDR_ERR("xdr_int()");
    if (xdr_int(&xdrs, &(P_node->l[2])) == FALSE)
      XDR_ERR("xdr_int()");

    /**
      * compute atomic positions
      */
    Lattice *latticeClass = Lattice::getInstance();

    latticeClass->getSiteInitialPosition(P_node->initial_position, P_node->l,
                                         P_lattice, quasicontinuum_id);

    P_node->position[0] = P_node->initial_position[0];
    P_node->position[1] = P_node->initial_position[1];
    P_node->position[2] = P_node->initial_position[2];

    /**
      * read in node mask
      */

    if (xdr_u_int(&xdrs, &(P_node->fix_mask)) == FALSE)
      XDR_ERR("xdr_int()");

    /**
      * save node number
      */

    P_node->number = i_node;
  }

  /*
   * close
    */

  xdr_destroy(&xdrs);
  fclose(input_file);

  return;
} // end of ReadQuasiData()

//
//  display error in command line input usage
//
void Input::ErrUsage(char *pname) {

  fprintf(stderr, "usage: %s [-i<init file>] [-d<data file>] "
                  "[-m<materials file>] [-n<number threads>] "
                  "-r<restart file>\n",
          pname);
  exit(EXIT_FAILURE);
  return;
}

//
//
//
void Input::systemInfoOutput() {
  //
  //  open the file
  //
  FILE *file;
  char *filename;

  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);

  filename = (char *)malloc(path_max * sizeof(char));
  sprintf(filename, "%s/system.info", d_outputDirectory);

  file = fopen(filename, "w");

  //
  //  completing few remaining data of Input
  //
  if (strcmp(d_potentialName, "MishinNiAl") == 0) {
    sprintf(d_material, "NiAl");
    sprintf(d_potentialType, "eam");
    sprintf(d_latticeType, "BCC");
  }

  if (strcmp(d_potentialName, "Ceria") == 0) {
    sprintf(d_material, "GaN");
    sprintf(d_potentialType, "pairwise");
    sprintf(d_latticeType, "BCC");
  }

  if (strcmp(d_potentialName, "NiMn_LJ") == 0) {
    sprintf(d_material, "NiMn");
    sprintf(d_potentialType, "pairwise");
    sprintf(d_latticeType, "BCC");
  }

  if (strcmp(d_potentialName, "GaNCoreShell") == 0) {
    sprintf(d_material, "GaN");
    sprintf(d_potentialType, "pairwise");
    sprintf(d_latticeType, "BCC");
  }

  if (strcmp(d_potentialName, "GaNCoreShellCorrect") == 0) {
    sprintf(d_material, "GaN");
    sprintf(d_potentialType, "pairwise");
    sprintf(d_latticeType, "BCC");
  }

  if (strcmp(d_potentialName, "PTO_Shimada") == 0) {
    sprintf(d_material, "PbTiO3");
    sprintf(d_potentialType, "pairwise");
    sprintf(d_latticeType, "BCC");
  }

  if (strcmp(d_potentialName, "ArLJ") == 0) {
    sprintf(d_material, "Ar");
    sprintf(d_potentialType, "pairwise");
    sprintf(d_latticeType, "BCC");
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      d_configurationMinMax.first[i] +=
          d_latticeMinMax.first[j] * d_latticeConstants[i][j];

      d_configurationMinMax.second[i] +=
          d_latticeMinMax.second[j] * d_latticeConstants[i][j];
    }
  }

  //
  //  writing the data to file
  //

  //
  //  summmary
  //
  fprintf(file, "Summary: ");

  fprintf(file, "%s ", d_material);

  fprintf(file, "%d lattices ", d_numQuasi);

  fprintf(file, "%s ", d_potentialType);

  fprintf(file, "%s ", d_potentialName);

  if (d_tempFlag == 0)
    fprintf(file, "zero temperature ");
  else
    fprintf(file, "finite temperature ");

  if (d_electrostaticFlag == 0)
    fprintf(file, "no electrostatics ");
  else
    fprintf(file, "electrostatics ");

  if (d_electrostaticFlag == 1) {
    int i1, i2, i3;
    i1 = 0;
    i2 = 0;
    i3 = 0;

    for (int j = 0; j < 6; j++) {
      if (d_electroBCs[j].first == 0)
        i1 = 1;
      else if (d_electroBCs[j].first == 1)
        i2 = 1;
      else if (d_electroBCs[j].first == 2)
        i3 = 1;
    }

    if (i1 == 1 && i2 == 0 && i3 == 0)
      fprintf(file, "charge neutral bc ");

    if (i2 == 1 && i1 == 0 && i3 == 0)
      fprintf(file, "free surface bc ");

    if (i3 == 1 && i1 == 0 && i1 == 0)
      fprintf(file, "specified charge density bc ");

    if (i1 == 1 && i2 == 1 && i3 == 0)
      fprintf(file, "charge neutral and free surface mixed bc ");

    if (i2 == 1 && i3 == 1 && i2 == 0)
      fprintf(file, "charge neutral and specified charge density mixed bc ");

    if (i1 == 1 && i3 == 1 && i2 == 0)
      fprintf(file, "free surface and specified charge density mixed bc ");

    if (i1 == 1 && i2 == 1 && i3 == 1)
      fprintf(file, "charge neutral, free surface and specified charge density "
                    "mixed bc ");

    fprintf(file, "with ");

    if (d_externalCharges.size() == 0)
      fprintf(file, "no external charge ");
    else
      fprintf(file, "external charge ");
  }

  fprintf(file, "subjected to ");

  if (d_externalLoadFlag == -1)
    fprintf(file, "no external deformation ");
  else if (d_externalLoadFlag == 0)
    fprintf(file, "uniform external deformation to all atoms ");
  else
    fprintf(file, "external deformation at surface ");

  if (d_voidFlag == 0)
    fprintf(file, "no void defect ");
  else
    fprintf(file, "with void defect ");

  if (d_indentFlag == 0)
    fprintf(file, "subject to no indentation ");
  else
    fprintf(file, "subject to indentation ");

  //
  //  Material Info
  //
  int counter = 1;

  fprintf(file, "\n");
  fprintf(file, "\n");
  fprintf(file, "%d. Material Info \n", counter);

  fprintf(file, " Material                        %s\n", d_material);
  fprintf(file, " Potential                       %s\n", d_potentialName);
  fprintf(file, " Potential Type                  %s\n", d_potentialType);

  counter++;

  //
  //  restart info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Restart Info \n", counter);
  if (d_qcOption.restart == OFF) {
    fprintf(file, " Status                          OFF\n");
  } else {
    fprintf(file, " Status                          ON\n");
    fprintf(file, " Restart File                    %s\n", d_restartFile);
  }

  counter++;

  //
  //  Minimization Info
  //
  // printf("d_minMethod == %d\n", d_minMethod);
  fprintf(file, "\n");
  fprintf(file, "%d. Minimization Info \n", counter);
  if (d_minMethod == 0) {
    fprintf(file, " Minimization method              Minimization wrt u and w "
                  "simultaneously\n");
  } else if (d_minMethod == 1) {
    fprintf(file, " Minimization method                 Minimization wrt u and "
                  "w alternatively\n");
    fprintf(file, " Minimization Method MaxIter                 %d\n",
            d_minMethodMaxIter);
    fprintf(file, " Minimization Tolerance                      %f\n",
            d_minMethodTol);
  } else if (d_minMethod == 2) {
    fprintf(file, " Minimization method                 Minimization wrt w "
                  "keeping u fixed at initial state\n");
  } else if (d_minMethod == 3) {
    fprintf(file, " Minimization method                 Minimization wrt u "
                  "keeping w fixed at initial state\n");
  } else {
    fprintf(file, "Check minimization method flag = %d \n", d_minMethod);
    exit(EXIT_FAILURE);
  }

  counter++;

  //
  //  Geometrical Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Geometrical Info \n", counter);

  fprintf(file, " Number of Lattices                  %d\n", d_numQuasi);

  fprintf(file, " Lattice Type                        %s\n", d_latticeType);

  fprintf(file, " Lattice Constant \n");
  fprintf(file, "   %lf     %lf     %lf\n", d_latticeConstants[0][0],
          d_latticeConstants[0][1], d_latticeConstants[0][2]);
  fprintf(file, "   %lf     %lf     %lf\n", d_latticeConstants[1][0],
          d_latticeConstants[1][1], d_latticeConstants[1][2]);
  fprintf(file, "   %lf     %lf     %lf\n", d_latticeConstants[2][0],
          d_latticeConstants[2][1], d_latticeConstants[2][2]);

  fprintf(file, " Lattice Size \n");
  fprintf(file, "   l_min =  %d    %d    %d\n", d_latticeMinMax.first[0],
          d_latticeMinMax.first[1], d_latticeMinMax.first[2]);
  fprintf(file, "   l_max = %d    %d    %d\n", d_latticeMinMax.second[0],
          d_latticeMinMax.second[1], d_latticeMinMax.second[2]);

  fprintf(file, " Initial Configuration Size\n");
  fprintf(file, "   x_min = %f        %f        %f\n",
          d_configurationMinMax.first[0], d_configurationMinMax.first[1],
          d_configurationMinMax.first[2]);
  fprintf(file, "   x_max = %f        %f        %f\n",
          d_configurationMinMax.second[0], d_configurationMinMax.second[1],
          d_configurationMinMax.second[2]);

  fprintf(file, " Lattice Shift\n");
  fprintf(file, "     Lattice Number      shift_x       shift_y       shift_z  "
                "     core       atomic mass\n");
  for (int i = 0; i < d_shift.size(); i++)
    fprintf(
        file,
        "     %d                  %f        %f       %f        %d        %lf\n",
        i + 1, d_shift[i][0], d_shift[i][1], d_shift[i][2],
        d_latticeInfo.second[i].second.first,
        d_latticeInfo.second[i].second.second);

  counter++;

  //
  //  Geometrical Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Configuration Info \n", counter);
  if (d_initialConfigurationFlag == 0)
    fprintf(file, " Initial Configuration                 Crystal lattice "
                  "configuration \n");
  else if (d_initialConfigurationFlag == 1)
    fprintf(file, " Initial Configuration                 Equilibrated from "
                  "crystal lattice configuration \n");
  else
    fprintf(
        file,
        " Initial Configuration                 Specified by input file \n");

  counter++;

  //
  //  universal constants
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Universal Constants \n", counter);
  fprintf(file, " Boltzman Constant               %f\n",
          d_universalConst.second[0]);
  fprintf(file, " MaxPlanck Constant              %f\n",
          d_universalConst.second[1]);
  fprintf(file, " Electric Constant               %f\n",
          d_universalConst.second[2]);

  counter++;

  //
  //  Temperature Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Temperature Info \n", counter);
  if (d_tempFlag == 0)
    fprintf(file, " Status                          zero temperature\n");
  else {
    fprintf(file, " Status                          finite temperature\n");
    fprintf(file, " temperature                     %f\n", d_temperature);
    fprintf(file, " tau                             %f\n", d_tau);

    fprintf(file, " Initial Frequency               ");

    for (int i = 0; i < d_numQuasi; i++)
      fprintf(file, "%f   ", d_initialFrequency[i]);
    fprintf(file, "\n");

    fprintf(file, " Quadrature Method               %d\n", d_quadratureMethod);
  }

  counter++;

  //
  //  Electrostatics Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Electrostatics Info \n", counter);
  if (d_electrostaticFlag == 0)
    fprintf(file, " Status                            no\n");
  else {
    fprintf(file, " Status                                    yes\n");

    if (d_externalCharges.size() != 0) {
      fprintf(file, " External Charges                          yes\n");
      fprintf(file, "     Charge Number       increment         x           y  "
                    "          z\n");
      for (int i = 0; i < d_externalCharges.size(); i++)
        fprintf(file, "     %d                    %f    %f   %f      %f\n",
                i + 1, d_externalCharges[i].first,
                d_externalCharges[i].second[0], d_externalCharges[i].second[1],
                d_externalCharges[i].second[2]);
    }
    fprintf(file, " Boundary Condition\n");
    fprintf(file, "     xBegin Surface              %d      %f\n",
            d_electroBCs[0].first, d_electroBCs[0].second);

    fprintf(file, "     xEnd Surface                %d      %f\n",
            d_electroBCs[1].first, d_electroBCs[1].second);

    fprintf(file, "     yBegin Surface              %d      %f\n",
            d_electroBCs[2].first, d_electroBCs[2].second);

    fprintf(file, "     yEnd Surface                %d      %f\n",
            d_electroBCs[3].first, d_electroBCs[3].second);

    fprintf(file, "     zBegin Surface              %d      %f\n",
            d_electroBCs[4].first, d_electroBCs[4].second);

    fprintf(file, "     zEnd Surface                %d      %f\n",
            d_electroBCs[5].first, d_electroBCs[5].second);
  }

  counter++;

  //
  //  Loading Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Loading Info \n", counter);
  if (d_externalLoadFlag == -1)
    fprintf(file, " Status                            no\n");
  else {
    fprintf(file, " Status                            yes\n");

    fprintf(file, " Number of Load Steps              %d\n", d_numLoads);

    if (d_externalLoadFlag == 0)
      fprintf(file, " Bulk or Surface                   bulk\n");
    else
      fprintf(file, " Bulk or Surface                   surface\n");

    fprintf(file, " Deformation Matrix \n");
    for (int i = 0; i < 3; i++)
      fprintf(file, "   %lf     %lf     %lf\n", d_deformationMatrix[i][0],
              d_deformationMatrix[i][1], d_deformationMatrix[i][2]);
  }

  counter++;

  //
  //  Void Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Void Info \n", counter);
  if (d_voidFlag == 0)
    fprintf(file, " Status                      no\n");
  else {
    fprintf(file, " Status                      yes\n");

    if (d_voidType[0] == 'C')
      fprintf(file, " Void Type                 Circle\n");
    else
      fprintf(file, " Void Type                 VPITT\n");

    fprintf(file, " Location                    %f  %f  %f\n", d_voidCenter[0],
            d_voidCenter[1], d_voidCenter[2]);
  }

  counter++;

  //
  //  Indent Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Indent Info \n", counter);
  if (d_voidFlag == 0)
    fprintf(file, " Status                      no\n");
  else {
    fprintf(file, " Status                      yes\n");
    fprintf(file, " Location                    %f %f %f \n",
            d_indentInitialLocation[0], d_indentInitialLocation[1],
            d_indentInitialLocation[2]);
    fprintf(file, " Location                    %f %f %f \n",
            d_indentIncrement[0], d_indentIncrement[1], d_indentIncrement[2]);
    fprintf(file, " Indent Radius               %f\n", d_indentRadius);
  }

  counter++;

  //
  //  Residual Force Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Residual Force Remove Info \n", counter);
  fprintf(file,
          " All     xBegin     xEnd     yBegin     yEnd     zBegin     zEnd\n");
  if (d_residualForceRemoveFlag[0] == 0)
    fprintf(file, " no");
  else
    fprintf(file, " yes");

  if (d_residualForceRemoveFlag[1] == 0)
    fprintf(file, "       no");
  else
    fprintf(file, "       yes");

  if (d_residualForceRemoveFlag[2] == 0)
    fprintf(file, "         no");
  else
    fprintf(file, "         yes");

  if (d_residualForceRemoveFlag[3] == 0)
    fprintf(file, "        no");
  else
    fprintf(file, "        yes");

  if (d_residualForceRemoveFlag[4] == 0)
    fprintf(file, "        no");
  else
    fprintf(file, "        yes");

  if (d_residualForceRemoveFlag[5] == 0)
    fprintf(file, "        no");
  else
    fprintf(file, "        yes");

  if (d_residualForceRemoveFlag[6] == 0)
    fprintf(file, "         no\n");
  else
    fprintf(file, "         yes\n");

  counter++;

  //
  //  CG Info
  //
  fprintf(file, "\n");
  fprintf(file, "%d. CG Info \n", counter);
  fprintf(file, " Tolerance                             %lf\n", d_cgTolerance);
  fprintf(file, " Max Iterations                        %d\n",
          d_cgMaxIterations);
  fprintf(file, " Line Search Tolerance                 %lf\n",
          d_cgLineSearchTolerance);
  fprintf(file, " Line Search Max Iterations            %d\n",
          d_cgLineSearchMaxIterations);
  fprintf(file, " Remesh Tolerance                      %lf\n",
          d_remeshTolerance);

  counter++;

  //
  //  number of threads
  //
  fprintf(file, "\n");
  fprintf(file, "%d. Thread Info \n", counter);
  fprintf(file, " Number of Threads                     %d", d_numberThreads);

  counter++;

  //
  //  done
  //
  fclose(file);
  return;
}

//
//  ReadRestartFile()
//
void Input::ReadRestartFile(int iQuasi, const char *filename,
                            struct all_node_list_t *P_all_node_list,
                            struct element_list_t *P_element_list,
                            struct indentor_t *P_indentor,
                            struct lattice_t *P_lattice, char *materialName) {
  FILE *file;
  XDR xdrs;

  struct node_list_t *P_node_list = &(P_all_node_list->node_list);

  /**
   * prepare xdr stream for reading
   */

  // d_print("restart filename = %s\n", filename);
  if ((file = open_data_file(filename, "r")) == NULL) {
    D_ERROR(filename);
    exit(EXIT_FAILURE);
  }

  xdrstdio_create(&xdrs, file, XDR_DECODE);

  //
  //  read material name
  //
  if (xdr_string(&xdrs, &materialName, _POSIX_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  /**
   * read lattice data
   */
  char *name = P_lattice->name;

  if (xdr_string(&xdrs, &name, LATTICE_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  if (xdr_enum(&xdrs, (enum_t *)&(P_lattice->type)) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_lattice->a1[0]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a1[1]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a1[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_lattice->a2[0]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a2[1]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a2[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_lattice->a3[0]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a3[1]) == FALSE)
    XDR_ERR("xdr_double()");
  if (xdr_double(&xdrs, &P_lattice->a3[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_int(&xdrs, &P_lattice->l_start[0]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_start[1]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_start[2]) == FALSE)
    XDR_ERR("xdr_int()");

  if (xdr_int(&xdrs, &P_lattice->l_end[0]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_end[1]) == FALSE)
    XDR_ERR("xdr_int()");
  if (xdr_int(&xdrs, &P_lattice->l_end[2]) == FALSE)
    XDR_ERR("xdr_int()");

  //
  //  update local data of Input class
  //
  d_latticeConstants[0][0] = P_lattice->a1[0];
  d_latticeConstants[0][1] = P_lattice->a1[1];
  d_latticeConstants[0][2] = P_lattice->a1[2];

  d_latticeConstants[1][0] = P_lattice->a2[0];
  d_latticeConstants[1][1] = P_lattice->a2[1];
  d_latticeConstants[1][2] = P_lattice->a2[2];

  d_latticeConstants[2][0] = P_lattice->a3[0];
  d_latticeConstants[2][1] = P_lattice->a3[1];
  d_latticeConstants[2][2] = P_lattice->a3[2];

  d_latticeMinMax.first[0] = P_lattice->l_start[0];
  d_latticeMinMax.first[1] = P_lattice->l_start[1];
  d_latticeMinMax.first[2] = P_lattice->l_start[2];

  d_latticeMinMax.second[0] = P_lattice->l_end[0];
  d_latticeMinMax.second[1] = P_lattice->l_end[1];
  d_latticeMinMax.second[2] = P_lattice->l_end[2];

  //
  //  read node data
  //
  int i_node;

#if defined(__QC_SGI) || defined(__QC_HPUX)
  size_t allocation_size;
  char *allocation;
#endif /* sgi */

  /**
   * initialize node list data
   */
  memset((void *)&P_all_node_list->new_node_list, 0,
         sizeof(struct node_list_t));

  pthread_mutex_init(&P_all_node_list->node_list.lock, NULL);
  pthread_mutex_init(&P_all_node_list->new_node_list.lock, NULL);

  /**
   * read number of nodes
   */
  if (xdr_int(&xdrs, &P_node_list->number_nodes) == 0)
    XDR_ERR("xdr_int");

#if defined(__QC_SGI) || defined(__QC_HPUX)

/**
  * allocate enough pages for both pointer to nodes and nodes; place
  * them and free
  */
//    allocation_size =
//   P_node_list->number_nodes*(sizeof(struct node_t*)+sizeof(struct node_t));

//  allocation=calloc_threaded(allocation_size, sizeof(char));
//  if (allocation == NULL)
//  {
//    D_ERROR("malloc()");
//    exit(EXIT_FAILURE);
//  }

// free(allocation);

#endif /* sgi */

  //
  //  loop over all nodes
  //

  // get Node instance
  Node *nodeC = Node::getInstance();
  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {

    //
    //  create a node
    //
    struct node_t *P_node;

    nodeC->makeActiveNode(P_node_list, i_node);

    P_node = P_node_list->nodes[i_node];

    P_node_list->nodes[i_node]->number = i_node;

    /**
     * initial position, initial frequency and initial temperature
     */
    if (xdr_double(&xdrs, &P_node->initial_position[0]) == FALSE ||
        xdr_double(&xdrs, &P_node->initial_position[1]) == FALSE ||
        xdr_double(&xdrs, &P_node->initial_position[2]) == FALSE ||
        xdr_double(&xdrs, &P_node->initial_frequency) == FALSE ||
        xdr_double(&xdrs, &P_node->initial_temperature) == FALSE ||
        xdr_double(&xdrs, &P_node->initial_tau) == FALSE)
      XDR_ERR("xdr_double");

    /**
     * position, frequency and temperature
     */
    if (xdr_double(&xdrs, &P_node->position[0]) == FALSE ||
        xdr_double(&xdrs, &P_node->position[1]) == FALSE ||
        xdr_double(&xdrs, &P_node->position[2]) == FALSE ||
        xdr_double(&xdrs, &P_node->frequency) == FALSE ||
        xdr_double(&xdrs, &P_node->temperature) == FALSE ||
        xdr_double(&xdrs, &P_node->tau) == FALSE)
      XDR_ERR("xdr_double");

    /**
     * fixity mask and fixity_w_mask
     */
    if (xdr_u_int(&xdrs, &P_node->fix_mask) == FALSE ||
        xdr_u_int(&xdrs, &P_node->fix_w_mask) == FALSE)
      XDR_ERR("xdr_u_int");

    /**
     * lattice site
     */
    if (xdr_int(&xdrs, &(P_node->l[0])) == FALSE)
      XDR_ERR("xdr_int()");
    if (xdr_int(&xdrs, &(P_node->l[1])) == FALSE)
      XDR_ERR("xdr_int()");
    if (xdr_int(&xdrs, &(P_node->l[2])) == FALSE)
      XDR_ERR("xdr_int()");

    //
    // depending on boundary_fixity_flag modify the fixity
    //
    if (boundary_fixity_flag == 0) {
      // be default set if free and modify it only if
      // node satisfies any condition below
      P_node->fix_mask = FREE_MASK;

      // fix XYZ at z=0 plane
      if (P_node->l[2] == 0)
        P_node->fix_mask = FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
    } else if (boundary_fixity_flag == 1) {
      // be default set if free and modify it only if
      // node satisfies any condition below
      P_node->fix_mask = FREE_MASK;

      // fix z at z=0 plane, y at y=0 plane, x at x=0 plane
      if (P_node->l[0] == 0 && P_node->l[1] != 0 && P_node->l[2] != 0) {
        // x=0 plane
        P_node->fix_mask |= FIX_X_MASK;
      }

      if (P_node->l[1] == 0 && P_node->l[0] != 0 && P_node->l[2] != 0) {
        // y=0 plane
        P_node->fix_mask |= FIX_Y_MASK;
      }

      if (P_node->l[2] == 0 && P_node->l[0] != 0 && P_node->l[1] != 0) {
        // z=0 plane
        P_node->fix_mask |= FIX_Z_MASK;
      }

      if (P_node->l[0] == 0 && P_node->l[1] == 0 && P_node->l[2] != 0) {
        // intersection of x=0 and y=0 plane
        P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK;
      }

      if (P_node->l[0] == 0 && P_node->l[2] == 0 && P_node->l[1] != 0) {
        // intersection of x=0 and z=0 plane
        P_node->fix_mask |= FIX_X_MASK | FIX_Z_MASK;
      }

      if (P_node->l[1] == 0 && P_node->l[2] == 0 && P_node->l[0] != 0) {
        // intersection of z=0 and y=0 plane
        P_node->fix_mask |= FIX_Z_MASK | FIX_Y_MASK;
      }

      if (P_node->l[0] == 0 && P_node->l[1] == 0 && P_node->l[2] == 0) {
        // intersection of x=0, y=0, z=0 plane
        P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
      }
    }
  }

  //
  //  dump element data
  //
  int i_elem;

  //
  //  read number of elements
  //
  if (xdr_int(&xdrs, &P_element_list->number_elements) == FALSE)
    XDR_ERR("xdr_int()");

#if defined(__QC_SGI) || defined(__QC_HPUX)

// allocation_size=
//   P_element_list->number_elements*(sizeof(struct element_t *)+
//         sizeof(struct element_t));
// allocation=calloc_threaded(allocation_size, sizeof(char));

// free(allocation);
#endif /* sgi */

  /**
   * loop over all elements
   */
  Element *elementC = Element::getInstance();
  MiscFunctions *miscC = MiscFunctions::getInstance();

  for (i_elem = 0; i_elem < P_element_list->number_elements; i_elem++) {
    struct element_t *P_element;

    double center[3];
    int tet[4];

    elementC->makeActiveElement(P_element_list, i_elem);

    P_element = P_element_list->elements[i_elem];

    /**
     * element_number
     */
    if (xdr_int(&xdrs, &P_element->number) == FALSE)
      XDR_ERR("xdr_int()");

    /**
     * element node numbers
     */
    if (xdr_int(&xdrs, &tet[0]) == FALSE || xdr_int(&xdrs, &tet[1]) == FALSE ||
        xdr_int(&xdrs, &tet[2]) == FALSE || xdr_int(&xdrs, &tet[3]) == FALSE)
      XDR_ERR("xdr_int()");

    /**
     * initialize element
     */
    P_element_list->elements[i_elem]->mask = ELEM_CLEAN_MASK;
    P_element_list->elements[i_elem]->number = i_elem;

    P_element_list->elements[i_elem]->node[0] = P_node_list->nodes[tet[0]];
    P_element_list->elements[i_elem]->node[1] = P_node_list->nodes[tet[1]];
    P_element_list->elements[i_elem]->node[2] = P_node_list->nodes[tet[2]];
    P_element_list->elements[i_elem]->node[3] = P_node_list->nodes[tet[3]];

    P_element_list->elements[i_elem]->cir_radius = miscC->findTetraCenter(
        P_node_list->nodes[tet[0]]->initial_position,
        P_node_list->nodes[tet[1]]->initial_position,
        P_node_list->nodes[tet[2]]->initial_position,
        P_node_list->nodes[tet[3]]->initial_position, center);

    P_element_list->elements[i_elem]->center[0] = center[0];
    P_element_list->elements[i_elem]->center[1] = center[1];
    P_element_list->elements[i_elem]->center[2] = center[2];
  }

  //
  //  read indentor data
  //
  if (xdr_double(&xdrs, &P_indentor->radius) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_indentor->constant) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_indentor->displacement_incr[0]) == FALSE ||
      xdr_double(&xdrs, &P_indentor->displacement_incr[1]) == FALSE ||
      xdr_double(&xdrs, &P_indentor->displacement_incr[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_indentor->initial_position[0]) == FALSE ||
      xdr_double(&xdrs, &P_indentor->initial_position[1]) == FALSE ||
      xdr_double(&xdrs, &P_indentor->initial_position[2]) == FALSE)
    XDR_ERR("xdr_double()");

  if (xdr_double(&xdrs, &P_indentor->position[0]) == FALSE ||
      xdr_double(&xdrs, &P_indentor->position[1]) == FALSE ||
      xdr_double(&xdrs, &P_indentor->position[2]) == FALSE)
    XDR_ERR("xdr_double()");

  //
  //  done
  //
  xdr_destroy(&xdrs);
  if (fclose(file)) {
    ERROR(filename);
    exit(EXIT_SUCCESS);
  }

  return;
}

//
//  ReadAtomicMass()
//
void Input::ReadAtomicMass(int iQuasi, double &mass) {
  // read from Input data d_mass
  if (iQuasi < 0 || iQuasi >= d_latticeInfo.second.size()) {
    d_print("Error in ReadAtomicMass()\n");
    exit(EXIT_FAILURE);
  }

  int found = -1;
  for (int i = 0; i < d_latticeInfo.second.size(); i++) {
    if (d_latticeInfo.second[i].first == iQuasi) {
      // found iQuasi
      mass = d_latticeInfo.second[i].second.second;
      found = 0;
    }
  }

  if (found == -1) {
    d_print("iQuasi = %d not found in d_mass. Function ReadAtomicMass().\n",
            iQuasi);
    exit(EXIT_FAILURE);
  }

  return;
}

//
//  IsQuasiCore()
//
int Input::IsQuasiCore(int iQuasi) {
  // loop over d_latticeInfo
  for (int i = 0; i < d_latticeInfo.second.size(); i++) {
    if (d_latticeInfo.second[i].first == iQuasi) {
      // iQuasi found
      if (d_latticeInfo.second[i].second.first == -1) {
        // it's a shell quasi
        return -1;
      } else {
        // it's a core quasi
        return 0;
      }
    }
  }

  // if reached here, problem in the code
  d_print("problem in IsQuasiCore()\n");
  D_ERROR("IsQuasiCore()");
  exit(EXIT_FAILURE);
}
}
