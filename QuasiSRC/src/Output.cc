//
// Qutput.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

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
#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found
#endif /* HAVE UNISTD_H */

#ifdef HAVE_ERRNO_H
#include <errno.h>
#else
#error errno.h not found.
#endif /* HAVE_ERRNO_H */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#error limits.h not found
#endif /* HAVE_LIMITS_H */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found
#endif /* HAVE_SYS_TYPES_H */

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#else
#error sys/stat.h not found.
#endif /* HAVE_SYS_STAT_H */

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

#include <time.h>

// c++ files
#include "CrossNeighborList.h"
#include "DataTypes.h"
#include "Electrostatics.h"
#include "Element.h"
#include "Error.h"
#include "ForceEnergyCalculation.h"
#include "Indent.h"
#include "Input.h"
#include "Lattice.h"
#include "Output.h"
#include "PairPotentials.h"
#include "Quasicontinua.h"
#include "RunData.h"
#include "Void.h"

#include "read_pipe.h"
#include "threads.h"

static const char *tecplot_header_fields[] = {
#define TITLE 0
    "TITLE",
#define VARIABLES 1
    "VARIABLES",
#define ZONE 2
    "ZONE",
#define FEPOINT 3
    "FEPOINT",
#define TETRAHEDRON 4
    "TETRAHEDRON"};

// from tecplot.c
#define DISPLACEMENT_SCALE 1.0
#define DEFAULT_RESTART_FILE "restart"

// from output.c
#define MINIMAL_PATH_MAX 14

#define LATTICE_FILENAME_PREFIX "atom_"
#define LATTICE_FILENAME_SUFFIX ".plt.gz"
#define NODE_FILENAME_PREFIX "node_"
#define NODE_FILENAME_SUFFIX ".plt.gz"
#define RESTART_PREFIX "restart_"
#define RESTART_SUFFIX ".gz"
#define INPUT_PREFIX "input_"
#define INPUT_SUFFIX ".inp.gz"
#define CLUSTER_PREFIX "cluster_"
#define CLUSTER_SUFFIX ".plt.gz"
#define INDENT_PREFIX "indenter_"
#define INDENT_SUFFIX ".plt"

// for crossneighborlist data output
#define CROSSNEIGHBORLIST_CLUSTER_PREFIX "cluster_data_"
#define CROSSNEIGHBORLIST_NEIGHBOR_PREFIX "neighbor_data_"

#define XDR_ERR(msg)                                                           \
  do {                                                                         \
    fprintf(stderr, "[%s:%d] %s\n", __FILE__, __LINE__, msg);                  \
    exit(EXIT_FAILURE);                                                        \
    /*CONSTCOND*/                                                              \
  } while (0)

static double factor_ke = 1.036427059e-4;

//
//  converts integer to binary
//
//  Source: Stackoverflow website
//
const char *byte_to_binary(int x) {
  static char b[9];
  b[0] = '\0';

  int z;
  for (z = 128; z > 0; z >>= 1) {
    strcat(b, ((x & z) == z) ? "1" : "0");
  }

  return b;
}

//
//  converts integer to binary xyz mask
//
char *xyz_fixity(unsigned int x) {
  char *data = NULL;

  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);
  data = (char *)malloc(path_max * sizeof(char));

  // X fixity
  if ((x & FIX_X_MASK) == FREE_MASK) {
    // X is free

    if ((x & FIX_Y_MASK) == FREE_MASK) {
      // Y is free

      if ((x & FIX_Z_MASK) == FREE_MASK) {
        // Z is free
        sprintf(data, "000");
        // printf("here 000 x = %d\n", x);
      } else {
        // Z is fixed
        sprintf(data, "100");
        // printf("here 100 x = %d\n", x);
      }
    } else {
      // Y is fixed

      if ((x & FIX_Z_MASK) == FREE_MASK) {
        // Z is free
        sprintf(data, "010");
        // printf("here 010 x = %d\n", x);
      } else {
        // Z is fixed
        sprintf(data, "110");
        // printf("here 110 x = %d\n", x);
      }
    }
  } else {
    if ((x & FIX_Y_MASK) == FREE_MASK) {
      // Y is free

      if ((x & FIX_Z_MASK) == FREE_MASK) {
        // Z is free
        sprintf(data, "001");
        // printf("here 001 x = %d\n", x);
      } else {
        // Z is fixed
        sprintf(data, "101");
        // printf("here 101 x = %d\n", x);
      }
    } else {
      // Y is fixed

      if ((x & FIX_Z_MASK) == FREE_MASK) {
        // Z is free
        sprintf(data, "011");
        // printf("here 011 x = %d\n", x);
      } else {
        // Z is fixed
        sprintf(data, "111");
        // printf("here 111 x = %d\n", x);
      }
    }
  }

  return data;
}

namespace quasicontinuum {
//
//
//

Output *Output::_instance = NULL;

FILE *Output::d_debugFile = NULL;

//
// constructor
//

Output::Output() {

  //
  //
  //
  return;
}

//
// destructor
//

Output::~Output() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Output *Output::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Output();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Output::destroyInstance() {

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
//  createOutputDirectory()
//
void Output::createOutputDirectory(int outputFlag, char *outputDirectory) {
  //
  //  perform all this only once
  //
  static int init = 0;

  if (init == 0) {
    long int path_max;
    path_max = pathconf(".", _PC_PATH_MAX);

    //
    //  allocate space to d_outputDirectory
    //
    d_outputDirectory = AllocateFilename();

    //
    //  for error checking
    //
    struct stat buf;

    // extern int errno;
    int saved_errno;
    int stat_err;

    //
    //  create local folder string for output
    //
    time_t d_t = time(NULL);
    struct tm d_tm = *localtime(&d_t);
    char *dirname;

    dirname = AllocateFilename();

    sprintf(&(dirname[0]), "%d-%d_%d:%d:%d", d_tm.tm_mon + 1, d_tm.tm_mday,
            d_tm.tm_hour, d_tm.tm_min, d_tm.tm_sec);

    //
    //  now find whether to create directory at default place or at place given
    //  by quasi.ini
    //

    //
    //  local_flag = 0 :  default output directory assuming outputDirectory
    //  variable is
    //                    already allocated
    //             = 1 :  create output folder inside outputDirectory
    //

    int local_flag;
    local_flag = 0;

    if (outputFlag == 1) {
      //
      //  for this case we need to create output directory at the address given
      //  by outputDirectory
      //

      //
      //  check 1 : whether outputDirectory has been allocated space or not
      //
      if (outputDirectory == NULL) {
        local_flag = 0;
        // printf("creating output folder at default location bcause
        // outputDirectory is not allocated in Input::mainInput() while flag is
        // 1.\n");
      }

      if (outputDirectory[0] == '\0') {
        local_flag = 0;
        // printf("creating output folder at default location because
        // output_directory is empty in quasi.ini\n");
      }

      local_flag = 1;
    } else {
      // printf("output_flag is 0 thus creating output directory at default
      // place\n");

      local_flag = 0;
    }

    //
    // depending on local_flag create directory
    //
    if (local_flag == 1) {
      // add new folder to outputDirectory
      strcat(outputDirectory, dirname);

      // call mkdir function

      // make sure that directory does not exist and create it
      saved_errno = errno;

      stat_err = stat(outputDirectory, &buf);

      switch (stat_err) {
      case 0:
        errno = EEXIST;
        ERROR(outputDirectory);
        exit(EXIT_FAILURE);
        break;

      case -1:
        if (errno != ENOENT) {
          ERROR(outputDirectory);
          exit(EXIT_FAILURE);
        }
        break;
      }

      errno = saved_errno;

      // create directory
      if (mkdir(outputDirectory, S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
        ERROR("mkdir()");
        exit(EXIT_FAILURE);
      }

      //
      //  save the directory at d_outputDirectory of Output class
      //
      strcat(d_outputDirectory, outputDirectory);
    } else {
      //
      //  use getcwd() to get current working directory
      //
      if (outputDirectory == NULL)
        outputDirectory = (char *)malloc(path_max * sizeof(char));
      else
        outputDirectory =
            (char *)realloc(outputDirectory, path_max * sizeof(char));

      // get current directory
      if (getcwd(outputDirectory, path_max) == NULL)
        ERROR("getcwd() error");

      //
      //  copy dirname at the end of outputDirectory
      //
      char *local_dirname;
      local_dirname = AllocateFilename();

      strcat(outputDirectory, "/");
      strcat(outputDirectory, dirname);

      //
      //  create output directory at default place
      //
      // make sure that directory does not exist and create it
      saved_errno = errno;

      stat_err = stat(outputDirectory, &buf);

      switch (stat_err) {
      case 0:
        errno = EEXIST;
        ERROR(outputDirectory);
        exit(EXIT_FAILURE);
        break;

      case -1:
        if (errno != ENOENT) {
          ERROR(outputDirectory);
          exit(EXIT_FAILURE);
        }
        break;
      }

      errno = saved_errno;

      // create directory
      if (mkdir(outputDirectory, S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
        ERROR("mkdir()");
        exit(EXIT_FAILURE);
      }

      //
      //  save the directory at d_outputDirectory of Output class
      //
      strcat(d_outputDirectory, outputDirectory);
    }

    //
    //  output the info
    //
    // printf("********        Output Directory        *******\n");
    d_print("Output directory = %s \n", outputDirectory);
    // printf("********            Done        *******\n");

    //
    //  create d_degugFilename and d_debugFile
    //
    d_debugFilename = (char *)malloc(path_max * sizeof(char));
    sprintf(d_debugFilename, "%s/debug_data.out", outputDirectory);

    // create file and open it
    d_debugFile = fopen(d_debugFilename, "w");

    // update the init
    init = 1;
  }

  //
  //
  return;
}

//
//  performOutputCentral()
//
void Output::performOutputCentral(unsigned int flags,
                                  std::vector<int> output_flag_vector) {
  //
  //  check if d_outputDirectory is empty, if yes then call
  //  createOutputDirectory()
  //
  if (d_outputDirectory == NULL) {
    char *dummy;

    createOutputDirectory(0, dummy);
  }

  // loop over iQuasi and perform output
  int numQuasi = Quasicontinua::getInstance()->size();

  d_print("Output : Performing Output...");
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++)
    performOutput(iQuasi, flags, output_flag_vector);

  d_print("done\n");

  //
  //  perform output of time data if needed
  //
  if (RunData::getInstance()->d_timeOut == true) {
    FILE *fileOut;
    char *filename = AllocateFilename();

    int numThreads = get_max_number_threads();

    d_print("******** here\n");
    sprintf(filename, "%s/time_data_%d", d_outputDirectory, numThreads);

    // open file
    if ((fileOut = fopen(filename, "a")) == NULL) {
      ERROR(filename);
      exit(EXIT_FAILURE);
    }

    // get time data from RunData class
    std::vector<std::vector<double>> data;

    std::vector<int> clusterNumbers = RunData::getInstance()->writeData(data);

    // first dump data which are independent of iQuasi counter
    fprintf(fileOut, "number of threads                   %d\n", numThreads);
    fprintf(fileOut, "total time                          %lf\n", data[9][0]);
    fprintf(fileOut, "total time residual force           %lf\n", data[10][0]);

    for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
      fprintf(fileOut, "Time data for Quasi = %d\n", iQuasi);
      fprintf(fileOut, "Number of cluster sites              %d\n",
              clusterNumbers[iQuasi]);
      fprintf(fileOut, "time cluster and neighbor build      %lf\n",
              data[0][iQuasi]);
      fprintf(fileOut, "time neighbor build                  %lf\n",
              data[1][iQuasi]);
      fprintf(fileOut, "time update cluster and neighbor     %lf\n",
              data[2][iQuasi]);
      fprintf(fileOut, "time pairwise calculation            %lf\n",
              data[3][iQuasi]);
      fprintf(fileOut, "time eam calculation                 %lf\n",
              data[4][iQuasi]);
      fprintf(fileOut, "time electric field calculation      %lf\n",
              data[5][iQuasi]);
      fprintf(fileOut, "time electrostatics calculation      %lf\n",
              data[6][iQuasi]);
      fprintf(fileOut, "time entropy calculation             %lf\n",
              data[7][iQuasi]);
      fprintf(fileOut, "time cluster force to nodes          %lf\n",
              data[8][iQuasi]);
    }

    if (fclose(fileOut) == EOF) {
      ERROR(filename);
      exit(EXIT_FAILURE);
    }

    free(filename);
  }

  return;
} // end of performOutputCentral()

//
// define performOutput()
//
void Output::performOutput(const int iQuasi, unsigned int flags,
                           std::vector<int> output_flag_vector) {
  // get iQuasi instance
  Quasicontinuum iQuasicontinuum =
      Quasicontinua::getInstance()->getQuasi(iQuasi);

  // get all data needed
  struct all_node_list_t all_node_list = iQuasicontinuum.getNodeList();
  struct element_list_t element_list = iQuasicontinuum.getElementList();
  struct lattice_t lattice = iQuasicontinuum.getLattice();
  struct indentor_t indentor = iQuasicontinuum.getIndentor();
  char *materialName = iQuasicontinuum.getMaterialName();

  //
  // perform output
  //
  // Disabling following output option so that code
  // does not perform it for large system.
  // 1. PerformLatticeOutput()
  // 2. PerformInputFormat()
  // 3. PerformClusterOutput()
  // 4. PerformIndentOutput()
  // 5. PerformCrossNeighborListClusterOutput()
  // 6. PerformCrossNeighborListNeighborOutput()
  // 7. PerformClusterForcesOutput()
  //

  // if (flags & FULL_LATTICE_OUTPUT_FLAG)
  //   PerformLatticeOutput(iQuasi,
  //     d_outputDirectory,
  //     &lattice,
  //     mark,
  //     cg_info_flag,
  //     iteration_number);

  if (flags & NODE_OUTPUT_FLAG) {
    PerformNodeOutput(iQuasi, d_outputDirectory, &(all_node_list.node_list),
                      &(element_list), output_flag_vector);
  }

  // if (flags & INPUT_FLAG)
  //   PerformInputFormat(iQuasi,
  //       d_outputDirectory,
  //       &(all_node_list.node_list),
  //       &lattice,
  //       materialName,
  //       mark,
  //       cg_info_flag,
  //       iteration_number);

  // if (flags & CLUSTER_FLAG)
  //   PerformClusterOutput(iQuasi,
  //     d_outputDirectory,
  //     &(all_node_list.node_list),
  //     &lattice,
  //     mark,
  //     cg_info_flag,
  //     iteration_number);

  // if (flags & INDENT_FLAG)
  //   PerformIndentOutput(iQuasi,
  //     d_outputDirectory,
  //     &indentor,
  //     mark,
  //     cg_info_flag,
  //     iteration_number);

  // if (flags & CROSSNEIGHBORLIST_CLUSTER_FLAG)
  //   PerformCrossNeighborListClusterOutput(iQuasi,
  //     d_outputDirectory,
  //     mark,
  //     cg_info_flag,
  //     iteration_number);

  // if (flags & CROSSNEIGHBORLIST_NEIGHBOR_FLAG)
  //   PerformCrossNeighborListNeighborOutput(iQuasi,
  //     d_outputDirectory,
  //     mark,
  //     cg_info_flag,
  //     iteration_number);

  // if(flags & CLUSTER_FORCE_OUTPUT)
  // {
  //   PerformClusterForcesOutput(iQuasi,
  //     d_outputDirectory,
  //     mark,
  //     cg_info_flag,
  //     iteration_number);
  // }

  if (flags & RESTART_FLAG)
    PerformRestartOutput(iQuasi, d_outputDirectory, &all_node_list,
                         &element_list, &lattice, &indentor, materialName,
                         output_flag_vector);

  return;
}

//
// define AllocateFilename()
//
char *Output::AllocateFilename(void) {
  char *filename;

  static long path_max = -1;

  /**
  * on init call pathconf and get PATH_MAX
  */

  if (path_max < 0) {

    // extern int errno;
    int saved_errno;

    saved_errno = errno;

    if ((path_max = pathconf(".", _PC_PATH_MAX)) == -1) {
      if (errno != saved_errno) {
        ERROR("pathconf()");
        exit(EXIT_FAILURE);
      }

      errno = saved_errno;

#ifdef PATH_MAX
      path_max = PATH_MAX;
#else
      path_max = MINIMAL_PATH_MAX;
#endif /* PATH_MAX */
    }
  }

  /**
    * allocate space
    */

  if ((filename = (char *)malloc(path_max * sizeof(char))) == NULL) {
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  return (filename);
}

//
// perform lattice output
//

void Output::PerformLatticeOutput(const int iQuasi, const char *dirname,
                                  struct lattice_t *P_lattice, int mark,
                                  int cg_info_flag, int iteration_number) {
  FILE *file;
  char *filename = AllocateFilename();
  int l[3];
  int n_sites = 0;

  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  // get shift
  std::vector<double> shift = quasicontinua->getShift(iQuasi);

  // fill in filename
  if (cg_info_flag == CONVERGED_OUTPUT)
    sprintf(filename, "%s/%squasi_%i_load_number_%05d%s", dirname,
            LATTICE_FILENAME_PREFIX, iQuasi + 1, mark, LATTICE_FILENAME_SUFFIX);

  if (cg_info_flag == ITERATION_COMBINED_OUTPUT)
    sprintf(filename, "%s/"
                      "iteration_combined_%squasi_%i_load_number_%05d_"
                      "iteration_number_%05d%s",
            dirname, LATTICE_FILENAME_PREFIX, iQuasi + 1, mark,
            iteration_number, LATTICE_FILENAME_SUFFIX);

  if (cg_info_flag == ITERATION_FREQUENCY_OUTPUT)
    sprintf(filename, "%s/"
                      "iteration_frequency_combined_%squasi_%i_load_number_%"
                      "05d_iteration_number_%05d%s",
            dirname, LATTICE_FILENAME_PREFIX, iQuasi + 1, mark,
            iteration_number, LATTICE_FILENAME_SUFFIX);

  if (cg_info_flag == ITERATION_POSITION_OUTPUT)
    sprintf(filename, "%s/"
                      "iteration_position_combined_%squasi_%i_load_number_%05d_"
                      "iteration_number_%05d%s",
            dirname, LATTICE_FILENAME_PREFIX, iQuasi + 1, mark,
            iteration_number, LATTICE_FILENAME_SUFFIX);

  if (cg_info_flag == CONVERGED_INTERMEDIATE_COMBINED)
    sprintf(filename, "%s/"
                      "converged_intermediate_combined_%squasi_%i_load_number_%"
                      "05d_iteration_number_%05d%s",
            dirname, LATTICE_FILENAME_PREFIX, iQuasi + 1, mark,
            iteration_number, LATTICE_FILENAME_SUFFIX);

  if (cg_info_flag == CONVERGED_INTERMEDIATE_FREQUENCY)
    sprintf(filename, "%s/"
                      "converged_intermediate_frequency_%squasi_%i_load_number_"
                      "%05d_iteration_number_%05d%s",
            dirname, LATTICE_FILENAME_PREFIX, iQuasi + 1, mark,
            iteration_number, LATTICE_FILENAME_SUFFIX);

  if (cg_info_flag == CONVERGED_INTERMEDIATE_POSITION)
    sprintf(filename, "%s/"
                      "converged_intermediate_position_%squasi_%i_load_number_%"
                      "05d_iteration_number_%05d%s",
            dirname, LATTICE_FILENAME_PREFIX, iQuasi + 1, mark,
            iteration_number, LATTICE_FILENAME_SUFFIX);

  // open compressed file
  file = open_write_pipe(filename);

  // compute number of sites

  Lattice *latticeClass = Lattice::getInstance();

  for (l[0] = P_lattice->l_start[0]; l[0] <= P_lattice->l_end[0]; l[0]++)
    for (l[1] = P_lattice->l_start[1]; l[1] <= P_lattice->l_end[1]; l[1]++)
      for (l[2] = P_lattice->l_start[2]; l[2] <= P_lattice->l_end[2]; l[2]++)
        if (latticeClass->isLatticeSite(l, P_lattice, iQuasi) == RETURN_SUCCESS)
          n_sites++;

  // write header
  for (l[0] = P_lattice->l_start[0]; l[0] <= P_lattice->l_end[0]; l[0]++)
    for (l[1] = P_lattice->l_start[1]; l[1] <= P_lattice->l_end[1]; l[1]++)
      for (l[2] = P_lattice->l_start[2]; l[2] <= P_lattice->l_end[2]; l[2]++) {

        struct element_t *P_element;
        double S[4];
        double s[4];
        double t;

        if (latticeClass->isLatticeSite(l, P_lattice, iQuasi) == RETURN_FAILURE)
          continue;

        /**
          * find element the site is in
          */

        P_element = Element::getInstance()->locateSiteElement(
            NULL, l, P_lattice, iQuasi);

        /**
          * find material coordinates
          */

        latticeClass->getSiteInitialState(S, l, P_lattice, iQuasi);

        /**
          * find spatial coordinates
          */

        latticeClass->getSiteCurrentState(s, l, P_lattice, iQuasi);

        //
        t = Quasicontinua::getInstance()->getTemperature();

        /**
          * add global shift vector
          */
        S[0] += shift[0];
        S[1] += shift[1];
        S[2] += shift[2];
        s[0] += shift[0];
        s[1] += shift[1];
        s[2] += shift[2];

        /**
          * print
          */

        fprintf(file, "% 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e "
                      "% 8.6e % 8.6e %d %d\n",
                S[0], S[1], S[2], S[3], t, s[0], s[1], s[2], s[3], t,
                P_element->number, 0);
      }

  /**
   * clean-up
   */

  fclose(file);
  free(filename);

  return;
} // end of PerformLatticeOutput()

//
// PerformNodeOutput()
//
void Output::PerformNodeOutput(const int iQuasi, const char *dirname,
                               const struct node_list_t *P_node_list,
                               const struct element_list_t *P_element_list,
                               std::vector<int> output_flag_vector) {
  FILE *file;
  char *filename_temp = AllocateFilename();
  char *filename = AllocateFilename();

  //
  // get the basic filename depending on the values of flags in
  // output_flag_vector
  //
  OutputDataFilename(filename_temp, iQuasi, output_flag_vector);

  // d_print("%s", filename_temp);

  // add the node_data to the filename
  sprintf(filename, "%s/node_data_%s%s", dirname, filename_temp,
          NODE_FILENAME_SUFFIX);

  int i_node;
  int i_elem;

  /**
   * vector to hold relative shift
   */
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  /**
   * open pipe to gzip
   */
  file = open_write_pipe(filename);

  /**
   * write header
   */
  //
  //  check if we are minimizing frequency alone, and using quasi harmonic
  //  approximation. For this case, we know minimizing freq analytically
  //
  int minMethod = PairPotentials::getInstance()->d_minMethod;
  int quadMethod = PairPotentials::getInstance()->d_quadMethod;
  int electro_enable = Electrostatics::getInstance()->isElectrostaticEnable();

  d_print("here in PerformNodeOutput()\n");

  // discard electrostatics for quasi harmonic
  //  approximation when all nodes are fixed, i.e. for frequency
  // minimization
  if (minMethod == 2 && quadMethod == 1)
    electro_enable = 0;

  if (electro_enable == 1) {
    if (quadMethod == 1) {
      std::vector<double> traceOfK =
          Quasicontinua::getInstance()->getQuasi(iQuasi).getTraceOfKData();
      double mass =
          Quasicontinua::getInstance()->getQuasi(iQuasi).getAtomicMass();

      fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
      fprintf(file,
              " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\" \"Tau\""
              " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\" \"tau\" \"ww_analytical\""
              " \"ux\" \"uy\" \"uz\" \"uw\" \"utau\""
              " \"vx\" \"vy\" \"vz\""
              " \"fx\" \"fy\" \"fz\" \"fw\""
              " \"ext_fx\" \"ext_fy\" \"ext_fz\" \"ext_fw\""
              " \"Ex\" \"Ey\" \"Ez\""
              " \"energy\" \"fixity\" \"fixity_xyz_bin\" \"fixityw\"\n",
              tecplot_header_fields[VARIABLES]);

      fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
              tecplot_header_fields[ZONE], P_node_list->number_nodes,
              P_element_list->number_elements, tecplot_header_fields[FEPOINT],
              tecplot_header_fields[TETRAHEDRON]);

      /**
       * print data
       *
       * once Electrostatics is fully prepared, replace
       * dummy by ElectricField
       */
      const std::vector<std::vector<double>> ElectricField =
          Electrostatics::getInstance()->getElectricField(iQuasi);

      for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
        fprintf(file, "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e"
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "% 6.4e % 6.4e "
                      "%d %4s %d\n",
                P_node_list->nodes[i_node]->initial_position[0] + shift[0],
                P_node_list->nodes[i_node]->initial_position[1] + shift[1],
                P_node_list->nodes[i_node]->initial_position[2] + shift[2],
                P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->initial_temperature,
                P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->position[0] + shift[0],
                P_node_list->nodes[i_node]->position[1] + shift[1],
                P_node_list->nodes[i_node]->position[2] + shift[2],
                P_node_list->nodes[i_node]->frequency,
                P_node_list->nodes[i_node]->temperature,
                P_node_list->nodes[i_node]->tau,
                sqrt(mass * traceOfK[i_node] / (3.0 * factor_ke)),
                P_node_list->nodes[i_node]->position[0] -
                    P_node_list->nodes[i_node]->initial_position[0],
                P_node_list->nodes[i_node]->position[1] -
                    P_node_list->nodes[i_node]->initial_position[1],
                P_node_list->nodes[i_node]->position[2] -
                    P_node_list->nodes[i_node]->initial_position[2],
                P_node_list->nodes[i_node]->frequency -
                    P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->tau -
                    P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->velocity[0],
                P_node_list->nodes[i_node]->velocity[1],
                P_node_list->nodes[i_node]->velocity[2],
                P_node_list->nodes[i_node]->acceleration[0],
                P_node_list->nodes[i_node]->acceleration[1],
                P_node_list->nodes[i_node]->acceleration[2],
                P_node_list->nodes[i_node]->force_frequency,
                P_node_list->nodes[i_node]->external_force[0],
                P_node_list->nodes[i_node]->external_force[1],
                P_node_list->nodes[i_node]->external_force[2],
                P_node_list->nodes[i_node]->external_freq_force,
                ElectricField[i_node][0], ElectricField[i_node][1],
                ElectricField[i_node][2],
                P_node_list->nodes[i_node]->energy.potential,
                P_node_list->nodes[i_node]->fix_mask,
                xyz_fixity(P_node_list->nodes[i_node]->fix_mask),
                P_node_list->nodes[i_node]->fix_w_mask);
    } else {
      fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
      fprintf(file, " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\" \"Tau\""
                    " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\" \"tau\""
                    " \"ux\" \"uy\" \"uz\" \"uw\" \"utau\""
                    " \"vx\" \"vy\" \"vz\""
                    " \"fx\" \"fy\" \"fz\" \"fw\""
                    " \"ext_fx\" \"ext_fy\" \"ext_fz\" \"ext_fw\""
                    " \"Ex\" \"Ey\" \"Ez\""
                    " \"energy\" \"fixity\" \"fixity_xyz_bin\" \"fixityw\"\n",
              tecplot_header_fields[VARIABLES]);

      fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
              tecplot_header_fields[ZONE], P_node_list->number_nodes,
              P_element_list->number_elements, tecplot_header_fields[FEPOINT],
              tecplot_header_fields[TETRAHEDRON]);

      /**
       * print data
       *
       * once Electrostatics is fully prepared, replace
       * dummy by ElectricField
       */
      const std::vector<std::vector<double>> ElectricField =
          Electrostatics::getInstance()->getElectricField(iQuasi);

      for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
        fprintf(file,
                "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % "
                "6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e"
                "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % "
                "6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                "%d %4s %d\n",
                P_node_list->nodes[i_node]->initial_position[0] + shift[0],
                P_node_list->nodes[i_node]->initial_position[1] + shift[1],
                P_node_list->nodes[i_node]->initial_position[2] + shift[2],
                P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->initial_temperature,
                P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->position[0] + shift[0],
                P_node_list->nodes[i_node]->position[1] + shift[1],
                P_node_list->nodes[i_node]->position[2] + shift[2],
                P_node_list->nodes[i_node]->frequency,
                P_node_list->nodes[i_node]->temperature,
                P_node_list->nodes[i_node]->tau,
                P_node_list->nodes[i_node]->position[0] -
                    P_node_list->nodes[i_node]->initial_position[0],
                P_node_list->nodes[i_node]->position[1] -
                    P_node_list->nodes[i_node]->initial_position[1],
                P_node_list->nodes[i_node]->position[2] -
                    P_node_list->nodes[i_node]->initial_position[2],
                P_node_list->nodes[i_node]->frequency -
                    P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->tau -
                    P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->velocity[0],
                P_node_list->nodes[i_node]->velocity[1],
                P_node_list->nodes[i_node]->velocity[2],
                P_node_list->nodes[i_node]->acceleration[0],
                P_node_list->nodes[i_node]->acceleration[1],
                P_node_list->nodes[i_node]->acceleration[2],
                P_node_list->nodes[i_node]->force_frequency,
                P_node_list->nodes[i_node]->external_force[0],
                P_node_list->nodes[i_node]->external_force[1],
                P_node_list->nodes[i_node]->external_force[2],
                P_node_list->nodes[i_node]->external_freq_force,
                ElectricField[i_node][0], ElectricField[i_node][1],
                ElectricField[i_node][2],
                P_node_list->nodes[i_node]->energy.potential,
                P_node_list->nodes[i_node]->fix_mask,
                xyz_fixity(P_node_list->nodes[i_node]->fix_mask),
                P_node_list->nodes[i_node]->fix_w_mask);
    }
  } else {
    if (quadMethod == 1) {
      d_print("output here\n");

      std::vector<double> traceOfK =
          Quasicontinua::getInstance()->getQuasi(iQuasi).getTraceOfKData();

      double mass =
          Quasicontinua::getInstance()->getQuasi(iQuasi).getAtomicMass();

      fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
      fprintf(file,
              " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\" \"Tau\""
              " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\" \"tau\" \"ww_analytical\""
              " \"ux\" \"uy\" \"uz\" \"uw\" \"utau\""
              " \"vx\" \"vy\" \"vz\""
              " \"fx\" \"fy\" \"fz\" \"fw\""
              " \"ext_fx\" \"ext_fy\" \"ext_fz\" \"ext_fw\" \"weight\""
              " \"energy\" \"fixity\" \"fixity_xyz_bin\" \"fixityw\"\n",
              tecplot_header_fields[VARIABLES]);
      fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
              tecplot_header_fields[ZONE], P_node_list->number_nodes,
              P_element_list->number_elements, tecplot_header_fields[FEPOINT],
              tecplot_header_fields[TETRAHEDRON]);

      /**
       * print data
       */
      // d_print("********** analytical freq for quasi harmonic **********\n");
      for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
        fprintf(file, "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e"
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "%d %4s %d\n",
                P_node_list->nodes[i_node]->initial_position[0] + shift[0],
                P_node_list->nodes[i_node]->initial_position[1] + shift[1],
                P_node_list->nodes[i_node]->initial_position[2] + shift[2],
                P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->initial_temperature,
                P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->position[0] + shift[0],
                P_node_list->nodes[i_node]->position[1] + shift[1],
                P_node_list->nodes[i_node]->position[2] + shift[2],
                P_node_list->nodes[i_node]->frequency,
                P_node_list->nodes[i_node]->temperature,
                P_node_list->nodes[i_node]->tau,
                sqrt(mass * traceOfK[i_node] / (3.0 * factor_ke)),
                P_node_list->nodes[i_node]->position[0] -
                    P_node_list->nodes[i_node]->initial_position[0],
                P_node_list->nodes[i_node]->position[1] -
                    P_node_list->nodes[i_node]->initial_position[1],
                P_node_list->nodes[i_node]->position[2] -
                    P_node_list->nodes[i_node]->initial_position[2],
                P_node_list->nodes[i_node]->frequency -
                    P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->tau -
                    P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->velocity[0],
                P_node_list->nodes[i_node]->velocity[1],
                P_node_list->nodes[i_node]->velocity[2],
                P_node_list->nodes[i_node]->acceleration[0],
                P_node_list->nodes[i_node]->acceleration[1],
                P_node_list->nodes[i_node]->acceleration[2],
                P_node_list->nodes[i_node]->force_frequency,
                P_node_list->nodes[i_node]->external_force[0],
                P_node_list->nodes[i_node]->external_force[1],
                P_node_list->nodes[i_node]->external_force[2],
                P_node_list->nodes[i_node]->external_freq_force,
                P_node_list->nodes[i_node]->weight,
                P_node_list->nodes[i_node]->energy.potential,
                P_node_list->nodes[i_node]->fix_mask,
                xyz_fixity(P_node_list->nodes[i_node]->fix_mask),
                P_node_list->nodes[i_node]->fix_w_mask);

        // d_print("%f\n", sqrt(traceOfK[i_node]/3.0));
      }

    } else {
      fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
      fprintf(file, " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\" \"Tau\""
                    " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\" \"tau\""
                    // " \"ux\" \"uy\" \"uz\" \"uw\""
                    // " \"vx\" \"vy\" \"vz\""
                    " \"fx\" \"fy\" \"fz\" \"fw\""
                    " \"ext_fx\" \"ext_fy\" \"ext_fz\" \"ext_fw\" \"weight\""
                    " \"energy\" \"fixity\" \"fixity_xyz_bin\" \"fixityw\"\n",
              tecplot_header_fields[VARIABLES]);
      fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
              tecplot_header_fields[ZONE], P_node_list->number_nodes,
              P_element_list->number_elements, tecplot_header_fields[FEPOINT],
              tecplot_header_fields[TETRAHEDRON]);

      /**
       * print data
       */

      for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
        fprintf(file, "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e"
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                      "%d %4s %d\n",
                P_node_list->nodes[i_node]->initial_position[0] + shift[0],
                P_node_list->nodes[i_node]->initial_position[1] + shift[1],
                P_node_list->nodes[i_node]->initial_position[2] + shift[2],
                P_node_list->nodes[i_node]->initial_frequency,
                P_node_list->nodes[i_node]->initial_temperature,
                P_node_list->nodes[i_node]->initial_tau,
                P_node_list->nodes[i_node]->position[0] + shift[0],
                P_node_list->nodes[i_node]->position[1] + shift[1],
                P_node_list->nodes[i_node]->position[2] + shift[2],
                P_node_list->nodes[i_node]->frequency,
                P_node_list->nodes[i_node]->temperature,
                P_node_list->nodes[i_node]->tau,
                // P_node_list->nodes[i_node]->position[0]-
                //   P_node_list->nodes[i_node]->initial_position[0],
                // P_node_list->nodes[i_node]->position[1]-
                //   P_node_list->nodes[i_node]->initial_position[1],
                // P_node_list->nodes[i_node]->position[2]-
                //   P_node_list->nodes[i_node]->initial_position[2],
                // P_node_list->nodes[i_node]->frequency -
                //   P_node_list->nodes[i_node]->initial_frequency,
                // P_node_list->nodes[i_node]->velocity[0],
                // P_node_list->nodes[i_node]->velocity[1],
                // P_node_list->nodes[i_node]->velocity[2],
                P_node_list->nodes[i_node]->acceleration[0],
                P_node_list->nodes[i_node]->acceleration[1],
                P_node_list->nodes[i_node]->acceleration[2],
                P_node_list->nodes[i_node]->force_frequency,
                P_node_list->nodes[i_node]->external_force[0],
                P_node_list->nodes[i_node]->external_force[1],
                P_node_list->nodes[i_node]->external_force[2],
                P_node_list->nodes[i_node]->external_freq_force,
                P_node_list->nodes[i_node]->weight,
                P_node_list->nodes[i_node]->energy.potential,
                P_node_list->nodes[i_node]->fix_mask,
                xyz_fixity(P_node_list->nodes[i_node]->fix_mask),
                P_node_list->nodes[i_node]->fix_w_mask);

      //
      //  For debug
      //
      // int debug = 1; // 0 - debug output, 1- no debug output
      // if(iQuasi == 0 && debug == 0)
      // {
      //   char *filename_1 = AllocateFilename();

      //   sprintf(filename_1, "%s/node_0_iteration_%d.gz",
      //     dirname,
      //     iteration_number);

      //   FILE * file_1;

      //   file_1 = open_write_pipe(filename_1);

      //   // WORKING
      //   // for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
      //   //   fprintf(file_1, "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e %
      //   6.4e % 6.4e % 6.4e % 6.4e %d %d\n",
      //   //     P_node_list->nodes[i_node]->initial_position[0]+shift[0],
      //   //     P_node_list->nodes[i_node]->initial_position[1]+shift[1],
      //   //     P_node_list->nodes[i_node]->initial_position[2]+shift[2],
      //   //     P_node_list->nodes[i_node]->initial_frequency,
      //   //     P_node_list->nodes[i_node]->initial_temperature,
      //   //     P_node_list->nodes[i_node]->position[0]+shift[0],
      //   //     P_node_list->nodes[i_node]->position[1]+shift[1],
      //   //     P_node_list->nodes[i_node]->position[2]+shift[2],
      //   //     P_node_list->nodes[i_node]->frequency,
      //   //     P_node_list->nodes[i_node]->temperature,
      //   //     P_node_list->nodes[i_node]->fix_mask,
      //   //     P_node_list->nodes[i_node]->fix_w_mask);

      //   // TESTING
      //   // Conclusion: some data creates problem
      //   // when included in output.
      //   for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
      //     fprintf(file_1, "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e %
      //     6.4e % 6.4e"
      //     "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e %
      //     6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
      //     "%d %d\n",
      //     P_node_list->nodes[i_node]->initial_position[0]+shift[0],
      //     P_node_list->nodes[i_node]->initial_position[1]+shift[1],
      //     P_node_list->nodes[i_node]->initial_position[2]+shift[2],
      //     P_node_list->nodes[i_node]->initial_frequency,
      //     P_node_list->nodes[i_node]->initial_temperature,
      //     P_node_list->nodes[i_node]->position[0]+shift[0],
      //     P_node_list->nodes[i_node]->position[1]+shift[1],
      //     P_node_list->nodes[i_node]->position[2]+shift[2],
      //     P_node_list->nodes[i_node]->frequency,
      //     P_node_list->nodes[i_node]->temperature,
      //     P_node_list->nodes[i_node]->position[0]-
      //     P_node_list->nodes[i_node]->initial_position[0],
      //     P_node_list->nodes[i_node]->position[1]-
      //     P_node_list->nodes[i_node]->initial_position[1],
      //     P_node_list->nodes[i_node]->position[2]-
      //     P_node_list->nodes[i_node]->initial_position[2],
      //     P_node_list->nodes[i_node]->frequency -
      //     P_node_list->nodes[i_node]->initial_frequency,
      //     P_node_list->nodes[i_node]->acceleration[0],
      //     P_node_list->nodes[i_node]->acceleration[1],
      //     P_node_list->nodes[i_node]->acceleration[2],
      //     P_node_list->nodes[i_node]->force_frequency,
      //     P_node_list->nodes[i_node]->external_force[0],
      //     P_node_list->nodes[i_node]->external_force[1],
      //     P_node_list->nodes[i_node]->external_force[2],
      //     P_node_list->nodes[i_node]->external_freq_force,
      //     P_node_list->nodes[i_node]->weight,
      //     P_node_list->nodes[i_node]->energy.potential,
      //     P_node_list->nodes[i_node]->fix_mask,
      //     P_node_list->nodes[i_node]->fix_w_mask);

      //   // count = count+1;

      //   if (fclose(file_1) == EOF)
      //   {
      //     ERROR("fclose()");
      //     exit(EXIT_FAILURE);
      //   }

      //   free(filename_1);
      // } // debug block
    }
  }

  /**
   * write connectivity
   */

  fprintf(file, "\n");

  // printf("element-node connectivity : iElem  node1  node2  node3  node4\n");

  for (i_elem = 0; i_elem < P_element_list->number_elements; i_elem++) {
    fprintf(file, "%d %d %d %d\n",
            P_element_list->elements[i_elem]->node[0]->number + 1,
            P_element_list->elements[i_elem]->node[1]->number + 1,
            P_element_list->elements[i_elem]->node[2]->number + 1,
            P_element_list->elements[i_elem]->node[3]->number + 1);

    // printf("%d %d %d %d\n",
    //         P_element_list->elements[i_elem]->node[0]->number+1,
    //         P_element_list->elements[i_elem]->node[1]->number+1,
    //         P_element_list->elements[i_elem]->node[2]->number+1,
    //         P_element_list->elements[i_elem]->node[3]->number+1);
  }

  /**
   * cleanup
   */

  if (fclose(file) == EOF) {
    ERROR("fclose()");
    exit(EXIT_FAILURE);
  }

  free(filename);

  return;
} // end of PerformNodeOutput()

//
// PerformInputFormat()
//
void Output::PerformInputFormat(const int iQuasi, const char *dirname,
                                struct node_list_t *P_node_list,
                                struct lattice_t *P_lattice, char *materialName,
                                int mark, int cg_info_flag,
                                int iteration_number) {
  FILE *file;
  XDR xdrs;
  char *filename = AllocateFilename();

  char *name;
  int i_node;

  /**
   * fill in filename
   */
  if (cg_info_flag == CONVERGED_OUTPUT) {
    // d_print("CONVERGED\n");
    sprintf(filename, "%s/%squasi_%i_load_number_%05d%s", dirname, INPUT_PREFIX,
            iQuasi + 1, mark, INPUT_SUFFIX);
  }

  if (cg_info_flag == ITERATION_OUTPUT) {
    // d_print("ITERATION\n");
    sprintf(filename,
            "%s/iteration_%squasi_%i_load_number_%05d_iteration_number_%05d%s",
            dirname, INPUT_PREFIX, iQuasi + 1, mark, iteration_number,
            INPUT_SUFFIX);
  }

  if (cg_info_flag == CONVERGED_INTERMEDIATE) {
    // d_print("CONVERGED INTERMEDIATE\n");
    sprintf(filename, "%s/"
                      "converged_intermediate_%squasi_%i_load_number_%05d_"
                      "iteration_number_%05d%s",
            dirname, INPUT_PREFIX, iQuasi + 1, mark, iteration_number,
            INPUT_SUFFIX);
  }

  /**
   * open compressed file
   */

  file = open_write_pipe(filename);

  /**
   * create XDR stream
   */

  xdrstdio_create(&xdrs, file, XDR_ENCODE);

  /**
   * write data; should look like read_data
   */

  //  create name just to write it in file, though not be using it
  //  name = material_name_is_not_required;

  if (xdr_string(&xdrs, &materialName, _POSIX_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  name = P_lattice->name;
  if (xdr_string(&xdrs, &name, LATTICE_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  if (xdr_enum(&xdrs, (enum_t *)&(P_lattice->type)) == FALSE)
    XDR_ERR("xdr_enum()");

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

  if (xdr_int(&xdrs, &P_node_list->number_nodes) == FALSE)
    XDR_ERR("xdr_int()");

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {

    struct node_t *P_node = P_node_list->nodes[i_node];

    if (xdr_int(&xdrs, &(P_node->l[0])) == FALSE)
      XDR_ERR("xdr_int()");
    if (xdr_int(&xdrs, &(P_node->l[1])) == FALSE)
      XDR_ERR("xdr_int()");
    if (xdr_int(&xdrs, &(P_node->l[2])) == FALSE)
      XDR_ERR("xdr_int()");

    if (xdr_u_int(&xdrs, &(P_node->fix_mask)) == FALSE)
      XDR_ERR("xdr_int()");
  }

  /**
   * close XDR stream
   */

  xdr_destroy(&xdrs);

  /**
   * close file
   */

  if (fclose(file) == EOF) {
    ERROR(filename);
    exit(EXIT_FAILURE);
  }

  /**
   * cleanup
   */

  free(filename);

  return;
} // end of PerformInputFormat()

//
//  PerformClusterOutput()
//
void Output::PerformClusterOutput(const int iQuasi, const char *dirname,
                                  struct node_list_t *P_node_list,
                                  struct lattice_t *P_lattice, int mark,
                                  int cg_info_flag, int iteration_number) {

  FILE *file;

  char *filename = AllocateFilename();

  int i_node;
  int number_sites = 0;

  /**
   * vector to hold relative shift
   */
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  /**
   * fill in filename
   */

  if (cg_info_flag == CONVERGED_OUTPUT) {
    sprintf(filename, "%s/%squasi_%i_load_number_%05d%s", dirname,
            CLUSTER_PREFIX, iQuasi + 1, mark, CLUSTER_SUFFIX);
  }

  if (cg_info_flag == ITERATION_OUTPUT) {
    sprintf(filename,
            "%s/iteration_%squasi_%i_load_number_%05d_iteration_number_%05d%s",
            dirname, CLUSTER_PREFIX, iQuasi + 1, mark, iteration_number,
            CLUSTER_SUFFIX);
  }

  if (cg_info_flag == CONVERGED_INTERMEDIATE) {
    sprintf(filename, "%s/"
                      "converged_intermediate_%squasi_%i_load_number_%05d_"
                      "iteration_number_%05d%s",
            dirname, CLUSTER_PREFIX, iQuasi + 1, mark, iteration_number,
            CLUSTER_SUFFIX);
  }

  /**
   * open compressed file
   */

  file = open_write_pipe(filename);

  /**
   * loop over all nodes and compute number of sites in all clusters
   */

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
    number_sites += P_node_list->nodes[i_node]->site_cluster.number_sites;

  /**
   * compose a header
   */

  fprintf(file, " %s = \"%s\"\n", "TITLE", "Atoms");
  fprintf(
      file,
      " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\" \"xx\" \"yy\" \"zz\" \"ww\" \"tt\" "
      "\"Element No.\" \"fixity\" \"fixityw\" \n",
      "VARIABLES");
  fprintf(file, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

  /**
   * plot sites in all clusters
   */
  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {

    struct node_t *P_node = P_node_list->nodes[i_node];
    int i_site;

    /**
     * loop over all sites in a cluster
     */

    for (i_site = 0; i_site < P_node->site_cluster.number_sites; i_site++) {

      double S[4];
      double s[4];
      double t;
      int *l = P_node->site_cluster.sites[i_site];

      /**
       * compute coordinates
       */

      Lattice::getInstance()->getSiteInitialState(S, l, P_lattice, iQuasi);
      Lattice::getInstance()->getSiteCurrentState(S, l, P_lattice, iQuasi);

      t = Quasicontinua::getInstance()->getTemperature();

      /**
       * add global shift vector
       */
      S[0] += shift[0];
      S[1] += shift[1];
      S[2] += shift[2];
      s[0] += shift[0];
      s[1] += shift[1];
      s[2] += shift[2];

      /**
       * do the output
       */

      fprintf(file, "% 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % "
                    "8.6e % 8.6e %d %d %d\n",
              S[0], S[1], S[2], S[3], t, s[0], s[1], s[2], s[3], t, 0, 0, 0);
    }
  }

  /**
   * close file
   */

  if (fclose(file) == EOF) {
    ERROR(filename);
    exit(EXIT_FAILURE);
  }

  /**
   * cleanup
   */

  free(filename);

#if 0
  //
  // for debug purpose
  //
  //  Plotting cluster sites both from node list and CrossNeighborList 
  //  withing box, to see if both have equal number of points
  //  
  //  writing it in format readable by matlab
  //
  if(iQuasi == 0)
  {
    FILE *file_1;
    FILE *file_2;

    file_1 = fopen("node_cluster_matlab.quasi","w");

    // 
    //  plot points within box
    //
    //  x = [0, 10.0000]
    //  y = [0, 10.0000]
    //  z = [0, 10.0000]
    //
    double x_min = 0.0000;
    double x_max = 10.0000;

    for(int iNode = 0; iNode < P_node_list->number_nodes; iNode++)
    {
      struct node_t *P_node = P_node_list->nodes[iNode];

      int iC;
      for(iC=0; iC < P_node->site_cluster.number_sites; iC++)
      {
        double x[3];

        int *l = P_node->site_cluster.sites[iC];

        // get current position
        Lattice::getInstance()->getSiteCurrentPosition(x, l, P_lattice, iQuasi);

        if(x[0] >= x_min && x[0] <= x_max)
        {
          if(x[1] >= x_min && x[1] <= x_max)
          {
            if(x[2] >= x_min && x[2] <= x_max)
              fprintf(file_1, "%f   %f  %f\n", x[0], x[1], x[2]);
          }
        }
      }
    }

    fclose(file_1);

    // get cluster data from crossneighbor list
    std::vector< std::vector< cluster_site_data_t > > iClusterData
      = CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

    // open file
    file_2 = fopen("crossneighborlist_cluster_matlab.quasi","w");

    // loop over sites
    for(int iBucket =0; iBucket<iClusterData.size(); iBucket++)
    {
      if(iClusterData[iBucket].size() != 0)
      {
        for(int iC=0; iC < iClusterData[iBucket].size(); iC++)
        {
          // check if it is inside box
          if(  iClusterData[iBucket][iC].first.second[0] >= x_min 
            && iClusterData[iBucket][iC].first.second[0] <= x_max)
          {
            if(  iClusterData[iBucket][iC].first.second[1] >= x_min 
              && iClusterData[iBucket][iC].first.second[1] <= x_max)
            {
              if(  iClusterData[iBucket][iC].first.second[2] >= x_min 
                && iClusterData[iBucket][iC].first.second[2] <= x_max)
              {
                fprintf(file_2, "%f   %f  %f\n",
                  iClusterData[iBucket][iC].first.second[0],
                  iClusterData[iBucket][iC].first.second[1],
                  iClusterData[iBucket][iC].first.second[2]);
              }
            }
          }
        } // loop over cluster sites
      }// if bucket not empty
    } // loop over buckets  

    // close files
    fclose(file_2);
  }

#endif
  //  ******* end of debug block *******

  return;
} // end of PerformClusterOutput()

//
// PerformIndentOutput()
//
void Output::PerformIndentOutput(const int iQuasi, const char *dirname,
                                 const struct indentor_t *P_indentor, int mark,
                                 int cg_info_flag, int iteration_number) {
  FILE *file;

  char *filename = AllocateFilename();

  /**
   * fill in filename
   */
  if (cg_info_flag == CONVERGED_OUTPUT) {
    sprintf(filename, "%s/%squasi_%i_load_number_%05d%s", dirname,
            INDENT_PREFIX, iQuasi + 1, mark, INDENT_SUFFIX);
  }

  if (cg_info_flag == ITERATION_OUTPUT) {
    sprintf(filename, "%s/iteration_%squasi_%i_load_number_%05d%s", dirname,
            INDENT_PREFIX, iQuasi + 1, mark, INDENT_SUFFIX);
  }

  if (cg_info_flag == CONVERGED_INTERMEDIATE) {
    sprintf(filename, "%s/converged_intermediate_%squasi_%i_load_number_%05d%s",
            dirname, INDENT_PREFIX, iQuasi + 1, mark, INDENT_SUFFIX);
  }

  /**
   * open file
   */

  if ((file = fopen(filename, "a")) == NULL) {
    ERROR(filename);
    exit(EXIT_FAILURE);
  }

  /**
   * write indenter data
   */

  fprintf(file, "% 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e % 8.6e\n",
          P_indentor->position[0], P_indentor->position[1],
          P_indentor->position[2], P_indentor->force[0], P_indentor->force[1],
          P_indentor->force[2], P_indentor->energy.potential);

  /**
   * close file
   */

  if (fclose(file) == EOF) {
    ERROR(filename);
    exit(EXIT_FAILURE);
  }

  /**
   * cleanup
   */

  free(filename);

  return;
} // end of PerformIndentOutput()

//
//  PerformCrossNeighborListClusterOutput()
//
void Output::PerformCrossNeighborListClusterOutput(const int iQuasi,
                                                   const char *dirname,
                                                   int mark, int cg_info_flag,
                                                   int iteration_number) {
  FILE *file;

  char *filename = AllocateFilename();

  int number_sites = 0;

  /**
   * vector to hold relative shift
   */
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  //
  //  First write CrossNeighborList cluster data
  //

  /**
   * fill in filename
   */
  if (cg_info_flag == CONVERGED_OUTPUT) {
    sprintf(filename, "%s/%squasi_%i_load_number_%05d%s", dirname,
            CROSSNEIGHBORLIST_CLUSTER_PREFIX, iQuasi + 1, mark, CLUSTER_SUFFIX);
  }

  if (cg_info_flag == ITERATION_OUTPUT) {
    sprintf(filename,
            "%s/iteration_%squasi_%i_load_number_%05d_iteration_number_%05d%s",
            dirname, CROSSNEIGHBORLIST_CLUSTER_PREFIX, iQuasi + 1, mark,
            iteration_number, CLUSTER_SUFFIX);
  }

  if (cg_info_flag == CONVERGED_INTERMEDIATE) {
    sprintf(filename, "%s/"
                      "converged_intermediate_%squasi_%i_load_number_%05d_"
                      "iteration_number_%05d%s",
            dirname, CROSSNEIGHBORLIST_CLUSTER_PREFIX, iQuasi + 1, mark,
            iteration_number, CLUSTER_SUFFIX);
  }

  // open compressed file
  file = open_write_pipe(filename);

  // get cluster data from CrossNeighborList()
  std::vector<std::vector<cluster_site_data_t>> iClusterData =
      CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  // get total number os cluster sites
  for (int iBucket = 0; iBucket < iClusterData.size(); iBucket++)
    number_sites += iClusterData[iBucket].size();

  // compose a header
  fprintf(file, " %s = \"%s\"\n", "TITLE", "Atoms");
  fprintf(file, " %s = \"l\" \"m\" \"n\" \"X\" \"Y\" \"Z\" \"W\" \n",
          "VARIABLES");
  fprintf(file, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

  // put cluster sites into file
  for (int iBucket = 0; iBucket < iClusterData.size(); iBucket++) {
    if (iClusterData[iBucket].size() != 0) {
      for (int iC = 0; iC < iClusterData[iBucket].size(); iC++) {
        fprintf(file, "% i % i % i % 8.6e % 8.6e % 8.6e % 8.6e \n",
                iClusterData[iBucket][iC].first.first[0],
                iClusterData[iBucket][iC].first.first[1],
                iClusterData[iBucket][iC].first.first[2],
                iClusterData[iBucket][iC].first.second[0],
                iClusterData[iBucket][iC].first.second[1],
                iClusterData[iBucket][iC].first.second[2],
                iClusterData[iBucket][iC].first.second[3]);
      } // loop over cluster sites
    }   // if bucket not empty
  }     // loop over buckets

  // // compose a header
  // fprintf(file, " %s = \"%s\"\n", "TITLE", "Atoms");
  // fprintf(file, " %s = \"X\" \"Y\" \"Z\" \n",
  //         "VARIABLES");
  // fprintf(file, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

  // // put cluster sites into file
  // for(int iBucket =0; iBucket<iClusterData.size(); iBucket++)
  // {
  //   if(iClusterData[iBucket].size() != 0)
  //   {
  //     for(int iC=0; iC < iClusterData[iBucket].size(); iC++)
  //     {
  //       fprintf(file, "% 8.6e % 8.6e % 8.6e \n",
  //       iClusterData[iBucket][iC].first.second[0],
  //       iClusterData[iBucket][iC].first.second[1],
  //       iClusterData[iBucket][iC].first.second[2]);
  //     } // loop over cluster sites
  //   }// if bucket not empty
  // } // loop over buckets

  if (fclose(file) == EOF) {
    ERROR(filename);
    exit(EXIT_FAILURE);
  }

  // cleanup
  free(filename);

  return;
}

//
//  PerformCrossNeighborListNeighborOutput()
//
void Output::PerformCrossNeighborListNeighborOutput(const int iQuasi,
                                                    const char *dirname,
                                                    int mark, int cg_info_flag,
                                                    int iteration_number) {
  //
  //  proceed only if EAM is present
  //
  if (PairPotentials::getInstance()->d_EAMFlag == true) {
    FILE *file;

    char *filename = AllocateFilename();

    int number_sites = 0;

    /**
     * vector to hold relative shift
     */
    std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

    //
    //  First write CrossNeighborList cluster data
    //

    /**
     * fill in filename
     */
    if (cg_info_flag == CONVERGED_OUTPUT) {
      sprintf(filename, "%s/%squasi_%i_load_number_%05d%s", dirname,
              CROSSNEIGHBORLIST_NEIGHBOR_PREFIX, iQuasi + 1, mark,
              CLUSTER_SUFFIX);
    }

    if (cg_info_flag == ITERATION_OUTPUT) {
      sprintf(
          filename,
          "%s/iteration_%squasi_%i_load_number_%05d_iteration_number_%05d%s",
          dirname, CROSSNEIGHBORLIST_NEIGHBOR_PREFIX, iQuasi + 1, mark,
          iteration_number, CLUSTER_SUFFIX);
    }

    if (cg_info_flag == CONVERGED_INTERMEDIATE) {
      sprintf(filename, "%s/"
                        "converged_intermediate_%squasi_%i_load_number_%05d_"
                        "iteration_number_%05d%s",
              dirname, CROSSNEIGHBORLIST_NEIGHBOR_PREFIX, iQuasi + 1, mark,
              iteration_number, CLUSTER_SUFFIX);
    }

    // open compressed file
    file = open_write_pipe(filename);

    // get cluster data from CrossNeighborList()
    std::vector<std::vector<neigh_site_data_t>> iNeighborData =
        CrossNeighborList::getInstance()->getQuasiNeighborData(iQuasi);

    // get total number os cluster sites
    for (int iBucket = 0; iBucket < iNeighborData.size(); iBucket++)
      number_sites += iNeighborData[iBucket].size();

    // compose a header
    fprintf(file, " %s = \"%s\"\n", "TITLE", "Atoms");
    fprintf(file, " %s = \"l\" \"m\" \"n\" \"X\" \"Y\" \"Z\" \"W\" \n",
            "VARIABLES");
    fprintf(file, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

    // put cluster sites into file
    for (int iBucket = 0; iBucket < iNeighborData.size(); iBucket++) {
      if (iNeighborData[iBucket].size() != 0) {
        for (int iC = 0; iC < iNeighborData[iBucket].size(); iC++) {
          fprintf(file, "% i % i % i % 8.6e % 8.6e % 8.6e % 8.6e \n",
                  iNeighborData[iBucket][iC].first.first[0],
                  iNeighborData[iBucket][iC].first.first[1],
                  iNeighborData[iBucket][iC].first.first[2],
                  iNeighborData[iBucket][iC].first.second[0],
                  iNeighborData[iBucket][iC].first.second[1],
                  iNeighborData[iBucket][iC].first.second[2],
                  iNeighborData[iBucket][iC].first.second[3]);
        } // loop over cluster sites
      }   // if bucket not empty
    }     // loop over buckets

    if (fclose(file) == EOF) {
      ERROR(filename);
      exit(EXIT_FAILURE);
    }

    // cleanup
    free(filename);
  }

  return;
}

//
// PerformClusterForcesOutput()
//
void Output::PerformClusterForcesOutput(const int iQuasi, const char *dirname,
                                        int mark, int cg_info_flag,
                                        int iteration_number) {

  FILE *file;
  char *filename = AllocateFilename();

  int i_node;
  int i_elem;

  /**
   * vector to hold relative shift
   */
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  /**
   * fill in filename
   */
  sprintf(filename, "%s/cluster_forces_%i%s", dirname, iQuasi + 1,
          NODE_FILENAME_SUFFIX);

  /**
   * fill in filename
   */
  if (cg_info_flag == CONVERGED_OUTPUT) {
    sprintf(filename, "%s/cluster_site_forces_quasi_%i_load_number_%05d%s",
            dirname, iQuasi + 1, mark, NODE_FILENAME_SUFFIX);
  }

  if (cg_info_flag == ITERATION_OUTPUT) {
    sprintf(filename, "%s/"
                      "iteration_cluster_site_forces_quasi_%i_load_number_%05d_"
                      "iteration_number_%05d%s",
            dirname, iQuasi + 1, mark, iteration_number, NODE_FILENAME_SUFFIX);
  }

  if (cg_info_flag == CONVERGED_INTERMEDIATE) {
    sprintf(filename, "%s/"
                      "converged_intermediate_cluster_site_forces_quasi_%i_"
                      "load_number_%05d_iteration_number_%05d%s",
            dirname, iQuasi + 1, mark, iteration_number, NODE_FILENAME_SUFFIX);
  }

  // get neighbor data and cluster data
  CrossNeighborList *crossC = CrossNeighborList::getInstance();

  const std::vector<std::vector<cluster_site_data_t>> iQuasi_clusterData =
      crossC->getQuasiClusterData(iQuasi);

  const std::vector<std::vector<std::pair<std::vector<double>, double>>>
      iQuasi_clusterForceEnergy =
          ForceEnergyCalculation::getInstance()->getQuasiClusterSiteForceData(
              iQuasi);

  int number_sites = 0;
  for (int i = 0; i < iQuasi_clusterForceEnergy.size(); i++) {
    number_sites = number_sites + iQuasi_clusterForceEnergy[i].size();
  }

  /**
   * open pipe to gzip
   */
  file = open_write_pipe(filename);

  fprintf(file, " %s = \"%s\"\n", "TITLE", "Atoms");

  fprintf(file, " %s = \"X\" \"Y\" \"Z\" \"W\""
                " \"fx\" \"fy\" \"fz\" \"fw\""
                " \"energy\" \n",
          tecplot_header_fields[VARIABLES]);
  fprintf(file, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

  /**
   * print data
   */

  for (int b = 0; b < iQuasi_clusterForceEnergy.size(); b++)
    for (int c = 0; c < iQuasi_clusterForceEnergy[b].size(); c++)
      fprintf(file,
              "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e\n",
              iQuasi_clusterData[b][c].first.second[0],
              iQuasi_clusterData[b][c].first.second[1],
              iQuasi_clusterData[b][c].first.second[2],
              iQuasi_clusterData[b][c].first.second[3],
              iQuasi_clusterForceEnergy[b][c].first[0],
              iQuasi_clusterForceEnergy[b][c].first[1],
              iQuasi_clusterForceEnergy[b][c].first[2],
              iQuasi_clusterForceEnergy[b][c].first[3],
              iQuasi_clusterForceEnergy[b][c].second);

  /**
   * cleanup
   */

  if (fclose(file) == EOF) {
    ERROR("fclose()");
    exit(EXIT_FAILURE);
  }

  free(filename);

  return;
} // end of PerformClusterForcesOutput()

//
//  PerformRestartOutput()
//
void Output::PerformRestartOutput(const int iQuasi, const char *dirname,
                                  struct all_node_list_t *P_all_node_list,
                                  struct element_list_t *P_element_list,
                                  struct lattice_t *P_lattice,
                                  struct indentor_t *P_indentor,
                                  char *materialName,
                                  std::vector<int> output_flag_vector) {
  char *filename = AllocateFilename();
  char *filename_temp = AllocateFilename();

  //
  // get the basic filename depending on the values of flags in
  // output_flag_vector
  //
  OutputDataFilename(filename_temp, iQuasi, output_flag_vector);

  // add the node_data to the filename
  sprintf(filename, "%s/restart_%s%s", dirname, filename_temp, RESTART_SUFFIX);

  //
  //  re insert void nodes
  if (Void::getInstance()->isVoidEnable() != 0)
    Void::getInstance()->insertVoidNodes(&(P_all_node_list->node_list));

  //
  //  dump data
  //
  FILE *file;
  XDR xdrs;

  //
  //  open file
  //
  if ((file = open_data_file(filename, "w")) == NULL) {
    D_ERROR(filename);
    exit(EXIT_FAILURE);
  }

  xdrstdio_create(&xdrs, file, XDR_ENCODE);

  //
  //  write material name
  //
  if (xdr_string(&xdrs, &materialName, _POSIX_NAME_MAX) == FALSE)
    XDR_ERR("xdr_string()");

  //
  //  dump lattice data
  //
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
  //  dump nodes
  //
  struct node_list_t *P_node_list = &(P_all_node_list->node_list);

  int i_node;

  //
  //  write number of nodes
  //
  if (xdr_int(&xdrs, &P_node_list->number_nodes) == 0)
    XDR_ERR("xdr_int");

  //
  //  loop over all nodes
  //
  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

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
  }

  //
  //  dump element data
  //
  int i_elem;

  //
  //  write number of elements
  //
  if (xdr_int(&xdrs, &P_element_list->number_elements) == FALSE)
    XDR_ERR("xdr_int()");

  //
  //  loop over all elements
  //
  for (i_elem = 0; i_elem < P_element_list->number_elements; i_elem++) {
    struct element_t *P_element = P_element_list->elements[i_elem];

    /**
     * element number
     */
    if (xdr_int(&xdrs, &P_element->number) == FALSE)
      XDR_ERR("xdr_int()");

    /*
     * element node numbers
     */
    if (xdr_int(&xdrs, &P_element->node[0]->number) == FALSE ||
        xdr_int(&xdrs, &P_element->node[1]->number) == FALSE ||
        xdr_int(&xdrs, &P_element->node[2]->number) == FALSE ||
        xdr_int(&xdrs, &P_element->node[3]->number) == FALSE)
      XDR_ERR("xdr_int()");
  }

  //
  //  dump indentor data
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
    exit(EXIT_FAILURE);
  }

  return;
}

//
//  testEnergyOutput()
//
void Output::testEnergyOutput(
    std::vector<std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
        *test_cluster_force_energy,
    std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>
        *test_node_force_energy,
    std::vector<
        std::vector<std::vector<std::pair<std::vector<double>, double>>>>
        *test_total_cluster_force_energy,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *test_total_node_force_energy) {
  static int count = 0;
  //
  //  check if d_outputDirectory is empty, if yes then call
  //  createOutputDirectory()
  //
  if (d_outputDirectory == NULL) {
    char *dummy;

    createOutputDirectory(0, dummy);
  }

  // loop over iQuasi and perform output
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  int numQuasi = quasicontinua->size();

  d_print("Output : Performing Output...");
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    // get iQuasi instance
    Quasicontinuum iQuasicontinuum =
        Quasicontinua::getInstance()->getQuasi(iQuasi);

    // get all data needed
    struct node_list_t *P_node_list =
        &(iQuasicontinuum.getNodeList().node_list);
    struct element_list_t *P_element_list = &(iQuasicontinuum.getElementList());
    struct lattice_t *P_lattice = &(iQuasicontinuum.getLattice());
    struct indentor_t *P_indentor = &(iQuasicontinuum.getIndentor());
    char *materialName = iQuasicontinuum.getMaterialName();

    /**
     * vector to hold relative shift
     */
    std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

    //
    //  open file
    //
    FILE *file;
    char *filename = AllocateFilename();

    sprintf(filename, "%s/test_node_output_quasi_%i_count_%i%s",
            d_outputDirectory, iQuasi, count, NODE_FILENAME_SUFFIX);

    file = open_write_pipe(filename);

    //
    //  write header
    //
    if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
      if (PairPotentials::getInstance()->d_EAMFlag == true) {
        fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
        fprintf(
            file,
            " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\""
            " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\""
            " \"ux\" \"uy\" \"uz\" \"uw\""
            " \"vx\" \"vy\" \"vz\""
            " \"fx_eam\" \"fy_eam\" \"fz_eam\" \"fw_eam\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\""
            " \"fx_electrostatics\" \"fy_electrostatics\" "
            "\"fz_electrostatics\" \"fw__electrostatics\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_eam\" \"energy_pairwise\" \"energy_entropy\""
            " \"energy_electrostatics\" \"energy_total\""
            " \"Ex\" \"Ey\" \"Ez\""
            " \"weight\" \"fixity\" \"fixityw\"\n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
                tecplot_header_fields[ZONE], P_node_list->number_nodes,
                P_element_list->number_elements, tecplot_header_fields[FEPOINT],
                tecplot_header_fields[TETRAHEDRON]);

        /**
         * print data
         *
         * once Electrostatics is fully prepared, replace
         * dummy by ElectricField
         */
        const std::vector<std::vector<double>> ElectricField =
            Electrostatics::getInstance()->getElectricField(iQuasi);

        for (int i_node = 0; i_node < P_node_list->number_nodes; i_node++)
          fprintf(
              file,
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % d % d \n",
              P_node_list->nodes[i_node]->initial_position[0] + shift[0],
              P_node_list->nodes[i_node]->initial_position[1] + shift[1],
              P_node_list->nodes[i_node]->initial_position[2] + shift[2],
              P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->initial_temperature,
              P_node_list->nodes[i_node]->position[0] + shift[0],
              P_node_list->nodes[i_node]->position[1] + shift[1],
              P_node_list->nodes[i_node]->position[2] + shift[2],
              P_node_list->nodes[i_node]->frequency,
              P_node_list->nodes[i_node]->temperature,
              P_node_list->nodes[i_node]->position[0] -
                  P_node_list->nodes[i_node]->initial_position[0],
              P_node_list->nodes[i_node]->position[1] -
                  P_node_list->nodes[i_node]->initial_position[1],
              P_node_list->nodes[i_node]->position[2] -
                  P_node_list->nodes[i_node]->initial_position[2],
              P_node_list->nodes[i_node]->frequency -
                  P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->velocity[0],
              P_node_list->nodes[i_node]->velocity[1],
              P_node_list->nodes[i_node]->velocity[2],
              (*test_node_force_energy)[iQuasi][i_node].first[0],
              (*test_node_force_energy)[iQuasi][i_node].first[1],
              (*test_node_force_energy)[iQuasi][i_node].first[2],
              (*test_node_force_energy)[iQuasi][i_node].first[3],
              (*test_node_force_energy)[iQuasi][i_node].first[4],
              (*test_node_force_energy)[iQuasi][i_node].first[5],
              (*test_node_force_energy)[iQuasi][i_node].first[6],
              (*test_node_force_energy)[iQuasi][i_node].first[7],
              (*test_node_force_energy)[iQuasi][i_node].first[8],
              (*test_node_force_energy)[iQuasi][i_node].first[9],
              (*test_node_force_energy)[iQuasi][i_node].first[10],
              (*test_node_force_energy)[iQuasi][i_node].first[11],
              (*test_node_force_energy)[iQuasi][i_node].first[12],
              (*test_total_node_force_energy)[iQuasi][i_node].first[0],
              (*test_total_node_force_energy)[iQuasi][i_node].first[1],
              (*test_total_node_force_energy)[iQuasi][i_node].first[2],
              (*test_total_node_force_energy)[iQuasi][i_node].first[3],
              (*test_node_force_energy)[iQuasi][i_node].second[0],
              (*test_node_force_energy)[iQuasi][i_node].second[1],
              (*test_node_force_energy)[iQuasi][i_node].second[2],
              (*test_node_force_energy)[iQuasi][i_node].second[3],
              (*test_total_node_force_energy)[iQuasi][i_node].second,
              ElectricField[i_node][0], ElectricField[i_node][1],
              ElectricField[i_node][2], P_node_list->nodes[i_node]->weight,
              P_node_list->nodes[i_node]->fix_mask,
              P_node_list->nodes[i_node]->fix_w_mask);
      } else {
        fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
        fprintf(
            file,
            " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\""
            " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\""
            " \"ux\" \"uy\" \"uz\" \"uw\""
            " \"vx\" \"vy\" \"vz\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\""
            " \"fx_electrostatics\" \"fy_electrostatics\" "
            "\"fz_electrostatics\" \"fw__electrostatics\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_pairwise\" \"energy_entropy\""
            " \"energy_electrostatics\" \"energy_total\""
            " \"Ex\" \"Ey\" \"Ez\""
            " \"weight\" \"fixity\" \"fixityw\"\n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
                tecplot_header_fields[ZONE], P_node_list->number_nodes,
                P_element_list->number_elements, tecplot_header_fields[FEPOINT],
                tecplot_header_fields[TETRAHEDRON]);

        /**
         * print data
         *
         * once Electrostatics is fully prepared, replace
         * dummy by ElectricField
         */
        const std::vector<std::vector<double>> ElectricField =
            Electrostatics::getInstance()->getElectricField(iQuasi);

        for (int i_node = 0; i_node < P_node_list->number_nodes; i_node++)
          fprintf(
              file,
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % d % d \n",
              P_node_list->nodes[i_node]->initial_position[0] + shift[0],
              P_node_list->nodes[i_node]->initial_position[1] + shift[1],
              P_node_list->nodes[i_node]->initial_position[2] + shift[2],
              P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->initial_temperature,
              P_node_list->nodes[i_node]->position[0] + shift[0],
              P_node_list->nodes[i_node]->position[1] + shift[1],
              P_node_list->nodes[i_node]->position[2] + shift[2],
              P_node_list->nodes[i_node]->frequency,
              P_node_list->nodes[i_node]->temperature,
              P_node_list->nodes[i_node]->position[0] -
                  P_node_list->nodes[i_node]->initial_position[0],
              P_node_list->nodes[i_node]->position[1] -
                  P_node_list->nodes[i_node]->initial_position[1],
              P_node_list->nodes[i_node]->position[2] -
                  P_node_list->nodes[i_node]->initial_position[2],
              P_node_list->nodes[i_node]->frequency -
                  P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->velocity[0],
              P_node_list->nodes[i_node]->velocity[1],
              P_node_list->nodes[i_node]->velocity[2],
              (*test_node_force_energy)[iQuasi][i_node].first[4],
              (*test_node_force_energy)[iQuasi][i_node].first[5],
              (*test_node_force_energy)[iQuasi][i_node].first[6],
              (*test_node_force_energy)[iQuasi][i_node].first[7],
              (*test_node_force_energy)[iQuasi][i_node].first[8],
              (*test_node_force_energy)[iQuasi][i_node].first[9],
              (*test_node_force_energy)[iQuasi][i_node].first[10],
              (*test_node_force_energy)[iQuasi][i_node].first[11],
              (*test_node_force_energy)[iQuasi][i_node].first[12],
              (*test_total_node_force_energy)[iQuasi][i_node].first[0],
              (*test_total_node_force_energy)[iQuasi][i_node].first[1],
              (*test_total_node_force_energy)[iQuasi][i_node].first[2],
              (*test_total_node_force_energy)[iQuasi][i_node].first[3],
              (*test_node_force_energy)[iQuasi][i_node].second[1],
              (*test_node_force_energy)[iQuasi][i_node].second[2],
              (*test_node_force_energy)[iQuasi][i_node].second[3],
              (*test_total_node_force_energy)[iQuasi][i_node].second,
              ElectricField[i_node][0], ElectricField[i_node][1],
              ElectricField[i_node][2], P_node_list->nodes[i_node]->weight,
              P_node_list->nodes[i_node]->fix_mask,
              P_node_list->nodes[i_node]->fix_w_mask);
      }
    } else {
      if (PairPotentials::getInstance()->d_EAMFlag == true) {
        fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
        fprintf(
            file,
            " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\""
            " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\""
            " \"ux\" \"uy\" \"uz\" \"uw\""
            " \"vx\" \"vy\" \"vz\""
            " \"fx_eam\" \"fy_eam\" \"fz_eam\" \"fw_eam\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_eam\" \"energy_pairwise\" \"energy_entropy\""
            " \"energy_total\""
            " \"weight\" \"fixity\" \"fixityw\"\n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
                tecplot_header_fields[ZONE], P_node_list->number_nodes,
                P_element_list->number_elements, tecplot_header_fields[FEPOINT],
                tecplot_header_fields[TETRAHEDRON]);

        /**
         * print data
         *
         * once Electrostatics is fully prepared, replace
         * dummy by ElectricField
         */
        for (int i_node = 0; i_node < P_node_list->number_nodes; i_node++)
          fprintf(
              file,
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % d % d \n",
              P_node_list->nodes[i_node]->initial_position[0] + shift[0],
              P_node_list->nodes[i_node]->initial_position[1] + shift[1],
              P_node_list->nodes[i_node]->initial_position[2] + shift[2],
              P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->initial_temperature,
              P_node_list->nodes[i_node]->position[0] + shift[0],
              P_node_list->nodes[i_node]->position[1] + shift[1],
              P_node_list->nodes[i_node]->position[2] + shift[2],
              P_node_list->nodes[i_node]->frequency,
              P_node_list->nodes[i_node]->temperature,
              P_node_list->nodes[i_node]->position[0] -
                  P_node_list->nodes[i_node]->initial_position[0],
              P_node_list->nodes[i_node]->position[1] -
                  P_node_list->nodes[i_node]->initial_position[1],
              P_node_list->nodes[i_node]->position[2] -
                  P_node_list->nodes[i_node]->initial_position[2],
              P_node_list->nodes[i_node]->frequency -
                  P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->velocity[0],
              P_node_list->nodes[i_node]->velocity[1],
              P_node_list->nodes[i_node]->velocity[2],
              (*test_node_force_energy)[iQuasi][i_node].first[0],
              (*test_node_force_energy)[iQuasi][i_node].first[1],
              (*test_node_force_energy)[iQuasi][i_node].first[2],
              (*test_node_force_energy)[iQuasi][i_node].first[3],
              (*test_node_force_energy)[iQuasi][i_node].first[4],
              (*test_node_force_energy)[iQuasi][i_node].first[5],
              (*test_node_force_energy)[iQuasi][i_node].first[6],
              (*test_node_force_energy)[iQuasi][i_node].first[7],
              (*test_node_force_energy)[iQuasi][i_node].first[8],
              (*test_total_node_force_energy)[iQuasi][i_node].first[0],
              (*test_total_node_force_energy)[iQuasi][i_node].first[1],
              (*test_total_node_force_energy)[iQuasi][i_node].first[2],
              (*test_total_node_force_energy)[iQuasi][i_node].first[3],
              (*test_node_force_energy)[iQuasi][i_node].second[0],
              (*test_node_force_energy)[iQuasi][i_node].second[1],
              (*test_node_force_energy)[iQuasi][i_node].second[2],
              (*test_total_node_force_energy)[iQuasi][i_node].second,
              P_node_list->nodes[i_node]->weight,
              P_node_list->nodes[i_node]->fix_mask,
              P_node_list->nodes[i_node]->fix_w_mask);
      } else {
        fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
        fprintf(
            file,
            " %s = \"X\" \"Y\" \"Z\" \"W\" \"T\""
            " \"xx\" \"yy\" \"zz\" \"ww\" \"tt\""
            " \"ux\" \"uy\" \"uz\" \"uw\""
            " \"vx\" \"vy\" \"vz\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_pairwise\" \"energy_entropy\""
            " \"energy_total\""
            " \"weight\" \"fixity\" \"fixityw\"\n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n",
                tecplot_header_fields[ZONE], P_node_list->number_nodes,
                P_element_list->number_elements, tecplot_header_fields[FEPOINT],
                tecplot_header_fields[TETRAHEDRON]);

        /**
         * print data
         *
         * once Electrostatics is fully prepared, replace
         * dummy by ElectricField
         */
        for (int i_node = 0; i_node < P_node_list->number_nodes; i_node++)
          fprintf(
              file,
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
              "% 6.4e % 6.4e % 6.4e % d % d \n",
              P_node_list->nodes[i_node]->initial_position[0] + shift[0],
              P_node_list->nodes[i_node]->initial_position[1] + shift[1],
              P_node_list->nodes[i_node]->initial_position[2] + shift[2],
              P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->initial_temperature,
              P_node_list->nodes[i_node]->position[0] + shift[0],
              P_node_list->nodes[i_node]->position[1] + shift[1],
              P_node_list->nodes[i_node]->position[2] + shift[2],
              P_node_list->nodes[i_node]->frequency,
              P_node_list->nodes[i_node]->temperature,
              P_node_list->nodes[i_node]->position[0] -
                  P_node_list->nodes[i_node]->initial_position[0],
              P_node_list->nodes[i_node]->position[1] -
                  P_node_list->nodes[i_node]->initial_position[1],
              P_node_list->nodes[i_node]->position[2] -
                  P_node_list->nodes[i_node]->initial_position[2],
              P_node_list->nodes[i_node]->frequency -
                  P_node_list->nodes[i_node]->initial_frequency,
              P_node_list->nodes[i_node]->velocity[0],
              P_node_list->nodes[i_node]->velocity[1],
              P_node_list->nodes[i_node]->velocity[2],
              (*test_node_force_energy)[iQuasi][i_node].first[4],
              (*test_node_force_energy)[iQuasi][i_node].first[5],
              (*test_node_force_energy)[iQuasi][i_node].first[6],
              (*test_node_force_energy)[iQuasi][i_node].first[7],
              (*test_node_force_energy)[iQuasi][i_node].first[8],
              (*test_total_node_force_energy)[iQuasi][i_node].first[0],
              (*test_total_node_force_energy)[iQuasi][i_node].first[1],
              (*test_total_node_force_energy)[iQuasi][i_node].first[2],
              (*test_total_node_force_energy)[iQuasi][i_node].first[3],
              (*test_node_force_energy)[iQuasi][i_node].second[1],
              (*test_node_force_energy)[iQuasi][i_node].second[2],
              (*test_total_node_force_energy)[iQuasi][i_node].second,
              P_node_list->nodes[i_node]->weight,
              P_node_list->nodes[i_node]->fix_mask,
              P_node_list->nodes[i_node]->fix_w_mask);
      }
    }

    /**
     * write connectivity
     */

    fprintf(file, "\n");

    // printf("element-node connectivity : iElem  node1  node2  node3 node4\n");

    for (int i_elem = 0; i_elem < P_element_list->number_elements; i_elem++) {
      fprintf(file, "%d %d %d %d\n",
              P_element_list->elements[i_elem]->node[0]->number + 1,
              P_element_list->elements[i_elem]->node[1]->number + 1,
              P_element_list->elements[i_elem]->node[2]->number + 1,
              P_element_list->elements[i_elem]->node[3]->number + 1);
    }

    /**
     * cleanup
     */
    if (fclose(file) == EOF) {
      ERROR("fclose()");
      exit(EXIT_FAILURE);
    }

    //
    //  write cluster data to file
    //
    FILE *file_1;

    sprintf(filename, "%s/test_cluster_output_quasi_%i_count_%i%s",
            d_outputDirectory, iQuasi, count, NODE_FILENAME_SUFFIX);

    file_1 = open_write_pipe(filename);

    //
    //  get the cluster site data from cross neighbor list
    //
    const std::vector<std::vector<cluster_site_data_t>> iQuasi_clusterData =
        CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

    int number_sites = 0;
    for (int i = 0; i < iQuasi_clusterData.size(); i++) {
      number_sites = number_sites + iQuasi_clusterData[i].size();
    }

    //
    //  title
    //
    fprintf(file_1, " %s = \"%s\"\n", "TITLE", "Atoms");

    if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
      if (PairPotentials::getInstance()->d_EAMFlag == true) {
        fprintf(
            file_1,
            " %s = \"X\" \"Y\" \"Z\" \"W\""
            " \"fx_eam\" \"fy_eam\" \"fz_eam\" \"fw_eam\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\" \"fx_electrostatics\""
            " \"fy_electrostatics\" \"fz_electrostatics\" \"fw_electrostatics\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_eam\" \"energy_pairwise\" \"energy_entropy\""
            " \"energy_electrostatics\" \"energy_total\" \n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file_1, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

        /**
         * print data
         */
        for (int b = 0; b < iQuasi_clusterData.size(); b++)
          for (int c = 0; c < iQuasi_clusterData[b].size(); c++) {
            fprintf(file_1,
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e \n",
                    iQuasi_clusterData[b][c].first.second[0],
                    iQuasi_clusterData[b][c].first.second[1],
                    iQuasi_clusterData[b][c].first.second[2],
                    iQuasi_clusterData[b][c].first.second[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[0],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[1],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[2],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[4],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[5],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[6],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[7],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[8],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[9],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[10],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[11],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[12],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[0],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[1],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[2],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[0],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[1],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[2],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[3],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].second);
          }
      } else {
        fprintf(
            file_1,
            " %s = \"X\" \"Y\" \"Z\" \"W\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\" \"fx_electrostatics\""
            " \"fy_electrostatics\" \"fz_electrostatics\" \"fw_electrostatics\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_pairwise\" \"energy_entropy\""
            " \"energy_electrostatics\" \"energy_total\" \n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file_1, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

        /**
         * print data
         */
        for (int b = 0; b < iQuasi_clusterData.size(); b++)
          for (int c = 0; c < iQuasi_clusterData[b].size(); c++) {
            fprintf(file_1,
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e \n",
                    iQuasi_clusterData[b][c].first.second[0],
                    iQuasi_clusterData[b][c].first.second[1],
                    iQuasi_clusterData[b][c].first.second[2],
                    iQuasi_clusterData[b][c].first.second[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[4],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[5],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[6],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[7],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[8],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[9],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[10],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[11],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[12],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[0],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[1],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[2],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[1],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[2],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[3],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].second);
          }
      }
    } else {
      if (PairPotentials::getInstance()->d_EAMFlag == true) {
        fprintf(
            file_1,
            " %s = \"X\" \"Y\" \"Z\" \"W\""
            " \"fx_eam\" \"fy_eam\" \"fz_eam\" \"fw_eam\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\""
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_eam\" \"energy_pairwise\" \"energy_entropy\""
            " \"energy_total\" \n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file_1, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

        /**
         * print data
         */
        for (int b = 0; b < iQuasi_clusterData.size(); b++)
          for (int c = 0; c < iQuasi_clusterData[b].size(); c++) {
            fprintf(file_1,
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e \n",
                    iQuasi_clusterData[b][c].first.second[0],
                    iQuasi_clusterData[b][c].first.second[1],
                    iQuasi_clusterData[b][c].first.second[2],
                    iQuasi_clusterData[b][c].first.second[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[0],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[1],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[2],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[4],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[5],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[6],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[7],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[8],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[0],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[1],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[2],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[0],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[1],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[2],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].second);
          }
      } else {
        fprintf(
            file_1,
            " %s = \"X\" \"Y\" \"Z\" \"W\""
            " \"fx_pairwise\" \"fy_pairwise\" \"fz_pairwise\" \"fw_pairwise\""
            " \"fw_entropy\" "
            " \"fx_total\" \"fy_total\" \"fz_total\" \"fw_total\""
            " \"energy_pairwise\" \"energy_entropy\""
            " \"energy_total\" \n",
            tecplot_header_fields[VARIABLES]);

        fprintf(file_1, " %s I= %d, F=%s\n", "ZONE", number_sites, "POINT");

        /**
         * print data
         */
        for (int b = 0; b < iQuasi_clusterData.size(); b++)
          for (int c = 0; c < iQuasi_clusterData[b].size(); c++) {
            fprintf(file_1,
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e "
                    "% 6.4e %6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e % 6.4e \n",
                    iQuasi_clusterData[b][c].first.second[0],
                    iQuasi_clusterData[b][c].first.second[1],
                    iQuasi_clusterData[b][c].first.second[2],
                    iQuasi_clusterData[b][c].first.second[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[4],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[5],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[6],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[7],
                    (*test_cluster_force_energy)[iQuasi][b][c].first[8],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[0],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[1],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[2],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].first[3],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[1],
                    (*test_cluster_force_energy)[iQuasi][b][c].second[2],
                    (*test_total_cluster_force_energy)[iQuasi][b][c].second);
          }
      }
    }

    /**
     * cleanup
     */
    if (fclose(file_1) == EOF) {
      ERROR("fclose()");
      exit(EXIT_FAILURE);
    }

    free(filename);
    // if( free(filename) )
    //   d_print("error in freeing filename\n");
  }

  count += 1;

  d_print("clearing the vector data...");
  test_cluster_force_energy->clear();
  test_node_force_energy->clear();
  test_total_cluster_force_energy->clear();
  test_total_node_force_energy->clear();
  d_print("done\n");

  d_print("done\n");

  return;
}

/**
 * @brief performNodeListOutput()
 */
void Output::performNodeListOutput(const char *file_name,
                                   const struct node_list_t *P_node_list) {
  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);

  char *local_name;
  local_name = (char *)malloc(path_max * sizeof(char));
  //
  // create and open file
  //
  sprintf(local_name, "%s/%s%s", d_outputDirectory, file_name,
          NODE_FILENAME_SUFFIX);

  /**
   * open pipe to gzip
   */
  FILE *file;
  file = open_write_pipe(local_name);

  /**
   * write header
   */
  fprintf(file, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
  fprintf(file, " %s = \"l\" \"m\" \"n\""
                " \"X\" \"Y\" \"Z\" \"W\" \"T\""
                " \"fixity\" \"fixityw\"\n",
          tecplot_header_fields[VARIABLES]);

  fprintf(file, " %s N= %d, E= %d, F=%s, ET=%s\n", tecplot_header_fields[ZONE],
          P_node_list->number_nodes, 0, tecplot_header_fields[FEPOINT],
          tecplot_header_fields[TETRAHEDRON]);

  /**
   * print data
   */
  double xmax[3];
  double xmin[3];

  for (int dof = 0; dof < 3; dof++) {
    xmax[dof] = 0.0;
    xmin[dof] = 0.0;
  }

  for (int i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    for (int dof = 0; dof < 3; dof++) {
      if (xmax[dof] < P_node_list->nodes[i_node]->initial_position[dof])
        xmax[dof] = P_node_list->nodes[i_node]->initial_position[dof];
      if (xmin[dof] > P_node_list->nodes[i_node]->initial_position[dof])
        xmin[dof] = P_node_list->nodes[i_node]->initial_position[dof];
    }

    fprintf(file, " %d %d %d %6.4e %6.4e %6.4e %6.4e %6.4e %d %d\n",
            P_node_list->nodes[i_node]->l[0], P_node_list->nodes[i_node]->l[1],
            P_node_list->nodes[i_node]->l[2],
            P_node_list->nodes[i_node]->initial_position[0],
            P_node_list->nodes[i_node]->initial_position[1],
            P_node_list->nodes[i_node]->initial_position[2],
            P_node_list->nodes[i_node]->initial_frequency,
            P_node_list->nodes[i_node]->initial_temperature,
            P_node_list->nodes[i_node]->fix_mask,
            P_node_list->nodes[i_node]->fix_w_mask);
  }
  printf("xmax = (%f, %f, %f), xmin = (%f, %f, %f)\n", xmax[0], xmax[1],
         xmax[2], xmin[0], xmin[1], xmin[2]);

  /**
   * cleanup
   */

  if (fclose(file) == EOF) {
    ERROR("fclose()");
    exit(EXIT_FAILURE);
  }

  free(local_name);

  return;
}

//
// Depending on values of various flags, it creates the suitable
// filename to write the data
//
// std::vector<int> flag_data
// flag_data[0] --> minimization_method
// flag_data[1] --> freq_or_position_minimization
//                  -1 for all method except for
//                  minimization_method 1
// flag_data[2] --> cg_iteration_number
// flag_data[3] --> alternate_iteration_number
// flag_data[4] --> 0 - ITERATION
//                  1 - CONVERGED
//                  2 - CONVERGED_ALTERNATE
//                  3 - CONVERGED_FINAL
// flag_data[4] --> load_number
//
void Output::OutputDataFilename(char *filename, int iQuasi,
                                std::vector<int> flags) {
  // check size of vector. If it is 1 then create default filename
  // as it is not called from inside the conjugate gradient method
  if (flags.size() == 1) {
    sprintf(filename, "output_quasi_%d", iQuasi);
    return;
  } else if (flags.size() == 0) {
    d_print("Check the flag vector supplied to performOutput()\n");
    exit(EXIT_FAILURE);
    return;
  }

  if (flags.size() < 5) {
    d_print("Check the flag vector supplied to performOutput()\n");
    exit(EXIT_FAILURE);
    return;
  }

  int minimization_method = flags[0];
  int freq_or_position_min = flags[1];
  int cg_iteration_number = flags[2];
  int alternate_iteration_number = flags[3];
  int cg_status = flags[4];
  int load_number = flags[5];

  switch (minimization_method) {
  case 0: {
    // miniming energy wrt u and w simultaneously
    if (cg_status == 0) {
      // inside cg iteration
      sprintf(filename,
              "cgmethod_0_load_%d_iteration_data_cg_iter_%d_iQuasi_%d",
              load_number, cg_iteration_number, iQuasi);

    } else if (cg_status == 1) {
      // converged cg
      sprintf(filename, "cgmethod_0_load_%d_converged_cg_quasi_%d", load_number,
              iQuasi);
    } else {
      // converged final
      sprintf(filename, "cgmethod_0_load_%d_converged_cg_final_quasi_%d",
              load_number, iQuasi);
    }

    return;
  } break;

  case 1: {
    // miniming energy wrt u and w alternatively
    if (cg_status == 0) {
      // inside cg iteration
      if (freq_or_position_min == 0) {
        // freq minimization in alternate method
        sprintf(filename, "cgmethod_1_load_%d_iteration_data_cg_iter_%d_"
                          "alternate_iter_%d_w_min_quasi_%d",
                load_number, cg_iteration_number, alternate_iteration_number,
                iQuasi);
      } else {
        // position minimization in alternate method
        sprintf(filename, "cgmethod_1_load_%d_iteration_data_cg_iter_%d_"
                          "alternate_iter_%d_u_min_quasi_%d",
                load_number, cg_iteration_number, alternate_iteration_number,
                iQuasi);
      }
    } else if (cg_status == 1) {
      // converged cg
      if (freq_or_position_min == 0) {
        // freq minimization in alternate method
        sprintf(
            filename,
            "cgmethod_1_load_%d_converged_cg_alternate_iter_%d_w_min_quasi_%d",
            load_number, alternate_iteration_number, iQuasi);
      } else {
        // position minimization in alternate method
        sprintf(
            filename,
            "cgmethod_1_load_%d_converged_cg_alternate_iter_%d_u_min_quasi_%d",
            load_number, alternate_iteration_number, iQuasi);
      }
    } else if (cg_status == 2) {
      // converged alternate
      sprintf(filename,
              "cgmethod_1_load_%d_converged_cg_alternate_iter_%d_quasi_%d",
              load_number, alternate_iteration_number, iQuasi);
    } else {
      // converged final
      sprintf(filename, "cgmethod_1_load_%d_converged_cg_final_quasi_%d",
              load_number, iQuasi);
    }

    return;
  } break;

  case 2: {
    // miniming energy wrt w
    if (cg_status == 0) {
      // inside cg iteration
      sprintf(filename, "cgmethod_2_load_%d_iteration_data_cg_iter_%d_quasi_%d",
              load_number, cg_iteration_number, iQuasi);

    } else if (cg_status == 1) {
      // converged cg
      sprintf(filename, "cgmethod_2_load_%d_converged_cg_quasi_%d", load_number,
              iQuasi);
    } else {
      // converged final
      sprintf(filename, "cgmethod_2_load_%d_converged_cg_final_quasi_%d",
              load_number, iQuasi);
    }

    return;
  } break;

  case 3: {
    // miniming energy wrt u
    if (cg_status == 0) {
      // inside cg iteration
      sprintf(filename, "cgmethod_3_load_%d_iteration_data_cg_iter_%d_quasi_%d",
              load_number, cg_iteration_number, iQuasi);

    } else if (cg_status == 1) {
      // converged cg
      sprintf(filename, "cgmethod_3_load_%d_converged_cg_quasi_%d", load_number,
              iQuasi);
    } else {
      // converged final
      sprintf(filename, "cgmethod_3_load_%d_converged_cg_final_quasi_%d",
              load_number, iQuasi);
    }

    return;
  } break;

  default: {
    d_print("Check output_flag_vector supplied to performOutput()\n");
    exit(EXIT_FAILURE);
  } break;
  }

  return;
}
}