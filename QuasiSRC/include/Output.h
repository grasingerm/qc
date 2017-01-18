//
// Output.h
//

#if !defined(OUTPUT_H)
#define OUTPUT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <vector>
#include "DataTypes.h"

// definitions
#define NO_OUTPUT_FLAG 0x00
#define FULL_LATTICE_OUTPUT_FLAG 0x01
#define NODE_OUTPUT_FLAG 0x02
#define RESTART_FLAG 0x04
#define INPUT_FLAG 0x08
#define CLUSTER_FLAG 0x10
#define INDENT_FLAG 0x20
#define CROSSNEIGHBORLIST_CLUSTER_FLAG 0x40
#define CROSSNEIGHBORLIST_NEIGHBOR_FLAG 0x80
#define CLUSTER_FORCE_OUTPUT 0x100

#define CONVERGED_OUTPUT 0
#define ITERATION_OUTPUT 1
#define CONVERGED_INTERMEDIATE 2
#define ITERATION_COMBINED_OUTPUT 3
#define ITERATION_FREQUENCY_OUTPUT 4
#define ITERATION_POSITION_OUTPUT 5
#define CONVERGED_INTERMEDIATE_COMBINED 6
#define CONVERGED_INTERMEDIATE_FREQUENCY 7
#define CONVERGED_INTERMEDIATE_POSITION 8

// tecplot parameters
#define OUTPUT_RET_OK 0
#define OUTPUT_RET_NOK -1

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Output {

  //
  // public methods
  //

public:
  static FILE *d_debugFile;
  char *d_debugFilename;

  /**
   * @brief getInstance.
   */
  static Output *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  //
  //  createOutputDirectory()
  //
  void createOutputDirectory(int outputFlag, char *outputDirectory);

  //
  //  performOutputCentral()
  //
  void performOutputCentral(unsigned int flags,
                            std::vector<int> output_flag_vector);

  /**
   * @brief performOutput()
   */
  void performOutput(const int iQuasi, unsigned int flags,
                     std::vector<int> output_flag_vector);

  /**
   * @brief performNodeListOutput()
   */
  void performNodeListOutput(const char *file_name,
                             const struct node_list_t *P_node_list);

  /**
   * @brief performNodeListOutput()
   */
  int tecplotOutput(struct all_node_list_t *node_list,
                    struct atom_list_t *atom_list,
                    struct element_list_t *element_list,
                    struct indentor_t *P_indentor, const char *title,
                    const char *fem_filename, const char *atom_filename,
                    const char *elem_filename, const char *energy_filename,
                    const double time, int i_load);

  void printCluster(const struct cluster_atom_list_t cluster,
                    const char *file_name, const int i_cluster);

  void printNodeNeighbors(const struct node_list_t node_list,
                          const struct atom_list_t atom_list,
                          const char *file_name, const int i_node);

  //
  //  testEnergyOutput()
  //
  void testEnergyOutput(
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
          *test_total_node_force_energy);

  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  Output();

  /**
   * @brief Copy constructor.
   */
  Output(Output const &);

  /**
   * @brief Assignment operator.
   */
  const Output &operator=(const Output &);

  /**
   * @brief Destructor.
   */
  ~Output();

  //
  // private methods
  //

private:
  char *AllocateFilename(void);

  void PerformLatticeOutput(const int iQuasi, const char *dirname,
                            struct lattice_t *P_lattice, int mark,
                            int cg_info_flag, int iteration_number);

  void PerformNodeOutput(const int iQuasi, const char *dirname,
                         const struct node_list_t *P_node_list,
                         const struct element_list_t *P_element_list,
                         std::vector<int> output_flag_vector);

  void PerformInputFormat(const int iQuasi, const char *dirname,
                          struct node_list_t *P_node_list,
                          struct lattice_t *P_lattice, char *materialName,
                          int mark, int cg_info_flag, int iteration_number);

  void PerformClusterOutput(const int iQuasi, const char *dirname,
                            struct node_list_t *P_node_list,
                            struct lattice_t *P_lattice, int mark,
                            int cg_info_flag, int iteration_number);

  void PerformIndentOutput(const int iQuasi, const char *dirname,
                           const struct indentor_t *P_indentor, int mark,
                           int cg_info_flag, int iteration_number);

  void PerformCrossNeighborListClusterOutput(const int iQuasi,
                                             const char *dirname, int mark,
                                             int cg_info_flag,
                                             int iteration_number);

  void PerformCrossNeighborListNeighborOutput(const int iQuasi,
                                              const char *dirname, int mark,
                                              int cg_info_flag,
                                              int iteration_number);

  void PerformClusterForcesOutput(const int iQuasi, const char *dirname,
                                  int mark, int cg_info_flag,
                                  int iteration_number);

  void PerformRestartOutput(const int iQuasi, const char *dirname,
                            struct all_node_list_t *P_all_node_list,
                            struct element_list_t *P_element_list,
                            struct lattice_t *P_lattice,
                            struct indentor_t *P_indentor, char *materialName,
                            std::vector<int> output_flag_vector);

  void OutputDataFilename(char *filename, int iQuasi,
                          std::vector<int> output_flag_vector);

  //
  // private data types
  //
private:
  static Output *_instance;
  char *d_outputDirectory;
};
}

#endif // OUTPUT_H
