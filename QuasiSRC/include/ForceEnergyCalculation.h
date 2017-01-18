#if !defined(FORCEENERGYCALCULATION_H)
#define FORCEENERGYCALCULATION_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <vector>

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class ForceEnergyCalculation {

public:
  /**
   * @brief getInstance.
   */
  static ForceEnergyCalculation *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  void setAtomisticDeadLoads(std::vector<std::vector<double>> forces,
                             const double element_size);

  void noAtomisticLoads(void);

  void setRemoveResidualFlags(const std::vector<int> flags);

  //
  //  removeResidualForces()
  //  This computes the force in reference condition for current mesh
  //  and depending on remove flags, it puts the forces in external force
  //  variable of node.
  //
  void removeResidualForces(void);

  void computeForceEnergy(int rebuild_neighbor_flag, bool compute_in_reference);

  //
  //  iQuasiForceEnergy()
  //
  //  Depending on potential and Electrostatics enabled in system,
  //  this function calls function corresponding to EAM, Pairwise interaction
  //  and Electrostatics.
  //
  void iQuasiForceEnergy(int iQuasi, bool compute_in_reference);

  void iQuasiForceEnergyHarmonic(int iQuasi, bool compute_in_reference);

  void iQuasiAddForceAndEnergyToNodes(int iQuasi);

  //
  //  iQuasiEAMForceEnergy()
  //
  //  Computes force and energy due to EAM potential and adds the
  //  contribution to relevant nodes.
  //
  void iQuasiEAMForceEnergy(int iQuasi);

  //
  //  iQuasiPairwiseForceEnergy()
  //
  //  Computes force and energy due to pairwise interaction and adds the
  //  contribution to relevant nodes.
  //
  void iQuasiPairwiseForceEnergy(int iQuasi);

  void iQuasiPairwiseForceEnergyHarmonic(int iQuasi);

  void iQuasiEntropyForceEnergy(int iQuasi);

  void iQuasiEntropyForceEnergyHarmonic(int iQuasi);

  //
  //  iQuasiElectrostaticsForceEnergy()
  //
  //  Computes force due to electrostatics and add the
  //  contribution to relevant nodes.
  //
  void iQuasiElectrostaticsForceEnergy(int iQuasi, bool compute_in_reference);

  void iQuasiAddExternalForces(int iQuasi);

  //
  //  iQuasiComputeTotalEnergy()
  //
  //  adds the energy of all the nodes to all_node_list data of iQuasi
  //
  void iQuasiComputeTotalEnergy(int iQuasi);

  //
  //  iQuasiComputeExternalEnergy()
  //
  //  Computes energy due to external forces and adds it to all_node_list data
  //  of iQuasi
  //
  void iQuasiComputeExternalEnergy(int iQuasi);

  void iQuasiAddEnergy(int iQuasi);

  void copyEnergyToInitialEnergy(void);

  void testForceEnergyCalculation(
      int rebuild_neighbor_flag, int rebuild_cluster_flag,
      std::vector<std::vector<
          std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
          &test_cluster_force_energy,
      std::vector<
          std::vector<std::pair<std::vector<double>, std::vector<double>>>>
          &test_node_force_energy,
      std::vector<
          std::vector<std::vector<std::pair<std::vector<double>, double>>>>
          &test_total_cluster_force_energy,
      std::vector<std::vector<std::pair<std::vector<double>, double>>>
          &test_total_node_force_energy);

  //
  //  testAddClusterForceEnergy()
  //
  //  flag : 0 - for eam
  //         1 - for pairwise
  //         2 - for entropy
  //         3 - electrostatics
  //
  void testAddForceEnergy(
      std::vector<std::vector<
          std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
          &test_cluster_force_energy,
      std::vector<
          std::vector<std::pair<std::vector<double>, std::vector<double>>>>
          &test_node_force_energy,
      std::vector<
          std::vector<std::vector<std::pair<std::vector<double>, double>>>>
          &test_total_cluster_force_energy,
      std::vector<std::vector<std::pair<std::vector<double>, double>>>
          &test_total_node_force_energy,
      int flag);

  std::vector<std::vector<std::pair<std::vector<double>, double>>> &
  getQuasiClusterSiteForceData(int iQuasi);

  std::vector<std::vector<double>> &getQuasiClusterTraceOfK(int iQuasi);

  std::vector<std::vector<std::pair<std::vector<double>, double>>>
  iQuasiTestPairwiseForceEnegry(int iQuasi);

  std::vector<std::vector<std::pair<std::vector<double>, double>>>
  iQuasiTestEntropyForceEnegry(int iQuasi);

  std::vector<std::vector<std::pair<std::vector<double>, double>>>
  iQuasiTestElectrostaticsForceEnegry(int iQuasi);

private:
  /**
   * @brief Constructor.
   */
  ForceEnergyCalculation();

  /**
   * @brief Copy constructor.
   */
  ForceEnergyCalculation(ForceEnergyCalculation const &);

  /**
   * @brief Assignment operator.
   */
  const ForceEnergyCalculation &operator=(const ForceEnergyCalculation &);

  /**
   * @brief Destructor.
   */
  ~ForceEnergyCalculation();

  void SetResidualFlagsForAdjacentNodes(const struct node_t *P_node,
                                        std::vector<bool> &node_residual_flags);

  void InitializeClusterForceEnergyData();

  void InitializeNewBucketLock();

  void InitializeTraceOfKData();
  
  void InitializeTestForceEnergyData(
      std::vector<std::vector<
          std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
          &test_cluster_force_energy,
      std::vector<
          std::vector<std::pair<std::vector<double>, std::vector<double>>>>
          &test_node_force_energy,
      std::vector<
          std::vector<std::vector<std::pair<std::vector<double>, double>>>>
          &test_total_cluster_force_energy,
      std::vector<std::vector<std::pair<std::vector<double>, double>>>
          &test_total_node_force_energy);

  void AddClusterSiteForcesToNodes(int iQuasi);

  void ComputeTraceOfKAtNodes(int iQuasi);

private:
  static ForceEnergyCalculation *_instance;
  int d_atomisticFlag;
  std::vector<std::vector<double>> d_atomisticDeadLoads;
  double d_atomisticElementSize;
  std::vector<bool> d_residualForcesFlags;
  std::vector<std::vector<std::vector<std::pair<std::vector<double>, double>>>>
      d_clusterForceEnergy;
  std::vector<std::vector<std::vector<double>>> d_traceOfK;
};
}

#endif // FORCEENERGYCALCULATION_H
