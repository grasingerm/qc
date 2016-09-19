//
// File:     Quasicontinuum.h
//
#if !defined(QUASICONTINUUM_H)
#define QUASICONTINUUM_H

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

#include <string>

#include "DataTypes.h"

//
//
//
namespace quasicontinuum {

/**
 * @brief Class implementing Quasicontinuum data structures and methods.
 */
class Quasicontinuum {
public:
  /**
   * @brief Constructor.
   *
   * @param temperature of quasi.
   * @param ac Command line options for input.
   * @param av Command line options for input.
   */
  Quasicontinuum(const int id, int ac, char *av[]);

  /**
   * @brief Destructor.
   */
  ~Quasicontinuum();

  /**
   * @brief Accessors for various data structures.
   */
  inline const lattice_t &getLattice() const { return d_lattice; }
  inline lattice_t &getLattice() { return d_lattice; }

  inline const element_list_t &getElementList() const { return d_element_list; }
  inline element_list_t &getElementList() { return d_element_list; }

  inline const all_node_list_t &getNodeList() const { return d_node_list; }
  inline all_node_list_t &getNodeList() { return d_node_list; }

  inline const indentor_t &getIndentor() const { return d_indentor; }
  inline indentor_t &getIndentor() { return d_indentor; }

  inline const char *getMaterialName() const { return d_materialName; }
  inline char *getMaterialName() { return d_materialName; }

  int isRestartOn() const;
  int isRestartOn();

  inline double getAtomicMass() const { return d_mass; }
  inline double getAtomicMass() { return d_mass; }

  /**
   * @brief Get forces on nodes.
   */
  const std::vector<double> getForces();

  /**
   * @brief Get energy.
   */
  inline double getEnergy() { return d_node_list.energy.total; }

  /**
   * @brief Update state of nodes.
   *
   * @param update_direction Direction to move for each dof.
   */
  void updateState(const std::vector<double> &update);

  /**
   * @brief Apply deformation to body.
   *
   * @param flag 0 = whole domain
   *             1 = fixed x positions only
   * @param def Deformation matrix.
   * @param shiftVector Shift vector for quasicontinuum
   */
  void applyDeformation(const int flag,
                        const std::vector< std::vector<double> > &def,
                        const std::vector<double> &shiftVector);

  /**
   * @brief Return number of nodes.
   *
   * @return Number of nodes in lattice.
   */
  inline int getNumberNodes() { return d_node_list.node_list.number_nodes; }

  /**
   * @brief Return current solution.
   *
   * @return Current solution.
   */
  const std::vector<double> getPosition();

  /**
   * @brief Return current frequency.
   *
   * @return Current solution.
   */
  const std::vector<double> getFrequency();

  /**
   * @brief Return current state.
   *
   * @return Current solution.
   */
  const std::vector<double> getState();

  /**
   * @brief Return precondition values.
   *
   * @return Precondition values.
   */
  const std::vector<double> getPreconditionValues();

  /**
    * @brief Get number of sites in node cluster.
    *
    * @param iNode Node to get number of sites for.
    *
    * @return Number of sites in iNode cluster.
    */
  inline int getNumberSites(int iNode) { 
    return d_node_list.node_list.nodes[iNode]->site_cluster.number_sites; 
  }

  /**
    * @brief Get location of iSite in iNode cluster.
    *
    * @param iNode Node which site is in.
    * @param iSite Site to get location of.
    *
    * @return Number of sites in iNode cluster.
    */
  std::pair< std::vector<int>, std::vector<double> >
  getClusterSiteLocationInNodeCluster(int iNode, int iSite);

  /**
    * @brief Get location of iSite in iNode cluster.
    *
    * @param iNode Node which site is in.
    * @param iSite Site to get location of.
    *
    * @return Number of sites in iNode cluster.
    */
  std::pair< std::vector<int>, std::vector<double> >
  getClusterSiteStateInNodeCluster(int iNode, int iSite);

  /**
    * @brief Get location of lattice site.
    *
    * @param site Site to get location of.
    *
    * @return location of site.
    */
  std::vector<double> getLatticeLocation(std::vector<int> site);

  /**
    * @brief Get location of lattice site.
    *
    * @param site Site to get location of.
    *
    * @return location of site.
    */
  std::vector<double> getLatticeState(std::vector<int> site);

  /**
    * @brief Get weights of nodes.
    *
    * @return Vector of node weights.
    */
  const std::vector<double> getWeights(void);

  /**
    * @brief Reset current position to initial position and return
    *        vector of current positions to be added back later.
    *
    * @return Vector of node positions.
    */
  const std::vector< std::vector<double> > resetToInitialPosition(void);

  /**
    * @brief Reset current state to initial state and return
    *        vector of current state to be added back later.
    *
    * @return Vector of node state.
    */
  const std::vector< std::vector<double> > resetToInitialState(void);

  /**
    * @brief Reset current position.
    *
    * @param currentPositions Vector of node positions.
    */
  void resetToCurrentPosition(
      const std::vector< std::vector<double> > currentPositions);

  /**
   * @brief Reset current state.
   *
   * @param currentState Vector of node state.
   */
  void resetToCurrentState(const std::vector< std::vector<double> > currentState);

  /**
    * @brief Get vector of sites within radius of point.
    *
    * @param siteLocation Location to find sites around.
    * @param cutoffRadius Radius to search around point.
    *
    * @return Vector of global state of sites
    */
  std::vector< std::pair< std::vector<int>, std::vector<double> > >
  getNeighborSitesAndStateAroundPoint(
      const std::pair< std::vector<int>, std::vector<double> > siteState,
      const double cutoffRadius);

  /**
    * @brief Get vector of sites within radius of point.
    *
    * @param siteLocation Location to find sites around.
    * @param cutoffRadius Radius to search around point.
    *
    * @return Vector of global state of sites
    */
  std::vector< std::vector<double> > getNeighborStateAroundPoint(
      const std::pair< std::vector<int>, std::vector<double> > siteState,
      const double cutoffRadius);

  //
  //  initializeExternalForces()
  //
  void initializeExternalForces();

  //
  //  initializeNodeEnergy()
  //
  void initializeNodeEnergy(void);

  //
  //  initializeNodeForces(void)
  //
  void initializeNodeForces(void);

  /**
   * @brief Get vector of sites and locations in cluster.
   *
   * @param site_index_guess Site guess for something close to site_location.
   * @param site_location Location to find cluster around.
   * @param cutoff_radius_squared Radius squared of cluster.
   * @param vector_return_sites Site indexes of atoms in cluster.
   * @param vector_return_locations Locations of atoms in cluster.
   */
  void
  getSitesInCluster(const std::vector<int> &site_index_guess,
                    const std::vector<double> &site_location,
                    const double &cutoff_radius_squared,
                    const std::vector<double> global_shift_vector,
                    std::vector< std::vector<int> > &vector_return_sites,
                    std::vector< std::vector<double> > &vector_return_locations);

  //
  //  modifyInitialConfiguration()
  //
  void modifyInitialConfiguration();

  //
  //  set or reset position fixity mask of all nodes
  //
  //  fixity_flag = 0 - reset the fixity mask
  //                1 - set fixity mask so that all
  //                    position dof are fixed
  //                2 - set fixity mask so that all
  //                    position dof are free
  //
  void setPositionFixity(const int fixity_flag);

  //
  //  set or reset frequency fixity mask of all nodes
  //
  //  fixity_flag = 0 - reset the fixity mask
  //                1 - set fixity mask so that all
  //                    frequency dof are fixed
  //                2 - set fixity mask so that all
  //                    frequency dof are free
  //
  void setFrequencyFixity(const int fixity_flag);

  //
  //  getTraceOfKData()
  //
  const std::vector<double> &getTraceOfKData() const;
  std::vector<double> &getTraceOfKData();

  //
  //  initializeTraceOfKNode()
  //
  void initializeTraceOfKData();

  //
  // private data
  //
private:
  int d_id;
  element_list_t d_element_list;
  all_node_list_t d_node_list;
  lattice_t d_lattice;
  indentor_t d_indentor;
  qc_options_t d_qc_options;
  // loads_t          d_loads;
  int d_force_flag;
  int d_energy_flag;
  int d_output_flag;
  mt_version_t d_mt_version;
  char d_materialName[100];
  double d_mass;
  std::vector<unsigned int> d_positionFixityOriginal;
  std::vector<unsigned int> d_frequencyFixityOriginal;
  std::vector<double> d_traceOfKNode;

  //
  // private methods
  //
private:
  /**
   * @brief Get site that is near a site
   *
   * @param site_index_guess Initial site guess.
   * @param return_site Site near initial guess that is in lattice.
   */
  void FindSiteInLatticeNearPoint(const std::vector<int> &site_index_guess,
                                  std::vector<int> &return_site) const;
};
}

#endif // QUASICONTINUUM_H
