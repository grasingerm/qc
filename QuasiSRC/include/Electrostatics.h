//
// Electrostatics.h
//

#if !defined(ELECTROSTATICS_H)
#define ELECTROSTATICS_H

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

#include "DataTypes.h"
#include "CrossNeighborList.h"

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class Electrostatics {

    //
    // public methods
    //

  public:
 
    /**
     * @brief getInstance.
     */
    static Electrostatics * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    //
    // isElectrostaticEnable()
    //
    int isElectrostaticEnable(void);

    //
    //  computeGeaometries()
    //
    void computeGeometries(void);

    //
    //  computeElectricField()
    //
    void computeElectricField(bool compute_in_reference);

    //
    //  computeForceOnNodeDueToSite()
    //  
    void computeForceOnNodeDueToSite(int iQuasi,
        cluster_site_data_t             iCSite_data,
        bool                            compute_in_reference);

    //
    //  computeForceOnSite()
    //
    void computeForceOnSite(int iQuasi,
      cluster_site_data_t                       iCSite_data,
      std::vector<double>                       &force,
      double                                    &energy,
      bool                                      compute_in_reference);

    /**
     * @brief Insert fixed charge for all Quasicontinuum and unit cell radius
     *        for exact calculation.
     *
     * @param charge Charge of atoms in Quasicontinua.
     * @param unitCellRadius Number of unit cells to calculate over for exact
     *        calculation.
     * @param elementSize Atomistic element size cutoff.
     * @param electricConstant Electric constant that depends on unit system.
     * @param faceMethod Method for calculating face charges.
     * @param boundaries Electrostatic boundary conditions, xBegin through
     *        zEnd.
     * @param integration_error Integration error tolerance.
     * @param num_shells Number of shells to include in cluster.
     */
    void insertInputData(const std::vector<double> charges,
            const double              atomisticRadius,
            const double              elementSize,
            const double              electricConstant,
            const int                 faceMethod,
            const std::vector<std::pair<int,double> > boundaries,
            const double              integration_error,
            const int                 num_shells); 

    /**
     * @brief Get electric field.
     *
     * @param iQuasi Quasicontinuum Id.
     */
    const std::vector<std::vector<double> > & getElectricField(const int i_quasi) const;

    /**
     * @brief Clear external point charges.
     */
    void clearExternalPointCharges(void);

   /**
     * @brief Clear geometries.
     *
     */
    void clearGeometries(void);     

    /**
     * @brief Compute electric field in reference.
     */
    void computeElectricFieldReference(void);

    /**
     * @brief Insert external point charges.
     *
     * @param externalCharge Vector of external charge locations and strength.
     */
    void insertExternalPointCharges(std::vector<std::pair<double,std::vector<double> > > externalCharge);    

    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    Electrostatics();

    /**
     * @brief Copy constructor.
     */
    Electrostatics(Electrostatics const&);

    /**
     * @brief Assignment operator.
     */
    const Electrostatics & operator=(const Electrostatics &);

    /**
     * @brief Destructor.
     */
    ~Electrostatics();

    //
    //  overloading the constructor
    //
    Electrostatics(int          default_electrostatics_disabled);

    /**
     * @brief Compute faces.
     *
     * @param P_element_list Pointer to element list.
     */
    void computeFaces(struct element_list_t *P_element_list);    

    /**
     * @brief Compute polarization sites.
     *
     * @param P_element_list Pointer to element list.
     */
    void computePolarizationSites(struct element_list_t *P_element_list,
        struct lattice_t                                *P_lattice);

    /**
     * @brief Compute atomistic elements.
     *
     * @param P_element_list Pointer to element list.
     */
    void computeAtomisticElements(struct element_list_t *P_element_list);

    /**
     * @brief Compute atomistic element locations.
     *
     * @param P_element_list Pointer to element list.
     */
    void computeAtomisticElementLocations(struct element_list_t *P_element_list);

    /**
     * @brief Compute atomistic nodes.
     *
     * @param P_node_list Pointer to node list.
     */
    void computeAtomisticNodes(struct node_list_t *P_node_list);

    /**
     * @brief Compute atomistic charges.
     */
    void computeAtomisticCharges(struct lattice_t *P_lattice,
                 struct element_list_t *P_element_list);

    /**
     * @brief Compute element polarization locations.
     */
    void computeElementPolarizationLocations();

    /**
     * @brief Compute element polarization locations.
     *
     * @param P_element_list Pointer to element list.
     */
    void computeElementPolarization(struct element_list_t *P_element_list);

    /**
     * @brief Compute polarization.
     *
     * @param P_element_list Pointer to base lattice element list.
     */
    void computePolarization(struct element_list_t *P_element_list);

    /**
     * @brief Compute charge densities.
     *
     * @param P_node_list Pointer to base lattice node list.
     */
    void computeChargeDensities(struct node_list_t *P_node_list,
                struct lattice_t   *P_lattice);

    /**
     * @brief Compute atomistic charge locations.
     */
    void computeAtomisticChargeState();

    /**
     * @brief Update cluster atom locations.
     */
    void UpdateClusterAtomLocations();

    /**
     * @brief Compute exact nodal electric field.
     */
    void computeElectricFieldDirectNodal();

    /**
     * @brief Compute contribution to field and potential from cluster
     *
     * @param compute_in_reference Flag to determine whether to compute in reference or not
     */
    void ComputeContributionFromClusters(const bool &compute_in_reference);

    /**
     * @brief Compute contribution to field and potential from cluster surfaces
     *
     * @param compute_in_reference Flag to determine whether to compute in reference or not
     */
    void ComputeContributionFromClusterSurface(const bool &compute_in_reference);

    /**
     * @brief Compute base node list hash table.
     *
     * @param P_node_list Pointer to base node list.
     */
    void ComputeBaseNodeHashList(struct node_list_t *P_node_list);

    /**
     * @brief Compute contribution to field and potential from charged faces.
     *
     * @param compute_in_reference Flag to determine whether to compute in reference or not
     */
    void ComputeContributionFromChargedFaces(const bool &compute_in_reference);

    /**
     * @brief Compute contribution to field and potential from fully atomistic atoms.
     *
     * @param compute_in_reference Flag to determine whether to compute in reference or not
     */
    void ComputeContributionFromAtomistic(const bool &compute_in_reference);

    /**
     * @brief Compute contribution to field and potential from fully external loads.
     *
     * @param compute_in_reference Flag to determine whether to compute in reference or not
     */
    void ComputeContributionFromExternal(const bool &compute_in_reference);    

    //
    // private data types
    //
  private:

    static Electrostatics*                        _instance;
    std::vector<double>                             d_fixedCharge;
    std::vector<std::vector<std::vector<int> > >    d_faces;
    std::vector<std::vector<std::pair<int,int> > >  d_faceElementList;
    std::vector<std::vector<int> >                  d_elementPolarizationSites;
    std::vector<std::vector<std::pair<std::vector<double>,double> > > d_elementPolarizationLocations;
    std::vector<std::vector<std::pair<double,std::vector<std::vector<double> > > > > d_chargeDensity;
    std::vector<std::vector<double> >               d_elementPolarization;
    std::vector<std::pair<double,std::vector<double> > >  d_externalPointCharges;
    std::vector<std::vector<std::vector<int> > >    d_atomisticCharges;
    std::vector<std::vector<std::vector<double> > > d_atomisticChargeState;
    std::vector<int>                                d_atomisticElements;
    std::vector<std::pair<int,std::vector<int> > >  d_atomisticNodes;
    std::vector<std::vector<std::vector<double> > > d_atomisticElementLocations;
    double                                          d_atomisticRadiusSquared;
    double                                          d_atomisticRadius;
    double                                          d_atomisticElementSize;
    std::vector<std::vector<std::vector<double> > > d_electricField;
    std::vector<std::vector<std::vector<double> > > d_electricFieldReference;
    std::vector<std::vector<double> >               d_electricPotential;
    std::vector<std::vector<double> >               d_electricPotentialReference;
    std::vector<std::vector<int> >                  d_baseElementList;
    //std::vector<std::vector<int> >                  d_quasiElementList;
    //std::vector<std::vector<int> >                  d_quasiElementListReference;
    double                                          d_electricConstant;
    int                                             d_faceMethod;
    std::vector<std::vector<std::vector<double> > > d_exactField;
    std::vector<std::vector<double> >               d_exactPotential;
    std::vector<std::pair<int,double> >             d_electroBoundaries;
    int                                             d_num_shells;
    double                                          d_integration_error;
    std::vector<std::vector<std::pair<std::vector<int>,int > > > d_base_node_hash_list;
    int                                             d_electrostaticsEnabled; 
    int                                             d_baseQuasi;
  };

}

#endif // ELECTROSTATICS_H
