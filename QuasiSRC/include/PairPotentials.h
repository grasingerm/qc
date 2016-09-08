//
// PairPotentials.h
//

#if !defined(PAIRPOTENTIALS_H)
#define PAIRPOTENTIALS_H

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
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for PairPotentials instance
   */
  class PairPotentials {

    //
    // public data types
    //
  public:

    bool d_EAMFlag;
    int  d_minMethod;
    int  d_quadMethod;
    int  d_statistics;    

    //
    // public methods
    //

  public:
 
    /**
     * @brief getInstance.
     */
    static PairPotentials * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    /**
     * @brief Insert a Potential.
     *
     * @param quasicontinuumPair Pair of Quasicontinuum that interact.
     * @param potentialType Type of potential interaction.
     * @param potentialParameters Parameters for potential.
     * @param EAMFlag True/false flag for EAM potential.
     */
    void insertNonEAMPotential(std::pair<int,int>  quasicontinuumPair,
		int                 potentialType,
		std::vector<double> potentialParameters,
		const bool          EAMFlag);

    /**
     * @brief Insert a EAMData Potential.
     *
     * @param quasicontinuumPair Pair of Quasicontinuum that interact.
     * @param potentialType Type of potential interaction.
     * @param potentialParameters Parameters for potential.
     * @param dataFile1 Datafile 1 for potential.
     * @param NGridPoints Number of grid point of data.
     * @param EAMFlag True/false flag for EAM potential.
     */
    void insertEAMPotential(std::pair<int,int>  quasicontinuumPair,
			 int                 potentialType,
			 std::vector<double> potentialParameters,
			 const char*         dataFile1,
			 const int           NGridPoints,
			 const bool          EAMFlag);

   /**
     * @brief Get force between pair of atoms.
     *
     * @param potentialNumber Potential number of interaction.
     * @param separation Distance between atoms.
     *
     * @return Force and energy on an atom.
     */
    std::pair<double,double> getPairForceAndEnergy(int    potentialNumber,
					  double separation);

    //
    //  getHarmonicApproxData()
    //  return energy, first derivative, second derivative for pairwise
    //  potential
    //
    void getHarmonicApproxData(double       &phi,
            double                          &dphi,
            double                          &d_dphi,
            int                             potentialNumber,
            double                          separation);

    /**
     * @brief Find potentialNumber for pair of Quasicontinuum, if they interact.
     *
     * @param quasicontinuumIdOne First Quasicontinuum Id.
     * @param quasicontinuumIdTwo Second Quasicontinuum Id.
     *
     * @return potentialNumber if pairs interact, -1 if pairs don't interact.
     */
    int doQuasicontinuumInteract(int quasicontinuumIdOne,
				 int quasicontinuumIdTwo);

    /**
     * @brief Get cutoff radius for potential.
     *
     * @param potentialNumber Potential number for Quasicontinuum
     *
     * @return cutoff radius.
     */
    double getCutoffRadius(int potentialNumber);

    //
    //  getCutoffRadiusNeighborList()
    //
    double getCutoffRadiusNeighborList(const int potentialNumber);

    //
    //  getCutoffRadiusClusterList()
    //
    double getCutoffRadiusClusterList(const int potentialNumber);

    /**
     * @brief Density from r separation distance.
     *
     * @param potentialNumber Potential number for Quasicontinuum interacting
     * @param r Separation distance.
     *
     * @return density and its derivative.
     */
    std::pair<double,double> getDensityAndDerivative(int    potentialNumber,
					double r);

    /**
     * @brief Embedding function from site density.
     *
     * @param potentialNumber Potential number for Quasicontinuum interacting
     * @param density Density at site.
     *
     * @return embedding function and its derivative.
     */
    std::pair<double,double> getEmbeddingFunctionAndDerivative(int    potentialNumber,
						  double density);

    //
    //  setUniversalConstants()
    //  constants : (max planck, boltzman, dielectric)
    void setUniversalConstants(std::vector<double>  constants);

    //
    //  getUniversalConstants()  
    //
    std::vector<double> getUniversalConstants();

    //
    //  getBoltzmanConstant()
    //
    double getBoltzmanConstant();

    //
    //  getBoltzmanConstant()
    //
    double getElectricConstant();    

      //
      //  getPotentialType()
      //
    int getPotentialType(int potentialNumber);    

    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    PairPotentials();

    /**
     * @brief Copy constructor.
     */
    PairPotentials(PairPotentials const&);

    /**
     * @brief Assignment operator.
     */
    const PairPotentials & operator=(const PairPotentials &);

    /**
     * @brief Destructor.
     */
    ~PairPotentials();

    /**
     * @brief Find force and energy for LJ interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> LJ(int    potentialNumber,
				double r);

   /**
     * @brief Find force and energy for EAMJohnson interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> EAMJohnson(int    potentialNumber,
					double r);

   /**
     * @brief Find force and energy for EAMData interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> EAMData(int    potentialNumber,
				     double r);

    /**
     * @brief Find force and energy for LJNiMn interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> LJNiMn(int    potentialNumber,
				    double r);

    /**
     * @brief Find force and energy for Buckingham interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> Buckingham(int    potentialNumber,
					double r);

    /**
     * @brief Find force and energy for Harmonic interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> Harmonic(int    potentialNumber,
				      double r);

    /**
     * @brief Find force and energy for Anharmonic interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    std::pair<double,double> Anharmonic(int    potentialNumber,
                      double r);    

    /**
     * @brief Find energy, first derivative and second derivative
     * for LJ interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     */
    void LJ_Harmonic(double       &phi,
            double                          &dphi,
            double                          &d_dphi,
            int                             potentialNumber,
            double                          separation);

    /**
     * @brief Find energy, first derivative and second derivative
     * LJNiMn interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    void LJNiMn_Harmonic(double       &phi,
            double                          &dphi,
            double                          &d_dphi,
            int                             potentialNumber,
            double                          separation);

    /**
     * @brief Find energy, first derivative and second derivative
     * Buckingham interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    void Buckingham_Harmonic(double       &phi,
            double                          &dphi,
            double                          &d_dphi,
            int                             potentialNumber,
            double                          separation);

    /**
     * @brief Find energy, first derivative and second derivative
     * Harmonic interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    void Harmonic_Harmonic(double       &phi,
            double                          &dphi,
            double                          &d_dphi,
            int                             potentialNumber,
            double                          separation);

    /**
     * @brief Find energy, first derivative and second derivative
     * Anharmonic interaction.
     *
     * @param potentialNumber Potential number to get parameters.
     * @param r Distance between atoms.
     *
     * @return force (first) and energy (second) values.
     */
    void Anharmonic_Harmonic(double       &phi,
            double                          &dphi,
            double                          &d_dphi,
            int                             potentialNumber,
            double                          separation);

    //
    // private data types
    //
  private:

    static PairPotentials*              _instance;
    std::vector< std::pair <int,int> >  d_PotentialInteractions;
    std::vector<int>                    d_PotentialType;
    std::vector< std::vector <double> > d_PotentialParameters;
    std::vector< double >               d_boltzmanConstant;
    std::vector< std::vector <double> > d_RhoX;
    std::vector< std::vector <double> > d_RhoY;
    std::vector< std::vector <double> > d_RhoB;
    std::vector< std::vector <double> > d_RhoC;
    std::vector< std::vector <double> > d_RhoD;
    std::vector< std::vector <double> > d_EmbedX;
    std::vector< std::vector <double> > d_EmbedY;
    std::vector< std::vector <double> > d_EmbedB;
    std::vector< std::vector <double> > d_EmbedC;
    std::vector< std::vector <double> > d_EmbedD;
    std::vector< std::vector <double> > d_PairX;
    std::vector< std::vector <double> > d_PairY;
    std::vector< std::vector <double> > d_PairB;
    std::vector< std::vector <double> > d_PairC;
    std::vector< std::vector <double> > d_PairD;

    // universal constants
    std::vector<double>                    c_universalConst;
  };

}

#endif // PAIRPOTENTIALS_H
