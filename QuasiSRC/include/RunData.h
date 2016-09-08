//
// RunData.h
//

#if !defined(RUNDATA_H)
#define RUNDATA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#include <iostream>
#include <vector>

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class RunData {

    //
    // public methods
    //

  public:

    bool d_timeOut;
 
    /**
     * @brief getInstance.
     */
    static RunData * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    //
    //  initializeData()
    //
    //  Only initialize if we want to compute the time required by
    //  varioous elements of code
    //
    //  Call this only after initializing the Quasicontinua class
    //
    void initializeData();

    //
    //  addComputeClusterAndNeighborListTime()
    //
    void addComputeClusterAndNeighborListTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addBuildNeighborListTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addUpdateClusterAndNeighborListTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addComputePairwiseTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addComputeEAMTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addComputeElectricFieldTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addComputeElectrostaticsTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addComputeEntropyTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addClusterForceToNodeForceTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addInitialToFirstForceCalculationTime(int iQuasi, double time);

    //
    //  addComputeResidualForceTime()
    //
    void addComputeResidualForceTime(int iQuasi, double time);

    //
    //  addBuildNeighborListTime()
    //
    void addNumberOfClusterSites(int iQuasi, int n);  

    //
    //  writeData()
    //
    std::vector<int> writeData(std::vector<std::vector<double>> & data);


    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    RunData();

    /**
     * @brief Copy constructor.
     */
    RunData(RunData const&);

    /**
     * @brief Assignment operator.
     */
    const RunData & operator=(const RunData &);

    /**
     * @brief Destructor.
     */
    ~RunData();

    //
    // private data types
    //
    //  Following data will be sum of time for all quasi. This sum of time
    //  will be put in zeroth element of vector.
    //  d_computeElectricFieldTime
    //  d_initialToFirstForceCalculationTime
    //  d_computeResidualForceTime
    //  
  private:

    static RunData*                        _instance;
    int                                    d_numQuasi;
    int                                    d_numberThreads;
    std::vector<double>                    d_computeClusterAndNeighborListTime;
    std::vector<double>                    d_buildNeighborListTime;
    std::vector<double>                    d_updateClusterAndNeighborListTime;
    std::vector<double>                    d_computePairwiseTime;
    std::vector<double>                    d_computeEAMTime;
    std::vector<double>                    d_computeElectricFieldTime;
    std::vector<double>                    d_computeElectrostaticsTime;
    std::vector<double>                    d_computeEntropyTime;
    std::vector<double>                    d_clusterForceToNodeForceTime;
    std::vector<double>                    d_initialToFirstForceCalculationTime;
    std::vector<double>                    d_computeResidualForceTime;
    std::vector<int>                    d_numberOfClusterSites;
  };

}

#endif // RUNDATA_H
