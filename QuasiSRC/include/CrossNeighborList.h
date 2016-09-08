//
// CrossNeighborList.h
//

#if !defined(CROSSNEIGHBORLIST_H)
#define CROSSNEIGHBORLIST_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

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

// data definition
typedef std::vector<int>  vec_int_t;
typedef std::vector<double>  vec_dbl_t;

typedef std::pair<vec_int_t, vec_dbl_t> pair_state_t;
typedef std::vector<pair_state_t> vec_state_t;

typedef std::pair<pair_state_t, vec_int_t > cluster_site_data_t;

typedef std::pair<pair_state_t, std::vector< std::vector< vec_int_t > > > neigh_site_data_t;

// global variables specific to CrossNeighborList

//
//
//
namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class CrossNeighborList {

    //
    // public methods
    //

  public:
 
    /**
     * @brief getInstance.
     */
    static CrossNeighborList * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    //
    //  computeCrossNeighborList()
    //
    //  rebuild_neighbor_flag = 0 : don't rebuild neighbors, just update state
    //                        = 1 : rebuild new neighbors
    //  rebuild_cluster_flag  = 0 : use old cluster data, update state
    //                        = 1 : built new cluster data from node data
    //
    void computeCrossNeighborList(const int rebuild_neighbor_flag,
              const int           rebuild_cluster_flag);

    //
    //  printCrossNeighborList()
    //
    void printCrossNeighborList(const int thread_flag);

    //
    //  getQuasiNeighborData()
    //
    std::vector< std::vector< neigh_site_data_t > > & 
    getQuasiNeighborData(const int iQuasi);

    //
    //  getClusterData()
    //
    std::vector< std::vector< std::vector< cluster_site_data_t > > > &
    getClusterData(void);

    //
    //  getQuasiClusterData()  
    //
    std::vector< std::vector< cluster_site_data_t > >  &
    getQuasiClusterData(const int iQuasi);    

    //
    //  getNodeInfoOfClusterSite()
    //
    std::vector<int> getNodeInfoOfClusterSite(int iQuasi, std::vector<int> clusterSite);

    //
    //  getNumberBuckets()
    //
    int getNumberBuckets();

    //
    //  getSizeOfClusterData()
    //
    int getSizeOfClusterData();

    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    CrossNeighborList();

    /**
     * @brief Copy constructor.
     */
    CrossNeighborList(CrossNeighborList const&);

    /**
     * @brief Assignment operator.
     */
    const CrossNeighborList & operator=(const CrossNeighborList &);

    /**
     * @brief Destructor.
     */
    ~CrossNeighborList();

    //
    //  UpdateLocation()
    //
    void UpdateLocation(enum mt_version_t mt_version);

    //
    //  ComputeNeighborList()
    //
    void ComputeNeighborList(enum mt_version_t  mt_version,
            const int         rebuild_cluster_flag,
            const int         rebuild_neighbor_flag);

    //
    //  BuildClusterAndComputeNeighborList()
    //
    void BuildClusterAndComputeNeighborList(enum mt_version_t mt_version,
            int         iQuasi);

    //
    //  ComputeNeighborListUsingCurrentClusterData()
    //
    void ComputeNeighborListUsingCurrentClusterData(enum mt_version_t mt_version,
            int         iQuasi);

    //
    // private data types
    //
  private:

    static CrossNeighborList*                        _instance;
    //  vector to hold all neighborlist 
    //
    //  d_neighborList[iQuasi[iBucket][iN].pair(pair_state_t, 
    //                                          quasi_cluster_data_t)
    //
    std::vector< std::vector< std::vector< neigh_site_data_t > > > d_neighborList;  

    std::vector< std::vector< std::vector< cluster_site_data_t > > >  d_clusterData;
  };

}

#endif // CROSSNEIGHBORLIST_H