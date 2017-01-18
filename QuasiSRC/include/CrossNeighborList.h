//
// CrossNeighborList.h
//

#if !defined(CROSSNEIGHBORLIST_H)
#define CROSSNEIGHBORLIST_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <vector>
#include "DataTypes.h"

// data definition
typedef std::vector<int> vec_int_t;
typedef std::vector<double> vec_dbl_t;

typedef std::pair<vec_int_t, vec_dbl_t> pair_state_t;
typedef std::vector<pair_state_t> vec_state_t;

typedef std::pair<pair_state_t, vec_int_t> cluster_site_data_t;

typedef std::pair<pair_state_t, std::vector<std::vector<vec_int_t>>>
    neigh_site_data_t;

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class CrossNeighborList {

public:
  /**
   * @brief getInstance.
   */
  static CrossNeighborList *getInstance();

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
                                const int rebuild_cluster_flag);

  void printCrossNeighborList(const int thread_flag);

  std::vector<std::vector<neigh_site_data_t>> &
  getQuasiNeighborData(const int iQuasi);

  std::vector<std::vector<std::vector<cluster_site_data_t>>> &
  getClusterData(void);

  std::vector<std::vector<cluster_site_data_t>> &
  getQuasiClusterData(const int iQuasi);

  std::vector<int> getNodeInfoOfClusterSite(int iQuasi,
                                            std::vector<int> clusterSite);

  int getNumberBuckets();
  int getSizeOfClusterData();

private:
  CrossNeighborList();
  CrossNeighborList(CrossNeighborList const &);
  const CrossNeighborList &operator=(const CrossNeighborList &);
  ~CrossNeighborList();
  void UpdateLocation(enum mt_version_t mt_version);
  void ComputeNeighborList(enum mt_version_t mt_version,
                           const int rebuild_cluster_flag,
                           const int rebuild_neighbor_flag);
  void BuildClusterAndComputeNeighborList(enum mt_version_t mt_version,
                                          int iQuasi);
  void ComputeNeighborListUsingCurrentClusterData(enum mt_version_t mt_version,
                                                  int iQuasi);

private:
  static CrossNeighborList *_instance;
  //  vector to hold all neighborlist
  //
  //  d_neighborList[iQuasi[iBucket][iN].pair(pair_state_t,
  //                                          quasi_cluster_data_t)
  //
  std::vector<std::vector<std::vector<neigh_site_data_t>>> d_neighborList;

  std::vector<std::vector<std::vector<cluster_site_data_t>>> d_clusterData;
};
}

#endif // CROSSNEIGHBORLIST_H
