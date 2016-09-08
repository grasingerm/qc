//
// CrossNeighborList.cc
//

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

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#include <iostream>
#include <vector>

// to get amount of time required in computing neighbor list
#include <time.h>
#include <sys/time.h>

#include "CrossNeighborList.h"
#include "DataTypes.h"
#include "Quasicontinua.h"
#include "Quasicontinuum.h"
#include "PairPotentials.h"
#include "Node.h"
#include "Element.h"
#include "Lattice.h"
#include "Error.h"
#include "Output.h"
#include "RunData.h"

#include "threads.h"
#include "monitor.h"

// static int l_check[3] = {11, 19, 17};

// 1  complete GetBucketNumberForNeighborListLocal()
// 2  complete GetCounterIndexOfNeighborSite_InNeighborListBucket_Local()
//
//

namespace quasicontinuum {
  //
  //  namespace for neighbor list computing
  //
  namespace
  {
    // cluster site count
    // static int cluster_sites_count = 0;
    // static pthread_mutex_t cluster_count_lock = PTHREAD_MUTEX_INITIALIZER;
    static std::vector<pthread_mutex_t> bucket_lock;
    static int numBuckets = 65535;
    static pthread_mutex_t lock_pint = PTHREAD_MUTEX_INITIALIZER;

    //
    //
    //
    void bucketLockInitLocal(void)
    {
      if(bucket_lock.size() == 0 || bucket_lock.size() < numBuckets)
      {
        bucket_lock.resize(numBuckets);
        for (int i = 0; i < numBuckets; ++i)
          bucket_lock[i] = PTHREAD_MUTEX_INITIALIZER;
      }

      return;
    }

    //
    //  CantorPairLocal()
    //
    int CantorPairLocal(int a, int b)
    {
      return ((a+b)*(a+b+1)/2+b);
    }

    //
    //
    //
    int GetBucketNumberLocal(vec_int_t siteLattice)
    {
    int a = siteLattice[0];
    int b = siteLattice[1];
    int c = siteLattice[2];

    // check for -ve values
    if(a < 0) a = -a;
    if(a < 0) b = -b;
    if(a < 0) c = -c;

    // get key for a and b
    int d = CantorPairLocal(a,b);

    // get key for c and d
    int key = CantorPairLocal(c,d) % numBuckets;

    // check for -ve values if key overruns (unlikely)
    if(key < 0) key = -key;

    // return
    return key;      
    }

    //
    //  SortNodeInfoLocal()
    //
    void SortNodeInfoLocal(vec_int_t & nodeInfo, const int iNode)
    {
      // counter
      int i;
      int found = -1;

      // find i such that nodeInfo[i] = iNode
      for(int j=0; j<4; j++)
        if(nodeInfo[j] == iNode)
        { 
          found =0;
          i = j;
          break;
        }

      // swap between 0 and i
      if(found == 0)
      {
        int temp = nodeInfo[0];
        nodeInfo[0] = nodeInfo[i];   // thus nodeInfo[0] = iNode
        nodeInfo[i] = temp;
      }
      else
      {
        nodeInfo.push_back(iNode);
      }

      return;
    } // end of SortNodeInfoLocal()  

    struct ProcessAllClusterOfNodeData_t{
      int             iQuasi;
      int             numQuasi;
      int             iNode;
      int             coreShellType;
      struct site_list_t   *P_iCSites;
      struct lattice_t     *P_iQuasi_lattice;
      std::vector< std::vector< std::vector< neigh_site_data_t > > >    *P_neighborList;
      std::vector< std::vector< std::vector< cluster_site_data_t > > >  *P_clusterData;
    };

    struct ComputeNeighborFromCurrentClusterData_t{
      int               iQuasi;
      int               numQuasi;
      int               coreShellType;
      std::vector< std::vector< std::vector< neigh_site_data_t > > >    *P_neighborList;
      std::vector< std::vector< std::vector< cluster_site_data_t > > >  *P_clusterData;  
    };

    //
    //
    //
    int FindSiteInNeighBucketLocal(std::vector< neigh_site_data_t > neighBucket,
        vec_int_t siteLattice)
    {
      int index = -1;
      if(neighBucket.size() != 0)
      {
        for(int i =0; i< neighBucket.size(); i++)
        {
          if(  neighBucket[i].first.first[0] == siteLattice[0] 
            && neighBucket[i].first.first[1] == siteLattice[1]
            && neighBucket[i].first.first[2] == siteLattice[2])
          {  
            index = i;
            return index;
          }
        }
      }

      // siteLattice not found
      return index;   
    } // end of FindSiteInNeighBucketLocal()

    //
    //
    //
    int FindSiteInClusterBucketLocal(std::vector< cluster_site_data_t > clusterBucket,
        vec_int_t siteLattice)
    {
      int index = -1;
      if(clusterBucket.size() != 0)
      {
        for(int i =0; i< clusterBucket.size(); i++)
        {
          if(  clusterBucket[i].first.first[0] == siteLattice[0] 
            && clusterBucket[i].first.first[1] == siteLattice[1]
            && clusterBucket[i].first.first[2] == siteLattice[2])
          {  
            index = i;
            return index;
          }
        }
      }

      // siteLattice not found
      return index;   
    } // end of FindSiteInNeighBucketLocal()    

    //
    //  ProcessClusterSiteLocal()
    //
    //  thread_flag = 0 called from single thread subroutine
    //              = 1 called from multi thread subroutine
    //
    void ProcessClusterSiteLocal(pair_state_t iC_state,
          vec_int_t       iC_Info,
          const int       iQuasi,
          const int       numQuasi,
          std::vector< std::vector< std::vector< neigh_site_data_t > > >    *P_neighborList,
          const int     thread_flag)
    {
      int jQuasi;
      int jNSite;

      // get necessary class instance
      Quasicontinua * quasicontinua = Quasicontinua::getInstance();
      PairPotentials * pairPotentials = PairPotentials::getInstance();

      // loop over all quasi and find the neighbors aroung iC_state
      for(jQuasi =0; jQuasi < numQuasi; jQuasi++)
      {
        // get potentialNumber for interaction between iQuasi and jQuasi
        const int potentialNumber = 
          pairPotentials->doQuasicontinuumInteract(iQuasi, jQuasi);

        // cut off radius to build neighbor list
        double cutoff_cluster_radius;
        if(potentialNumber != -1)
        {
          // quasicontinuum interacts and thus, call PairPotentials to get cutoff radius
          cutoff_cluster_radius = 
            pairPotentials->getCutoffRadiusNeighborList(potentialNumber);
        }
        else
        {
          // quasicontinuum don't interact
          cutoff_cluster_radius = 0.00001;
        }

        if( potentialNumber != -1)
        {
          // iQuasi and jQuasi interact

          // get list of neigbors around iC_state at jQuasi
          vec_state_t jNeighSitesAndStates = 
            quasicontinua->getQuasi(jQuasi).getNeighborSitesAndStateAroundPoint(iC_state, cutoff_cluster_radius);

          // loop over neighbor sites and add it to neighbor list
          for(jNSite=0; jNSite < jNeighSitesAndStates.size(); jNSite++)
          {
            // get bucket
            int jBucket = 
              GetBucketNumberLocal(jNeighSitesAndStates[jNSite].first);

            // before we search for neighbor site in bucket, we must lock it
            // as other thread might be writing it to this bucket, when
            // this thread is searching for the same neighbor site.
            pthread_mutex_lock(&(bucket_lock[jBucket]));

            // check if jNSite is already assigned to d_neighborList
            int Index_jNSite = 
              FindSiteInNeighBucketLocal((*P_neighborList)[jQuasi][jBucket],
                jNeighSitesAndStates[jNSite].first);

            if(Index_jNSite == -1)
            {
              // site not found on bucket
              neigh_site_data_t data;
              data.first = jNeighSitesAndStates[jNSite];

              data.second.resize(numQuasi);

              data.second[iQuasi].push_back(iC_Info);

              // push back the data to neighbor list
              (*P_neighborList)[jQuasi][jBucket].push_back(data);
            } // site not found on bucket
            else
            {
              // site is found
              (*P_neighborList)[jQuasi][jBucket][Index_jNSite].second[iQuasi].push_back(iC_Info);
            } // site found

            pthread_mutex_unlock(&(bucket_lock[jBucket]));
          } // loop over neigh sites
        } // if quasis interact
      } // jQuasi

      return;
    } // end of  ProcessClusterSiteLocal()

    //
    //
    //
    void * ProcessAllClusterOfNodeWorker(void * arg)
    {
      const int   iQuasi = 
        ((struct ProcessAllClusterOfNodeData_t *)arg)->iQuasi;
      const int numQuasi =
        ((struct ProcessAllClusterOfNodeData_t *)arg)->numQuasi;
      const int iNode = 
        ((struct ProcessAllClusterOfNodeData_t *)arg)->iNode;
      const int coreShellType = 
        ((struct ProcessAllClusterOfNodeData_t *)arg)->coreShellType;
      struct site_list_t  *P_iCSites = 
        ((struct ProcessAllClusterOfNodeData_t *)arg)->P_iCSites;
      struct lattice_t     *P_iQuasi_lattice = 
        ((struct ProcessAllClusterOfNodeData_t *)arg)->P_iQuasi_lattice;
      std::vector< std::vector< std::vector< neigh_site_data_t > > >    *P_neighborList =
        ((struct ProcessAllClusterOfNodeData_t *)arg)->P_neighborList;
      std::vector< std::vector< std::vector< cluster_site_data_t > > >  *P_clusterData =
        ((struct ProcessAllClusterOfNodeData_t *)arg)->P_clusterData;

      int number_threads=get_number_threads();
      int my_id=get_my_tid();
      int number_thread;
      int cluster_site_start;
      int cluster_site_end;

      get_share(my_id, number_threads, P_iCSites->number_sites,
          &number_thread, &cluster_site_start, &cluster_site_end);

      int i_cluster;

      Lattice * latticeC = Lattice::getInstance();

      //  PairPotentials
      PairPotentials * pairC = PairPotentials::getInstance();

      // get shift
      vec_dbl_t iShift = Quasicontinua::getInstance()->getShift(iQuasi);

      int thread_flag = 1; 

      for(i_cluster = cluster_site_start; i_cluster <= cluster_site_end; i_cluster++)
      {
        // create cluster site data to insert it into d_clusterData
        cluster_site_data_t iC_data;

        for(int dof=0; dof< 3; dof++)
          iC_data.first.first.push_back(P_iCSites->sites[i_cluster][dof]);



        // read data from Lattice class function
        latticeC->getSiteCurrentStateAndNodeInfo(iC_data, 
          P_iQuasi_lattice, iQuasi);

        // account for global shift
        for(int dof =0; dof<3; dof++)
          iC_data.first.second[dof] = iC_data.first.second[dof] + iShift[dof];

        // get i_bucket for cluster site
        int i_bucket = 
          GetBucketNumberLocal( iC_data.first.first );

        // check if cluster site is already in cluster data bucket
        pthread_mutex_lock(&(bucket_lock[i_bucket]));
        int index_iC = FindSiteInClusterBucketLocal((*P_clusterData)[iQuasi][i_bucket], iC_data.first.first);
        pthread_mutex_unlock(&(bucket_lock[i_bucket]));

        // index_iC must be -1 because by design no node should have common cluster site
        if(index_iC != -1)
        {
          d_print("Error in cluster data and node's cluster data\n");
          D_ERROR("ProcessAllClusterOfNodeWorker()");
          exit(EXIT_FAILURE);
        }

        // sort node info
        SortNodeInfoLocal(iC_data.second, iNode);

        // lock bucket and insert the data into cluster data
        pthread_mutex_lock(&(bucket_lock[i_bucket]));
        (*P_clusterData)[iQuasi][i_bucket].push_back(iC_data);
        index_iC = (*P_clusterData)[iQuasi][i_bucket].size() -1;
        pthread_mutex_unlock(&(bucket_lock[i_bucket]));

        // create cluster info vector to be put into d_neighborList
        vec_int_t iC_Info;
        iC_Info.resize(2);

        iC_Info[0] = i_bucket;
        iC_Info[1] = index_iC;

        //
        // now compute all neighbor sites of i_cluster
        //  Note : call this only if EAM is enabled and if quasi is of shell type
        //
        if(pairC->d_EAMFlag == true && coreShellType == -1)
          ProcessClusterSiteLocal(iC_data.first, 
            iC_Info, 
            iQuasi, 
            numQuasi, 
            P_neighborList,
            thread_flag);
      } // loop over cluster sites

      return((void *)NULL);

    } // end of ProcessAllClusterOfNodeWorker()

    //
    //
    //
    void * ComputeNeighborFromCurrentClusterWorker(void * arg)
    {
      const int  iQuasi =
        ((struct ComputeNeighborFromCurrentClusterData_t *)arg)->iQuasi;
      const int numQuasi =
        ((struct ComputeNeighborFromCurrentClusterData_t *)arg)->numQuasi;
      const int coreShellType = 
        ((struct ComputeNeighborFromCurrentClusterData_t *)arg)->coreShellType;
      std::vector< std::vector< std::vector< neigh_site_data_t > > >    *P_neighborList = 
        ((struct ComputeNeighborFromCurrentClusterData_t *)arg)->P_neighborList;
      std::vector< std::vector< std::vector< cluster_site_data_t > > >  *P_clusterData = 
        ((struct ComputeNeighborFromCurrentClusterData_t *)arg)->P_clusterData;

      int number_threads=get_number_threads();
      int my_id=get_my_tid();
      int number_thread;
      int bucket_start;
      int bucket_end;

      get_share(my_id, number_threads, numBuckets,
          &number_thread, &bucket_start, &bucket_end);

      int iCluster;

      vec_dbl_t iShift = 
        Quasicontinua::getInstance()->getShift(iQuasi);

      Lattice * latticeC = Lattice::getInstance();

      PairPotentials * pairC = PairPotentials::getInstance();

      int thread_flag = 1; 
      int iBucket;

      for(iBucket = bucket_start; iBucket <= bucket_end; iBucket++)
      {
        for(iCluster = 0; iCluster < (*P_clusterData)[iQuasi][iBucket].size();
          iCluster++ )
        {
          // update the location of iCluster
          std::vector<double> iC_state;
          std::vector<int> iC_site;

          for(int dof =0; dof < 3; dof++)
            iC_site.push_back((*P_clusterData)[iQuasi][iBucket][iCluster].first.first[dof]);

          // pthread_mutex_lock(&lock_pint);
          // std::cout<<"processing cluster site "<<iBucket<<"   "<<iCluster<<std::endl;
          // pthread_mutex_unlock(&lock_pint);

          latticeC->getSiteCurrentStateNew(iC_state,
              iC_site,
              iQuasi);

          // account for global shift
          for(int dof =0; dof<3; dof++)
            iC_state[dof] += iShift[dof];

          pthread_mutex_lock(&(bucket_lock[iBucket]));
          (*P_clusterData)[iQuasi][iBucket][iCluster].first.second[0] = iC_state[0];
          (*P_clusterData)[iQuasi][iBucket][iCluster].first.second[1] = iC_state[1];
          (*P_clusterData)[iQuasi][iBucket][iCluster].first.second[2] = iC_state[2];
          (*P_clusterData)[iQuasi][iBucket][iCluster].first.second[3] = iC_state[3];
          pthread_mutex_unlock(&(bucket_lock[iBucket]));

          vec_int_t iC_Info;
          iC_Info.resize(2);

          iC_Info[0]= iBucket;
          iC_Info[1] = iCluster;

          //
          // use update cluster site to find it's neighbors
          //  Call this if EAM is enabled and if quasi is of "shell" type
          //
          if(pairC->d_EAMFlag == true && coreShellType == -1)
            ProcessClusterSiteLocal((*P_clusterData)[iQuasi][iBucket][iCluster].first, 
              iC_Info,
              iQuasi,
              numQuasi,
              P_neighborList,
              thread_flag);       
        } // loop over cluster site
      } // loop over bucket

      return((void *) NULL);

    } // end of ComputeNeighborFromCurrentClusterWorker()    
  } // end of namespace for neighbor list computing

  //
  //  namespace for updating locations
  //
  namespace
  {
    struct UpdateClusterLocationsData_t{
      int       iQuasi;
      vec_dbl_t       iShift;
      std::vector< std::vector< cluster_site_data_t > > *P_clusterBucket;
    };

    struct UpdateNeighborLocationsData_t{
      int       iQuasi;
      vec_dbl_t       iShift;
      std::vector< std::vector< neigh_site_data_t > > *P_neighborBucket;
    };

    //
    //
    //
    void UpdateClusterLocationsLocal(const int iQuasi,
          vec_dbl_t       iShift,
          const int       bucket_start,
          const int       bucket_end,
          std::vector< std::vector< cluster_site_data_t > > *P_clusterBucket)
    {
      int iBucket;
      int iCluster;

      Lattice * latticeC = Lattice::getInstance();

      for(iBucket = bucket_start; iBucket <= bucket_end; iBucket++)
      {
        for(iCluster =0; iCluster < (*P_clusterBucket)[iBucket].size(); 
          iCluster++)
        {
          std::vector<double> iC_state;
          std::vector<int> iC_site;
          for(int dof =0; dof < 3; dof++)
            iC_site.push_back((*P_clusterBucket)[iBucket][iCluster].first.first[dof]);

          latticeC->getSiteCurrentStateNew(iC_state,
            iC_site,
            iQuasi);

          // account for global shift
          for(int dof =0; dof<3; dof++)
            iC_state[dof] += iShift[dof];

          // lock the bucket and update state of cluster site
          pthread_mutex_lock(&(bucket_lock[iBucket]));
          (*P_clusterBucket)[iBucket][iCluster].first.second[0] = iC_state[0];
          (*P_clusterBucket)[iBucket][iCluster].first.second[1] = iC_state[1];
          (*P_clusterBucket)[iBucket][iCluster].first.second[2] = iC_state[2];
          (*P_clusterBucket)[iBucket][iCluster].first.second[3] = iC_state[3];
          pthread_mutex_unlock(&(bucket_lock[iBucket]));
        } // loop over cluster
      } // loop over buckets

      return;
    } // end of UpdateClusterLocationsLocal()

    //
    //
    //
    void UpdateNeighborLocationsLocal(const int iQuasi,
          vec_dbl_t       iShift,
          const int       bucket_start,
          const int       bucket_end,
          std::vector< std::vector< neigh_site_data_t > > *P_neighborBucket)
    {
      int iBucket;
      int iNeighbor;

      Lattice * latticeC = Lattice::getInstance();

      for(iBucket = bucket_start; iBucket <= bucket_end; iBucket++)
      {
        for(iNeighbor =0; iNeighbor < (*P_neighborBucket)[iBucket].size(); 
          iNeighbor++)
        {
          std::vector<double> iN_state;
          std::vector<int> iN_site;
          for(int dof=0; dof < 3; dof++)
            iN_site.push_back((*P_neighborBucket)[iBucket][iNeighbor].first.first[dof]);

          latticeC->getSiteCurrentStateNew(iN_state,
            iN_site,
            iQuasi);

          // account for global shift
          for(int dof =0; dof<3; dof++)
            iN_state[dof] += iShift[dof];

          // lock the bucket and update state of cluster site
          pthread_mutex_lock(&(bucket_lock[iBucket]));
          (*P_neighborBucket)[iBucket][iNeighbor].first.second[0] = iN_state[0];
          (*P_neighborBucket)[iBucket][iNeighbor].first.second[1] = iN_state[1];
          (*P_neighborBucket)[iBucket][iNeighbor].first.second[2] = iN_state[2];
          (*P_neighborBucket)[iBucket][iNeighbor].first.second[3] = iN_state[3];
          pthread_mutex_unlock(&(bucket_lock[iBucket]));

          // for debug
          iN_state.clear();
          iN_site.clear();
        } // loop over cluster
      } // loop over buckets

      return;
    } // end of UpdateNeighborLocationsLocal()

    //
    //
    //
    void * UpdateClusterLocationsWorker(void * arg)
    {
      const int iQuasi =
        ((struct UpdateClusterLocationsData_t *)arg)->iQuasi;
      vec_dbl_t iShift = 
        ((struct UpdateClusterLocationsData_t *)arg)->iShift;
      std::vector< std::vector< cluster_site_data_t > > *P_clusterBucket =
        ((struct UpdateClusterLocationsData_t *)arg)->P_clusterBucket;

      int number_threads=get_number_threads();
      int my_id=get_my_tid();
      int number_thread;
      int bucket_start;
      int bucket_end;

      get_share(my_id, number_threads, numBuckets,
          &number_thread, &bucket_start, &bucket_end);

      // call function to process range of buckets assigned to current thread
      UpdateClusterLocationsLocal(iQuasi,
        iShift,
        bucket_start,
        bucket_end,
        P_clusterBucket);

      return((void *) NULL);
    } // end of UpdateClusterLocationsWorker()

    //
    //
    //
    void * UpdateNeighborLocationsWorker(void * arg)
    {
      const int iQuasi =
        ((struct UpdateNeighborLocationsData_t *)arg)->iQuasi;
      vec_dbl_t iShift = 
        ((struct UpdateNeighborLocationsData_t *)arg)->iShift;
      std::vector< std::vector< neigh_site_data_t > > *P_neighborBucket =
        ((struct UpdateNeighborLocationsData_t *)arg)->P_neighborBucket;

      int number_threads=get_number_threads();
      int my_id=get_my_tid();
      int number_thread;
      int bucket_start;
      int bucket_end;

      get_share(my_id, number_threads, numBuckets,
          &number_thread, &bucket_start, &bucket_end);

      // call function to process range of buckets assigned to current thread
      UpdateNeighborLocationsLocal(iQuasi,
        iShift,
        bucket_start,
        bucket_end,
        P_neighborBucket);

      return((void *) NULL);
    } // end of UpdateClusterLocationsWorker()   

  } // end of namespace for updating locations

  //
  //
  //

  CrossNeighborList* CrossNeighborList::_instance = NULL;

  //
  // constructor
  //

  CrossNeighborList::CrossNeighborList()
  {

    //
    //
    //
    return;

  }

  //
  // destructor
  //

  CrossNeighborList::~CrossNeighborList()
  {

    //
    //
    //
    return;

  }

  //
  // getInstance method
  //

  CrossNeighborList*
  CrossNeighborList::getInstance()
  {

    //
    // if not created, create
    //
    if(_instance == NULL){
      _instance = new CrossNeighborList();
    }

    //
    // return instance
    //
    return _instance;

  }

  //
  // destroy instance method
  //

  void
  CrossNeighborList::destroyInstance()
  {

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
  //  computeCrossNeighborList()
  //
  void CrossNeighborList::computeCrossNeighborList(const int rebuild_neighbor_flag,
          const int           rebuild_cluster_flag)
  {
    // check if it is multi threading or single threading
    int number_threads = get_max_number_threads();
    enum mt_version_t mt_version = MULTI_THREADED;

    int numQuasi = Quasicontinua::getInstance()->size();

    if(number_threads == 1)
    {
      mt_version = SINGLE_THREADED;
    }

    // we use bucket_lock in local function, irrespective of single/multi 
    // thread
    // if(mt_version == MULTI_THREADED)
    bucketLockInitLocal();

    // ComputeNeighborList(mt_version, 1, rebuild_neighbor_flag);

    // return;

    if(rebuild_cluster_flag == 1 )
    {
      // prepare new cluster data ---> then prepare new neighbor data
      d_print("Build Neighbor List...");
      ComputeNeighborList(mt_version, rebuild_cluster_flag, rebuild_neighbor_flag);

      return;
    }
    else
    { 
      // check for rebuild_neighbor_flag
      if(rebuild_neighbor_flag == 1)
      {
        // prepare new neighbor list using old cluster data
        d_print("Neighbor List...");
        ComputeNeighborList(mt_version, rebuild_cluster_flag, rebuild_neighbor_flag);

        return;
      }
      else
      {
        // update location of all sites
        d_print("Update Neighbor List...");
        UpdateLocation(mt_version);

        return;
      }
    }

    // if reached here, then there is problem in the code
    d_print("check flags passed to computeCrossNeighborList()\n");
    D_ERROR("computeCrossNeighborList()");
    exit(EXIT_FAILURE);
  } // end of computeCrossNeighborList()

  //
  //  printCrossNeighborList()
  //
  void CrossNeighborList::printCrossNeighborList(const int thread_flag)
  {
    FILE * file_1;
    FILE * file_2;
    // read file based on thread flag
    if(thread_flag == 0)
    {
      file_1 = fopen("cluster_data_single_thread_print.quasi", "w");
      file_2 = fopen("neighbor_data_single_thread_print.quasi", "w");
    }
    else
    {
      file_1 = fopen("cluster_data_multi_thread_print.quasi", "w");
      file_2 = fopen("neighbor_data_multi_thread_print.quasi", "w");
    }

    // print cluster and neighbor data for iQuasi =0
    for(int iBucket =0; iBucket < numBuckets; iBucket++)
    {
      if(d_clusterData[0][iBucket].size() != 0)
      {
        for(int i = 0; i< d_clusterData[0][iBucket].size(); i++)
        {
          //fprintf(file_1, "( iBucket, iC) = (%i, %i) \n", iBucket, i);
          fprintf(file_1,"(%i, %i, %i) , (%f, %f, %f) \n", 
            d_clusterData[0][iBucket][i].first.first[0],
            d_clusterData[0][iBucket][i].first.first[1],
            d_clusterData[0][iBucket][i].first.first[2],
            d_clusterData[0][iBucket][i].first.second[0],
            d_clusterData[0][iBucket][i].first.second[1],
            d_clusterData[0][iBucket][i].first.second[2]);
          // printf("%i    %i    %i    %f    %f    %f\n",
          //   d_clusterData[0][iBucket][i].first.first[0],
          //   d_clusterData[0][iBucket][i].first.first[1],
          //   d_clusterData[0][iBucket][i].first.first[2],
          //   d_clusterData[0][iBucket][i].first.second[0],
          //   d_clusterData[0][iBucket][i].first.second[1],
          //   d_clusterData[0][iBucket][i].first.second[2]);      
        }
      }

      if(d_neighborList[0][iBucket].size() != 0)
      {
        for(int i = 0; i< d_neighborList[0][iBucket].size(); i++)
        {
          //fprintf(file_2, "( iBucket, iN) = (%i, %i) \n", iBucket, i);
          fprintf(file_2,"(%i, %i, %i) , (%f, %f, %f) \n", 
            d_neighborList[0][iBucket][i].first.first[0],
            d_neighborList[0][iBucket][i].first.first[1],
            d_neighborList[0][iBucket][i].first.first[2],
            d_neighborList[0][iBucket][i].first.second[0],
            d_neighborList[0][iBucket][i].first.second[1],
            d_neighborList[0][iBucket][i].first.second[2]);
        }
      }      
    }

    fclose(file_1);
    fclose(file_2);

    return;
  }

  //
  //  ComputeNeighborList()
  //
  void CrossNeighborList::ComputeNeighborList(enum mt_version_t mt_version,
        const int         rebuild_cluster_flag,
        const int         rebuild_neighbor_flag)
  {
    // resize the neighbor list data
    Quasicontinua * quasicontinua = Quasicontinua::getInstance();

    int numQuasi = quasicontinua->size();
    int iQuasi;

    //
    // time related data
    //
    int t_size = 2 * numQuasi;

    timeval buildTime_t[t_size];
    timeval neighborList_t[t_size];

    double elapsedTime;

    //
    //  get pairPotentials instance
    //
    PairPotentials * pairC = PairPotentials::getInstance();

    // building cluster list and neighbor list every time, instead of subjecting it to flags
    if(pairC->d_EAMFlag == true)
    {
      d_neighborList.clear();
      d_neighborList.resize(numQuasi);
      for(int iQ =0; iQ < numQuasi; iQ++)
        d_neighborList[iQ].resize(numBuckets);
    }

    if(rebuild_cluster_flag == 1)
    {
      // we need to recreate cluster data in this case. Thus resetting the data
      d_clusterData.clear();
      d_clusterData.resize(numQuasi);
      for(int iQuasi = 0; iQuasi < numQuasi; iQuasi ++ )
        d_clusterData[iQuasi].resize(numBuckets);              
    }

    // loop over iQuasi and all iCluster site
    for(iQuasi =0; iQuasi < numQuasi; iQuasi++)
    {
      // depending on cluster flag call function
      if(rebuild_cluster_flag == 1)
      {
        int t_1 = iQuasi*2 + 0;
        int t_2 = iQuasi*2 + 1;
        
        gettimeofday(&buildTime_t[t_1],NULL);      

        d_print("cluster list...");
        BuildClusterAndComputeNeighborList(mt_version, iQuasi);

        gettimeofday(&buildTime_t[t_2],NULL);      
      }
      else
      {
        int t_1 = iQuasi*2 + 0;
        int t_2 = iQuasi*2 + 1;

        gettimeofday(&neighborList_t[t_1],NULL);      

        d_print("neighbor list...");
        ComputeNeighborListUsingCurrentClusterData(mt_version, iQuasi);

        gettimeofday(&neighborList_t[t_2],NULL);      
      }
    }

    // compute time if required
    static int clusterNumFlag =0;
    if(RunData::getInstance()->d_timeOut == true)
    {
      for(iQuasi =0; iQuasi < numQuasi; iQuasi++)
      {
        // time required 
        if(rebuild_cluster_flag == 1)
        {
          int t_1 = iQuasi*2 + 0;
          int t_2 = iQuasi*2 + 1;        

          elapsedTime = (buildTime_t[t_2].tv_sec - buildTime_t[t_1].tv_sec);
          elapsedTime += (buildTime_t[t_2].tv_usec - buildTime_t[t_1].tv_usec) / 1000000.0;
          RunData::getInstance()->addComputeClusterAndNeighborListTime(iQuasi,
            elapsedTime);
        }
        else
        {
          int t_1 = iQuasi*2 + 0;
          int t_2 = iQuasi*2 + 1;        

          elapsedTime = (neighborList_t[t_2].tv_sec - neighborList_t[t_1].tv_sec);
          elapsedTime += (neighborList_t[t_2].tv_usec - neighborList_t[t_1].tv_usec) / 1000000.0;

          RunData::getInstance()->addBuildNeighborListTime(iQuasi,
            elapsedTime);
        }

        // compute number of cluster sites
        if(clusterNumFlag < numQuasi)
        {
          int n =0;
          for(int iB=0; iB < d_clusterData[iQuasi].size(); iB++)
            n += d_clusterData[iQuasi][iB].size();

          RunData::getInstance()->addNumberOfClusterSites(iQuasi, n);

          // set flag to 1
          clusterNumFlag += 1;
        }
      }
    }
    return;
  } // end of ComputeNeighborList()

  //
  //  BuildClusterAndComputeNeighborList()
  //
  void CrossNeighborList::BuildClusterAndComputeNeighborList(enum mt_version_t mt_version,
          int         iQuasi)
  {
    // get Quasicontinua instance
    Quasicontinua * quasicontinua = Quasicontinua::getInstance();

    // PairPotentials
    PairPotentials * pairC = PairPotentials::getInstance();

    // get Lattice instance
    Lattice * latticeC = Lattice::getInstance();

    // get shift vector
    vec_dbl_t iShift = 
      quasicontinua->getShift(iQuasi);    

    // get number of quasi
    int numQuasi = quasicontinua->size();       

    switch(mt_version)
    {
      case SINGLE_THREADED:
      {
        // std::cout<<"CrossNeighborlist : Single Threaded"<<std::endl;

        int iNode;
        int i_cluster;

        struct all_node_list_t all_node_list = 
          quasicontinua->getQuasi(iQuasi).getNodeList();

        struct lattice_t iQuasi_lattice = 
          quasicontinua->getQuasi(iQuasi).getLattice();   

        for(iNode = 0; iNode < all_node_list.node_list.number_nodes; iNode++)
        {
          // get handle of iNode
          struct node_t * P_node = all_node_list.node_list.nodes[iNode];
          struct site_list_t iCSites = P_node->site_cluster;

          // loop over cluster sites
          for(i_cluster=0; i_cluster < iCSites.number_sites; i_cluster++)
          {
            //std::cout<<"cluster site = "<<i_cluster<<std::endl;
            // create cluster site data to insert it into d_clusterData
            cluster_site_data_t iC_data;

            for(int dof=0; dof< 3; dof++)
              iC_data.first.first.push_back(iCSites.sites[i_cluster][dof]);

            // read data from Lattice class function
            latticeC->getSiteCurrentStateAndNodeInfo(iC_data, 
              &iQuasi_lattice, iQuasi);

            // account for global shift
            for(int dof =0; dof<3; dof++)
              iC_data.first.second[dof] = iC_data.first.second[dof] + iShift[dof];

            // get i_bucket for cluster site
            int i_bucket = 
              GetBucketNumberLocal( iC_data.first.first );

            // check if cluster site is already in cluster data bucket
            int index_iC = FindSiteInClusterBucketLocal(d_clusterData[iQuasi][i_bucket], iC_data.first.first);

            // index_iC must be -1 because by design no node should have common cluster site
            if(index_iC != -1)
            {
              d_print("Error in cluster data and node's cluster data\n");
              D_ERROR("BuildClusterAndComputeNeighborList()");
              exit(EXIT_FAILURE);
            }

            // sort node info
            SortNodeInfoLocal(iC_data.second, iNode);

            // insert the data into cluster data
            d_clusterData[iQuasi][i_bucket].push_back(iC_data);
            
            // location of current cluster site in bucket
            index_iC = d_clusterData[iQuasi][i_bucket].size() -1;

            // create cluster info vector to be put into d_neighborList
            vec_int_t iC_Info;
            iC_Info.resize(2);

            iC_Info[0] = i_bucket;
            iC_Info[1] = index_iC;

            //
            // now compute all neighbor sites of i_cluster
            //  Note : do this only if EAM is enabled and if quasi is "shell" type
            //
            int thread_flag = 0;
            int coreShellType = quasicontinua->getCoreShell(iQuasi);

            if(pairC->d_EAMFlag == true && coreShellType == -1)
              ProcessClusterSiteLocal(iC_data.first, 
                iC_Info, 
                iQuasi, 
                numQuasi, 
                &d_neighborList,
                thread_flag);
          } // loop over cluster sites
        } // iNode

        return;
      } // single threaded

      break;

      case MULTI_THREADED:
      {
        // std::cout<<"CrossNeighborlist : Multi Threaded"<<std::endl;
        int iNode;
        int i_cluster;

        struct all_node_list_t all_node_list = 
          quasicontinua->getQuasi(iQuasi).getNodeList();

        struct lattice_t iQuasi_lattice = 
          quasicontinua->getQuasi(iQuasi).getLattice();

        for(iNode = 0; iNode < all_node_list.node_list.number_nodes; iNode++)
        {
          // get handle of iNode
          struct node_t * P_node = all_node_list.node_list.nodes[iNode];
          struct site_list_t iCSites = P_node->site_cluster;

          // call multi thread only if number of cluster site is more than
          // number of max threads.
          int number_threads = get_max_number_threads();
          if(iCSites.number_sites >= number_threads)
          {
            // call multi thread
            struct ProcessAllClusterOfNodeData_t data;

            int coreShellType = quasicontinua->getCoreShell(iQuasi);

            data.iQuasi    = iQuasi;
            data.numQuasi  = numQuasi;
            data.iNode     = iNode;
            data.coreShellType = coreShellType;
            data.P_iCSites = &iCSites;
            data.P_iQuasi_lattice = &iQuasi_lattice;
            data.P_neighborList = &d_neighborList;
            data.P_clusterData = &d_clusterData;

            thread_monitor(ProcessAllClusterOfNodeWorker, (void *)&data,
              get_max_number_threads());
          }
          else
          {
            // loop over cluster sites
            for(i_cluster=0; i_cluster < iCSites.number_sites; i_cluster++)
            {
              // create cluster site data to insert it into d_clusterData
              cluster_site_data_t iC_data;

              for(int dof=0; dof< 3; dof++)
                iC_data.first.first.push_back(iCSites.sites[i_cluster][dof]);

              // read data from Lattice class function
              latticeC->getSiteCurrentStateAndNodeInfo(iC_data, 
                &iQuasi_lattice, iQuasi);

              // account for global shift
              for(int dof =0; dof<3; dof++)
                iC_data.first.second[dof] = iC_data.first.second[dof] + iShift[dof];

              // get i_bucket for cluster site
              int i_bucket = 
                GetBucketNumberLocal( iC_data.first.first );

              // check if cluster site is already in cluster data bucket
              int index_iC = FindSiteInClusterBucketLocal(d_clusterData[iQuasi][i_bucket], iC_data.first.first);

              // index_iC must be -1 because by design no node should have common cluster site
              if(index_iC != -1)
              {
                d_print("Error in cluster data and node's cluster data\n");
                D_ERROR("BuildClusterAndComputeNeighborList()");
                exit(EXIT_FAILURE);
              }

              // sort node info
              SortNodeInfoLocal(iC_data.second, iNode);

              // insert the data into cluster data
              d_clusterData[iQuasi][i_bucket].push_back(iC_data);
              
              index_iC = d_clusterData[iQuasi][i_bucket].size() -1;

              // create cluster info vector to be put into d_neighborList
              vec_int_t iC_Info;
              iC_Info.resize(2);

              iC_Info[0] = i_bucket;
              iC_Info[1] = index_iC;

              //
              // now compute all neighbor sites of i_cluster
              //  Note : call this only if EAM is enabled and if quasi is of shell type
              //
              int thread_flag = 0;  // 0: single, 1: multi
              int coreShellType = quasicontinua->getCoreShell(iQuasi);              

              if(pairC->d_EAMFlag == true && coreShellType == -1)
                ProcessClusterSiteLocal(iC_data.first, 
                  iC_Info, 
                  iQuasi, 
                  numQuasi, 
                  &d_neighborList,
                  thread_flag);
            } // loop over cluster sites   

          } // end of if-else
        } // iNode

        return;
      } // multi threaded

      break;
    } // switch

    // std::cout<<"CrossNeighborList : Done"<<std::endl;

    return;
  } // end of BuildClusterAndComputeNeighborList()

  //
  //  ComputeNeighborListUsingCurrentClusterData()
  //
  void CrossNeighborList::ComputeNeighborListUsingCurrentClusterData(
            enum mt_version_t mt_version,
            int         iQuasi)
  {
    // d_print("Neighbor List ...");

    // get Quasicontinua instance
    Quasicontinua * quasicontinua = Quasicontinua::getInstance();
    int numQuasi = quasicontinua->size();
    Lattice * latticeC = Lattice::getInstance();

    PairPotentials * pairC = PairPotentials::getInstance();

    switch(mt_version)
    {
      case SINGLE_THREADED:
      {
        // std::cout<<"CrossNeighborList : Single Threaded"<<std::endl;
        int iBucket;
        int iCluster;
        int thread_flag = 0;  // single thread

        vec_dbl_t iShift = quasicontinua->getShift(iQuasi);

        // loop over all cluster sites in buckets
        for(iBucket =0; iBucket< numBuckets; iBucket++)
        {
          for(iCluster=0; iCluster < d_clusterData[iQuasi][iBucket].size();
            iCluster++)
          {
            // update the location of i_cluster
            std::vector<int> iC_site;
            std::vector<double> iC_state;

            for(int dof =0; dof < 3; dof++)
              iC_site.push_back(d_clusterData[iQuasi][iBucket][iCluster].first.first[dof]);

            latticeC->getSiteCurrentStateNew(iC_state,
              iC_site,
              iQuasi);

            // account for global shift
            for(int dof =0; dof<3; dof++)
              iC_state[dof] += iShift[dof];

            d_clusterData[iQuasi][iBucket][iCluster].first.second[0] = iC_state[0];
            d_clusterData[iQuasi][iBucket][iCluster].first.second[1] = iC_state[1];
            d_clusterData[iQuasi][iBucket][iCluster].first.second[2] = iC_state[2];
            d_clusterData[iQuasi][iBucket][iCluster].first.second[3] = iC_state[3];

            vec_int_t iC_Info;
            iC_Info.resize(2);

            iC_Info[0] = iBucket;
            iC_Info[1] = iCluster;

            // 
            // use update cluster site to find it's neighbors
            // Note : do this only if EAM is enabled and if quasi is of "shell" type
            //
            int coreShellType = quasicontinua->getCoreShell(iQuasi);
            if(pairC->d_EAMFlag == true && coreShellType == -1)
              ProcessClusterSiteLocal(d_clusterData[iQuasi][iBucket][iCluster].first, 
                iC_Info,
                iQuasi,
                numQuasi,
                &d_neighborList,
                thread_flag);
          }// loop over cluster sites
        } // loop over buckets
      } // single threaded

      break;

      case MULTI_THREADED:
      {
        // std::cout<<"CrossNeighborList : Multi Threaded"<<std::endl;
        struct ComputeNeighborFromCurrentClusterData_t data;
        int coreShellType = quasicontinua->getCoreShell(iQuasi);

        data.iQuasi = iQuasi;
        data.numQuasi = numQuasi;
        data.coreShellType = coreShellType;
        data.P_neighborList = &d_neighborList;
        data.P_clusterData = &d_clusterData;

        thread_monitor(ComputeNeighborFromCurrentClusterWorker, (void*)&data,
          get_max_number_threads());
      } // multi threaded

      break;
    } // switch

    // std::cout<<"CrossNeighborList : Done"<<std::endl;

    return;
  } // end of ComputeNeighborListUsingCurrentClusterData()

  //
  //  UpdateLocation()
  //
  void CrossNeighborList::UpdateLocation(enum mt_version_t mt_version)
  {
    // d_print("Update Neighbor List ...");

    int numQuasi = Quasicontinua::getInstance()->size();

    PairPotentials * pairC = PairPotentials::getInstance();

    // time data
    timeval t1, t2;
    double elapsedTime;  
    
    switch(mt_version)
    {
      case SINGLE_THREADED:
      {
        // std::cout<<"CrossNeighborList : Sinlge Threaded"<<std::endl;

        int iQuasi;
        // first update location of cluster list
        for(int iQuasi = 0; iQuasi< numQuasi; iQuasi++)
        {
          // get time
          gettimeofday(&t1,NULL); 

          // get global shift
          vec_dbl_t iShift = Quasicontinua::getInstance()->getShift(iQuasi);

          // first update location of clusters
          int bucket_start = 0;
          int bucket_end = d_clusterData[iQuasi].size();

          UpdateClusterLocationsLocal(iQuasi, 
            iShift,
            bucket_start, 
            bucket_end - 1,
            &(d_clusterData[iQuasi]));

          // Update neighbor locations only if EAM is enabled
          if(pairC->d_EAMFlag == true)
          {
            bucket_end = d_neighborList[iQuasi].size();

            UpdateNeighborLocationsLocal(iQuasi,
              iShift,
              bucket_start,
              bucket_end - 1,
              &(d_neighborList[iQuasi]));
          }

          // get time
          gettimeofday(&t2,NULL); 
          elapsedTime = (t2.tv_sec - t1.tv_sec);
          elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0; 

          if(RunData::getInstance()->d_timeOut == true)
            RunData::getInstance()->addUpdateClusterAndNeighborListTime(iQuasi,
              elapsedTime);          
        } // loop over iQuasi
      } // single threaded

      break;

      case MULTI_THREADED:
      {
        // std::cout<<"CrossNeighborList : Multi Threaded"<<std::endl;
        int iQuasi;

        // loop over quasis
        for(iQuasi =0; iQuasi < numQuasi; iQuasi++)
        {
          // get time
          gettimeofday(&t1,NULL);           
          
          // get global shift
          vec_dbl_t iShift = Quasicontinua::getInstance()->getShift(iQuasi);

          // first update cluster data by calling thread function
          struct UpdateClusterLocationsData_t data_1;
          data_1.iQuasi = iQuasi;
          data_1.iShift = iShift;
          data_1.P_clusterBucket = &(d_clusterData[iQuasi]);

          thread_monitor(UpdateClusterLocationsWorker, (void *)&data_1,
            get_max_number_threads());

          // call htread function to update neighbors
          if(pairC->d_EAMFlag == true)
          {
            struct UpdateNeighborLocationsData_t data_2;
            data_2.iQuasi = iQuasi;
            data_2.iShift = iShift;
            data_2.P_neighborBucket = &(d_neighborList[iQuasi]);

            thread_monitor(UpdateNeighborLocationsWorker, (void *)&data_2,
              get_max_number_threads());          
          }

          // get time
          gettimeofday(&t2,NULL); 
          elapsedTime = (t2.tv_sec - t1.tv_sec);
          elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0; 

          if(RunData::getInstance()->d_timeOut == true)
            RunData::getInstance()->addUpdateClusterAndNeighborListTime(iQuasi,
              elapsedTime); 
        } // loop over iQuasi
      } // multi threaded

      break;


    } // switch

    return;
  } // end of UpdateLocation()

  //
  //  getQuasiNeighborData()
  //
  std::vector< std::vector< neigh_site_data_t > > & 
  CrossNeighborList::getQuasiNeighborData(const int iQuasi)  
  {
    if(iQuasi == -1 || iQuasi > d_neighborList.size())
    {
      ERROR("getCrossNeighborListNeighborData()");
      exit(EXIT_FAILURE);
    }

    return d_neighborList[iQuasi];
  }

  //
  //  getQuasiClusterData
  //
  std::vector< std::vector< cluster_site_data_t > > &
  CrossNeighborList::getQuasiClusterData(const int iQuasi)
  {
    if(iQuasi == -1 || iQuasi > d_clusterData.size())
    {
      ERROR("getCrossNeighborListClusterData()");
      exit(EXIT_FAILURE);
    }

    return d_clusterData[iQuasi];
  }    

  //
  //  getClusterData
  //
  std::vector< std::vector< std::vector< cluster_site_data_t > > > &
  CrossNeighborList::getClusterData(void)
  {
    return d_clusterData;
  }      

  //
  //  getNodeInfoOfClusterSite()
  //
  std::vector<int> CrossNeighborList::getNodeInfoOfClusterSite(int iQuasi, 
    std::vector<int> clusterSite)  
  {
    if(iQuasi > d_clusterData.size() || iQuasi < 0)
    {
      d_print("wrong iQuasi passed to CrossNeighborlist::getNodeInfoOfClusterSite\n");
      D_ERROR("getNodeInfoOfClusterSite()");
      exit(EXIT_FAILURE);
    }

    // 
    std::vector<int> returnVec;

    // get bucket number
    int bucket = GetBucketNumberLocal(clusterSite);

    // find site in bucket
    int index = -1;

    for(int i =0; i < d_clusterData[iQuasi][bucket].size(); i++)
    {
      if(   d_clusterData[iQuasi][bucket][i].first.first[0] == clusterSite[0]
        &&  d_clusterData[iQuasi][bucket][i].first.first[1] == clusterSite[1]
        &&  d_clusterData[iQuasi][bucket][i].first.first[2] == clusterSite[2])
        index = i;
    }

    if(index == -1)
    {
      d_print("site (%d, %d, %d) of quasi %d not found in cluster data\n",
        clusterSite[0],
        clusterSite[1],
        clusterSite[2],
        iQuasi);

      return returnVec;
    }

    returnVec = d_clusterData[iQuasi][bucket][index].second;

    return returnVec;
  }

  //
  //  getNumberBuckets()
  //
  int CrossNeighborList::getNumberBuckets()
  {
    return numBuckets;
  }

  //
  //  getSizeOfClusterData
  //
  int CrossNeighborList::getSizeOfClusterData(void)
  {
    return d_clusterData.size();
  }       

} // end of quasicontinuum namespace

