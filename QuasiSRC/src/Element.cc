//
// Element.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#if !defined(_REENTRANT)
#  define _REENTRANT
#endif
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#else
#error No standard libraries available.
#endif /* STDC_HEADERS */

#ifdef HAVE_SCHED_H
#include <sched.h>
#else
#error sched.h not available.
#endif /* HAVE_SCHED_H */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h header not available
#endif /* HAVE_SYS_TYPES_H */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error No math.h found.
#endif /* HAVE_MATH_H */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not available.
#endif /* HAVE_UNISTD_H */

#ifdef HAVE_SEARCH_H
#include <search.h>
#else
#error search.h not found
#endif

#include "Element.h"
#include "DataTypes.h"
#include "Error.h"
#include "MiscFunctions.h"
#include "Shape.h"

#include "threads.h"
#include "monitor.h"
#include "numa.h"
#include "site_element_map_cache.h"
#include "site_element_map.h"

// global variables
#define NODE_ALLOC_CHUNK 16   // used in addNodeLocal()
#define ELEMENT_ALLOC_CHUNK 16
#define N_ALLOC_PAGES 2   // used in make_active_element()
// static variables for makeActiveElements()
static int      elem_bucket;
static int      init_active_function = 0;
static long int page_size;

// ****hash table related definitions**** //
#define WIDTH_0 5
#define WIDTH_1 5
#define WIDTH_2 5

#define MODULO_BITS 21 
#define MODULO_MASK ~(~0 << MODULO_BITS)
/* #define NUMBER_BUCKETS ((unsigned) 1 << MODULO_BITS) */
#define NUMBER_BUCKETS 64997
/* #define NUMBER_BUCKETS ((unsigned) 1 << (WIDTH_0+WIDTH_1+WIDTH_2)) */

//
//
//

namespace quasicontinuum {

  // splitting unnamed namespace for different functions.

  //
  // namespace for binElements() and cleanElements()
  //
  namespace
  {
    // static variables for makeActiveElements()
    static int      elem_bucket;
    static int      init_active_function = 0;

    // ***binElements related stuff*** //
    struct binElementsData_t{
      struct element_list_t *P_element_list;
      struct node_list_t    *P_node_list;
    };

    //
    // addElementLocal()
    //
    void addElementLocal(struct element_t   *P_element,
                struct node_t               *P_node)
    {
      const size_t P_elem_len = sizeof(struct element_t *);
      int          i_element = P_node->element_list.number_elements + 1;

      // check for allocated space
      if( P_node->element_list.alloc_space == 0 )
      {
        if( (P_node->element_list.elements = 
          (struct element_t **) malloc(ELEMENT_ALLOC_CHUNK*P_elem_len))
          == NULL )
        {
          PERROR_MALLOC();
          ERROR("xmalloc()");
          exit(EXIT_FAILURE);
        }

        P_node->element_list.alloc_space = ELEMENT_ALLOC_CHUNK;
      }
      else if( i_element > P_node->element_list.alloc_space )
      {
        P_node->element_list.alloc_space += ELEMENT_ALLOC_CHUNK;

        // realloc the space
        P_node->element_list.elements = 
          (struct element_t **) realloc(P_node->element_list.elements,
              P_node->element_list.alloc_space*P_elem_len);

        if( P_node->element_list.elements == NULL )
        {
          PERROR_MALLOC();
          ERROR("realloc()");
          exit(EXIT_FAILURE);
        }
      }
      
      #if defined(__QC_SGI)
        change_memory_placement(P_node->element_list.elements,
          P_node->element_list.alloc_space*P_elem_len, 0);
      #endif /* sgi */

      // insert element
      P_node->element_list.number_elements++;
      i_element--;
      P_node->element_list.elements[i_element]=P_element;

      return;
    } // end of addElementLocal()

    //
    // addNodeLocal()
    //
    void addNodeLocal(struct node_t     *P_node,
            struct node_t               *P_node_add)
    {
      const size_t P_node_len=sizeof(struct node_t *);
      int j_node;

      // scan the list of nodes associated to P_node, if found, we
      // don't need to add the node to P_node's neighbor node list.
      for(j_node=0; j_node < P_node->node_list.number_nodes; j_node++)
        if( P_node->node_list.nodes[j_node] == P_node_add ) 
          return;

      // P_node_add has not been found on the list; add it
      // check if there is enough space first
      j_node = P_node->node_list.number_nodes+1;

      if( P_node->node_list.alloc_space == 0 )
      {
        if( (P_node->node_list.nodes = 
            (struct node_t **) malloc(NODE_ALLOC_CHUNK*P_node_len)) 
            == NULL )
        {
          PERROR_MALLOC();
          ERROR("malloc()");
          exit(EXIT_FAILURE);
        }

        P_node->node_list.alloc_space = NODE_ALLOC_CHUNK;
      }
      else if( j_node > P_node->node_list.alloc_space )
      {
        P_node->node_list.alloc_space += NODE_ALLOC_CHUNK;

        P_node->node_list.nodes = 
          (struct node_t **) realloc(P_node->node_list.nodes, 
                P_node->node_list.alloc_space*P_node_len);

        if( P_node->node_list.nodes == NULL )
        {
          PERROR_MALLOC();
          ERROR("realloc()");
          exit(EXIT_FAILURE);
        }
      }

      #if defined(__QC_SGI)
        change_memory_placement(P_node->node_list.nodes, 
              P_node->node_list.alloc_space*P_node_len, 0);
      #endif /* sgi */

      // add node
      P_node->node_list.number_nodes++;
      j_node--;
      P_node->node_list.nodes[j_node] = P_node_add;

      return;      
    } // end of addNodeLocal()

    //
    // binElementsProcessingLocal()
    //
    void binElementsProcessingLocal(struct element_list_t    *P_element_list,
              struct node_list_t                        *P_node_list,
              const int                                 i_start,
              const int                                 i_end)
    {
      struct element_t *P_element;
      struct node_t    *P_node; 

      int j_node;
      int i_elem;

      // loop over nodes
      for(j_node=i_start; j_node <= i_end; j_node++)
      {
        P_node = P_node_list->nodes[j_node];

        // initialize
        P_node->element_list.number_elements = 0;

        // scan all the elements
        for(i_elem=0; i_elem < P_element_list->number_elements; i_elem++)
        {
          P_element = P_element_list->elements[i_elem];

          // compare all four vertices
          if( P_element->node[0] == P_node ) 
            addElementLocal( P_element, P_node );
          if( P_element->node[1] == P_node ) 
            addElementLocal( P_element, P_node );
          if( P_element->node[2] == P_node ) 
            addElementLocal( P_element, P_node );
          if( P_element->node[3] == P_node ) 
            addElementLocal( P_element, P_node );
        }

        // loop over all adjacent elements and pick up *nodes skipping the
        // * current node
        for(i_elem=0; i_elem < P_node->element_list.number_elements; i_elem++)
        {
          P_element = P_node->element_list.elements[i_elem];

          // add node to P_node's neighbor node list
          if( P_element->node[0] != P_node ) 
            addNodeLocal(P_node, P_element->node[0]);

          if( P_element->node[1] != P_node ) 
            addNodeLocal(P_node, P_element->node[1]);

          if( P_element->node[2] != P_node ) 
            addNodeLocal(P_node, P_element->node[2]);

          if( P_element->node[3] != P_node ) 
            addNodeLocal(P_node, P_element->node[3]);
        }
      }

      return;
    } // end of binElementsProcessingLocal()

    //
    //  binElementsWorker()
    //
    void * binElementsWorker(void * arg)
    {
      struct element_list_t *P_element_list =
        ((struct binElementsData_t *) arg)->P_element_list;
      struct node_list_t *P_node_list = 
        ((struct binElementsData_t *) arg)->P_node_list;

      int number_threads=get_number_threads();
      int my_id=get_my_tid();
      int number_node_thread;
      int number_node_start;
      int number_node_end;

      // get share
      get_share(my_id, number_threads, P_node_list->number_nodes, 
            &number_node_thread, &number_node_start, &number_node_end);
  
      // call binElementsProcessing to do the work
      binElementsProcessingLocal(P_element_list, P_node_list, 
              number_node_start, number_node_end);

      return((void *) NULL);
    } // end of binElementsWorker()

    // ***end of binElements related stuff ***//

    // ***cleamElements related stuff*** //  

    //
    //  cleanElemRefsLocal()
    //
    void cleanElemRefsLocal(struct node_t  *P_node)
    {
      pthread_mutex_lock( &P_node->lock );
      
      P_node->weight=0;
      P_node->site_cluster.number_sites=0;
      P_node->node_list.number_nodes=0;
      P_node->element_list.number_elements=0;

      
      if( P_node->mass != (double *) NULL )
        free(P_node->mass);
        P_node->mass=(double *) NULL;
      
      pthread_mutex_unlock( &P_node->lock );

      return;
    } //end of cleanElemRefsLocal()

    //
    //  zeroElementLocal()
    //
    void zeroElementLocal(struct element_t  *P_element)
    {
      P_element->mask=ELEM_CLEAN_MASK;

      // look at all nodes of this element; if nodes has not been already
      // zeroed do it
      cleanElemRefsLocal( P_element->node[0] );
      cleanElemRefsLocal( P_element->node[1] );
      cleanElemRefsLocal( P_element->node[2] );
      cleanElemRefsLocal( P_element->node[3] );

      return;      
    } // end of zeroElementLocal()

    //
    //  cleanELementsProcessing()
    //
    void cleanElementsProcessingLocal(struct element_list_t *P_element_list,
              const int              i_start,
              const int              i_end)
    {
      int i_tetra;

      for(i_tetra=i_start; i_tetra <= i_end; i_tetra++)
        zeroElementLocal( P_element_list->elements[i_tetra] );
      return;
    } // end of cleanElementsProcessing()

    //
    //  cleanElementsWorker()
    //
    void * cleanElementsWorker(void * arg)
    {
      struct element_list_t *P_element_list=(struct element_list_t *) arg;
  
      int number_threads=get_number_threads();
      int my_id=get_my_tid();
      int number_tetra_thread;
      int number_tetra_start;
      int number_tetra_end;

      // get share

      get_share(my_id, number_threads, P_element_list->number_elements, 
            &number_tetra_thread, &number_tetra_start, &number_tetra_end);

  
      cleanElementsProcessingLocal(P_element_list, number_tetra_start,
          number_tetra_end);

      return((void *) NULL);
    } // end of cleanElementsWorker()

    // ***end of cleanElements related stuff*** //


  } // end of namespace for binElements() and cleanElements()

  //
  //
  //

  Element* Element::_instance = NULL;

  //
  // constructor
  //

  Element::Element()
  {

    //
    //
    //
    return;

  }

  //
  // destructor
  //

  Element::~Element()
  {

    //
    //
    //
    return;

  }

  //
  // getInstance method
  //

  Element*
  Element::getInstance()
  {

    //
    // if not created, create
    //
    if(_instance == NULL){
      _instance = new Element();
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
  Element::destroyInstance()
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

// +++++++++++ Public functions +++++++++++++++++++++//

  // ***************************************//
  // I = binelements.h

  //
  //  binElements()
  //
  void Element::binElements(struct element_list_t      *P_element_list,
             struct node_list_t                 *P_node_list,
             enum mt_version_t                  mt_version)
  {
    struct binElementsData_t binElementsData;
  
    switch( mt_version )
    {
      case MULTI_THREADED:
        binElementsData.P_element_list = P_element_list;
        binElementsData.P_node_list    = P_node_list;

        thread_monitor(binElementsWorker, (void *) &binElementsData,
              get_max_number_threads());
      break;

      case SINGLE_THREADED:
        binElementsProcessingLocal(P_element_list, P_node_list, 0,
              P_node_list->number_nodes-1);
      break;
    }

    return;
  } // end of binElements()


  // ***************************************//
  // I = element.h   

  //
  // makeActiveElement()
  //
  void Element::makeActiveElement(struct element_list_t    *P_element_list,
                  const int                        i_elem)
  {
    init_active_function = 0;

    int j_elem;

    // initialize
    if (init_active_function == 0)
    {
      //extern int errno;
      int old_errno;

      old_errno=errno;

      // for NUMA we really want to allocate more then a page aligned on
      // page boundary
      page_size = sysconf(_SC_PAGESIZE);

      if(page_size < 0 && errno != old_errno)
      {
        ERROR("sysconf()");
        exit(EXIT_FAILURE);
      }

      // compute node bucket; we want it to be as close to 
      // page_size as possible
      elem_bucket = (int)(N_ALLOC_PAGES*page_size/sizeof(struct element_t));

      if(elem_bucket==0) elem_bucket=1;

      init_active_function = 1;
    }

    // check inactive queue
    if( P_element_list->number_inactive_elements > 0 )
    {
      j_elem =--P_element_list->number_inactive_elements;
      
      P_element_list->elements[i_elem] 
        = P_element_list->inactive_elements[j_elem];

      return;
    }

    // there are no inactive elements
    if( P_element_list->alloc_space < i_elem+1 )
    {
      // there is not enough space in the arena-add some more
      if( P_element_list->alloc_space == 0 )
        P_element_list->elements =
          (struct element_t **)malloc(elem_bucket*sizeof(struct element_t *));
      else
        P_element_list->elements =
          (struct element_t **) realloc(P_element_list->elements, 
              (P_element_list->alloc_space + elem_bucket)*
                sizeof(struct element_t *));

      #if defined(__QC_SGI)
          change_memory_placement(P_element_list->elements, 
              (P_element_list->alloc_space + elem_bucket)*
              sizeof(struct element_t *), 0);
      #endif /* sgi */

      if( P_element_list->elements == NULL )
      {
        PERROR_MALLOC();
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      P_element_list->alloc_space += elem_bucket;

      // allocate bucket of space
      #if defined(__QC_SGI)
        P_element_list->elements[i_elem] =
          (struct element_t *) memalign(page_size, 
            elem_bucket*sizeof(struct element_t));

        change_memory_placement(P_element_list->elements[i_elem], 
            elem_bucket*sizeof(struct element_t), 0); 
      #else
        P_element_list->elements[i_elem] = 
          (struct element_t *) calloc(elem_bucket, sizeof(struct element_t));
      #endif /* sgi */

      if( P_element_list->elements[i_elem] == NULL )
      {
        PERROR_MALLOC();
        exit(EXIT_FAILURE);
      }

      #if defined(__QC_SGI)
        memset(P_element_list->elements[i_elem], 0, 
          elem_bucket*sizeof(struct element_t));
      #endif /* sgi */

      pthread_mutex_init(&P_element_list->elements[i_elem]->lock, NULL);

      // index into the bucket
      for(j_elem=1; j_elem < elem_bucket; j_elem++)
      {
        P_element_list->elements[i_elem+j_elem] = 
          P_element_list->elements[i_elem+j_elem-1] + 1;

        pthread_mutex_init(&P_element_list->elements[i_elem+j_elem]->lock, 
          NULL);
      }
    }

    return;
  } // end of makeActiveElement()

  //
  // makeInactiveElement()
  //
  void Element::makeInactiveElement(struct element_list_t    *P_element_list,
                  const int                          i_elem)  
  {
    // check space on the inactive queue
    if( P_element_list->alloc_inactive_space < i_elem+1 )
    {
      if( P_element_list->alloc_inactive_space == 0 )
        P_element_list->inactive_elements =
        (struct element_t **) malloc(elem_bucket*sizeof(struct element_t *));
      else
        P_element_list->inactive_elements = 
        (struct element_t **) realloc(P_element_list->inactive_elements,
          (P_element_list->alloc_inactive_space+elem_bucket)*
          sizeof(struct element_t *));

      if( P_element_list->inactive_elements == NULL )
      {
        PERROR_MALLOC();
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      P_element_list->alloc_inactive_space += elem_bucket;
    }

    P_element_list->inactive_elements[P_element_list->number_inactive_elements]
      = P_element_list->elements[i_elem];

    P_element_list->number_inactive_elements++;

    P_element_list->elements[i_elem] = NULL;

    return;
  } // end of makeInactiveElement()


  // ***************************************//
  // I = free_element.h  

  //
  // cleanElements()
  //
  void Element::cleanElements(struct element_list_t   *P_element_list,
            enum mt_version_t        mt_version)
  {
    int i_elem = P_element_list->number_elements;
    int j_elem;

    // clean all entries
    switch(mt_version)
    {
      case MULTI_THREADED:  
        thread_monitor( cleanElementsWorker, (void *) P_element_list,
              get_max_number_threads() );
      break;

      case SINGLE_THREADED:
        cleanElementsProcessingLocal(P_element_list, 0,
            P_element_list->number_elements-1);
      break;
    }

    // clean P_element_list itself
    for(j_elem=0; j_elem < P_element_list->number_inactive_elements; j_elem++)
        P_element_list->elements[i_elem+j_elem] = 
          P_element_list->inactive_elements[j_elem];

    P_element_list->number_inactive_elements=0;

    P_element_list->number_elements=0;

    return;   
  }


  // ***************************************//
  // I = site_element_map_cache.h 

  //
  // checkHashTableAllocation()
  //
  void Element::checkHashTableAllocation(const int iQuasi)
  {
    // call function check_hash_table_allocation() of 
    // site_element_map_cache.c
    check_hash_table_allocation(iQuasi);

    return;
  }

  //
  // hashTableClean()
  //
  void Element::hashTableClean(const int iQuasi)  
  {
    // call hash_table_clean() of site_element_map_cache.c
    hash_table_clean(iQuasi);

    return;
  }


  // ***************************************//
  // I = site_element_map.h
  
  //
  //  initializeSiteElementMapData()
  //
  void 
  Element::initializeSiteElementMapData(const int iQuasi,
    struct element_list_t *P_element_list)
  {
    // call initialize_site_element_map_data() of site_element_map.c
    initialize_site_element_map_data(iQuasi, P_element_list);

    return;
  }

  //
  //  cleanSiteElementMapData()
  //
  void Element::cleanSiteElementMapData(const int iQuasi)
  {
    // call clean_site_element_map_data() of site_element_map.c
    clean_site_element_map_data(iQuasi);

    return;
  }  

  //
  //  locateSiteElement()
  //
  struct element_t * Element::locateSiteElement(double       *shape,
        const int               l[3],
        struct lattice_t *P_lattice,
        const int               iQuasi)
  {
    // call locate_site_element() of site_element_map.c
    return locate_site_element(shape, l, P_lattice, iQuasi);
  }


  // ***************************************//
  // I = delaunay_3d.h  

  //
  //  computeElementNumberSites()
  //
  void 
  Element::computeElementNumberSites(struct element_list_t  *P_element_list, 
              const struct lattice_t *P_lattice)
  {
    double atomic_volume;
    
    int i_elem;

    // compute attomic volume
    switch (P_lattice->type) 
    {
      case FCC:
        atomic_volume = 2.0*P_lattice->a1[0]*P_lattice->a1[0]*P_lattice->a1[0];
      break;

      case BCC:
        atomic_volume = P_lattice->a1[0]*P_lattice->a2[1]*P_lattice->a3[2];
      break;
    }

    // loop over all elements
    for (i_elem = 0; i_elem < P_element_list->number_elements; i_elem++) 
    {
      struct element_t *P_element = P_element_list->elements[i_elem];

      P_element->number_sites = 
        fabs(MiscFunctions::getInstance()->tetraVolume(P_element->node[0]->initial_position,
        P_element->node[1]->initial_position,
        P_element->node[2]->initial_position,
        P_element->node[3]->initial_position)
        )/atomic_volume;
      /* P_element->number_sites = get_tet_volume(P_element)/atomic_volume; */
    }

    /**
     * done
     */

    return;
  } // end of computeElementNumberSites()

  // ***************************************//
  // I = triangulation_client.h

  //
  //  createElementsFromMeshData()
  //
  //  P_mesh_data is a list of new tet elements from delaunay_3d
  //  In this function, we will create new elements from this tets.
  //  This will follow the procedure of making elements from 
  //  tets in triangulation_cleint.h
  //
  void 
  Element::createElementsFromMeshData(struct element_list_t *P_element_list,
            struct node_list_t      *P_node_list,
            struct lattice_t        *P_lattice,
            struct mesh_data_t      *P_mesh_data)  
  {
    // set number of elements equal to number of tets in P_mesh_data
    P_element_list->number_elements = P_mesh_data->number_tets;

    // loop over all the elements of P_mesh_data and create new element data
    // and fill it
    for(int i_elem = 0; i_elem < P_mesh_data->number_tets; i_elem++)
    {
      // get the node number which are at vertices of tets
      int node_0 = P_mesh_data->tets[i_elem][0];
      int node_1 = P_mesh_data->tets[i_elem][1];
      int node_2 = P_mesh_data->tets[i_elem][2];
      int node_3 = P_mesh_data->tets[i_elem][3];

      double center[3];
      // make the i_elem active
      makeActiveElement(P_element_list, i_elem);

      // fill in the data into element
      P_element_list->elements[i_elem]->mask=ELEM_CLEAN_MASK;
      P_element_list->elements[i_elem]->number=i_elem;
      
      P_element_list->elements[i_elem]->node[0] = P_node_list->nodes[node_0];
      P_element_list->elements[i_elem]->node[1] = P_node_list->nodes[node_1];
      P_element_list->elements[i_elem]->node[2] = P_node_list->nodes[node_2];
      P_element_list->elements[i_elem]->node[3] = P_node_list->nodes[node_3];
      
      P_element_list->elements[i_elem]->cir_radius 
        = MiscFunctions::getInstance()->
          findTetraCenter(P_node_list->nodes[node_0]->initial_position,
            P_node_list->nodes[node_1]->initial_position,
            P_node_list->nodes[node_2]->initial_position,
            P_node_list->nodes[node_3]->initial_position,
            center);
      
      P_element_list->elements[i_elem]->center[0]=center[0];
      P_element_list->elements[i_elem]->center[1]=center[1];
      P_element_list->elements[i_elem]->center[2]=center[2];
    }

    //  compute number of sites inside element 
    computeElementNumberSites(P_element_list, P_lattice);

    // investigate memory placement
    #if defined(__QC_SGI)
      // investigate memory placement
      investigate_element_memory_placement(P_element_list);
    #endif /* sgi */    

    return;
  } // end of createElementsFromMeshData()


  // ***************************************//
  // I = strain.c  

  //
  //  computeDeformGradInfinityNorm()
  //
  double Element::computeDeformGradInfinityNorm(double F[3][3])
  {
    double max=0.0;
    double sum;

    int i;

    for(i=0; i < 3; i++){
      sum=fabs(F[i][0])+fabs(F[i][1])+fabs(F[i][2]);
      max= ((max) > (sum) ? (max) : (sum));
    }

    return(max);    
  } // end of computeDeformGradInfinityNorm()

  //
  //  computeDeformGrad()
  //
  void Element::computeDeformGrad(struct element_t  *P_element)
  {
    struct node_t *P_node1=P_element->node[0];
    struct node_t *P_node2=P_element->node[1];
    struct node_t *P_node3=P_element->node[2];
    struct node_t *P_node4=P_element->node[3];

    double grad_1[3];
    double grad_2[3];
    double grad_3[3];
    double grad_4[3];

    int i_coor;
    int j_coor;

    // compute gradient of shape function for each node of element.
    Shape * shapeC = Shape::getInstance();

    shapeC->computeShapeFunctionGradient(P_node1->initial_position,
         P_node2->initial_position,
         P_node3->initial_position,
         P_node4->initial_position,
         grad_1);
   
    shapeC->computeShapeFunctionGradient(P_node2->initial_position,
         P_node1->initial_position,
         P_node3->initial_position,
         P_node4->initial_position,
         grad_2);

    shapeC->computeShapeFunctionGradient(P_node3->initial_position,
         P_node2->initial_position,
         P_node1->initial_position,
         P_node4->initial_position,
         grad_3);

    shapeC->computeShapeFunctionGradient(P_node4->initial_position,
         P_node2->initial_position,
         P_node3->initial_position,
         P_node1->initial_position,
         grad_4);

    // write the gradient to element data
    for(i_coor=0; i_coor < 3; i_coor++)
      for(j_coor=0; j_coor < 3; j_coor++)
        P_element->F[i_coor][j_coor]=
            P_node1->position[i_coor]*grad_1[j_coor]
          + P_node2->position[i_coor]*grad_2[j_coor]
          + P_node3->position[i_coor]*grad_3[j_coor]
          + P_node4->position[i_coor]*grad_4[j_coor];
  
    return;    
  } // end of computeDeformGrad()

  //
  //  computeElementStrain()
  //
  //  compute strain within an element; returns second strain invariant
  //
  double Element::computeElementStrain(struct element_t *P_element,
             int         i_load)
  {
    double (*F)[3];
    
    double p, deviator[3];
    int k_coor;
    int restart_load = 0;

    // deformation gradient
    computeDeformGrad(P_element);

    F=P_element->F;
    
    #if 0
      /**
       * Simple Shear Calculations:
       * Subtract 'homogeneous' deformation
       * gradient from F before calculating
       * the Lagrangian strain
       *
       *
       *        | 1   g   0 |
       * [Fh] = | 0   1   0 |
       *        | 0   0   1 |
       *
       *           | 1  -g   0 |
       * [Fh]^-1 = | 0   1   0 |
       *           | 0   0   1 |
       *
       *
       * Remeshing is based on non-homogeneous 
       * Lagrangian strain.
       * 
       *
       * [F] = [Fh].[Fnh]
       * [Fnh] = [Fh]^-1.[F]
       *
       * only elements on the top row change:
       */
      
      F[0][0] -= F[1][0]*(i_load + restart_load)*STRETCH;
      F[0][1] -= F[1][1]*(i_load + restart_load)*STRETCH;
      F[0][2] -= F[1][2]*(i_load + restart_load)*STRETCH;
    #endif
    
    // fill-in Lagrangian strain diagonal
    P_element->strain[0]=0.0;

    for(k_coor=0; k_coor < 3; k_coor++)
      P_element->strain[0] += F[k_coor][0]*F[k_coor][0];

    P_element->strain[0] = 0.5*(P_element->strain[0]-1.0);
    
    P_element->strain[1]=0.0;

    for(k_coor=0; k_coor < 3; k_coor++)
      P_element->strain[1] += F[k_coor][1]*F[k_coor][1];

    P_element->strain[1] = 0.5*(P_element->strain[1]-1.0);  

    
    P_element->strain[2]=0.0;

    for(k_coor=0; k_coor < 3; k_coor++)
      P_element->strain[2] += F[k_coor][2]*F[k_coor][2];

    P_element->strain[2] = 0.5*(P_element->strain[2]-1.0);

    /*
     * off-diagonal
     */
    
    P_element->strain[3]=0.0;

    for(k_coor=0; k_coor < 3; k_coor++)
      P_element->strain[3] += F[k_coor][0]*F[k_coor][1];

    P_element->strain[3] = 0.5*P_element->strain[3];

    P_element->strain[4]=0.0;

    for(k_coor=0; k_coor < 3; k_coor++)
      P_element->strain[4] += F[k_coor][0]*F[k_coor][2];

    P_element->strain[4] = 0.5*P_element->strain[4];  

    P_element->strain[5]=0.0;

    for(k_coor=0; k_coor < 3; k_coor++)
      P_element->strain[5] += F[k_coor][1]*F[k_coor][2];

    P_element->strain[5] = 0.5*P_element->strain[5];


    // return second invariant of strain tensor
    /*    return(P_element->strain[3]*P_element->strain[3]+ */
    /*     P_element->strain[4]*P_element->strain[4]+ */
    /*     P_element->strain[5]*P_element->strain[5]- */
    /*     P_element->strain[0]*P_element->strain[1]- */
    /*     P_element->strain[0]*P_element->strain[2]- */
    /*     P_element->strain[1]*P_element->strain[2]); */
    
    // compute deviatoric part
    p = (P_element->strain[0] + P_element->strain[1] 
        + P_element->strain[2]) / 3;

    deviator[0] = P_element->strain[0] - p;
    deviator[1] = P_element->strain[1] - p;
    deviator[2] = P_element->strain[2] - p;
      
    // return second invariant of deviator of strain tensor
    
    return(P_element->strain[3]*P_element->strain[3]+
     P_element->strain[4]*P_element->strain[4]+
     P_element->strain[5]*P_element->strain[5]-
     deviator[0]*deviator[1]-
     deviator[0]*deviator[2]-
     deviator[1]*deviator[2]);
  } // end of computeElementStrain()

}