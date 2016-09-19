/**
 * $Log: numa.c,v $
 * Revision 1.8  2002/12/04 17:45:43  fago
 * Ported to automake 1.7.1, autoconf 2.56. perror -> ERROR. Fixed build bugs,
 * especially on HPUX.
 *
 * Revision 1.7  2002/03/07 23:53:07  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.6  2001/10/03 22:29:52  knap
 * Checked in the latest changes.
 *
 * Revision 1.5  2000/07/25 17:40:53  knap
 * Added some ccNUMA specific functions for HP-UX.
 *
 * Revision 1.4  2000/05/10 18:42:49  knap
 * Linted.
 *
 * Revision 1.3  2000/03/17 18:58:46  knap
 * get_area_placement() changed to allow for faster allocation tracing.
 *
 * Revision 1.2  2000/02/24 01:12:19  knap
 * Various corrections aimed at improving scalability.
 *
 * Revision 1.1  2000/02/17 21:02:01  knap
 * Initial revs.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#if defined(__QC_HPUX)
#include <errno.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#endif /* hpux */

#if defined(__QC_SGI)
#include <alloca.h>
#include <errno.h>
#include <invent.h>
#include <limits.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/pmo.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/sysmp.h>
#include <sys/types.h>
#include <unistd.h>
#endif /* sgi */

#if defined(HAVE_NUMA_H)
#define _GNU_SOURCE
#include <numa.h>
#include <pthread.h>
#include <sched.h>
#include <stdlib.h>
#endif /* numa */

#include "Error.h"
#include "numa.h"

#if !defined(lint)
const static char rcsid[] = "$Id: numa.c,v 1.8 2002/12/04 17:45:43 fago Exp $";
#endif /* !lint */

#if defined(__QC_SGI)

/**
 * missing prototypes
 */

extern int __pm_get_page_info(void *base_addr, size_t length,
                              pm_pginfo_t *pginfo_buf, int buf_len);
extern void migr_policy_args_init(migr_policy_uparms_t *p);

/**
 * compare two device numbers
 */

static int compare_devs(const void *dev1, const void *dev2) {

  if (*((dev_t *)dev1) == *((dev_t *)dev2))
    return (0);

  return (1);
}

/**
 * given a base address, number of elements and element size display
 * nodes that this memory area has been placed on
 */

#define PGINFO_MAX_NUM_PAGES 1024

/* ARGSUSED4 */
void get_area_placement(caddr_t vaddr, size_t n_elems, size_t elem_size,
                        dev_t *dev_list, int dev_list_len, int *n_devs) {

  pm_pginfo_t *pginfo_buf;

  static int init = 0;
  static long int page_size;
  int i_page;
  int n_pages;
  int n_pages_read;

  /**
   * initialize
   */

  if (init == 0) {

    if ((page_size = sysconf(_SC_PAGESIZE)) < 0) {
      ERROR("sysconf()");
      exit(EXIT_FAILURE);
    }

    init = 1;
  }

  /**
   * get number of pages
   */

  n_pages = (int)((n_elems * elem_size) / page_size);

  n_pages += 2; /* worst case scenario infect two neighboring pages */

  /**
   * size of pginfo_buf is currently 1024
   */

  if (n_pages <= PGINFO_MAX_NUM_PAGES) {

    pginfo_buf = alloca(n_pages * sizeof(pm_pginfo_t));

    /**
     * get placement of all pages
     */

    if ((n_pages_read = __pm_get_page_info(vaddr, n_elems * elem_size,
                                           pginfo_buf, n_pages)) < 0) {
      ERROR("__pm_get_page_info()");
      exit(EXIT_FAILURE);
    }

    /**
     * loop over all pages
     */

    for (i_page = 0; i_page < n_pages_read; i_page++)
      lsearch(&pginfo_buf[i_page].node_dev, dev_list, (size_t *)n_devs,
              sizeof(dev_t), compare_devs);

  } else {

    caddr_t addr_start;
    caddr_t addr_end;

    n_pages = PGINFO_MAX_NUM_PAGES;
    pginfo_buf = alloca(PGINFO_MAX_NUM_PAGES * sizeof(pm_pginfo_t));

    addr_start = vaddr;
    addr_end = vaddr + PGINFO_MAX_NUM_PAGES * page_size;

    while (addr_end < (vaddr + n_elems * elem_size)) {

      if ((n_pages_read = __pm_get_page_info(addr_start, n_pages * page_size,
                                             pginfo_buf, n_pages)) < 0) {
        ERROR("__pm_get_page_info()");
        exit(EXIT_FAILURE);
      }

      /**
       * loop over all pages
       */

      for (i_page = 0; i_page < n_pages_read; i_page++)
        lsearch(&pginfo_buf[i_page].node_dev, dev_list, (size_t *)n_devs,
                sizeof(dev_t), compare_devs);

      addr_start = addr_end;
      addr_end += PGINFO_MAX_NUM_PAGES * page_size;
    }

    /**
     * process the reminder
     */

    if ((n_pages_read = __pm_get_page_info(
             addr_start, vaddr + n_elems * elem_size - addr_start, pginfo_buf,
             n_pages)) < 0) {
      ERROR("__pm_get_page_info()");
      exit(EXIT_FAILURE);
    }

    for (i_page = 0; i_page < n_pages_read; i_page++)
      lsearch(&pginfo_buf[i_page].node_dev, dev_list, (size_t *)n_devs,
              sizeof(dev_t), compare_devs);
  }

  return;
}

/**
 * change placement policy for memory area identified with address and
 * size
 */

#define NUMBER_CPUS_PER_NODE 2

static pmo_handle_t placement;

void create_locality_domain(int number_threads) {

  policy_set_t policy;
  pmo_handle_t mld;
  pmo_handle_t mldset;
  migr_policy_uparms_t migr_parms;

  int mld_radius = 0;

  /**
   * compute radius of MLD
   */

  do {
    mld_radius++;
  } while ((1 << mld_radius) < number_threads);

  /* mld_radius--;*/

  /* printf("MLD radius %d\n", mld_radius); */

  /**
   * initialize MLD
   */

  if ((mld = mld_create(mld_radius, (long)RLIM_INFINITY)) < 0) {

    /**
     * return if memory locality domains are not supoorted
     * or exit otherwise
     */

    if (errno == ENOTSUP) {
      fprintf(stderr, "Support for memory locality domains (MLD) disabled on"
                      " this system architecture.\n");
      fprintf(stderr, "MLD will not be used.\n");
      return;
    }

    ERROR("mld_create()");
    exit(EXIT_FAILURE);
  }

  /**
   * create MLD set
   */

  if ((mldset = mldset_create(&mld, 1)) < 0) {
    ERROR("mldset_create()");
    exit(EXIT_FAILURE);
  }

  /**
   * place MLD set
   */

  if (mldset_place(mldset, TOPOLOGY_FREE, NULL, 0, RQMODE_MANDATORY) < 0) {
    ERROR("mldset_place()");
    exit(EXIT_FAILURE);
  }

  /**
   * attach self to the MLD set
   */

  if (process_mldlink(0, mld, RQMODE_MANDATORY) < 0) {
    ERROR("process_mldlink()");
    exit(EXIT_FAILURE);
  }

  /**
   * fill policy
   */

  pm_filldefault(&policy);

  /**
   * Placement
   */

  /*
  policy.placement_policy_name = "PlacementRoundRobin";
  policy.placement_policy_args = (void *) mldset;
  */

  policy.placement_policy_name = "PlacementDefault";
  policy.placement_policy_args = (void *)number_threads;

  /* policy.placement_policy_name = "PlacementFirstTouch"; */

  /*
  policy.placement_policy_name = "PlacementFixed";
  policy.placement_policy_args = (void *) mld;
  */
  /* policy.placement_policy_name = "PlacementCacheColor";
     policy.placement_policy_args = (void *) mld; */

  /*
  policy.placement_policy_name = "PlacementThreadLocal";
  policy.placement_policy_args = (void *) mldset;
  */

  policy.policy_flags = POLICY_CACHE_COLOR_FIRST;

  /**
   * Migration
   */

  migr_policy_args_init(&migr_parms);
  migr_parms.migr_base_enabled = 1;
  migr_parms.migr_base_threshold = 90;
  migr_parms.migr_dampening_enabled = 0;

  policy.migration_policy_name = "MigrationControl";
  policy.migration_policy_args = (void *)&migr_parms;

  /**
   * create memory placement
   */

  if ((placement = pm_create(&policy)) < 0) {
    ERROR("pm_create()");
    exit(EXIT_FAILURE);
  }

  /**
   * set policy as a default policy for stack. text and data
   */

  if (pm_setdefault(placement, MEM_STACK) < 0) {
    ERROR("pm_setdefault()");
    exit(EXIT_FAILURE);
  }

  if (pm_setdefault(placement, MEM_TEXT) < 0) {
    ERROR("pm_setdefault()");
    exit(EXIT_FAILURE);
  }

  if (pm_setdefault(placement, MEM_DATA) < 0) {
    ERROR("pm_setdefault()");
    exit(EXIT_FAILURE);
  }

  return;
}

void change_memory_placement(caddr_t vaddr, size_t size, int number_threads) {

  /**
   * associate placement and address
   */

  if (pm_attach(placement, vaddr, size) < 0) {
    ERROR("pm_atach()");
    exit(EXIT_FAILURE);
  }

  return;
}

/**
 * attatch thread to CPU
 */

void bind_thread(int cpu) {

  int max_cpus;

  if ((max_cpus = sysmp(MP_NAPROCS)) == -1) {
    ERROR("sysmp(MP_NAPROCS)");
    exit(EXIT_FAILURE);
  }

  if (cpu >= max_cpus)
    cpu -= max_cpus;

  if (sysmp(MP_MUSTRUN, cpu) == -1) {
    ERROR("sysmp(MP_MUSTRUN,...)");
    exit(EXIT_SUCCESS);
  }

  printf("[tid:%d] bound to cpu %d\n", pthread_self(), sysmp(MP_GETMUSTRUN));

  return;
}

#endif /* sgi */

#if defined(HAVE_NUMA_H)

/**
 * attatch thread to CPU
 */

void bind_thread(int thread_id) {

  /**
   * get max number of memory nodes
   */
  const int NUM_NODES = numa_num_configured_nodes();

  /**
   * get max number of cpus
   */
  const int NUM_CPUs = numa_num_configured_cpus();

  /**
   * get max number of cpus
   */
  const int CPUs_PER_NODE = NUM_CPUs / NUM_NODES;

  /**
   * choose threading scheme
   */
  const int THREADING_SCHEME = 1;

  /**
   * create cpu affinity mask
   */
  cpu_set_t cpuset;

  /**
   * zero out affinity mask
   */
  CPU_ZERO(&cpuset);

  /**
   * variable to hold cpu id, initialized to -1 for error checking
   */
  int cpu_id = -1;

  /**
   * switch between thread binding schemes
   */
  switch (THREADING_SCHEME) {

  /**
   * bind thread using cpu id = thread id
   */
  case 0: {

    /**
     * change thread id to  beginning if greater than number of cpus
     * only needed for more threads than cpus
     */
    if (thread_id >= NUM_NODES * CPUs_PER_NODE)
      thread_id = thread_id % (NUM_NODES * CPUs_PER_NODE);

    /**
     * calculate cpu_id
     */
    cpu_id = thread_id;

  }

  break;

  /**
   * bind thread giving each core 1 thread before adding 2nd thread
   */
  case 1: {

    /**
     * change thread id to  beginning if greater than number of cpus
     * only needed for more threads than cpus
     */
    if (thread_id >= NUM_NODES * CPUs_PER_NODE)
      thread_id = thread_id % (NUM_NODES * CPUs_PER_NODE);

    /**
     * toggle between assigning first and second thread to a core
     */
    if (thread_id < (NUM_NODES * CPUs_PER_NODE) / 2) {

      /**
       * calculate cpu_id
       */
      cpu_id = 2 * thread_id;

    } else if (thread_id >= (NUM_NODES * CPUs_PER_NODE) / 2) {

      /**
       * calculate cpu_id
       */
      cpu_id = (2 * thread_id) % (NUM_NODES * CPUs_PER_NODE) + 1;
    }

  }

  break;
#if 0
    /**
     * fill all cpus on a node first
     */
  case 2:
    {

      /**
       * change thread id to  beginning if greater than number of cpus
       * only needed for more threads than cpus
       */
      if (thread_id >= NUM_NODES * CPUs_PER_NODE)
  	thread_id = thread_id % (NUM_NODES * CPUs_PER_NODE);

      /**
       * loop over all nodes
       */
      for (i_node = 1; i_node <= NUM_NODES; ++i_node){

  	/**
  	 * toggle between nodes
  	 */
  	if(thread_id < i_node * CPUs_PER_NODE){

  	  /**
  	   * toggle between cores
  	   */
  	  if(thread_id < ((i_node * CPUs_PER_NODE) - CPUs_PER_NODE / 2)){

  	    /**
  	     * calculate cpu_id
  	     */
  	    cpu_id = (thread_id % CPUs_PER_NODE) * 2);

  	  }

  	}

      }

    break;
#endif
  }

  /**
   * check to make sure cpu_id has been set
   */
  if (cpu_id == -1) {
    ERROR("cpu_id");
    exit(EXIT_FAILURE);
  }

  /**
   * set cpu affinity
   */
  CPU_SET(cpu_id, &cpuset);

  /**
   * set thread to run on specific cpu
   */
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

  /**
   * print thread binding
   */
  printf("[tid:%d] bound to cpu %d\n", thread_id, cpu_id);

  /**
   *
   */
  return;
}

#endif /* numa */

#if defined(__QC_HPUX)

/**
 * HP-UX specific functions to access numa
 */

static pthread_spu_t current_spu;

void bind_to_spu_init(void) {

  if (pthread_processor_id_np(PTHREAD_GETFIRSTSPU_NP, &current_spu,
                              current_spu)) {
    fprintf(stderr, "pthread_processor_id_np() failed\n");
    exit(EXIT_FAILURE);
  }

  return;
}

pthread_spu_t bind_to_spu(pthread_t tid) {

  pthread_spu_t answer;
  int err;

  /**
   * bind thread to the current spu
   */

  if (pthread_processor_bind_np(PTHREAD_BIND_FORCED_NP, &answer, current_spu,
                                tid)) {
    fprintf(stderr, "pthread_processor_bind_np() failed\n");
    exit(EXIT_FAILURE);
  }

  /**
   * get next spu
   */

  if (err = pthread_processor_id_np(PTHREAD_GETNEXTSPU_NP, &current_spu,
                                    current_spu)) {

    /**
     * pthread_processor_id_np() may fail if the last processor
     * has been used; reset current_spu back to the first one
     */

    switch (err) {

    case EINVAL:
      bind_to_spu_init();
      break;

    default:
      fprintf(stderr, "pthread_processor_id_np() failed\n");
      exit(EXIT_FAILURE);
      break;
    }
  }

  return (answer);
}

/**
 * add threads to the appropriate locality domain
 */

static pthread_ldom_t current_ldom;
static int current_ldom_spu;

void bind_to_ldom_init(void) {

  if (pthread_ldom_id_np(PTHREAD_GETFIRSTLDOM_NP, &current_ldom,
                         current_ldom)) {
    fprintf(stderr, "pthread_ldom_id_np() failed\n");
    exit(EXIT_FAILURE);
  }

  current_ldom_spu = 0;

  return;
}

pthread_ldom_t bind_to_ldom(pthread_t tid) {

  pthread_ldom_t answer;
  int n_ldom_spu;
  int err;

  /**
   * bind thread to a current LDOM
   */

  if (pthread_ldom_bind_np(&answer, current_ldom, tid)) {
    fprintf(stderr, "pthread_ldom_bind_np() failed\n");
    exit(EXIT_FAILURE);
  }

  /**
   * increase the number of threads placed in the current
   * LDOM
   */

  if (pthread_num_ldomprocs_np(&n_ldom_spu, current_ldom)) {
    fprintf(stderr, "pthread_num_ldomprocs_np() failed\n");
    exit(EXIT_FAILURE);
  }

  if (current_ldom_spu++ == n_ldom_spu) {

    /**
     * switch to the next LDOM
     */

    if (err = pthread_ldom_id_np(PTHREAD_GETNEXTLDOM_NP, &current_ldom,
                                 current_ldom)) {

      switch (err) {

      case EINVAL:
        bind_to_ldom_init();
        break;

      default:
        fprintf(stderr, "pthread_ldom_id_np() failed\n");
        exit(EXIT_FAILURE);
        break;
      }
    }
  }

  return (answer);
}

#endif /* hpux */
