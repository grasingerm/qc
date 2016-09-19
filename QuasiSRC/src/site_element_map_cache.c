/**
 * $Log: site_element_map_cache.c,v $
 * Revision 1.8  2002/03/07 23:53:09  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.7  2001/10/03 22:29:52  knap
 * Checked in the latest changes.
 *
 * Revision 1.6  2000/07/19 18:54:18  knap
 * Corrected rw-locks for Solaris.
 *
 * Revision 1.5  2000/05/24 18:35:39  knap
 * Ported to Linux and OSF.
 *
 * Revision 1.4  2000/05/10 18:42:57  knap
 * Linted.
 *
 * Revision 1.3  2000/05/09 20:05:44  knap
 * corrected cast to lfind().
 *
 * Revision 1.2  2000/04/20 01:09:11  knap
 * Ported to HP-UX.
 *
 * Revision 1.1  2000/04/12 22:07:21  knap
 * Removed references to atoms.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#if defined(__QC_LINUX)
#define _XOPEN_SOURCE 500
#endif /* __linux__ */

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
/* No _REENTRANT */
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#ifdef STDC_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_STRING_H
#include <search.h>
#else
#error search.h not found
#endif /* HAVE STRING_H */

#if defined(__QC_SUN)
#ifdef HAVE_SYNCH_H
#include <synch.h>
#endif /* HAVE_SYNCH_H */
#endif /* sun */

#include "C_Interface.h"
#include "DataTypes.h"
#include "Error.h"

#include "site_element_map_cache.h"

#if !defined(lint)
static const char rcsid[] =
    "$Id: site_element_map_cache.c,v 1.8 2002/03/07 23:53:09 knap Exp $";
#endif /* !lint */

/**
 * hash table defs
 */

#define WIDTH_0 5
#define WIDTH_1 5
#define WIDTH_2 5

#define MODULO_BITS 21
#define MODULO_MASK ~(~0 << MODULO_BITS)
/* #define NUMBER_BUCKETS ((unsigned) 1 << MODULO_BITS) */
#define NUMBER_BUCKETS 64997
/* #define NUMBER_BUCKETS ((unsigned) 1 << (WIDTH_0+WIDTH_1+WIDTH_2)) */

/**
 * local data defs
 */

struct datum_t {
  int l[3];
  struct element_t *P_element;
  double shape[4];
};

struct hash_table_t {
  struct hash_bucket_t {
    struct datum_t *data;
    int n_data;
    int alloc_space;
#define BUCKET_INCREMENT_ALLOC 20

/* Should check if types defined, but can't get it working reliably */
#if defined(__QC_SUN)
    rwlock_t lock;
#else
    pthread_rwlock_t lock;
#endif /* sun */

    /* pthread_mutex_t  lock; */
  } buckets[NUMBER_BUCKETS];
};

/**
 * local data
 */

static struct hash_table_t *hash_table;
/*
static unsigned int mask_0 = ~( ~0 << WIDTH_0);
static unsigned int mask_1 = ~( ~0 << WIDTH_1);
static unsigned int mask_2 = ~( ~0 << WIDTH_2);
*/
static pthread_rwlock_t rwlock =
    PTHREAD_RWLOCK_INITIALIZER; // FIXME: May not work on Sun
static pthread_mutex_t hash_table_alloc_lock = PTHREAD_MUTEX_INITIALIZER;
static int hash_table_allocated = 0;

/**
 * local protos
 */

static struct element_t *hash_table_search_element_unlocked(const int l[3],
                                                            const int iQuasi);

/**
 * hash function; take WIDTH_* bits from each site coordinate and
 * create (WIDTH_0+WIDTH_1+WIDTH_2) index
 */

/**
#define hash_function(l) ((l)[2] & mask_2) + (((l)[1] & mask_1) << WIDTH_2) + \
                         (((l)[0] & mask_0) << (WIDTH_1+WIDTH_2))
**/
/**
#define hash_function(l) ((((l)[0])+((l)[1] << 1)+((l)[2] << 3))%32768)
**/
/**
#define hash_function(l) ((((l)[0] << 1)+((l)[1] << 4)+((l)[2] << 8
))&MODULO_MASK)
**/

/*
#define WIDTH_BIT_MASK (~((~(unsigned) 0) << WIDTH_0))

#define BITS_L1_MASK (WIDTH_BIT_MASK << 1)
#define BITS_L2_MASK (WIDTH_BIT_MASK << 4)
#define BITS_L3_MASK (WIDTH_BIT_MASK << 7)

#define hash_function(l) (((l)[0] & BITS_L1_MASK) | ((l[1]) & BITS_L2_MASK) |
((l)[2] & BITS_L3_MASK))
*/

/**
 * hash function
 */

static int hash_function(const int l[3]) {

  int n_bytes = sizeof(int[3]);
  const unsigned char *bytes = (const unsigned char *)&(l[0]);
  int h;
  int i_byte;

  h = bytes[0];

  for (i_byte = 1; i_byte < n_bytes; i_byte++)
    h = ((h * 16) + bytes[i_byte]) % NUMBER_BUCKETS;

  return (h);
}

/**
 * init
 */

static void hash_table_init(const int iQuasi) {
  /**
   * check if allocation has occurred
   */
  if (iQuasi >= hash_table_allocated) {

    /**
     * check if this is first allocation
     */
    if (iQuasi == 0) {

      /**
       * malloc hash_table
       */
      hash_table = malloc(sizeof(struct hash_table_t));

    } else {

      /**
       * reallocate hash_table
       */
      hash_table =
          realloc(hash_table, (iQuasi + 1) * sizeof(struct hash_table_t));
    }

    /**
     * check to make sure hash_table allocated
     */
    if (hash_table == NULL) {
      PERROR_MALLOC();
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }
  }

  /**
   * define bucket
   */
  int j_bucket;

  /**
    * counter to go up to currentId
    */
  int j_Quasi;

  /**
   * set n_data and alloc_space in current Quasicontinuum to 0
   */
  for (j_bucket = 0; j_bucket < NUMBER_BUCKETS; j_bucket++) {
    hash_table[iQuasi].buckets[j_bucket].n_data = 0;
    hash_table[iQuasi].buckets[j_bucket].alloc_space = 0;
  }

  /**
   * loop over Quasicontinuum
   */
  for (j_Quasi = 0; j_Quasi < (iQuasi + 1); j_Quasi++) {

    /**
     * loop over each bucket
     */
    for (j_bucket = 0; j_bucket < NUMBER_BUCKETS; j_bucket++) {
/* if (pthread_mutex_init(&(hash_table.buckets[i_bucket].lock), NULL)) { */

#if defined(__QC_SUN) /* Sun has pthread_rwlock_init, but doesn't work right   \
                         */
      if (rwlock_init(&(hash_table[j_Quasi].buckets[j_bucket].lock),
                      USYNC_THREAD, NULL)) {
        ERROR("rwlock_init()");
#else
#ifdef HAVE_PTHREAD_RWLOCK_INIT
      if (pthread_rwlock_init(&(hash_table[j_Quasi].buckets[j_bucket].lock),
                              NULL)) {
        ERROR("pthread_mutex_init()");
#else
#error No pthread_rwlock_init() found
#endif /* HAVE_PTHREAD_RWLOCK_INIT */
#endif /* SUN */

        exit(EXIT_FAILURE);
      }
    }
  }

  /**
   * get write lock for hash_table_allocated
   */
  pthread_rwlock_wrlock(&rwlock);

  /**
   * increment hash table
   */
  ++hash_table_allocated;

  /**
   * unlock write lock for hash_table_allocated
   */
  pthread_rwlock_unlock(&rwlock);

  /**
   *
   */
  return;
}

/**
 * add data to a bucket
 */

static void add_data_to_bucket(struct hash_bucket_t *P_bucket,
                               const struct datum_t *P_datum) {

  int i_datum = P_bucket->n_data + 1;

  /**
   * check for allocated space
   */

  if (P_bucket->alloc_space == 0) {

    if ((P_bucket->data =
             malloc(BUCKET_INCREMENT_ALLOC * sizeof(struct datum_t))) == NULL) {
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    P_bucket->alloc_space = BUCKET_INCREMENT_ALLOC;

  } else if (i_datum > P_bucket->alloc_space) {

    P_bucket->alloc_space += BUCKET_INCREMENT_ALLOC;

    // if(P_bucket->alloc_space > 100)
    //   printf("P_bucket->alloc_space = %d\n", P_bucket->alloc_space);

    if ((P_bucket->data = realloc(P_bucket->data, sizeof(struct datum_t) *
                                                      P_bucket->alloc_space)) ==
        NULL) {
      ERROR("realloc()");
      exit(EXIT_FAILURE);
    }
  }

  P_bucket->n_data++;
  i_datum--;
  P_bucket->data[i_datum] = *P_datum;

  return;
}

/**
 * clear all hash table locks
 */
static void hash_table_clear_locks(const int iQuasi) {
  /**
   * define bucket
   */
  int j_bucket;

  /**
   * counter to go up to currentId
   */
  int j_Quasi;

  /**
   * loop over Quasicontinuum
   */
  for (j_Quasi = 0; j_Quasi < iQuasi; j_Quasi++) {

    /**
     * loop over all buckets in current Quasicontinuum
     */
    for (j_bucket = 0; j_bucket < NUMBER_BUCKETS; j_bucket++) {
/* if (pthread_mutex_init(&(hash_table.buckets[i_bucket].lock), NULL)) { */

#if defined(__QC_SUN) /* Sun has pthread_rwlock_init, but doesn't work right   \
                         */
      if (rwlock_destroy(&(hash_table[j_Quasi].buckets[j_bucket]).lock)) {
        ERROR("rwlock_init()");
#else
#ifdef HAVE_PTHREAD_RWLOCK_INIT
      if (pthread_rwlock_destroy(
              &(hash_table[j_Quasi].buckets[j_bucket].lock))) {
        ERROR("pthread_mutex_init()");
#else
#error No pthread_rwlock_init() found
#endif /* HAVE_PTHREAD_RWLOCK_INIT */
#endif /* SUN */

        exit(EXIT_FAILURE);
      }
    }
  }

  /**
   *
   */
  return;
}

/**
 * check allocation
 */
void check_hash_table_allocation(const int iQuasi) {

  /**
   * check if allocation has occurred
   */
  if (iQuasi >= hash_table_allocated) {

    /**
     * check if this is not first allocation
     */
    if (iQuasi > 0) {

      /**
       * clear all locks
       */
      hash_table_clear_locks(iQuasi);
    }

    /**
     * initialize hash table
     */
    printf("Initializing hash table...");
    hash_table_init(iQuasi);
    printf("Completed\n");
  }

  /**
   *
   */
  return;
}

/**
 * add site <-> element pair to the hash table
 */

void hash_table_add_site(double *shape, const int l[3],
                         struct element_t *P_element,
                         struct lattice_t *P_lattice, const int iQuasi) {

  double X[3];

  struct node_t *P_node_0 = P_element->node[0];
  struct node_t *P_node_1 = P_element->node[1];
  struct node_t *P_node_2 = P_element->node[2];
  struct node_t *P_node_3 = P_element->node[3];

  struct datum_t datum;
  int hash_table_index;

  /**
   * get read lock for hash_table_allocated
   */
  pthread_rwlock_rdlock(&rwlock);

  /**
   * check if allocation has occurred
   */
  if (iQuasi >= hash_table_allocated) {

    /**
     * unlock read lock for hash_table_allocated
     */
    pthread_rwlock_unlock(&rwlock);

    /**
     * lock allocation mutex
     */
    pthread_mutex_lock(&hash_table_alloc_lock);

    /**
     * check allocation and reallocate if necessary
     */
    check_hash_table_allocation(iQuasi);

    /**
     * lock allocation mutex
     */
    pthread_mutex_unlock(&hash_table_alloc_lock);

  } else {

    /**
     * unlock read lock for hash_table_allocated
     */
    pthread_rwlock_unlock(&rwlock);
  }

  /**
   * fill in datum
   */

  datum.l[0] = l[0];
  datum.l[1] = l[1];
  datum.l[2] = l[2];
  datum.P_element = P_element;

  /**
   * find site position
   */

  getSiteInitialPosition(X, l, P_lattice, iQuasi);

  /**
   * fill in shape funs
   */

  datum.shape[0] = computeShapeFunction(
      P_node_0->initial_position, P_node_1->initial_position,
      P_node_2->initial_position, P_node_3->initial_position, X);

  datum.shape[1] = computeShapeFunction(
      P_node_1->initial_position, P_node_0->initial_position,
      P_node_2->initial_position, P_node_3->initial_position, X);

  datum.shape[2] = computeShapeFunction(
      P_node_2->initial_position, P_node_1->initial_position,
      P_node_0->initial_position, P_node_3->initial_position, X);

  datum.shape[3] = computeShapeFunction(
      P_node_3->initial_position, P_node_1->initial_position,
      P_node_2->initial_position, P_node_0->initial_position, X);

  if (shape) {

    shape[0] = datum.shape[0];
    shape[1] = datum.shape[1];
    shape[2] = datum.shape[2];
    shape[3] = datum.shape[3];
  }

  /**
   * compute hash table index
   */

  hash_table_index = hash_function(l);

/**
 * get write mutex associated with the bucket
 */

#if defined(__QC_SUN) /* Sun has pthread_rwlock_wrlock, doesn't work */
  rw_wrlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#ifdef HAVE_PTHREAD_RWLOCK_WRLOCK
  pthread_rwlock_wrlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#error pthread_rwlock_wrlock() not found
#endif /* HAVE_PTHREAD_RWLOCK_WRLOCK */
#endif /* SUN */

  /**
   * check if l is in the hash_table (we do MT, so other thread
   * may have done it)
   */

  if (hash_table_search_element_unlocked(l, iQuasi) != NULL) {

#if defined(__QC_SUN) /* Sun has pthread_rwlock_unlock, but doesn't work */
    rw_unlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#ifdef HAVE_PTHREAD_RWLOCK_UNLOCK
    pthread_rwlock_unlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#error pthread_rwlock_unlock() not found
#endif /* HAVE_PTHREAD_RWLOCK_UNLOCK */
#endif /* SUN */

    return;
  }

  /**
   * add datum to the bucket
   */

  add_data_to_bucket(&(hash_table[iQuasi].buckets[hash_table_index]), &datum);

/**
 * unlock lock
 */

#if defined(__QC_SUN) /* Sun has pthread version, but doesn't work? */
  rw_unlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#ifdef HAVE_PTHREAD_RWLOCK_UNLOCK
  pthread_rwlock_unlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#error pthread_rwlock_unlock() not found
#endif /* HAVE_PTHREAD_RWLOCK_UNLOCK */
#endif /* SUN */

  return;
}

/**
 * function to compare data in a bucket
 */

static int compare_datum(const void *d0, const void *d1) {

  const int *l0 = ((const struct datum_t *)d0)->l;
  const int *l1 = ((const struct datum_t *)d1)->l;

  if ((l0[0] == l1[0]) && (l0[1] == l1[1]) && (l0[2] == l1[2]))
    return (0);
  /* return(memcmp(l0, l1, sizeof(int[3]))); */
  return (1);
}

/**
 * search hash table for particular site<->element pair
 */

static struct element_t *hash_table_search_element_unlocked(const int l[3],
                                                            const int iQuasi) {

  struct datum_t datum;
  struct datum_t *P_datum;
  size_t n_data;
  int hash_table_index;

  /**
   * get read lock for hash_table_allocated
   */
  pthread_rwlock_rdlock(&rwlock);

  /**
   * check if allocation has occurred
   */
  if (iQuasi >= hash_table_allocated) {

    /**
     * unlock read lock for hash_table_allocated
     */
    pthread_rwlock_unlock(&rwlock);

    /**
     * lock allocation mutex
     */
    pthread_mutex_lock(&hash_table_alloc_lock);

    /**
     * check allocation and reallocate if necessary
     */
    check_hash_table_allocation(iQuasi);

    /**
     * lock allocation mutex
     */
    pthread_mutex_unlock(&hash_table_alloc_lock);

  } else {

    /**
     * unlock read lock for hash_table_allocated
     */
    pthread_rwlock_unlock(&rwlock);
  }

  /**
   * get index to hash_table
   */

  hash_table_index = hash_function(l);

  /**
   * search bucket for site
   */

  datum.l[0] = l[0];
  datum.l[1] = l[1];
  datum.l[2] = l[2];

  n_data = hash_table[iQuasi].buckets[hash_table_index].n_data;
  P_datum =
      lfind(&datum, &(hash_table[iQuasi].buckets[hash_table_index].data[0]),
            (size_t *)&(n_data), sizeof(struct datum_t), compare_datum);

  if (P_datum == NULL)
    return (NULL);

  return (P_datum->P_element);
}

/**
 * MT-safe version of hash table search
 * new version does not depend on current id of quasi
 */

struct element_t *hash_table_search_element(double *shape, const int l[3],
                                            const int iQuasi) {

  struct datum_t datum;
  struct datum_t *P_datum;
  size_t n_data;
  int hash_table_index;

  /**
   * get read lock for hash_table_allocated
   */
  pthread_rwlock_rdlock(&rwlock);

  /**
   * check if allocation has occurred
   */
  if (iQuasi >= hash_table_allocated) {

    /**
     * unlock read lock for hash_table_allocated
     */
    pthread_rwlock_unlock(&rwlock);

    /**
     * lock allocation mutex
     */
    pthread_mutex_lock(&hash_table_alloc_lock);

    /**
     * check allocation and reallocate if necessary
     */
    check_hash_table_allocation(iQuasi);

    /**
     * lock allocation mutex
     */
    pthread_mutex_unlock(&hash_table_alloc_lock);

  } else {

    /**
     * unlock read lock for hash_table_allocated
     */
    pthread_rwlock_unlock(&rwlock);
  }

  /**
   * get index to hash_table
   */

  hash_table_index = hash_function(l);

  /**
   * search bucket for site
   */

  datum.l[0] = l[0];
  datum.l[1] = l[1];
  datum.l[2] = l[2];

/**
 * acquire lock
 */

#if defined(__QC_SUN) /* Sun pthread version doesn't work? */
  rw_rdlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#ifdef HAVE_PTHREAD_RWLOCK_RDLOCK
  pthread_rwlock_rdlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#error no pthread_rwlock_rdlock() found.
#endif /* HAVE_PTHREAD_RWLOCK_RDLOCK */
#endif /* SUN */

  /**
   * locate datum in the bucket
   */

  n_data = hash_table[iQuasi].buckets[hash_table_index].n_data;
  P_datum =
      lfind(&datum, &(hash_table[iQuasi].buckets[hash_table_index].data[0]),
            (size_t *)&(n_data), sizeof(struct datum_t), compare_datum);

/**
 * release lock
 */

#if defined(__QC_SUN) /* Sun's pthread doesn't work? */
  rw_unlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#ifdef HAVE_PTHREAD_RWLOCK_UNLOCK
  pthread_rwlock_unlock(&(hash_table[iQuasi].buckets[hash_table_index].lock));
#else
#error pthread_rwlock_unlock() not found
#endif /* HAVE_PTHREAD_RWLOCK_UNLOCK */
#endif /* SUN */

  if (P_datum == NULL)
    return (NULL);

  if (shape) {

    shape[0] = P_datum->shape[0];
    shape[1] = P_datum->shape[1];
    shape[2] = P_datum->shape[2];
    shape[3] = P_datum->shape[3];
  }

  return (P_datum->P_element);
}

/**
 * clean cache
 */

void hash_table_clean(const int iQuasi) {

  int i_bucket;

  for (i_bucket = 0; i_bucket < NUMBER_BUCKETS; i_bucket++)
    hash_table[iQuasi].buckets[i_bucket].n_data = 0;

  return;
}

/* /\** */
/*  * dump contents of the hash table */
/*  *\/ */

/* void */
/* hash_table_dump(void) */
/* { */

/*   int i_bucket; */

/*   for (i_bucket = 0; i_bucket < NUMBER_BUCKETS; i_bucket++) */
/*     printf("%d %d\n", i_bucket, hash_table.buckets[i_bucket].n_data); */

/*   return; */

/* } */

/* /\** */
/*  * compute hash function for all cluster sites in a lattice */
/*  *\/ */

/* void */
/* check_hash_table(struct node_list_t *P_node_list) */
/* { */

/*   int i_node; */
/*   int i_site; */

/*   /\** */
/*    * initialize table */
/*    *\/ */

/*    hash_table_init(); */

/*    /\** */
/*     * put all sites into the hash table */
/*     *\/ */

/*    for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) { */

/*      struct node_t *P_node = P_node_list->nodes[i_node]; */

/*      for (i_site = 0; i_site < P_node->site_cluster.number_sites; i_site++) {
 */

/*        struct datum_t datum; */
/*        int hash_table_index; */

/*        hash_table_index = hash_function(P_node->site_cluster.sites[i_site]);
 */

/*        datum.l[0] = P_node->site_cluster.sites[i_site][0]; */
/*        datum.l[1] = P_node->site_cluster.sites[i_site][1]; */
/*        datum.l[2] = P_node->site_cluster.sites[i_site][2]; */
/*        datum.P_element = NULL; */

/*        add_data_to_bucket(&(hash_table.buckets[hash_table_index]), &datum);
 */

/*      } */

/*    } */

/*    /\** */
/*     * dump contents of the hash table */
/*     *\/ */

/*    hash_table_dump(); */

/*   return; */

/* } */

/* /\** */
/*  * print the distribution of sites in the hash table for all  */
/*  * cluster atoms */
/*  *\/ */

/* void */
/* print_hash_distribution(const struct node_list_t *P_node_list) */
/* { */

/*   int distrib[NUMBER_BUCKETS]; */
/*   int i_node; */

/*   /\** */
/*    * initialize distrib table */
/*    *\/ */

/*   memset(distrib, '\0', sizeof(distrib)); */

/*   /\** */
/*    * fill distrib table */
/*    *\/ */

/*   for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) { */

/*     struct node_t *P_node = P_node_list->nodes[i_node]; */
/*     int i_site; */

/*     for (i_site = 0; i_site < P_node->site_cluster.number_sites; i_site++) */
/*       (distrib[hash_function(P_node->site_cluster.sites[i_site])])++; */

/*   } */

/*   /\** */
/*    * print distrib table */
/*    *\/ */

/*   for (i_node = 0; i_node < NUMBER_BUCKETS; i_node++) */
/*     printf("%d %d\n", i_node, distrib[i_node]); */

/*   return; */

/* } */
