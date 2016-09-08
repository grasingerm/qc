/**
 * $Id: numa.h,v 1.2 2003/08/14 20:59:28 knap Exp $
 *
 * $Log: numa.h,v $
 * Revision 1.2  2003/08/14 20:59:28  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.4  2002/03/08 00:06:50  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.3  2001/10/03 22:29:48  knap
 * Checked in the latest changes.
 *
 * Revision 1.2  2000/07/25 17:40:51  knap
 * Added some ccNUMA specific functions for HP-UX.
 *
 * Revision 1.1  2000/02/24 01:23:28  knap
 * Initial rev.
 *
 */

#ifndef _NUMA_H_
#define _NUMA_H_

#if defined(__cplusplus)
namespace quasicontinuum {
  extern "C" {
#endif /* __cplusplus */

/* IRIX specific functions */

#if defined(__QC_SGI)

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found
#endif /* HAVE_SYS_TYPES_H */

extern void create_locality_domain(int number_threads);

extern void get_area_placement(caddr_t vaddr,
			       size_t  n_elems,
			       size_t  elem_size,
			       dev_t  *dev_list,
			       int     dev_list_len,
			       int    *n_dev);
extern void change_memory_placement(caddr_t vaddr,
				    size_t  size,
				    int     number_threads);
extern void bind_thread(int cpu);
#endif /* sgi */

/*  LINUX specific functions */

#if defined(HAVE_NUMA_H)
extern void bind_thread(int thread_id);

#endif /* numa */

/* HP-UX specific functions */

#if defined(__QC_HPUX)

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

    /* typedef int pthread_ldom_t; */

extern void bind_to_spu_init(void);
extern pthread_spu_t bind_to_spu(pthread_t tid);
extern void bind_to_ldom_init(void);
extern pthread_ldom_t bind_to_ldom(pthread_t tid);

#endif /* hpux */

#if defined(__cplusplus)
  }
}
#endif /* __cplusplus */

#endif /* _NUMA_H_ */
