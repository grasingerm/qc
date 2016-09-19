/**
 * pthread related utility routines
 *
 * $Log: threads.c,v $
 * Revision 1.8  2002/12/09 18:24:03  fago
 * Another minor fix for AIX build.
 *
 * Revision 1.7  2002/12/06 01:49:45  fago
 * Significant fixes to automake build especially on AIX and HPUX.
 *
 * Revision 1.6  2002/03/07 23:53:09  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.5  2001/10/03 22:29:53  knap
 * Checked in the latest changes.
 *
 * Revision 1.4  2000/07/25 18:16:55  knap
 * Added launch policy for threads on HP-UX.
 *
 * Revision 1.3  2000/05/24 18:35:41  knap
 * Ported to Linux and OSF.
 *
 * Revision 1.2  2000/04/20 01:09:12  knap
 * Ported to HP-UX.
 *
 * Revision 1.1  2000/01/07 00:26:56  knap
 * Moved to src.
 *
 * Revision 1.10  1999/09/02 23:48:46  knap
 * Corrected all calls to thread_monitor() to get new interface.
 *
 * Revision 1.9  1999/07/23 22:40:26  knap
 * _REENTRANT caused some complaints on osf.
 *
 * Revision 1.8  1999/04/19 19:02:56  knap
 * Modified amount of worle per thread.
 *
 * Revision 1.7  1999/04/13 21:48:07  knap
 * Addedd missing pthread.h include.
 *
 * Revision 1.6  1999/04/11 19:33:12  knap
 * Added pthread support for OSF1.
 *
 * Revision 1.5  1999/04/07 22:50:11  knap
 * Linted.
 *
 * Revision 1.4  1999/04/07 19:25:48  knap
 * Minor corrections.
 *
 * Revision 1.3  1999/03/31 19:47:41  knap
 * Added:
 * get_share()
 *
 * Revision 1.2  1999/03/22 21:03:48  knap
 * Added get_share() to manage data distribution across multiple threads.
 *
 * Revision 1.1.1.1  1999/02/12 18:55:14  knap
 * Importing initial version
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#if !defined(_REENTRANT)
#define _REENTRANT
#endif
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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found.
#endif /* HAVE_UNISTD_H */

#ifdef HAVE_ASSERT_H
#include <assert.h>
#else
#error assert.h not found.
#endif /* HAVE_ASSERT_H */

#ifdef HAVE_SYS_MPCTL_H
#include <sys/mpctl.h>
#endif /*  HAVE_SYS_MPCTL_H */

#if defined(__QC_SGI)
#include "numa.h"
#endif /* sgi */

#if !defined(lint)
#include "threads.h"
#endif /* !lint */

#if !defined(lint)
static char rcsid[] = "$Id: threads.c,v 1.8 2002/12/09 18:24:03 fago Exp $";
#endif /* !lint */

#define N_THREADS_ENV "NUMBER_THREADS"

static int max_number_threads;
static int number_threads = 1;

/**
 *initialize;
 *
 * number of threads to use is supplied by:
 * 1) -n command line option
 * 2) NUMBER_THREADS enviroment variable
 * in case both are set, command line takes precedence
 *
 */

void qc_thread_init(int n_threads) {
  char *env;

  number_threads = 1;

#ifdef HAVE_PTHREAD_NUM_PROCESSORS_NP
#ifndef __QC_AIX /* Doesn't seem to work on AIX */
  number_threads = pthread_num_processors_np();
#endif /* !__QC_AIX */
#else  /* HAVE_PTHREAD_NUM_PROCESSORS_NP */

#ifdef HAVE__SC_NPROCESSORS_ONLN
  number_threads = sysconf(_SC_NPROCESSORS_ONLN);
#else /* HAVE__SC_NPROCESSORS_ONLN */

#ifdef HAVE__SC_NPROC_ONLN
  number_threads = sysconf(_SC_NPROC_ONLN);
#endif /* HAVE__SC_NPROC_ONLN */
#endif /* HAVE__SC_NPROCESSORS_ONLN */
#endif /* HAVE_PTHREAD_NUM_PROCESSORS_NP */

  if (n_threads >= 1)
    number_threads = n_threads;

  /**
   * get env
   */

  if ((env = getenv(N_THREADS_ENV)) != NULL)
    number_threads = atoi(env);

  /**
   * perform some checks on the number of threads
   */

  if (number_threads < 1) {
    fprintf(stderr, "Number of threads has to be >= 1."
                    "Resetting to %d\n",
            1);
    number_threads = 1;
  }

  /**
   * save max_number_threads
   */

  max_number_threads = number_threads;

  /**
   * some init required by different thread implementations
   */

  printf("configuring %d max number threads\n", max_number_threads);

#ifdef HAVE_PTHREAD_LAUNCH_POLICY_NP
#ifndef __QC_AIX /* Doesn't work on AIX now?! */

  /**
   * set default launch policy to PTHREAD_POLICY_FILL_NP
   */

  if (pthread_launch_policy_np(PTHREAD_POLICY_FILL_NP, NULL,
                               PTHREAD_SELFTID_NP))
    fprintf(stderr, "pthread_launch_policy_np() failed\n");

#endif /* !__QC_AIX */
#endif /* HAVE_PTHREAD_LAUNCH_POLICY_NP */

#if defined(__QC_SGI)
  create_locality_domain(number_threads);
#endif /* sgi */

#ifdef HAVE_PTHREAD_SETCONCURRENCY
  pthread_setconcurrency(number_threads);
#endif /* HAVE_PTHREAD_SETCONCURRENCY */

  return;
}

/**
 * get max_number_threads
 */

int get_max_number_threads(void) { return (max_number_threads); }

/**
 * get number of threads
 */

int get_number_threads(void) { return (number_threads); }

/**
 * set number threads
 * **** NOT MT SAFE ****
 */

void set_number_threads(const int n_number_threads) {

  number_threads = n_number_threads;

  return;
}

/**
 * get share of data
 */

void get_share(const int my_id, const int number_threads, const int number_data,
               int *number_data_thread, int *number_data_start,
               int *number_data_end) {
  /**
   * check if number_data >= number_threads; if not signal en error
   * and abort
   */

  assert(number_data >= number_threads);
  assert(my_id < number_threads);

  *number_data_thread = number_data / number_threads;

  *number_data_start = my_id * (*number_data_thread);

  /**
   * last thread
   */
  if (my_id == number_threads - 1) {
    *number_data_end = number_data - 1;
    *number_data_thread = *number_data_end - *number_data_start + 1;
  } else
    *number_data_end = *number_data_start + *number_data_thread - 1;

  return;
}
