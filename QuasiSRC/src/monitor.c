/**
 * $Log: monitor.c,v $
 * Revision 1.9  2002/03/07 23:53:06  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.8  2001/10/03 22:29:51  knap
 * Checked in the latest changes.
 *
 * Revision 1.7  2000/07/19 18:52:02  knap
 * Corrected thread stack size increase.
 *
 * Revision 1.6  2000/07/16 02:31:52  knap
 * Somewhat tweaked code for increasing threads stack.
 *
 * Revision 1.5  2000/05/24 18:35:37  knap
 * Ported to Linux and OSF.
 *
 * Revision 1.4  2000/04/20 01:09:02  knap
 * Ported to HP-UX.
 *
 * Revision 1.3  2000/02/24 01:12:17  knap
 * Various corrections aimed at improving scalability.
 *
 * Revision 1.2  2000/01/15 17:33:44  knap
 * Range search routines make lots of recursive function calling. This may
 * cause some problems for thread stacks (overruns etc.). Added calls
 * to obtain and change default thread stack size.
 *
 * Revision 1.1  2000/01/07 00:26:55  knap
 * Moved to src.
 *
 * Revision 1.10  1999/12/16 21:42:04  knap
 * Added second order quad and simple clusters.
 *
 * Revision 1.9  1999/09/02 23:48:37  knap
 * Corrected all calls to thread_monitor() to get new interface.
 *
 * Revision 1.8  1999/07/28 00:07:51  knap
 * Linted.
 *
 * Revision 1.7  1999/07/23 22:40:19  knap
 * _REENTRANT caused some complaints on osf.
 *
 * Revision 1.6  1999/07/22 20:58:10  knap
 * malloc.h include inserted.
 *
 * Revision 1.5  1999/04/11 19:26:26  knap
 * Added pthread support for OSF1.
 *
 * Revision 1.4  1999/04/07 22:50:02  knap
 * Linted.
 *
 * Revision 1.3  1999/04/01 21:34:32  knap
 * Calls to free() changed to free_e().
 *
 * Revision 1.2  1999/03/31 20:00:24  knap
 * Minor spelling errors.
 *
 * Revision 1.1  1999/03/30 19:54:01  knap
 * Initial rev.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

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
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_SCHED_H
#include <sched.h>
#endif /* HAVE_SCHED_H */

#if defined(__QC_SGI) || defined(_HAVE_NUMA_H)
#include <numa.h>
#include "numa.h"
#endif /* sgi || linux */

#include "threads.h"
#include "Error.h"
#include "monitor.h"
#include "DataTypes.h"
#include "C_Interface.h"   // for free_e()

#if !defined(lint)
static char rcsid[]="$Id: monitor.c,v 1.9 2002/03/07 23:53:06 knap Exp $";
#endif /* !lint */


/**
 * simple monitor for spreading loops over multiple threads;
 * takes a function and calls it by multiple threads; 
 */

struct thread_id_t {
  pthread_key_t tid;
  pthread_key_t id;
};

struct monitor_t{
  void             *(*function)(void *);
  void               *arg;
  pthread_key_t       id_key;
  enum {TID, ID}      id_type;
  pthread_mutex_t     start_lock;
  pthread_mutex_t     end_lock;
  pthread_cond_t      start_cond;
  pthread_cond_t      end_cond;
  int                 work;
  int                 queue;
  int                 runners;
  int                 init;
  int                *data_taken;
};

static struct monitor_t  monitor;

/**
 * get my_id
 */

int 
get_my_tid(void)
{
  struct thread_id_t *thread_id;

  thread_id = pthread_getspecific(monitor.id_key);
  return(thread_id->id);

}

/**
 * compute id of a thread; try to maintain close affinity
 * with thread ids; not MT-safe.
 */

static int
compute_thread_id(int tid,
		  int number_threads)
{

  int id   = tid;
  int flag = 1;
  int i_left;
  int i_right;

  /**
   * data_taken extends from 0 and ends at number_threads-1
   */
  
  if (id >= number_threads) id = number_threads - 1;

  i_left  = id;
  i_right = id;

  /**
   * shortcircuit here if data_taken is avaiable
   */

  if (!monitor.data_taken[id]) {
    monitor.data_taken[id] = 1;
    return(id);
  }

  /**
   * id data is gone try locate the closest one
   */

  while(flag) {

    i_left--;
    i_right++;

    /**
     * check i_left
     */

    if (i_left >= 0 && !monitor.data_taken[i_left]) {
      monitor.data_taken[i_left] = 1;
      return(i_left);
    }

    /*
     * check right
     */

    if (i_right < number_threads && !monitor.data_taken[i_right]) {
      monitor.data_taken[i_right] = 1;
      return(i_right);
    }
  
  }
  
  /* NOTREACHED */
  return(-1);
  
}


/**
 * monitor thread
 */

static void * /*ARGSUSED 0*/
monitor_thread(void *arg)
{

  static pthread_mutex_t tid_lock = PTHREAD_MUTEX_INITIALIZER;
  static int tid = 0;

  struct thread_id_t thread_id; 

  /**
   * set tid
   */

  pthread_mutex_lock(&tid_lock);

  thread_id.tid = tid++;
  pthread_setspecific(monitor.id_key, (void *) &thread_id);
  
  pthread_mutex_unlock(&tid_lock);

#if defined(__QC_SGI)
  /**
   * attatch thread to cpu;
   */

  /* bind_thread(thread_id.tid); */

#endif /* sgi */

#if defined(HAVE_NUMA_H)
  /**
   * attatch thread to cpu;
   */
  bind_thread(thread_id.tid); 

#endif /* numa */

  /**
   * loop while waiting for work
   */
  
  for(;;){
    /**
     * wait here for work
     */

    pthread_mutex_lock(&(monitor.start_lock));
    
    while (monitor.work == 0)
      pthread_cond_wait(&(monitor.start_cond), &(monitor.start_lock));

    monitor.work--;   
    thread_id.id = compute_thread_id(thread_id.tid, get_number_threads());
    /* thread_id.id = monitor.queue++; */

    pthread_mutex_unlock(&(monitor.start_lock));    
    
    /**
     * call per thread function
     */

    monitor.function(monitor.arg);

    /**
     * signal parent thread
     */

    pthread_mutex_lock(&(monitor.end_lock));

    monitor.runners--;
    pthread_cond_signal(&(monitor.end_cond));

    pthread_mutex_unlock(&(monitor.end_lock));

  }
  
  /*NOTREACHED*/
  return((void *) NULL);
}




void 
thread_monitor(void *(*function)(void *), 
	       void   *arg,
	       int     number_threads)
{

  static struct thread_id_t thread_id;

  /**
   * init
   */

  if( monitor.init == 0 ){

    pthread_t      *thread_ids;
    pthread_attr_t  attr;

    int             i_thread;
    int             max_number_threads=get_max_number_threads();

    /**
     * check availability of numa
     */
#if defined(HAVE_NUMA_H)
    if(numa_available() == -1)
      printf("NUMA functions not available");
#endif /* numa */

    /**
     * init data
     */

    monitor.function = function;
    monitor.arg      = arg;
    monitor.work     = 0;
    monitor.queue    = 0;
    monitor.runners  = 0;
    monitor.init     = 0;
    
    /**
     * allocate thread specific data
     */

    pthread_key_create(&monitor.id_key, NULL);
    monitor.id_type = TID;

    /**
     * if some functions that call thread functions themself (eg.
     * get_my_tid()) from the main thread, set the per thread data
     * for the main thread
     */

    thread_id.tid = 0;
    pthread_setspecific(monitor.id_key, (void *) &thread_id);

    /**
     * initialize locks and cond vars
     */

    pthread_mutex_init(&(monitor.start_lock), NULL);
    pthread_mutex_init(&(monitor.end_lock), NULL);

    pthread_cond_init(&(monitor.start_cond), NULL);
    pthread_cond_init(&(monitor.end_cond), NULL);

    /**
     * initialize space for data_taken
     */

    if ((monitor.data_taken = malloc(max_number_threads*sizeof(int))) == NULL){
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    /**
     * allocate space for thread ids
     */

    thread_ids=(pthread_t *) malloc(max_number_threads*sizeof(pthread_t));
    if( thread_ids == (pthread_t *) NULL ){
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    /*
     * initialize thread attributes to detached, system
     */

    pthread_attr_init( &attr );
    pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_DETACHED);

#ifdef HAVE_PTHREAD_SCOPE_SYSTEM
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );
#endif /* HAVE_PTHREAD_SCOPE_SYSTEM */

    /**
     * in case thread stack is too small use these to alter it
     */
#if 0
    {
      size_t size;

      pthread_attr_getstacksize(&attr, &size);
      pthread_attr_setstacksize(&attr, 10*size);
      pthread_attr_getstacksize(&attr, &size);
      printf("stack size is %d\n", size);

#ifdef HAVE_PTHREAD_ATTR_GETGUARDSIZE
      pthread_attr_getguardsize(&attr, &size);
      pthread_attr_setguardsize(&attr, 10*size);
      pthread_attr_getguardsize(&attr, &size);
      printf("stack guard size is %d\n", size);
#endif /* HAVE_PTHREAD_ATTR_GETGUARDSIZE */
    }
#endif
    /**/

    /**
     * spawn threads
     */

    for(i_thread=0; i_thread < max_number_threads; i_thread++)
      if( pthread_create(&(thread_ids[i_thread]), &attr, monitor_thread,
    			 (void *) NULL) ){
    	ERROR("pthread_create()");
    	exit(EXIT_FAILURE);
      }

    /* 
     * clean up
     */

    pthread_attr_destroy(&attr);
    (void) free_e(thread_ids);

    monitor.init=1;
  }

  /**
   * init data
   */

  set_number_threads( number_threads );
  
  monitor.function=function;
  monitor.arg=arg;

  pthread_mutex_lock( &(monitor.start_lock) );
  
  monitor.work=number_threads;
  monitor.runners=number_threads;
  monitor.queue=0;
  memset(monitor.data_taken, '\0', number_threads*sizeof(int));

  pthread_cond_broadcast( &(monitor.start_cond) );

  pthread_mutex_unlock( &(monitor.start_lock) );

#ifdef HAVE_SCHED_H
  sched_yield();
#endif /* HAVE_SCHED_H */

  /**
   * wait for threads to finish processing
   */

  pthread_mutex_lock( &(monitor.end_lock) );

  while( monitor.runners != 0 )
    pthread_cond_wait(&(monitor.end_cond), &(monitor.end_lock));

  pthread_mutex_unlock( &(monitor.end_lock) );  


  return;
}
