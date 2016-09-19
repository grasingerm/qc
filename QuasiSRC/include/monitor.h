/*
 * $Id: monitor.h,v 1.2 2003/08/14 20:59:28 knap Exp $
 *
 * $Log: monitor.h,v $
 * Revision 1.2  2003/08/14 20:59:28  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.4  2002/03/08 00:06:50  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.3  1999/09/02 23:15:41  knap
 * Corrected monitor to do MxN mapping.
 *
 * Revision 1.2  1999/07/22 18:52:14  knap
 * Added enum mt_version definition.
 *
 * Revision 1.1  1999/03/30 19:54:17  knap
 * Initial rev.
 *
 */
#ifndef _MONITOR_H_
#define _MONITOR_H_

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#if defined(__cplusplus)
namespace quasicontinuum {
extern "C" {
#endif /* __cplusplus */

extern int get_my_tid(void);
extern void thread_monitor(void *(*function)(void *), void *arg,
                           int number_threads);

#if defined(__cplusplus)
}
}
#endif /* __cplusplus */

#endif /* _MONITOR_H_ */
