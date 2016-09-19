/**
 * $Id: threads.h,v 1.2 2003/08/14 20:59:29 knap Exp $
 *
 * $Log: threads.h,v $
 * Revision 1.2  2003/08/14 20:59:29  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.5  2002/12/06 01:49:38  fago
 * Significant fixes to automake build especially on AIX and HPUX.
 *
 * Revision 1.4  2001/10/03 22:29:49  knap
 * Checked in the latest changes.
 *
 * Revision 1.3  1999/09/02 23:15:42  knap
 * Corrected monitor to do MxN mapping.
 *
 * Revision 1.2  1999/03/22 21:04:42  knap
 * Added get_share() proto.
 *
 * Revision 1.1.1.1  1999/02/12 18:55:15  knap
 * Importing initial version
 *
 */
#ifndef _THREADS_H_
#define _THREADS_H_

#if defined(__cplusplus)
namespace quasicontinuum {
extern "C" {
#endif /* __cplusplus */

extern void qc_thread_init(int n_threads);
extern int get_max_number_threads(void);
extern int get_number_threads(void);
extern void set_number_threads(const int n_number_threads);
extern void get_share(const int my_id, const int number_threads,
                      const int number_data, int *number_data_thread,
                      int *number_data_start, int *number_data_end);

#if defined(__cplusplus)
}
}
#endif /* __cplusplus */

#endif /* _THREADS_H_ */
