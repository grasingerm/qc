/**
 * $Id: util.h,v 1.1 1999/12/18 17:50:28 knap Exp $
 *
 * $Log: util.h,v $
 * Revision 1.1  1999/12/18 17:50:28  knap
 * Initial rev.
 *
 */
#ifndef _UTIL_H_
#define _UTIL_H_

extern ssize_t readn(int fildes, void *buf, size_t nbyte);
extern ssize_t writen(int fildes, const void *buf, size_t nbyte);

#endif /* _UTIL_H_ */

