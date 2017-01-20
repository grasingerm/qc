/*
 * $Log: util.c,v $
 * Revision 1.3  2002/12/06 01:50:06  fago
 * Significant fixes to automake build especially on AIX and HPUX.
 *
 * Revision 1.2  2002/03/07 23:25:38  knap
 * Changed the build procedure to automake/autoconf.
 *
 * Revision 1.1  2000/01/03 22:48:00  knap
 * Moved to ./src.
 *
 * Revision 1.1  1999/12/18 17:50:26  knap
 * Initial rev.
 *
 * Revision 1.1  1999/03/29 21:41:50  knap
 * Initial rev.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <unistd.h>
#include <errno.h>
#include "util.h"

#if !defined(lint)
static char rcsid[]="$Id: util.c,v 1.3 2002/12/06 01:50:06 fago Exp $";
#endif /* !lint */

/**
 * function top read/write exactly n bytes
 */


ssize_t   /* Read "n" bytes from a descriptor. */
readn(int fd, void *vptr, size_t n)
{
  size_t  nleft;
  ssize_t nread;
  char    *ptr;
  
  ptr = vptr;
  nleft = n;
  while (nleft > 0) {
    if ( (nread = read(fd, ptr, nleft)) < 0) {
      if (errno == EINTR)
	nread = 0;              /* and call read() again
				 */
      else
	return(-1);
    } else if (nread == 0)
      break;                          /* EOF */
    
    nleft -= nread;
    ptr   += nread;
  }
  return(n - nleft);              /* return >= 0 */
}

ssize_t  /* Write "n" bytes to a descriptor. */
writen(int fd, const void *vptr, size_t n)
{
  size_t          nleft;
  ssize_t         nwritten;
  const char      *ptr;
  
  ptr = vptr;
  nleft = n;
  while (nleft > 0) {
                if ( (nwritten = write(fd, ptr, nleft)) <= 0) {
		  if (errno == EINTR)
		    nwritten = 0;           /* and call write() agai
					       n */
		  else
		    return(-1);                     /* error */
                }
		
                nleft -= nwritten;
                ptr   += nwritten;
  }
  return(n);
}
/* end writen */


