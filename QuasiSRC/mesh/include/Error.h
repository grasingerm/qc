// Error.h

#ifndef ERROR_H
#define ERROR_H

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found.
#endif /* HAVE_SYS_TYPES_H */

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#else
#error signal.h not found.
#endif /* HAVE_SIGNAL_H */

#ifdef STDC_HEADERS
#include <stdio.h>
#include <string.h>
#else
#error Standard C library headers not found.
#endif /* STDC_HEADERS */

#ifdef HAVE_SCHED_H
#include <sched.h>
#else
#error sched.h not found.
#endif /* HAVE_SCHED_H */

static int ERROR_SPIN = 0;

#define ERROR(mesg) do{fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__);       \
                       perror(mesg); } while(ERROR_SPIN)


#define ERROR2(mesg,error) do{fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__); \
			      fprintf(stderr, "%s: %s\n", mesg,                \
		                      strerror(error)); } while(ERROR_SPIN)


#define PERROR_MALLOC() { (void) fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__);\
                         perror( "malloc()" );                                \
                      }
#endif /* ERROR_H */

