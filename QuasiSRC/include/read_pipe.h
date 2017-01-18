/**
 * $Id: read_pipe.h,v 1.2 2003/08/14 20:59:29 knap Exp $
 *
 * $Log: read_pipe.h,v $
 * Revision 1.2  2003/08/14 20:59:29  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.3  2002/03/08 00:06:50  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.2  1999/12/02 01:34:32  knap
 * Added read_pipe() and write_pipe() protos.
 *
 * Revision 1.1  1999/07/30 20:01:11  knap
 * Initial revision.
 *
 */
#ifndef _READ_PIPE_H_
#define _READ_PIPE_H_

#include <stdio.h>

#if defined(__cplusplus)
namespace quasicontinuum {
extern "C" {
#endif /* __cplusplus */

extern FILE *open_read_pipe(const char *data_file_name);
extern FILE *open_write_pipe(const char *data_file_name);
extern FILE *open_data_file(const char *data_file_name, const char *type);

#if defined(__cplusplus)
}
}
#endif /* __cplusplus */

#endif /* _READ_PIPE_H_ */
