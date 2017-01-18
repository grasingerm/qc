/**
 * $Log: read_pipe.c,v $
 * Revision 1.3  2002/03/07 23:53:08  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.2  2000/05/24 18:35:38  knap
 * Ported to Linux and OSF.
 *
 * Revision 1.1  2000/01/07 00:26:55  knap
 * Moved to src.
 *
 * Revision 1.4  2000/01/04 22:31:23  knap
 * Missing pthread_attr_destroy() added.
 *
 * Revision 1.3  1999/12/02 01:38:39  knap
 * O_TRUNC added to the file creation flags in open_write_pipe().
 *
 * Revision 1.2  1999/12/02 01:33:38  knap
 * Added write_pipe().
 *
 * Revision 1.1  1999/07/30 20:01:09  knap
 * Initial revision.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/wait.h>

#ifndef WEXITSTATUS
#define WEXITSTATUS(stat_val) ((unsigned)(stat_val) >> 8)
#endif
#ifndef WIFEXITED
#define WIFEXITED(stat_val) (((stat_val)&255) == 0)
#endif

#include "Error.h"
#if !defined(lint)
#include "read_pipe.h"
#endif /* !lint */

#if !defined(lint)
static const char rcsid[] =
    "$Id: read_pipe.c,v 1.3 2002/03/07 23:53:08 knap Exp $";
#endif /* !lint */

/**
 * worker detached thread that simply calls waitpid to reap a child and
 * exits
 */

static void * /*ARGSUSED0*/
    waitpid_thread(void *arg) {

  /*LINTED*/
  pid_t pid = -1;

  pid = waitpid(pid, NULL, 0);
  if ((pid == (pid_t)-1) && (errno != EINTR)) {
    ERROR("waitpid()");
    exit(EXIT_FAILURE);
  }

  return ((void *)NULL);
}

/**
 * start waitpid thread; we do not want to call waitpid outside this
 * function once data is read; we leave a detached thread that will
 * block in waitpid(). Once child exits waitpid() will return and
 * thread will exit
 */

static void start_waitpid_thread(void) {

  pthread_t thread_id;
  pthread_attr_t thread_attr;
  int pthread_error;
  /**
   * create detached thread to handle SIGCLD
   */

  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED);

  if (pthread_error = pthread_create(&thread_id, &thread_attr, waitpid_thread,
                                     (void *)NULL)) {
    ERROR2("pthread_create()", pthread_error);
    printf("%d\n", pthread_error);
    exit(EXIT_FAILURE);
  }

  pthread_attr_destroy(&thread_attr);

  return;
}

/**
 * given name of gzipped file set up a pipe to gunzip
 * input:
 *      data_file_name-name of data file
 * output:
 *      FILE *
 */

FILE *open_read_pipe(const char *data_file_name) {

  pid_t pid;

  int pipe_fd[2];

  /**
   * create pipe
   */

  if (pipe(pipe_fd) == -1) {
    ERROR("pipe()");
    exit(EXIT_FAILURE);
  }

  /**
   * fork a child
   */

  pid = fork();

  switch (pid) {

  case -1:

    ERROR("fork()");
    exit(EXIT_FAILURE);
    break;

  case 0: /* child */

#if defined(__QC_LINUX)
  {
    pid_t pid1;

    pid1 = fork();

    switch (pid1) {

    case -1:

      ERROR("fork()");
      exit(EXIT_FAILURE);
      break;

    default:

      _exit(EXIT_SUCCESS);
      break;

    case 0:

#endif /* __linux__ */

      /**
       * close reading side of the pipe;
       * some OSes use pipe_fd[0] for reading and pipe_fd[1] for writing;
       * some use bidirectional pipes; we want to be safe here
       */

      close(pipe_fd[0]);

      /**
       * attach writing side to stdout
       */

      close(STDOUT_FILENO);
      dup(pipe_fd[1]);
      close(pipe_fd[1]);

      /**
       * close stdin and stderr
       */

      close(STDIN_FILENO);
      /* close(STDERR_FILENO); */

      /**
       * exec gunzip
       */

      if (execlp("gunzip", "gunzip", "-c", data_file_name, (char *)NULL) ==
          -1) {
        ERROR("execl()");
        /* send SIGKILL to parent */
        kill(getppid(), SIGKILL);
        exit(EXIT_FAILURE);
      }

      break;
#if defined(__QC_LINUX)
    }
    break;
  }
#endif /* __linux__ */

  default: /* parent */

    /* close writing side of the pipe */

    close(pipe_fd[1]);
#if defined(__QC_LINUX)
    waitpid_thread(NULL);
#else
    start_waitpid_thread();
#endif /* __linux__ */
    return (fdopen(pipe_fd[0], "r"));
    /* NOTREACHED */
    break;
  }

  /* NOTREACHED */
  return ((FILE *)NULL);
}

/**
 * given name of gzipped file set up a pipe to gzip
 * input:
 *      data_file_name-name of data file
 * output:
 *      FILE *
 */

FILE *open_write_pipe(const char *data_file_name) {

  pid_t pid;

  int pipe_fd[2];
  int wfd;

  /**
   * create pipe
   */

  if (pipe(pipe_fd) == -1) {
    ERROR("pipe()");
    exit(EXIT_FAILURE);
  }

  /**
   * fork a child
   */

  pid = fork();

  switch (pid) {

  case -1:

    ERROR("fork()");
    exit(EXIT_FAILURE);
    break;

  case 0: /* child */

#if defined(__QC_LINUX)
  {
    pid_t pid1;

    pid1 = fork();

    switch (pid1) {

    case -1:

      ERROR("fork()");
      exit(EXIT_FAILURE);
      break;

    default:

      _exit(EXIT_SUCCESS);
      break;

    case 0:

#endif /* __linux__ */

      /**
       * close writing side of the pipe;
       * some OSes use pipe_fd[0] for reading and pipe_fd[1] for writing;
       * some use bidirectional pipes; we want to be safe here
       */

      close(pipe_fd[1]);

      /**
       * attach writing side to stdin
       */

      close(STDIN_FILENO);
      dup(pipe_fd[0]);
      close(pipe_fd[0]);

      /**
       * attach STDOUT to a file
       */

      if ((wfd = open(data_file_name, O_WRONLY | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH |
                          S_IWOTH)) == -1) {
        ERROR(data_file_name);
        exit(EXIT_FAILURE);
      }

      close(STDOUT_FILENO);
      dup(wfd);
      close(wfd);

      /**
       * exec gunzip
       */

      if (execlp("gzip", "gzip", "-c", (char *)NULL) == -1) {
        ERROR("execl()");
        /* send SIGKILL to parent */
        kill(getppid(), SIGKILL);
        exit(EXIT_FAILURE);
      }

      break;
#if defined(__QC_LINUX)
    }
    break;
  }
#endif /* __linux__ */

  default: /* parent */

    /* close reading side of the pipe */

    close(pipe_fd[0]);

#if defined(__QC_LINUX)
    waitpid_thread(NULL);
#else
    start_waitpid_thread();
#endif /* __linux__ */
    return (fdopen(pipe_fd[1], "w"));
    /* NOTREACHED */
    break;
  }

  /* NOTREACHED */
  return ((FILE *)NULL);
}

/**
 * get input stream; two things are done here:
 * 1) if file has extension gz gzipped data file is assumed
 * 2) if file does not have gz extension regular data file is assumed
 *
 * input:
 *       data_file_name-name of data file
 *       type          -same as in fopen()
 * output:
 *       FILE *
 */

FILE *open_data_file(const char *data_file_name, const char *type) {

  FILE *data_file = NULL;

  int len;

  enum { REGULAR, GZIPPED } data_file_type;

  char *str;

  /**
   * check if data_file_name ends with ".gz"
   */

  len = strlen(data_file_name);
  str = (char *)data_file_name + len - 3;

  if (strcmp(str, ".gz") == 0)
    data_file_type = GZIPPED;
  else
    data_file_type = REGULAR;

  switch (data_file_type) {

  case REGULAR:

    /**
     * check type
     */

    if (strcmp(type, "w") == 0)
      data_file = fopen(data_file_name, "w");
    if (strcmp(type, "r") == 0)
      data_file = fopen(data_file_name, "r");
    break;

  case GZIPPED:

    /**
     * check type
     */
    if (strcmp(type, "w") == 0)
      data_file = open_write_pipe(data_file_name);
    if (strcmp(type, "r") == 0)
      data_file = open_read_pipe(data_file_name);
    break;
  }

  return (data_file);
}
