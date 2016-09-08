/**
 * $Log: triangulation_client.c,v $
 * Revision 1.20  2002/12/09 17:10:41  fago
 * Minor build fix.
 *
 * Revision 1.19  2002/12/06 01:50:05  fago
 * Significant fixes to automake build especially on AIX and HPUX.
 *
 * Revision 1.18  2002/12/04 17:46:00  fago
 * Ported to automake 1.7.1, autoconf 2.56. perror -> ERROR. Fixed build bugs,
 * especially on HPUX.
 *
 * Revision 1.17  2002/03/07 23:25:38  knap
 * Changed the build procedure to automake/autoconf.
 *
 * Revision 1.16  2001/10/03 22:29:53  knap
 * Checked in the latest changes.
 *
 * Revision 1.15  2000/07/25 18:21:48  knap
 * Changed how memory is allocated on ccNUMA systems.
 *
 * Revision 1.14  2000/07/16 17:20:39  knap
 * Ported to Linux.
 *
 * Revision 1.13  2000/05/24 18:35:45  knap
 * Ported to Linux and OSF.
 *
 * Revision 1.12  2000/05/23 21:11:04  knap
 * Corrected some OSF1 errors.
 *
 * Revision 1.11  2000/05/18 18:32:58  knap
 * Some header added to make SUN compiler happy.
 *
 * Revision 1.10  2000/05/18 17:45:49  knap
 * Both client and server code changed to use a local socket instead of pipes.
 *
 * Revision 1.9  2000/05/10 04:41:12  knap
 * Removed ALL remaining references to atoms.
 *
 * Revision 1.8  2000/04/26 20:52:49  knap
 * Triangulation is based now on lattice coordinates (integers).
 *
 * Revision 1.7  2000/04/20 01:09:14  knap
 * Ported to HP-UX.
 *
 * Revision 1.6  2000/04/12 22:07:28  knap
 * Removed references to atoms.
 *
 * Revision 1.5  2000/02/17 20:53:53  knap
 * Corrected to place memeory more efficiently.
 *
 * Revision 1.4  2000/01/24 17:42:19  knap
 * Added exit handler for server shutdown.
 *
 * Revision 1.3  2000/01/12 22:15:49  knap
 * Client stub obtains path to the server from QUASI_HOME env.
 *
 * Revision 1.2  2000/01/07 00:05:41  knap
 * Corrected arguments to exec().
 *
 * Revision 1.1  2000/01/03 22:48:00  knap
 * Moved to ./src.
 *
 * Revision 1.3  1999/12/20 22:18:36  knap
 * Added output statement in check_for_duplicate_nodes().
 *
 * Revision 1.2  1999/12/20 22:15:46  knap
 * Added duplicate nodes check.
 *
 * Revision 1.1  1999/12/18 17:50:26  knap
 * Initial rev.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#if defined(__QC_SGI)
#include <pthread.h>
#include <invent.h>
#include <limits.h>
#endif /* sgi */

#if defined(__QC_HPUX)
#include <pthread.h>
#include <string.h>
#endif /* hpux */

#if defined(__QC_SUN) || defined(__QC_SGI) || defined(__QC_HPUX) || defined(__QC_LINUX) || defined(__QC_AIX)
#include <errno.h>
#endif /* sun || sgi || hpux || __linux__ || AIX */

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#else
#error signal.h not found.
#endif /* HAVE_SIGNAL_H */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#error limits.h not found
#endif /* HAVE_LIMITS_H */

#ifdef HAVE_UNISTD_H 
#include <unistd.h>
#else
#error unistd.h not found
#endif /* HAVE UNISTD_H */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found
#endif /* HAVE_SYS_TYPES_H */

#ifdef HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>
#else
#error sys/utsname.h not found.
#endif /* HAVE_SYS_UTSNAME_H */

#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#else
#error sys/socket.h not found.
#endif /* HAVE_SYS_SOCKET_H */

#ifdef HAVE_SYS_UN_H
#include <sys/un.h>
#else
#error sys/un.h not found.
#endif /* HAVE_SYS_UN_H */

#ifdef HAVE_ASSERT_H
#include <assert.h>
#else
#error assert.h not found.
#endif /* HAVE_ASSERT_H */

#include "Error.h"
#include "DataTypes.h"
#include "triangulation_client.h"
#include "triangulation_server.h"

#if !defined(lint)
static const char rcsid[]="$Id: triangulation_client.c,v 1.20 2002/12/09 17:10:41 fago Exp $";
#endif /* !lint */

#define TRIANGULATION_SERVER "triang_server"

static int server_pid;

/**
 * handler to kill server at exit()
 */

static void
kill_server_handler(void)
{

  kill(server_pid, SIGINT);

}

/**
 * create system extension to TRIANGULATION_SERVER
 */

static void
create_triang_server_file_name(char *filename)
{

  char   *triang_server_path;

  /**
   * get path to triang server from enviroment and if present add it;
   * otherwise just create it
   */
  
  if( (triang_server_path=getenv("QUASI_HOME")) != NULL )
    sprintf(filename, "%s/%s", triang_server_path, TRIANGULATION_SERVER);
  else  
    sprintf(filename, "%s", TRIANGULATION_SERVER);
  
  return;
}


/**
 * setup triangulation server returning two fds one for reading and
 * one for writing
 */

static pid_t
setup_server(void)
{

  char triang_server_filename[PATH_MAX+1];
  pid_t pid;

  /**
   * fork here
   */

  
  pid=fork();

  switch( pid ){

  case -1:
    ERROR("fork()");
    exit(EXIT_FAILURE);
    break;

    /*
     * child
     */
  case 0:
    
    /**
     * compose full path to the triangulation server
     */

    create_triang_server_file_name(triang_server_filename);

    /**
     * exec server
     */
     
    execl( triang_server_filename, triang_server_filename, NULL );
    ERROR(triang_server_filename);
    kill( getppid(), SIGKILL ); /* parent does not wait() */
    exit(EXIT_FAILURE);

    break;

    /**
     * parent
     */

  default:

    /**
     * install handler to kill server at exit
     */

    if(atexit(kill_server_handler)){
      ERROR("atexit()");
      exit(EXIT_FAILURE);
    }

    break;

  }
  

  return( pid );

}

/**
 * connect to triangulation server using local socket
 */

static int
connect_server(pid_t pid)
{
  
  struct sockaddr_un sin;
  
  extern int errno;
  int sd;

  char server_name[_POSIX_PATH_MAX+1];

  /**
   * create server_name
   */

  sprintf(server_name, "%s.%d", TRIANGULATION_SOCKET_PREFIX, pid);

  /**
   *  create socket
   */

#if defined(__QC_HPUX) || defined(__QC_SUN) || defined(__QC_OSF) || defined(__QC_AIX)
  sin.sun_family = AF_UNIX;
#else
  sin.sun_family = AF_LOCAL;
#endif /* hpux */

  strcpy(sin.sun_path, server_name);

  /**
   * loop waiting for the server to start
   */

  while(1){

#if defined(__QC_HPUX) || defined(__QC_SUN) || defined(__QC_OSF) || defined(__QC_AIX)
    if((sd=socket(AF_UNIX, SOCK_STREAM, 0)) < 0){
#else
    if((sd=socket(AF_LOCAL, SOCK_STREAM, 0)) < 0){
#endif /* hpux */

      ERROR(sin.sun_path);
      exit(EXIT_FAILURE);
    }
  
    /**
     * connect to server on a socket
     */

    if(connect(sd, (const struct sockaddr *) &sin, sizeof(sin)) < 0)
      switch(errno){
	
      case ENOENT:
	close(sd);
	break;

      case ECONNREFUSED:
	close(sd);
	break;
	
      default:
	fprintf(stderr, "connect():");
	ERROR(sin.sun_path);
	unlink(server_name);
	exit(EXIT_FAILURE);
	break;
	
      }
    else break;
    
  }

  /**
   * now that we are connected unlink the socket
   */

  if (unlink(server_name) == -1) {
    ERROR(server_name);
    exit(EXIT_FAILURE);
  }

  return(sd);

}


/**
 * Create elements
 */

void 
create_elements_external(struct node_list_t    *P_node_list, 
        struct mesh_data_t       *P_mesh_data)
{

  static pid_t pid;

  static int init_done=0;
  static int fd;
   
  int i_node;
  int i_tet;

  /**
   * initialize; start server and connect to it
   */

  if( init_done == 0 ){

    server_pid = pid = setup_server();
    fd = connect_server(pid);
    init_done = 1;

  }

  printf("child %d\n", pid);

  /**
   * write nodal data to w_file
   */

  if( writen( fd, &P_node_list->number_nodes, sizeof(int) ) < 0 ){
    ERROR( "writen()" );
    exit( EXIT_FAILURE );
  }

  for( i_node=0; i_node < P_node_list->number_nodes; i_node++ ) {
    struct node_t *P_node = P_node_list->nodes[i_node];
    double lattice_coord[3];

    lattice_coord[0] = P_node->l[0];
    lattice_coord[1] = P_node->l[1];
    lattice_coord[2] = P_node->l[2];
    
    /* if( writen( wfd, P_node_list->nodes[i_node]->initial_position, */
    if( writen( fd, &(lattice_coord[0]), sizeof(double[3]) ) < 0 ){
      ERROR( "writen()" );
      exit( EXIT_FAILURE);
    }
  }

  
  /*
   * read number of elements
   */

  if( readn( fd, &P_mesh_data->number_tets, sizeof(int) ) < 0 ){
    ERROR( "readn()" );
    exit( EXIT_SUCCESS );
  }

  // allocate space in P_mesh_data->tets so that it can hold above number
  // of tets
  P_mesh_data->tets = 
    (int (*)[4]) malloc(P_mesh_data->number_tets*sizeof(int[4]));

  for( i_tet=0; i_tet < P_mesh_data->number_tets; i_tet++ ){

    int    tet[4];

    if( readn( fd, tet, sizeof(int[4]) ) < 0 ){
      ERROR( "readn()" );
      exit( EXIT_FAILURE );
    }

    // substract 1 from each vertices of tetra so that we have currect node
    // number to write create elements.
    (tet[0])--;
    (tet[1])--;
    (tet[2])--;
    (tet[3])--;

    /**
     * fill in P_mesh_data
     */

    P_mesh_data->tets[i_tet][0] = tet[0];
    P_mesh_data->tets[i_tet][1] = tet[1];
    P_mesh_data->tets[i_tet][2] = tet[2];
    P_mesh_data->tets[i_tet][3] = tet[3];
  }

  /**
   * read mesh quality data
   */

  if( readn( fd, &P_mesh_data->eta_min, sizeof(double) ) < 0 ){
    ERROR( "readn()" );
    exit( EXIT_FAILURE );
  }

  if( readn( fd, &P_mesh_data->rho_min, sizeof(double) ) < 0 ){
    ERROR( "readn()" );
    exit( EXIT_FAILURE );
  }

  if( readn( fd, &P_mesh_data->sigma_min, sizeof(double) ) < 0 ){
    ERROR( "readn()" );
    exit( EXIT_FAILURE );
  }

  printf( "eta=%e rho=%e sigma=%e\n", P_mesh_data->eta_min, 
    P_mesh_data->rho_min, P_mesh_data->sigma_min );

  /* perform_output( (struct all_node_list_t *)P_node_list, P_element_list, */
  /* 		  NULL, NULL, NULL, -900,  */
  /* 		  NODE_OUTPUT_FLAG); */

  return;
}




