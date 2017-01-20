/**
 * $Log: triangulation_server.c,v $
 * Revision 1.16  2003/07/07 19:32:24  knap
 * The default sizes of fixed arrays used by GEOMPACK are determined by
 * compilation mode (32/64 bit).
 *
 * Revision 1.15  2002/12/06 01:50:06  fago
 * Significant fixes to automake build especially on AIX and HPUX.
 *
 * Revision 1.14  2002/12/04 17:46:00  fago
 * Ported to automake 1.7.1, autoconf 2.56. perror -> ERROR. Fixed build bugs,
 * especially on HPUX.
 *
 * Revision 1.13  2002/03/07 23:25:38  knap
 * Changed the build procedure to automake/autoconf.
 *
 * Revision 1.12  2001/10/03 22:29:53  knap
 * Checked in the latest changes.
 *
 * Revision 1.11  2000/07/19 18:54:59  knap
 * Corrected some compile errors on Solaris.
 *
 * Revision 1.10  2000/07/16 17:08:13  knap
 * Ported to 64-bit HP-UX.
 *
 * Revision 1.9  2000/05/30 23:19:20  knap
 * Added proper handling of error from GEOMPACK.
 *
 * Revision 1.8  2000/05/24 18:35:46  knap
 * Ported to Linux and OSF.
 *
 * Revision 1.7  2000/05/23 21:11:05  knap
 * Corrected some OSF1 errors.
 *
 * Revision 1.6  2000/05/18 17:45:50  knap
 * Both client and server code changed to use a local socket instead of pipes.
 *
 * Revision 1.5  2000/04/26 20:52:48  knap
 * Triangulation is based now on lattice coordinates (integers).
 *
 * Revision 1.4  2000/04/20 01:09:15  knap
 * Ported to HP-UX.
 *
 * Revision 1.3  2000/04/12 22:07:29  knap
 * Removed references to atoms.
 *
 * Revision 1.2  2000/01/07 00:04:53  knap
 * remived sleep(10) and corrected memeory leak.
 *
 * Revision 1.1  2000/01/03 22:48:00  knap
 * Moved to ./src.
 *
 * Revision 1.1  1999/12/18 17:50:26  knap
 * Initial rev.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

#include "triangulation_server.h"
#include "util.h"
#include "Error.h"

#if !defined(lint)
static const char rcsid[]="$Id: triangulation_server.c,v 1.16 2003/07/07 19:32:24 knap Exp $";
#endif /* !lint */

#ifdef __cplusplus
  extern "C" {
#endif

void GETERR_F77( int *auxerr );
double SANGMN_F77( double A[3], double B[3], double C[3], double D[3],
			    double SANG[4] );
double RADRTH_F77( double A[3], double B[3], double C[3], 
        double D[3] );
double EMNRTH_F77( double A[3], double B[3], double C[3], 
        double D[3] );
int PRIME_F77( int *K );
void INITCB_F77( double *tol );
void DTRIW3_F77( int *NPT, int *SIZHT, int *MAXBF, int *MAXFC,
			  double (*VCL)[3], int *VM, int *NBF, int *NFC,
			  int *NFACE, int *NTETRA, int (*BF)[3], 
			  int (*FC)[7], int *HT );
void TETLST_F77( int *nfc, int *vm, int (*fc)[7], int *number_tetra,
			  int (*tetra)[4] ); 
void IMPTR3_F77( long *BNDCON, long *POSTLT, int *CRIT, int *NPT, 
			  int *SIZHT, int *MAXFC, double (*VCL)[3], int *VM,
			  int *NFC, int *NTETRA, int (*BF)[3], int (*FC)[7],
			  int *HT, int *NFACE );
void IMPTRF_F77(long *BNDCON, int *CRIT, int *NPT, int *SIZHT,
			 int *MAXFC, double (*VCL)[3], int *VM, int *NFC,
			 int *NTETRA, int *HDAVFC, int (*BF)[3], int (*FC)[7],
			 int *HT);

#ifdef __cplusplus
  }
#endif

#if !defined(MAX) && !defined(MIN)
#  define MIN(a, b) ((a) < (b) ? (a) : (b))
#  define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif /* !MAX && !MIN */

/**
 * if the triangulation server calls exit() kill the parent since 
 * the server exits on error ONLY
 */

static void
kill_parent(void)
{

  pid_t ppid;

  if ((ppid = getppid()) == (pid_t ) -1) {
    ERROR("getppid()");
    return;
  }

  if (kill(ppid, SIGINT) == -1) {
    ERROR("kill()");
    return;
  }
    
  return;
}

/**
 * create IPC server
 */

static int
create_server(void)
{

  struct sockaddr_un sin;
  int sd;
  char server_name[_POSIX_PATH_MAX+1];

  /**
   * create server_name
   */

  sprintf(server_name, "%s.%d", TRIANGULATION_SOCKET_PREFIX, getpid());

  /**
   * clear stalled address
   */

  if(unlink(server_name) < 0)
    if(errno != ENOENT){
      ERROR("unlink()");
      exit(EXIT_FAILURE);
    }

  /**
   * create socket
   */

#if defined(__QC_HPUX) || defined(__QC_SUN) || defined(__QC_OSF) || defined(__QC_AIX)
  sin.sun_family = AF_UNIX;
#else
  sin.sun_family = AF_LOCAL;
#endif /* hpux */
  strcpy(sin.sun_path, server_name);

#if defined(__QC_HPUX) || defined(__QC_SUN) || defined(__QC_OSF) || defined(__QC_AIX)
  if((sd=socket(AF_UNIX, SOCK_STREAM, 0)) < 0){
#else
  if((sd=socket(AF_LOCAL, SOCK_STREAM, 0)) < 0){
#endif /* hpux */

    ERROR(sin.sun_path);
    exit(EXIT_FAILURE);
  }

  /**
   * bind socket to address
   */

  if(bind(sd, (const struct sockaddr *) &sin, sizeof(sin)) < 0){
    ERROR(sin.sun_path);
    exit(EXIT_SUCCESS);
  }

  /**
   * set up listen queue with backlog of 5
   */

  if(listen(sd, 5) < 0){
    ERROR("listen()");
    exit(EXIT_SUCCESS);
  }
  

  return(sd);

}


/**
 * extract error
 */

static int 
extract_error( void )
{
  extern void GETERR_F77( int *auxerr );

  int auxerr;

  GETERR_F77( &auxerr );

  return( auxerr );
}

/**
 * compute quality coeffs
 */

static void
check_mesh_quality( double *min_sigma, double *min_rho, double *min_eta,
                    double (*nodes)[3], int number_nodes, int (*tetra)[4],
                    int number_tetra )
{
  extern double SANGMN_F77( double A[3], double B[3], double C[3], double D[3],
			    double SANG[4] );
  extern double RADRTH_F77( double A[3], double B[3], double C[3], 
			    double D[3] );
  extern double EMNRTH_F77( double A[3], double B[3], double C[3], 
			    double D[3] );
  
  double sang[4];
  double sigma;
  double rho;
  double eta;

  int i_elem;

  /**
   * initialize
   */

  *min_sigma=*min_rho=*min_eta=2.0;

  /**
   * loop over all tetra
   */

  for( i_elem=0; i_elem < number_tetra; i_elem++ ){

    int a=tetra[i_elem][0]-1;
    int b=tetra[i_elem][1]-1;
    int c=tetra[i_elem][2]-1;
    int d=tetra[i_elem][3]-1;

    sigma = SANGMN_F77( nodes[a], nodes[b], nodes[c], nodes[d], sang );
    rho   = RADRTH_F77( nodes[a], nodes[b], nodes[c], nodes[d] );
    eta   = EMNRTH_F77( nodes[a], nodes[b], nodes[c], nodes[d] );

    *min_sigma=MIN( *min_sigma, sigma );
    *min_rho=MIN( *min_rho, rho );
    *min_eta=MIN( *min_eta, eta );
    
  }

  *min_sigma=*min_sigma*1.5*sqrt(6.0);

  return;

}

/**
 * use GEOMPACK to generate and improve triangulation
 */

#if defined(WITH_64_BIT_ENABLED)

#define MAX_BF 2000000
/*Jason: increase MAX_FC, don't know if this is an issue*/
/*#define MAX_FC 16000000*/
#define MAX_FC 1600000000
#define MAXHT  8011
#define MAXTH  8000000
#define MAXVC  2000000
#define TOL    1.0e-9

#else

#define MAX_BF 200000
/*Jason: increase MAX_FC, don't know if this is an issue*/
/*#define MAX_FC 1600000*/
#define MAX_FC 16000000
#define MAXHT  8011
#define MAXTH  800000
#define MAXVC  200000
#define TOL    1.0e-9

#endif /* WITH_64_BIT_ENABLED */


static int (*generate_tets( int     *number_tetra, 
			    double  *min_eta,
			    double  *min_rho,
			    double  *min_sigma,
			    double (*nodes)[3],
			    int      number_nodes ))[4]
{

  /**
   * prototypes
   */
  extern int PRIME_F77( int *K );
  extern void INITCB_F77( double *tol );
  extern void DTRIW3_F77( int *NPT, int *SIZHT, int *MAXBF, int *MAXFC,
			  double (*VCL)[3], int *VM, int *NBF, int *NFC,
			  int *NFACE, int *NTETRA, int (*BF)[3], 
			  int (*FC)[7], int *HT );
  extern void TETLST_F77( int *nfc, int *vm, int (*fc)[7], int *number_tetra,
			  int (*tetra)[4] ); 
  extern void IMPTR3_F77( long *BNDCON, long *POSTLT, int *CRIT, int *NPT, 
			  int *SIZHT, int *MAXFC, double (*VCL)[3], int *VM,
			  int *NFC, int *NTETRA, int (*BF)[3], int (*FC)[7],
			  int *HT, int *NFACE );
  extern void IMPTRF_F77(long *BNDCON, int *CRIT, int *NPT, int *SIZHT,
			 int *MAXFC, double (*VCL)[3], int *VM, int *NFC,
			 int *NTETRA, int *HDAVFC, int (*BF)[3], int (*FC)[7],
			 int *HT);
  /**
   * 
   */

  double tol=TOL;

  long bndcon;
  long postlt;

  int (*tetra)[4];
  int (*bf)[3];
  int (*fc)[7];
  int *ht;
  int *vm;
  int sizeht;
  int maxbf=MAX_BF;
  int maxfc=MAX_FC;
  int nface=0;
  int nfc=0;
  int nbf=0;
  int i_node;
  int err_code;
  int crit;

  /**
   * initialize
   */

  *number_tetra=0;

  /**
   * allocate auxiliary data
   */

  sizeht=number_nodes*3/2;

  sizeht = PRIME_F77( &sizeht );

  if( (ht=calloc(sizeht,sizeof(int))) == NULL ){
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  if( (bf=calloc(maxbf,sizeof(int[3]))) == NULL ){
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  if( (fc=calloc(maxfc,sizeof(int[7]))) == NULL ){
    ERROR("malloc()");
    exit(EXIT_SUCCESS);
  }

  if( (vm=calloc((number_nodes+1),sizeof(int))) == NULL ){
    ERROR("malloc()");
    exit(EXIT_SUCCESS);
  }
  
  /**
   * initialize vm 
   */
  
  for( i_node=0; i_node < number_nodes; i_node++ ){
    vm[i_node] = i_node+1;
  }

  /**
   * call DTRIW3
   */

  printf("Generating initial mesh..\n");

  INITCB_F77( &tol );
  DTRIW3_F77( &number_nodes, &sizeht, &maxbf, &maxfc, nodes, vm, &nbf,
	   &nfc, &nface, number_tetra, bf, fc, ht );

  if( (err_code=extract_error()) ){
    fprintf(stderr, "[dtriw3_] error no. %d\n", err_code);
    exit(EXIT_FAILURE);
  }

  printf("Initial mesh successfully generated.\n");

  /**
   * improve triangulation
   */

  bndcon=0L;
  postlt=1L;

  crit=3;
  printf("Improving mesh-method %d\n", crit);

  IMPTR3_F77( &bndcon, &postlt, &crit, &number_nodes, &sizeht, &maxfc, nodes,
 	   vm, &nfc, number_tetra, bf, fc, ht, &nface );
  IMPTRF_F77(&bndcon, &crit, &number_nodes, &sizeht, &maxfc, nodes, vm, &nfc,
 	     number_tetra, &(fc[1][6]), bf, fc, ht);

  if( (err_code=extract_error()) ){
    fprintf(stderr, "[imptr3_] error no. %d\n", err_code);
    exit(EXIT_FAILURE);
  }

  crit=2;
  printf("Improving mesh-method %d\n", crit);

  IMPTR3_F77( &bndcon, &postlt, &crit, &number_nodes, &sizeht, &maxfc, nodes,
 	   vm, &nfc, number_tetra, bf, fc, ht, &nface );
  IMPTRF_F77(&bndcon, &crit, &number_nodes, &sizeht, &maxfc, nodes, vm, &nfc,
 	  number_tetra, &(fc[1][6]), bf, fc, ht);

  if( (err_code=extract_error()) ){
    fprintf(stderr, "[imptr3_] error no. %d\n", err_code);
    exit(EXIT_FAILURE);
  }

  crit=3;
  printf("Improving mesh-method %d\n", crit);

  IMPTR3_F77( &bndcon, &postlt, &crit, &number_nodes, &sizeht, &maxfc, nodes,
 	   vm, &nfc, number_tetra, bf, fc, ht, &nface );
  IMPTRF_F77(&bndcon, &crit, &number_nodes, &sizeht, &maxfc, nodes, vm, &nfc,
 	 number_tetra, &(fc[1][6]), bf, fc, ht);

  if( (err_code=extract_error()) ){
    fprintf(stderr, "[imptr3_] error no. %d\n", err_code);
    exit(EXIT_FAILURE);
  }

  crit=2;
  printf("Improving mesh-method %d\n", crit);

  IMPTR3_F77( &bndcon, &postlt, &crit, &number_nodes, &sizeht, &maxfc, nodes,
 	      vm, &nfc, number_tetra, bf, fc, ht, &nface );
  IMPTRF_F77(&bndcon, &crit, &number_nodes, &sizeht, &maxfc, nodes, vm, &nfc,
 	     number_tetra, &(fc[1][6]), bf, fc, ht);

 if( (err_code=extract_error()) ){
    fprintf(stderr, "[imptr3_] error no. %d\n", err_code);
    exit(EXIT_FAILURE);
  }

  /**
   * allocate space for tetra
   */

  if( (tetra=malloc(*number_tetra*sizeof(int[4]))) == NULL ){
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  TETLST_F77( &nfc, vm, fc, number_tetra, tetra );

  if( (err_code=extract_error()) ){
    fprintf(stderr, "[tetlst_] error no. %d\n", err_code);
    exit(EXIT_FAILURE);
  }

  /**
   * obtain mesh quality indicators
   */
  
  check_mesh_quality( min_sigma, min_rho, min_eta, nodes, number_nodes, 
		      tetra, *number_tetra );


  /**
   * cleanup
   */

  free( ht );
  free( bf );
  free( fc );
  free( vm );

  return(tetra);

}



/**
 * this server starts with stdin and stdout connected to quasi;
 */ 


int 
main( int    ac,
      char **av )
{

  double (*position)[3];
  double   min_eta;
  double   min_rho;
  double   min_sigma;

  int (*tets)[4];
  int number_nodes;
  int number_tetra;
  int sd;
  int fd;
  /* int i_node; */

  /**
   * install exit handler that kills parent on exit
   */
 
  if(atexit(kill_parent)) {
    ERROR("atexit()");
    exit(EXIT_FAILURE);
  }
  
  /**
   * setup local socket
   */

  sd=create_server();

  /**
   * we connect once so there is a single accept() call only
   */
  
  if ((fd = accept(sd, NULL, NULL)) == -1) {
    ERROR("accept()");
    exit(EXIT_FAILURE);
  }

  /**
   * loop here waiting for data on stdin
   */

  for( ;; ){

    int n_bytes;

    /**
     * block in call to readn()
     */

    if( (n_bytes=readn(fd, &number_nodes, sizeof(int))) < 0 ){
      ERROR("readn()");
      exit(EXIT_FAILURE);
    }
    
    /**
     * allocate space for position
     */

    if( (position=malloc(number_nodes*sizeof(double[3]))) == NULL ){
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    /**
     * read in position in one chunk
     */

    if( readn(fd, position[0], number_nodes*sizeof(double[3]) )< 0 ){
      ERROR( "readn()" );
      exit( EXIT_FAILURE );
    }

    /**
     * process nodal data and generate tets
     */

    tets=generate_tets( &number_tetra, &min_eta, &min_rho, &min_sigma,
			position, number_nodes );

    /**
     * write element data to stdout
     */

    if( writen(fd, &number_tetra, sizeof(int) ) < 0){
      ERROR( "writen()" );
      exit( EXIT_FAILURE );
    }

    if( writen(fd, tets, number_tetra*sizeof(int[4]) ) < 0){
      ERROR( "writen()" );
      exit( EXIT_FAILURE );
    }

    /**
     * write mesh quality data
     */

    if( writen(fd, &min_eta, sizeof(double) ) < 0){
      ERROR( "writen()" );
      exit( EXIT_FAILURE );
    }

    if( writen(fd, &min_rho, sizeof(double) ) < 0){
      ERROR( "writen()" );
      exit( EXIT_FAILURE );
    }

    if( writen(fd, &min_sigma, sizeof(double) ) < 0){
      ERROR( "writen()" );
      exit( EXIT_FAILURE );
    }
     

    /**
     * cleanup
     */

    free(tets);
    free(position);


  }


  // this is weird...
  exit(EXIT_SUCCESS);
  return(EXIT_SUCCESS);

}

