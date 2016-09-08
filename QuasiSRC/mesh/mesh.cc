/**
 * $Log: indent_fcc.c,v $
 * Revision 1.11  2002/04/05 08:34:08  fago
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.10  2000/05/18 17:42:44  knap
 * Changed to match new format of the restart file.
 *
 * Revision 1.9  2000/04/20 01:08:58  knap
 * Ported to HP-UX.
 *
 * Revision 1.8  2000/04/12 22:06:49  knap
 * Removed references to atoms.
 *
 * Revision 1.3  2000/03/17 18:24:43  knap
 * Some minor correction.
 *
 * Revision 1.2  2000/02/24 01:12:00  knap
 * Various corrections aimed at improving scalability.
 *
 * Revision 1.1  2000/02/11 17:53:40  knap
 * Initial rev.
 *
 * Revision 1.7  1999/12/16 21:43:14  knap
 * Added second order quad and simple clusters.
 *
 * Revision 1.6  1999/12/01 23:14:14  knap
 * Minor corrections.
 *
 * Revision 1.5  1999/09/14 01:49:10  knap
 * *** empty log message ***
 *
 * Revision 1.4  1999/09/02 23:41:59  knap
 * Initial rev.
 *
 * Revision 1.3  1999/07/27 22:02:22  knap
 * Changed config file parser.
 *
 * Revision 1.2  1999/07/02 17:25:42  knap
 * Initial revision.
 *
 * Revision 1.1  1999/03/18 22:01:47  knap
 * Basics ready.
 *
 *
 * generate fcc lattice for nanoindentation
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#error limits.h not found.
#endif /* HAVE_LIMITS_H */

#ifdef HAVE_RPC_RPC_H
#include <rpc/rpc.h>
#else
#error rpc/rpc.h not found.
#endif /* HAVE_RPC_RPC_H */

#ifdef HAVE_RPC_XDR_H
#include <rpc/xdr.h>
#else
#error rpc/xdr.h not found.
#endif /* HAVE_RPC_XDR_H */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found.
#endif /* HAVE_UNISTD_H */

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#else
#error signal.h not found.
#endif /* HAVE_SIGNAL_H */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found.
#endif /* HAVE_SYS_TYPES_H */

#ifdef HAVE_SCHED_H
#include <sched.h>
#else
#error sched.h not found.
#endif /* HAVE_SCHED_H */

#include <vector>

#include "read_pipe.h"
#include "DataTypes.h"

static const char *tecplot_header_fields[] = 
{
#define TITLE       0
  "TITLE",
#define VARIABLES   1
  "VARIABLES",
#define ZONE        2
  "ZONE",
#define FEPOINT     3
  "FEPOINT",
#define TETRAHEDRON 4
  "TETRAHEDRON"
};

#ifndef lint
static char rcsid[]="$Id: indent_fcc.c,v 1.11 2002/04/05 08:34:08 fago Exp $";
#endif /* !lint */

static int ERROR_SPIN = 0;

#define ERROR(mesg) do{fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__);       \
                       perror(mesg); } while(ERROR_SPIN)

/**
 * read config file
 * input:
 *       output_file_name-pointer to space when the name of output file is
 *                        saved
 *       material        -space for material name
 *       config_file_name-name of config file
 *       nx, ny, nz      -sample size
 *       ax, ay, az      -atomistic region size
 *       ix, iy          -rectangular indentor size
 * output:
 *
 */

#define DEFAULT_INPUT_EXT ".ini"
#define DATA_FILENAME_SUFFIX ".inp.gz"
#define MAX_LATTICE_NUMBER 50
#define NODE_FILENAME_SUFFIX    ".plt.gz"

static const char *keys[]={
#define DATA_FILE        0
  "data_file",
#define MESH_FLAG        1
  "mesh_flag",
#define NUM_LATTICE      2
  "numLattice",
#define MATERIAL         3
  "material",
#define LATTICE_A        4
  "a1",
#define LATTICE_B        5
  "a2",
#define LATTICE_C        6
  "a3",
#define LATTICE_TYPE     7
  "lattice_type",
#define ATOMISTIC        8
  "atomistic",
#define FULL             9
  "full",
#define MESH_LEVELS      10
  "mesh_levels",
#define N_LEVEL          11
  "n_level",
#define PROBLEM_SIZE     12
  "problem_size",
#define ZPLANE           13
  "zplane",
#define SYSTEM_NAME      14
  "system_name",
#define POSITION_FLAG    15
  "position_flag",
#define SHIFT            16
  "shift",
#define DEBUGDATA_FLAG   17
  "debug_data_flag",
};

/**
 * XDR error function
 */

static void
xdr_err(const char *msg)
{
  fprintf(stderr, "%s\n", msg);
  exit(EXIT_FAILURE);
}

//
//  converts integer to binary xyz mask
//
char *xyz_fixity(unsigned int x)
{
    char * data = NULL;

    long int path_max;
    path_max = pathconf(".", _PC_PATH_MAX);
    data = (char *) malloc(path_max*sizeof(char));    

    // X fixity
    if((x & FIX_X_MASK) == FREE_MASK)
    {
      // X is free

      if((x & FIX_Y_MASK) == FREE_MASK)
      {
        // Y is free

        if((x & FIX_Z_MASK) == FREE_MASK)
        {
          // Z is free
          sprintf(data,"000");
          // printf("here 000 x = %d\n", x);
        }
        else
        {
          // Z is fixed
          sprintf(data,"100");
          // printf("here 100 x = %d\n", x);
        }        
      }
      else
      {
        // Y is fixed

        if((x & FIX_Z_MASK) == FREE_MASK)
        {
          // Z is free
          sprintf(data,"010");
          // printf("here 010 x = %d\n", x);
        }
        else
        {
          // Z is fixed
          sprintf(data,"110");
          // printf("here 110 x = %d\n", x);
        }                
      }
    }
    else
    {
      if((x & FIX_Y_MASK) == FREE_MASK)
      {
        // Y is free

        if((x & FIX_Z_MASK) == FREE_MASK)
        {
          // Z is free
          sprintf(data,"001");
          // printf("here 001 x = %d\n", x);
        }
        else
        {
          // Z is fixed
          sprintf(data,"101");
          // printf("here 101 x = %d\n", x);
        }        
      }
      else
      {
        // Y is fixed

        if((x & FIX_Z_MASK) == FREE_MASK)
        {
          // Z is free
          sprintf(data,"011");
          // printf("here 011 x = %d\n", x);
        }
        else
        {
          // Z is fixed
          sprintf(data,"111");
          //printf("here 111 x = %d\n", x);
        }                
      }
    }



    // if((x & FIX_X_MASK) != FREE_MASK && 
    //    (x & FIX_Y_MASK) != FREE_MASK &&
    //    (x & FIX_Z_MASK) != FREE_MASK)
    //   sprintf(data,"111");
    // else if((x & FIX_X_MASK) != FREE_MASK && 
    //         (x & FIX_Y_MASK) != FREE_MASK &&
    //         (x & FIX_Z_MASK) == FREE_MASK)
    //   sprintf(data,"011");
    // else if((x & FIX_X_MASK) != FREE_MASK && 
    //         (x & FIX_Y_MASK) == FREE_MASK &&
    //         (x & FIX_Z_MASK) != FREE_MASK )
    //   sprintf(data,"101");
    // else if((x & FIX_X_MASK) == FREE_MASK && 
    //         (x & FIX_Y_MASK) != FREE_MASK &&
    //         (x & FIX_Z_MASK) != FREE_MASK)
    //   sprintf(data,"110");    
    // else if((x & FIX_X_MASK) != FREE_MASK && 
    //         (x & FIX_Y_MASK) == FREE_MASK &&
    //         (x & FIX_Z_MASK) == FREE_MASK )
    //   sprintf(data,"001");
    // else if((x & FIX_X_MASK) == FREE_MASK && 
    //         (x & FIX_Y_MASK) != FREE_MASK &&
    //         (x & FIX_Z_MASK) == FREE_MASK )
    //   sprintf(data,"010");
    // else if((x & FIX_X_MASK) == FREE_MASK && 
    //         (x & FIX_Y_MASK) == FREE_MASK &&
    //         (x & FIX_Z_MASK) != FREE_MASK )
    //   sprintf(data,"100");    
    // else
    //   sprintf(data,"000");


    //   b[0] = "1";
    // else
    //   b[0] = "0";

    // if((x & FIX_Y_MASK) != FREE_MASK)
    //   b[1] = "1";
    // else
    //   b[1] = "0";
    
    // if((x & FIX_Z_MASK) != FREE_MASK)
    //   b[2] = "1";    
    // else
    //   b[2] = "0";

    // b[3] = '\0';

    return data;
}

//
//  converts integer to binary
//
//  Source: Stackoverflow website
//
const char *byte_to_binary(int x)
{
    static char b[9];
    b[0] = '\0';

    int z;
    for (z = 128; z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}

//
// read data from mesh.ini to generate initial configuration data
// /
void
readDataFile(const char             *input_file_name,
  char                              *data_file_name,
  char                              *system_name,
  int                               &mesh_flag,
  int                               &position_flag,
  int                               &numLattice,
  char **                           material_name,
  std::vector<std::vector<double>>  &lattice_constant_vec,
  int                               &lattice_type,
  std::vector<int>                  &sample_size,
  std::vector<int>                  &atomistic_size,
  int                               &n_levels,
  int                               **d_inc,
  int                               &problem_size,
  std::vector<std::vector<double>>  &shiftVector,
  int                               &zplaneOnly,
  int                               &debugDataFlag,
  int                               ac,
  char                              *av[])
{
  FILE *input_file;
  
  char local_input_file_name[PATH_MAX+1];
  char line[LINE_MAX+1];
  int levelMalloc = 0;

  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);    

  char *local_material_name;
  local_material_name = (char *) malloc(path_max*sizeof(char));

  //
  //  if input_file_name is empty, create default input_file_name
  //
  if( input_file_name == (char *) NULL )
  {
    sprintf(local_input_file_name, "%s%s", av[0], DEFAULT_INPUT_EXT);
    input_file_name=local_input_file_name;
  }
    
  /**
   * open file
   */
  if( (input_file=fopen(input_file_name, "r")) ==  NULL )
  {
    ERROR(input_file_name);
    exit(EXIT_FAILURE);
  }

  //
  //  some elementary operations on variables to be read
  //
  memset(material_name, 0, MAX_LATTICE_NUMBER * sizeof(char *));
  for(int i=0; i < MAX_LATTICE_NUMBER; i++)
  { 
    if((material_name[i] = (char *) malloc(100* sizeof(char))) == NULL)
    {
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }
  }

  lattice_constant_vec.clear();
  lattice_constant_vec.resize(3);
  for(int dof=0; dof < 3; dof++)
  {
    lattice_constant_vec[dof].push_back(0.0);
    lattice_constant_vec[dof].push_back(0.0);
    lattice_constant_vec[dof].push_back(0.0);
  }

  sample_size.clear();
  sample_size.push_back(0);
  sample_size.push_back(0);
  sample_size.push_back(0);
   
  atomistic_size.clear();
  atomistic_size.push_back(0);
  atomistic_size.push_back(0);
  atomistic_size.push_back(0);

  shiftVector.clear();

  /**
   * process lines
   */
  while( fgets(line, LINE_MAX, input_file) != NULL )
  {
    /**
     * skip comments and empty lines
     */
    if( line[0]=='#' || line[0]=='\n' ) 
      continue;

    /**
     *  read data
     */
    if( strstr(line, keys[DATA_FILE]) != NULL )
      sscanf(line, "%*s %s", data_file_name);

    if( strstr(line, keys[SYSTEM_NAME]) != NULL )
      sscanf(line, "%*s %s", system_name);

    if( strstr(line, keys[MESH_FLAG]) != NULL )
      sscanf(line, "%*s %d", &mesh_flag);

    if( strstr(line, keys[POSITION_FLAG]) != NULL )
    {
      sscanf(line, "%*s %d", &position_flag);    
    }    

    if( strstr(line, keys[NUM_LATTICE]) != NULL )
    {
      sscanf(line, "%*s %d", &numLattice);    
    }

    if( strstr(line, keys[LATTICE_A]) != NULL )
    {
      double a;
      double b;
      double c;
      sscanf(line, "%*s %lf %lf %lf", &a, &b, &c);

      lattice_constant_vec[0][0] = a;
      lattice_constant_vec[0][1] = b;
      lattice_constant_vec[0][2] = c;
    }

    if( strstr(line, keys[LATTICE_B]) != NULL )
    {
      double a;
      double b;
      double c;
      sscanf(line, "%*s %lf %lf %lf", &a, &b, &c);

      lattice_constant_vec[1][0] = a;
      lattice_constant_vec[1][1] = b;
      lattice_constant_vec[1][2] = c;
    }

    
    if( strstr(line, keys[LATTICE_C]) != NULL )
    {
      double a;
      double b;
      double c;
      sscanf(line, "%*s %lf %lf %lf", &a, &b, &c);

      lattice_constant_vec[2][0] = a;
      lattice_constant_vec[2][1] = b;
      lattice_constant_vec[2][2] = c;
    }


    if( strstr(line, keys[LATTICE_TYPE]) != NULL )
      sscanf(line, "%*s %d", &lattice_type);    

    if( strstr(line, keys[ATOMISTIC]) != NULL )
    {
      int l;
      int m;
      int n;

      sscanf(line, "%*s %d %d %d", &l, &m, &n);

      atomistic_size[0] = l;
      atomistic_size[1] = m;
      atomistic_size[2] = n;
    }

    if( strstr(line, keys[FULL]) != NULL )
    {
      int l;
      int m;
      int n;

      sscanf(line, "%*s %d %d %d", &l, &m, &n);

      sample_size[0] = l;
      sample_size[1] = m;
      sample_size[2] = n;
    }    

    if( strstr(line, keys[MESH_LEVELS]) != NULL )
    {
      int n;
      sscanf(line, "%*s %d", &n_levels);
    }

    if( strstr(line, keys[N_LEVEL]) != NULL )
    {
      int levelNum;
      int level;
      sscanf(line, "%*s %d %d",&levelNum,&level);
      
      if(levelMalloc == 0)
      {
        if((*d_inc = (int *) malloc(sizeof(int) * n_levels)) == NULL)
        {
          ERROR("malloc()");
          exit(EXIT_FAILURE);
        }
      
        levelMalloc++;
      }

      (*d_inc)[levelNum] = level;
    }

    if( strstr(line, keys[PROBLEM_SIZE]) != NULL )
      sscanf(line, "%*s %d", &problem_size);   

    if( strstr(line, keys[SHIFT]) != NULL )
    {
      std::vector<double> data;
      for(int i=0; i < 3; i++)
        data.push_back(0.0);

      sscanf(line, "%*s  %*d   %lf  %lf  %lf", 
        &data[0],
        &data[1],
        &data[2]);  

      shiftVector.push_back(data);
    }

    if( strstr(line, keys[ZPLANE]) != NULL )
      sscanf(line, "%*s %d", &zplaneOnly);    

    //
    //  debugDataFlag = 0 : no debug o/p
    if( strstr(line, keys[DEBUGDATA_FLAG]) != NULL )
      sscanf(line, "%*s %d", &debugDataFlag);    

    //
    //  read material name
    //
    if( strstr(line, keys[MATERIAL]) != NULL )
    {
      int n;
      sscanf(line, "%*s %d %s", &n, local_material_name);

      if(n < numLattice)
        strcpy(material_name[n], local_material_name);
    }
  }

  /**
   * close file
   */
  fclose(input_file);
  
  return;
} 

/**
 * check if this atom should be marked as node
 */
static int 
checkForNode(int n1, int n2, int n3, int nx, int ny, int nz,
         int ax, int ay, int az, int lattice_type)
{
  int center[3];
  int coor[3];
  int p=4;
  int l1,l2,l3;
  int shiftx = 0;
  int shifty = 0;


  if (nx%2==1){
    center[0]=(nx-1)/2;
    shiftx = 1;
  }
  else
    center[0]=nx/2;

  if (ny%2==1){
    center[1]=(ny-1)/2;
    shifty = 1;
  }
  else
    center[1]=ny/2;

  center[2]=nz;
            
  l1=-n1+n2+n3;
  l2= n1-n2+n3;
  l3= n1+n2-n3;

  /** 
   * atomistic region
   */

  if( n1 >= center[0] - ax && n1 <= center[0] + ax + shiftx &&
      n2 >= center[1] - ay && n2 <= center[1] + ay + shifty &&
      n3 >= center[2] - az && n3 <= center[2]) return(1);

  /**
   * free surface atoms
   */

  /* if( n1 <= 1 || n1 >= nx-1 ) return(1);
     if( n2 <= 1 || n2 >= ny-1 ) return(1);
     if( n3 >= nz-1)             return(1);
  */

  /**
   * corners of the sample are always nodes
   */

  if( n1 == 0  && n2 == 0  && n3 == 0 ) return(1);
  if( n1 == nx && n2 == 0  && n3 == 0 ) return(1);
  if( n1 == 0  && n2 == ny && n3 == 0 ) return(1);
  if( n1 == nx && n2 == ny && n3 == 0 ) return(1);

  if( n1 == 0  && n2 == 0  && n3 == nz ) return(1);
  if( n1 == nx && n2 == 0  && n3 == nz ) return(1);
  if( n1 == 0  && n2 == ny && n3 == nz ) return(1);
  if( n1 == nx && n2 == ny && n3 == nz ) return(1);

  /**
   *
   */
  if(lattice_type == 0)
  {
    coor[0]=ax;
    coor[1]=ay;
    coor[2]=az;

    while(1)
    {
      coor[0] += p;
      coor[1] += p;
      coor[2] += p;

      /* if( coor[0] > center[0] || coor[1] > center[1] || coor[2] > center[2] ) */
      if( coor[0] > center[0] && coor[1] > center[1] && coor[2] > center[2] )
        break;

      if( n1 >= center[0]-coor[0] && n1 <= center[0]+coor[0] &&
          n2 >= center[1]-coor[1] && n2 <= center[1]+coor[1] &&
          n3 >= center[2]-coor[2] && n3 <= center[2]+coor[2] )
        if( l1%(p) == 0 && l2%(p) == 0 && l3%(p) == 0 ) 
          return(1);

      p *= 2;
    }
  }

  return(0);
}

/**
 *
 */
#define NODES_ALLOC_BUCKET 1000

int
main(int argc, char **argv)
{
  using namespace quasicontinuum;
  //
  //  local variables
  //
   
  struct lattice_t lattice;

  struct node_t *nodes;
  int  n_nodes;
  int  new_n_nodes;
  int  n_alloc_nodes;
  int n_sites = 0;
  int fixity;
  int i_node;  

  char *name;

  int center[3];
  
  int n1;
  int n2;
  int n3;

  int l1;
  int l2;
  int l3;

  //
  // data to be read from input file
  //
  FILE *data_file;
  char data_filename[256];
  char *system_name;
  int  mesh_flag;
  int  position_flag;
  int numLattice;
  char *material_name[MAX_LATTICE_NUMBER];
  std::vector<std::vector<double>> lattice_constant_vec;
  int lattice_type;
  std::vector<int> sample_size;
  std::vector<int> atomistic_size;
  int n_levels;
  int *d_inc;
  int problem_size;
  std::vector<std::vector<double>>  shiftVector;
  int zplaneOnly = 0; // default values is 0
  int debugDataFlag = 0; // default : no debug data output

  /**
   * process input file
   */

  if( argc == 1 ) 
    readDataFile((char *) NULL,
      data_filename,
      system_name,
      mesh_flag,
      position_flag,
      numLattice,
      material_name,
      lattice_constant_vec,
      lattice_type,
      sample_size,
      atomistic_size,
      n_levels,
      &d_inc,
      problem_size,
      shiftVector,
      zplaneOnly,
      debugDataFlag,
      argc,
      argv);
  else
    readDataFile(argv[1],
      data_filename,
      system_name,
      mesh_flag,
      position_flag,
      numLattice,
      material_name,
      lattice_constant_vec,
      lattice_type,
      sample_size,
      atomistic_size,
      n_levels,
      &d_inc,
      problem_size,
      shiftVector,
      zplaneOnly,
      debugDataFlag,
      argc,
      argv);  

  // printf("problem_size = %d\n",problem_size);

  //
  //  writing data to local data
  //
  int nx = sample_size[0];
  int ny = sample_size[1];
  int nz = sample_size[2];

  int ax = atomistic_size[0];
  int ay = atomistic_size[1];
  int az = atomistic_size[2];

  /**
   * position of an indentor
   */
  if (nx%2==1)
    center[0]=(nx-1)/2;
  else
    center[0]=nx/2;

  if (ny%2==1)
    center[1]=(ny-1)/2;
  else
    center[1]=ny/2;

  center[2]=nz;  

  // copy materials.dat to material.data_file
  if (lattice_type == 0)
  {
    lattice.type = FCC;
    sprintf(&(lattice.name[0]), "%s", "fcc");
  }

  if (lattice_type == 1 || lattice_type == 2)
  {
    lattice.type = BCC;
    sprintf(&(lattice.name[0]), "%s", "bcc");
  }

  /**
   * get lattice parameter
   */
  lattice.a1[0]=lattice_constant_vec[0][0];
  lattice.a1[1]=lattice_constant_vec[0][1];
  lattice.a1[2]=lattice_constant_vec[0][2];

  lattice.a2[0]=lattice_constant_vec[1][0];
  lattice.a2[1]=lattice_constant_vec[1][1];
  lattice.a2[2]=lattice_constant_vec[1][2];
  
  lattice.a3[0]=lattice_constant_vec[2][0];
  lattice.a3[1]=lattice_constant_vec[2][1];
  lattice.a3[2]=lattice_constant_vec[2][2];

  lattice.l_start[0] = 0;
  lattice.l_start[1] = 0;
  lattice.l_start[2] = 0;

  lattice.l_end[0]=0;
  lattice.l_end[1]=0;
  lattice.l_end[2]=0;

  //
  //  generate node data
  //
  // allocate first bucket for nodes
  if((nodes =(struct node_t *) malloc(NODES_ALLOC_BUCKET*sizeof(struct node_t))) == NULL)
  {
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  n_nodes=0;
  n_alloc_nodes=NODES_ALLOC_BUCKET;

  //
  //  fixity method :
  //
  //  flag = 1 - Fixed at bottom, free at top, free z direction
  //  on sides.
  //
  //  flag = 2 - Came with the code. Not intuitive and clear how
  //  the fixity mask is defined
  //
  int flag_fixity = 1;

  switch (lattice.type)
  {
    case FCC:
    {
      for(n3=0; n3 <= nz; n3++)
        for(n2=0; n2 <= ny; n2++)
          for(n1=0; n1 <= nx; n1++)
          {
            l1=-n1+n2+n3;
            l2= n1-n2+n3;
            l3= n1+n2-n3;

            if( l1%2 == 0 && l2%2 == 0 && l3%2 == 0 )
            {
              n_sites++;

              /**
               * 
               */
              if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
              if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
              if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

              /**
               * fixity data
               */

              fixity = FREE_MASK;
          
              if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
           
              if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
              if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
          
              /**
               * surface data
               */

              if( n3 == 0 || n3 == nz)
                fixity |= SURFACE_MASK | SURFACE_Z_MASK;
              if( n2 == 0 || n2 == ny )
                fixity |= SURFACE_MASK | SURFACE_Y_MASK;
              if( n1 == 0 || n1 == nx )
                fixity |= SURFACE_MASK | SURFACE_X_MASK;


              if( checkForNode(n1, n2, n3, nx, ny, nz, ax, ay, az, lattice_type) )
              {
                /**
                 * check space for nodes
                 */

                if(n_nodes >= n_alloc_nodes)
                {
                  n_alloc_nodes += NODES_ALLOC_BUCKET;
                  if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                          sizeof(struct node_t))) == NULL)
                  {
                    ERROR("realloc()");
                    exit(EXIT_SUCCESS);
                  }
                }
      
                nodes[n_nodes].l[0] = n1;
                nodes[n_nodes].l[1] = n2;
                nodes[n_nodes].l[2] = n3;
                nodes[n_nodes].fix_mask = fixity;

                n_nodes++;
              }
            }
          }
    }
    break;

    case BCC:
    {
      // printf("nlevels = %d\n", n_levels);
      // for(int i=0; i < n_levels; i++)
      //   printf("%dth level = %d\n",i, d_inc[i]);
      if(mesh_flag == 0)
      {
        int level_count;
        int corner_check;
        int i_count;
        int j_count;

        int d_level[n_levels];

        for(i_count=0;i_count<n_levels;i_count++)
        { 
          d_level[i_count]=0;
          for(j_count=0;j_count<=i_count;j_count++)
          {
            d_level[i_count] += d_inc[j_count];
          }
        }
    
        int shiftx=0;
        int shifty=0;
        int shiftz=0;

        if (nx%2==1) shiftx=2;
        if (ny%2==1) shifty=2;
        if (nz%2==1) shiftz=1;

        int nzValue;
        if(problem_size == 0)
          nzValue = nz;
        else if(problem_size == 1)
        {
          nzValue = 2*nz;
        }

        for(n3=0; n3 <= nzValue; n3++)
          for(n2=0; n2 <= ny; n2++)
            for(n1=0; n1 <= nx; n1++)
            {
              n_sites++;

              /**
               * 
               */

              if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
              if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
              if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

              /**
               * fixity data
               */

              fixity = FREE_MASK;

              if(flag_fixity == 1)
              {              
                if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
              }
              else
              {
                if (   n1==0 || n1==nx || n1 == d_inc[n_levels-1] 
                    || n1 == nx - d_inc[n_levels-1] 
                    || n2==0 || n2==ny || n2 == d_inc[n_levels-1] 
                    || n2 == ny - d_inc[n_levels-1] 
                    || n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] 
                    || n3 == nzValue - d_inc[n_levels-1]
                   )
                  fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
              }

              /**
               * surface data
               */

              if( n3 == 0 || n3 == nzValue)
                fixity |= SURFACE_MASK | SURFACE_Z_MASK;
              if( n2 == 0 || n2 == ny )
                fixity |= SURFACE_MASK | SURFACE_Y_MASK;
              if( n1 == 0 || n1 == nx )
                fixity |= SURFACE_MASK | SURFACE_X_MASK;


              if(checkForNode(n1, n2, n3, nx, ny, nz, ax, ay, az, lattice_type))
              {
                /**
                 * check space for nodes
                 */

                if(n_nodes >= n_alloc_nodes)
                {
                  n_alloc_nodes += NODES_ALLOC_BUCKET;
                  if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                      sizeof(struct node_t))) == NULL)
                  {
                    ERROR("realloc()");
                    exit(EXIT_SUCCESS);
                  }
                }
        
                nodes[n_nodes].l[0] = n1;
                nodes[n_nodes].l[1] = n2;
                nodes[n_nodes].l[2] = n3;
                nodes[n_nodes].fix_mask = fixity;

                n_nodes++;
              }
            }

        //
        // check for z plane only
        //
        if(zplaneOnly == 0)
        {
          /**
           * generate nodes z plane
           */
          for (level_count=0;level_count<n_levels;level_count++)
          {
            n3=nz-az-d_level[level_count];
            for(n2=center[1]-ay-d_level[level_count]; n2 <= shifty+center[1]+ay+d_level[level_count]; n2=n2+d_inc[level_count])
            {
              if (n2==ny+1 && ny%2==1)
                n2--;
              for(n1=center[0]-ax-d_level[level_count]; n1 <= shiftx+center[0]+ax+d_level[level_count]; n1=n1+d_inc[level_count])
              {
                if (n1==nx+1 && nx%2==1)
                  n1--;

                corner_check=0;

                l1=-n1+n2+n3;
                l2= n1-n2+n3;
                l3= n1+n2-n3;

                if ((lattice.type == FCC && l1%2 == 0 && l2%2 == 0 && l3%2 == 0) || lattice.type == BCC)
                {
                  if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
                  if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
                  if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

                  /**
                   * fixity data
                   */

                  fixity = FREE_MASK;

                  if(flag_fixity == 1)
                  { 
                    if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                    if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                    if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
                  }
                  else
                  {
                    if ( n1==0 || n1==nx || n1 == d_inc[n_levels-1] 
                        || n1 == nx - d_inc[n_levels-1] 
                        || n2==0 || n2==ny || n2 == d_inc[n_levels-1] 
                        || n2 == ny - d_inc[n_levels-1] 
                        || n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] 
                        || n3 == nzValue - d_inc[n_levels-1]
                       ) 
                      fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
                  }

                  /**
                   * surface data
                   */

                  if( n3 == 0 || n3 == nzValue)
                    fixity |= SURFACE_MASK | SURFACE_Z_MASK;
                  if( n2 == 0 || n2 == ny )
                    fixity |= SURFACE_MASK | SURFACE_Y_MASK;
                  if( n1 == 0 || n1 == nx )
                    fixity |= SURFACE_MASK | SURFACE_X_MASK;


                  /**
                   * check space for nodes
                   */

                  if(n_nodes >= n_alloc_nodes)
                  {
                    n_alloc_nodes += NODES_ALLOC_BUCKET;

                    if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*sizeof(struct node_t))) == NULL)
                    {
                      ERROR("realloc()");
                      exit(EXIT_SUCCESS);
                    }
                  }

                  if( n1 == 0  && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == 0 ) corner_check=1;

                  if( n1 == 0  && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == nz ) corner_check=1;

                  if(corner_check==0)
                  {
                    nodes[n_nodes].l[0] = n1;
                    nodes[n_nodes].l[1] = n2;
                    nodes[n_nodes].l[2] = n3;
                    nodes[n_nodes].fix_mask = fixity;

                    n_nodes++;
                  }
                }
              }
            }
          }

          /*generate nodes +x plane*/
          for (level_count=0;level_count<n_levels;level_count++)
          {
            n1=shiftx+center[0]+ax+d_level[level_count];
            
            if (n1==nx+1 && nx%2==1) 
              n1--;
            
            for(n2=center[1]-ay-d_level[level_count]; n2 <= shifty+ center[1]+ay+d_level[level_count]; n2=n2+d_inc[level_count])
            {
              if (n2==ny+1 && ny%2==1)
                n2--;

              for(n3=nz-az-d_level[level_count]+d_inc[level_count]; n3 <= nz+shiftz; n3=n3+d_inc[level_count])
              {
                if (n3==nz+shiftz && nz%2==1)
                  n3--;
                
                corner_check=0;

                l1=-n1+n2+n3;
                l2= n1-n2+n3;
                l3= n1+n2-n3;

                if ((lattice.type == FCC && l1%2 == 0 && l2%2 == 0 && l3%2 == 0) || lattice.type == BCC)
                {
                  if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
                  if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
                  if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

                  /**
                   * fixity data
                   */

                  fixity = FREE_MASK;

                  if(flag_fixity == 1)
                  { 
                    if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                    if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                    if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
                  }
                  else
                  {                  
                    if (  n1==0 || n1==nx || n1 == d_inc[n_levels-1] 
                      ||  n1 == nx - d_inc[n_levels-1] 
                      ||  n2==0 || n2==ny || n2 == d_inc[n_levels-1] 
                      ||  n2 == ny - d_inc[n_levels-1] 
                      ||  n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] 
                      ||  n3 == nzValue - d_inc[n_levels-1]
                      )
                      fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
                  }

                  /**
                   * surface data
                   */

                  if( n3 == 0 || n3 == nzValue)
                    fixity |= SURFACE_MASK | SURFACE_Z_MASK;
                  if( n2 == 0 || n2 == ny )
                    fixity |= SURFACE_MASK | SURFACE_Y_MASK;
                  if( n1 == 0 || n1 == nx )
                    fixity |= SURFACE_MASK | SURFACE_X_MASK;

                  /**
                   * check space for nodes
                   */

                  if(n_nodes >= n_alloc_nodes)
                  {
                    n_alloc_nodes += NODES_ALLOC_BUCKET;
                    if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                        sizeof(struct node_t))) == NULL)
                    {
                      ERROR("realloc()");
                      exit(EXIT_SUCCESS);
                    }
                  }

                  if( n1 == 0  && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == 0 ) corner_check=1;

                  if( n1 == 0  && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == nz ) corner_check=1;
                    
                  if(corner_check==0)
                  {
                    nodes[n_nodes].l[0] = n1;
                    nodes[n_nodes].l[1] = n2;
                    nodes[n_nodes].l[2] = n3;
                    nodes[n_nodes].fix_mask = fixity;

                    n_nodes++;
                  }
                }
              }
            }
          }

          /*generate nodes -x plane*/
          for (level_count=0;level_count<n_levels;level_count++)
          {
            n1=center[0]-ax-d_level[level_count];
          
            for(n2=center[1]-ay-d_level[level_count]; n2 <= shifty+center[1]+ay+d_level[level_count]; n2=n2+d_inc[level_count])
            {
              if (n2==ny+1 && ny%2==1)
                n2--;
              for(n3=nz-az-d_level[level_count]+d_inc[level_count]; n3 <= shiftz+nz; n3=n3+d_inc[level_count])
              {
                if (n3==nz+shiftz && nz%2==1)
                  n3--;

                corner_check=0;

                l1=-n1+n2+n3;
                l2= n1-n2+n3;
                l3= n1+n2-n3;

                if ((lattice.type == FCC && l1%2 == 0 && l2%2 == 0 && l3%2 == 0) || lattice.type == BCC)
                {
                  if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
                  if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
                  if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

                  /**
                   * fixity data
                   */

                  fixity = FREE_MASK;

                  if(flag_fixity == 1)
                  { 
                    if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                    if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                    if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
                  }
                  else
                  {  
                    if (n1==0 || n1==nx || n1 == d_inc[n_levels-1] || n1 == nx - d_inc[n_levels-1] ||
                        n2==0 || n2==ny || n2 == d_inc[n_levels-1] || n2 == ny - d_inc[n_levels-1] ||
                        n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] || n3 == nzValue - d_inc[n_levels-1]) 
                      fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
                  }

                  /**
                   * surface data
                   */

                  if( n3 == 0 || n3 == nzValue)
                    fixity |= SURFACE_MASK | SURFACE_Z_MASK;
                  if( n2 == 0 || n2 == ny )
                    fixity |= SURFACE_MASK | SURFACE_Y_MASK;
                  if( n1 == 0 || n1 == nx )
                    fixity |= SURFACE_MASK | SURFACE_X_MASK;

                  /**
                   * check space for nodes
                   */

                  if(n_nodes >= n_alloc_nodes)
                  {
                    n_alloc_nodes += NODES_ALLOC_BUCKET;
                    if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                            sizeof(struct node_t))) == NULL)
                    {
                      ERROR("realloc()");
                      exit(EXIT_SUCCESS);
                    }
                  }

                  if( n1 == 0  && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == 0 ) corner_check=1;

                  if( n1 == 0  && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == nz ) corner_check=1;
        
                  if(corner_check==0)
                  {
                    nodes[n_nodes].l[0] = n1;
                    nodes[n_nodes].l[1] = n2;
                    nodes[n_nodes].l[2] = n3;
                    nodes[n_nodes].fix_mask = fixity;

                    n_nodes++;
                  }
                }
              }
            }
          }

          /*generate nodes -y plane*/
          for (level_count=0;level_count<n_levels;level_count++)
          {
            n2=center[1]-ay-d_level[level_count];
          
            for(n1=center[0]-ax-d_level[level_count]+d_inc[level_count]; n1 <= shiftx+center[0]+ax+d_level[level_count]-d_inc[level_count]; n1=n1+d_inc[level_count])
            {
              if (n1==nx+1 && nx%2==1)
                n1--;
              for(n3=nz-az-d_level[level_count]+d_inc[level_count]; n3 <= nz+shiftz; n3=n3+d_inc[level_count])
              {
                if (n3==nz+shiftz && nz%2==1)
                  n3--;

                corner_check=0;

                l1=-n1+n2+n3;
                l2= n1-n2+n3;
                l3= n1+n2-n3;

                if ((lattice.type == FCC && l1%2 == 0 && l2%2 == 0 && l3%2 == 0) || lattice.type == BCC)
                {
                  if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
                  if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
                  if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

                  /**
                   * fixity data
                   */

                  fixity = FREE_MASK;

                  if(flag_fixity == 1)
                  { 
                    if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                    if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                    if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
                  }
                  else
                  {                    
                    if (n1==0 || n1==nx || n1 == d_inc[n_levels-1] || n1 == nx - d_inc[n_levels-1] ||
                        n2==0 || n2==ny || n2 == d_inc[n_levels-1] || n2 == ny - d_inc[n_levels-1] ||
                        n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] || n3 == nzValue - d_inc[n_levels-1]) 
                      fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
                  }
                  
                  /**
                   * surface data
                   */

                  if( n3 == 0 || n3 == nzValue)
                    fixity |= SURFACE_MASK | SURFACE_Z_MASK;
                  if( n2 == 0 || n2 == ny )
                    fixity |= SURFACE_MASK | SURFACE_Y_MASK;
                  if( n1 == 0 || n1 == nx )
                    fixity |= SURFACE_MASK | SURFACE_X_MASK;


                  /**
                   * check space for nodes
                   */

                  if(n_nodes >= n_alloc_nodes)
                  {
                    n_alloc_nodes += NODES_ALLOC_BUCKET;
                    if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                            sizeof(struct node_t))) == NULL)
                    {
                      ERROR("realloc()");
                      exit(EXIT_SUCCESS);
                    }
                  }

                  if( n1 == 0  && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == 0 ) corner_check=1;

                  if( n1 == 0  && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == nz ) corner_check=1;
        
                  if(corner_check==0)
                  {
                    nodes[n_nodes].l[0] = n1;
                    nodes[n_nodes].l[1] = n2;
                    nodes[n_nodes].l[2] = n3;
                    nodes[n_nodes].fix_mask = fixity;

                    n_nodes++;
                  }
                }
              }
            }
          }

          /*generate nodes +y plane*/
          for (level_count=0;level_count<n_levels;level_count++)
          {
            n2=shifty+center[1]+ay+d_level[level_count];
            if (n2==ny+1 && ny%2==1) 
              n2--;
            
            for(n1=center[0]-ax-d_level[level_count]+d_inc[level_count]; n1 <= shiftx+center[0]+ax+d_level[level_count]-d_inc[level_count]; n1=n1+d_inc[level_count])
            {
              if (n1==nx+1 && nx%2==1)
                n1--;
              
              for(n3=nz-az-d_level[level_count]+d_inc[level_count]; n3 <= shiftz+nz; n3=n3+d_inc[level_count])
              {
                if (n3==nz+shiftz && nz%2==1)
                  n3--;

                corner_check=0;

                l1=-n1+n2+n3;
                l2= n1-n2+n3;
                l3= n1+n2-n3;

                if ((lattice.type == FCC && l1%2 == 0 && l2%2 == 0 && l3%2 == 0) || lattice.type == BCC)
                {
                  if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
                  if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
                  if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

                  /**
                   * fixity data
                   */

                  fixity = FREE_MASK;

                  if(flag_fixity == 1)
                  { 
                    if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                    if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                    if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
                  }
                  else
                  {                    
                    if (n1==0 || n1==nx || n1 == d_inc[n_levels-1] || n1 == nx - d_inc[n_levels-1] ||
                        n2==0 || n2==ny || n2 == d_inc[n_levels-1] || n2 == ny - d_inc[n_levels-1] ||
                        n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] || n3 == nzValue - d_inc[n_levels-1]) 
                      fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
                  }

                  /**
                   * surface data
                   */

                  if( n3 == 0 || n3 == nzValue)
                    fixity |= SURFACE_MASK | SURFACE_Z_MASK;
                  if( n2 == 0 || n2 == ny )
                    fixity |= SURFACE_MASK | SURFACE_Y_MASK;
                  if( n1 == 0 || n1 == nx )
                    fixity |= SURFACE_MASK | SURFACE_X_MASK;


                  /**
                   * check space for nodes
                   */

                  if(n_nodes >= n_alloc_nodes)
                  {
                    n_alloc_nodes += NODES_ALLOC_BUCKET;
                    if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                            sizeof(struct node_t))) == NULL)
                    {
                      ERROR("realloc()");
                      exit(EXIT_SUCCESS);
                    }
                  }
                    
                  if( n1 == 0  && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == 0 ) corner_check=1;

                  if( n1 == 0  && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == nz ) corner_check=1;
              
                  if(corner_check==0)
                  {
                    nodes[n_nodes].l[0] = n1;
                    nodes[n_nodes].l[1] = n2;
                    nodes[n_nodes].l[2] = n3;
                    nodes[n_nodes].fix_mask = fixity;

                    n_nodes++;
                  }
                }
              }
            }
          }
        }
        else
        {
          printf("here\n");
          /**
           * generate nodes z plane
           */
          for (level_count=0;level_count<n_levels;level_count++)
          {
            n3=nz-az-d_level[level_count];
            for(n2=center[1]-ay;
                n2 <= shifty+center[1]+ay;
                n2=n2+d_inc[level_count])
            {
              if (n2==ny+1 && ny%2==1)
                n2--;
              
              for(n1=center[0]-ax;
                  n1 <= shiftx+center[0]+ax;
                  n1=n1+d_inc[level_count])
              {
                if (n1==nx+1 && nx%2==1)
                  n1--;

                corner_check=0;

                l1=-n1+n2+n3;
                l2= n1-n2+n3;
                l3= n1+n2-n3;

                if ((lattice.type == FCC && l1%2 == 0 && l2%2 == 0 && l3%2 == 0) || lattice.type == BCC)
                {
                  if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
                  if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
                  if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

                  /**
                   * fixity data
                   */

                  fixity = FREE_MASK;

                  if(flag_fixity == 1)
                  { 
                    if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                    if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                    if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
                  }
                  else
                  {                  
                    if (n1==0 || n1==nx || n1 == d_inc[n_levels-1] || n1 == nx - d_inc[n_levels-1] ||
                        n2==0 || n2==ny || n2 == d_inc[n_levels-1] || n2 == ny - d_inc[n_levels-1] ||
                        n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] || n3 == nzValue - d_inc[n_levels-1]) 
                      fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
                  }
                  
                  /**
                   * surface data
                   */

                  if( n3 == 0 || n3 == nzValue)
                    fixity |= SURFACE_MASK | SURFACE_Z_MASK;
                  if( n2 == 0 || n2 == ny )
                    fixity |= SURFACE_MASK | SURFACE_Y_MASK;
                  if( n1 == 0 || n1 == nx )
                    fixity |= SURFACE_MASK | SURFACE_X_MASK;

                  /**
                   * check space for nodes
                   */

                  if(n_nodes >= n_alloc_nodes)
                  {
                    n_alloc_nodes += NODES_ALLOC_BUCKET;
                    if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                            sizeof(struct node_t))) == NULL)
                    {
                      ERROR("realloc()");
                      exit(EXIT_SUCCESS);
                    }
                  }

                  if( n1 == 0  && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == 0 ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == 0 ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == 0 ) corner_check=1;

                  if( n1 == 0  && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == 0  && n3 == nz ) corner_check=1;
                  if( n1 == 0  && n2 == ny && n3 == nz ) corner_check=1;
                  if( n1 == nx && n2 == ny && n3 == nz ) corner_check=1;
                    
                  if(corner_check==0)
                  {
                    nodes[n_nodes].l[0] = n1;
                    nodes[n_nodes].l[1] = n2;
                    nodes[n_nodes].l[2] = n3;
                    nodes[n_nodes].fix_mask = fixity;

                    n_nodes++;
                  }
                }
              }
            }
          }
        }

        /**
         * generate second half
         */
        if(problem_size == 1)
        {
          int halfNumNodes = n_nodes;
          int iNode = 0;
          for(iNode = 0; iNode < halfNumNodes; ++iNode)
          {
            /**
             * check to add
             */
            if(nodes[iNode].l[2] != nz)
            {
              /**
               * check space for nodes
               */
              if(n_nodes >= n_alloc_nodes)
              {
                n_alloc_nodes += NODES_ALLOC_BUCKET;
                if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                  sizeof(struct node_t))) == NULL)
                {
                  ERROR("realloc()");
                  exit(EXIT_SUCCESS);
                }
              }

              n1 = nodes[iNode].l[0];
              n2 = nodes[iNode].l[1];
              n3 = (nz - nodes[iNode].l[2]) + nz;

              if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
              if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
              if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

              /**
               * fixity data
               */
              fixity = FREE_MASK;

              if(flag_fixity == 1)
              { 
                if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

                if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
                if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
              }
              else
              {                
                if (n1==0 || n1==nx || n1 == d_inc[n_levels-1] || n1 == nx - d_inc[n_levels-1] ||
                    n2==0 || n2==ny || n2 == d_inc[n_levels-1] || n2 == ny - d_inc[n_levels-1] ||
                    n3==0 || n3==nzValue || n3 == d_inc[n_levels-1] || n3 == nzValue - d_inc[n_levels-1]) 
                  fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
              }

              /**
               * surface data
               */
              if( n3 == 0 || n3 == nzValue)
                fixity |= SURFACE_MASK | SURFACE_Z_MASK;
              if( n2 == 0 || n2 == ny )
                fixity |= SURFACE_MASK | SURFACE_Y_MASK;
              if( n1 == 0 || n1 == nx )
                fixity |= SURFACE_MASK | SURFACE_X_MASK;


              nodes[n_nodes].l[0] = n1;
              nodes[n_nodes].l[1] = n2;
              nodes[n_nodes].l[2] = n3;
              nodes[n_nodes].fix_mask = fixity;

              n_nodes++;
            }
          }
        }
      } // end of if
      else
      {
        for(n3=0; n3 <= nz; n3++)
          for(n2=0; n2 <= ny; n2++)
            for(n1=0; n1 <= nx; n1++)
            {
              n_sites++;

              /**
               * 
               */
              if (n1 > lattice.l_end[0]) lattice.l_end[0] = n1;
              if (n2 > lattice.l_end[1]) lattice.l_end[1] = n2;
              if (n3 > lattice.l_end[2]) lattice.l_end[2] = n3;

              /**
               * fixity data
               */

              fixity = FREE_MASK;
          
              if( n3==0 ) fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
           
              if (n1==0 || n1==nx) fixity |= FIX_X_MASK | FIX_Y_MASK;
              if (n2==0 || n2==ny) fixity |= FIX_X_MASK | FIX_Y_MASK;
          
              /**
               * surface data
               */

              if( n3 == 0 || n3 == nz)
                fixity |= SURFACE_MASK | SURFACE_Z_MASK;
              if( n2 == 0 || n2 == ny )
                fixity |= SURFACE_MASK | SURFACE_Y_MASK;
              if( n1 == 0 || n1 == nx )
                fixity |= SURFACE_MASK | SURFACE_X_MASK;

              if( checkForNode(n1, n2, n3, nx, ny, nz, ax, ay, az, lattice_type) )
              {
                /**
                 * check space for nodes
                 */

                if(n_nodes >= n_alloc_nodes)
                {
                  n_alloc_nodes += NODES_ALLOC_BUCKET;
                  if((nodes = (struct node_t *) realloc(nodes, n_alloc_nodes*
                          sizeof(struct node_t))) == NULL)
                  {
                    ERROR("realloc()");
                    exit(EXIT_SUCCESS);
                  }
                }
      
                nodes[n_nodes].l[0] = n1;
                nodes[n_nodes].l[1] = n2;
                nodes[n_nodes].l[2] = n3;
                nodes[n_nodes].fix_mask = fixity;

                n_nodes++;
              }
            }
      } // end of else
    } // end of case BBC
    break;
  } // end of switch

  /**
   * fix all atoms/nodes that are between surface and next plane 
   * of nodes
   */
  /*
  fix_extra_nodes(atoms, n_atoms, nodes, n_nodes, lattice_constant);
  */

  if(position_flag == 1)
  {
    fixity = FREE_MASK;
    fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
    // printf("fixity= %d",fixity);
  }

  //
  //  define path_max to allocate memory for filenames
  //
  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);   

  for(int i=0 ; i< numLattice; i++)
  { 
    XDR   xdrs;
    //
    // create filename for ith quasi and open the file
    //
    char *local_data_filename;
    local_data_filename = (char *) malloc(path_max*sizeof(char));

    if(i == 0)
    {
      // use the defualt filename given by readDataFile()
      sprintf(&(local_data_filename[0]),"%s%s", data_filename, DATA_FILENAME_SUFFIX);
    }
    else
    {
      sprintf(&(local_data_filename[0]),"%s_%i%s", 
        data_filename, 
        i+1,
        DATA_FILENAME_SUFFIX);
    }

    //
    //  open file
    //
    if( (data_file = open_write_pipe(local_data_filename)) == NULL )
    {
      ERROR(local_data_filename);
      exit(EXIT_FAILURE);
    }    

    xdrstdio_create(&xdrs, data_file, XDR_ENCODE);

    printf("data filename = %s\n",local_data_filename);
    printf("material name = %s\n",material_name[i]);

    /**
     * print material name
     */
     // char name_name[256];
    if( xdr_string(&xdrs, &(material_name[i]), _POSIX_NAME_MAX) == FALSE )
      xdr_err("xdr_string()");
     // if( xdr_string(&xdrs, &(name_name[0]), _POSIX_NAME_MAX) == FALSE )
      // xdr_err("xdr_string()");

    /**
     * write lattice
     */

    name = lattice.name;
    if (xdr_string(&xdrs, &name, LATTICE_NAME_MAX) == FALSE)
      xdr_err("xdr_string()");

    if (xdr_enum(&xdrs, (enum_t *) &(lattice.type)) == FALSE)
      xdr_err("xdr_enum()");

    if (xdr_double(&xdrs, &lattice.a1[0]) == FALSE) xdr_err("xdr_double()");
    if (xdr_double(&xdrs, &lattice.a1[1]) == FALSE) xdr_err("xdr_double()");
    if (xdr_double(&xdrs, &lattice.a1[2]) == FALSE) xdr_err("xdr_double()");

    if (xdr_double(&xdrs, &lattice.a2[0]) == FALSE) xdr_err("xdr_double()");
    if (xdr_double(&xdrs, &lattice.a2[1]) == FALSE) xdr_err("xdr_double()");
    if (xdr_double(&xdrs, &lattice.a2[2]) == FALSE) xdr_err("xdr_double()");
    
    if (xdr_double(&xdrs, &lattice.a3[0]) == FALSE) xdr_err("xdr_double()");
    if (xdr_double(&xdrs, &lattice.a3[1]) == FALSE) xdr_err("xdr_double()");
    if (xdr_double(&xdrs, &lattice.a3[2]) == FALSE) xdr_err("xdr_double()");

    if (xdr_int(&xdrs, &lattice.l_start[0]) == FALSE) xdr_err("xdr_int()");
    if (xdr_int(&xdrs, &lattice.l_start[1]) == FALSE) xdr_err("xdr_int()");
    if (xdr_int(&xdrs, &lattice.l_start[2]) == FALSE) xdr_err("xdr_int()");

    if (xdr_int(&xdrs, &lattice.l_end[0]) == FALSE) xdr_err("xdr_int()");
    if (xdr_int(&xdrs, &lattice.l_end[1]) == FALSE) xdr_err("xdr_int()");
    if (xdr_int(&xdrs, &lattice.l_end[2]) == FALSE) xdr_err("xdr_int()");

    /**
     * write nodes
     */

    if (xdr_int(&xdrs, &n_nodes) == FALSE)
      xdr_err("xdr_int()");


    for (i_node=0; i_node < n_nodes; i_node++)
    {
      if (xdr_int(&xdrs, &nodes[i_node].l[0]) == FALSE) xdr_err("xdr_int()");
      if (xdr_int(&xdrs, &nodes[i_node].l[1]) == FALSE) xdr_err("xdr_int()");
      if (xdr_int(&xdrs, &nodes[i_node].l[2]) == FALSE) xdr_err("xdr_int()");

      if(position_flag == 1)
      {
        nodes[i_node].fix_mask = fixity;
        // nodes[i_node].fix_mask |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;
				printf("here position_flag = 1\n");
      }

      if (xdr_u_int(&xdrs, &nodes[i_node].fix_mask) == FALSE)
        xdr_err("xdr_u_int()");
    }

    // close the files
    xdr_destroy(&xdrs);
    fclose(data_file);
  }

  if (debugDataFlag == 1)
  {
    //
    //  generate simple file to anylize the output of this program
    //
    FILE * data_file_1;
    data_file_1 = fopen("dump.out","w");

    char local_data_filename[256];

    sprintf(&(local_data_filename[0]),"%s%s", data_filename, DATA_FILENAME_SUFFIX);  

    fprintf(data_file_1, "data_file                                           %s\n",
      local_data_filename);
    fprintf(data_file_1, "system_name                                           %s\n",
      system_name);
    fprintf(data_file_1, "mesh_flag                                           %d\n",
      mesh_flag);
    fprintf(data_file_1, "numLattice                                           %d\n",
      numLattice);
    fprintf(data_file_1, "a1                                           %f   %f   %f\n",
      lattice_constant_vec[0][0],
      lattice_constant_vec[0][1],
      lattice_constant_vec[0][2]);
    fprintf(data_file_1, "a2                                           %f   %f   %f\n",
      lattice_constant_vec[1][0],
      lattice_constant_vec[1][1],
      lattice_constant_vec[1][2]);
    fprintf(data_file_1, "a3                                           %f   %f   %f\n",
      lattice_constant_vec[2][0],
      lattice_constant_vec[2][1],
      lattice_constant_vec[2][2]);  
    fprintf(data_file_1, "lattice_type                                           %d\n",
      lattice_type);
    fprintf(data_file_1, "atomistic                                           %d  %d  %d\n",
      atomistic_size[0],
      atomistic_size[1],
      atomistic_size[2]);
    fprintf(data_file_1, "full                                           %d  %d  %d\n",
      sample_size[0],
      sample_size[1],
      sample_size[2]);
    fprintf(data_file_1, "mesh_levels                                           %d\n",
      n_levels);
    for(int i=0; i < n_levels; i++)
    {
      fprintf(data_file_1, "n_level                                           %d  %d\n",
        i,
        d_inc[i]);
    }      
    fprintf(data_file_1, "problem_size                                           %d\n",
      problem_size);
    for(int i=0; i < numLattice; i++)
    {
      fprintf(data_file_1, "material                                           %d   %s\n",
        i,
        material_name[i]);
    }

    fprintf(data_file_1, "+++++++++++++++++++++++++++++++++++++++++++++++\n");
    fprintf(data_file_1, "++++++++++++++++ output data ++++++++++++++++++\n");
    fprintf(data_file_1, "+++++++++++++++++++++++++++++++++++++++++++++++\n");
    
    fprintf(data_file_1, "material    ");
    for(int i=0; i < numLattice; i++)
    {
      fprintf(data_file_1, "(%d, %s)   ",
        i,
        material_name[i]);
    }
    fprintf(data_file_1,"\n");

    name = lattice.name;
    fprintf(data_file_1, "lattice name    %s\n",
      name);
    if(lattice.type == FCC )
      fprintf(data_file_1, "lattice type    FCC\n");
    else
      fprintf(data_file_1, "lattice type    BCC\n");

    fprintf(data_file_1, "lattice constants : \n");
    fprintf(data_file_1, "a1       %f  %f  %f\n",
      lattice.a1[0],
      lattice.a1[1],
      lattice.a1[2]);
    fprintf(data_file_1, "a2       %f  %f  %f\n",
      lattice.a2[0],
      lattice.a2[1],
      lattice.a2[2]);
    fprintf(data_file_1, "a3       %f  %f  %f\n",
      lattice.a3[0],
      lattice.a3[1],
      lattice.a3[2]);  
    fprintf(data_file_1, "l_start       %d  %d  %d\n",
      lattice.l_start[0],
      lattice.l_start[1],
      lattice.l_start[2]);    
    fprintf(data_file_1, "l_end       %d  %d  %d\n",
      lattice.l_end[0],
      lattice.l_end[1],
      lattice.l_end[2]);

    /**
     * write nodes
     */
    fprintf(data_file_1,"+++++++++++++++++++  node data  ++++++++++++++++++\n");
    fprintf(data_file_1, "number nodes    %d\n", n_nodes);

    for (i_node=0; i_node < n_nodes; i_node++)
    {
      fprintf(data_file_1, "%d  %d  %d\n",
        nodes[i_node].l[0],
        nodes[i_node].l[1],
        nodes[i_node].l[2]);
    }

    // close the data_file
    fclose(data_file_1);

    //
    // write node data alone in other file for matlab/visit
    //
    char *data_filename_2;
    data_filename_2 = (char *) malloc(path_max*sizeof(char));

    sprintf(data_filename_2, "node%s", NODE_FILENAME_SUFFIX);

    FILE * data_file_2;
    // data_file_2 = fopen("node.plt","w");
    data_file_2 = open_write_pipe(data_filename_2);

    fprintf(data_file_2, " %s = \"%s\"\n", tecplot_header_fields[TITLE], "QC FEM");
    fprintf(data_file_2, " %s = \"X\" \"Y\" \"Z\" \"fixity\" \"fixity_bin\"\n", 
      tecplot_header_fields[VARIABLES]);

    fprintf(data_file_2, " %s N= %d, E= %d, F=%s, ET=%s\n", 
      tecplot_header_fields[ZONE],
      n_nodes, 
      n_nodes,
      tecplot_header_fields[FEPOINT], 
      tecplot_header_fields[TETRAHEDRON]);


    for (i_node=0; i_node < n_nodes; i_node++)
    {
      fprintf(data_file_2, "%f  %f  %f  %d  %4s\n",
        (float) nodes[i_node].l[0],
        (float) nodes[i_node].l[1],
        (float) nodes[i_node].l[2],
        nodes[i_node].fix_mask,
        xyz_fixity(nodes[i_node].fix_mask));
    }

    // create dummy connectivity.
    fprintf(data_file_2, "\n");
    for (int i_elem = 0; i_elem < n_nodes; i_elem++)
    {
      fprintf(data_file_2, "%d %d %d %d\n", 
              i_elem+0,
              i_elem+1,
              i_elem+2,
              i_elem+3);
    }    

    fclose(data_file_2);  
    free(data_filename_2);  

    //
    //  dump the crystal structure information
    //
    FILE * data_file_3;
    data_file_3 = fopen("crystal_structure.out","w");

    int atom_count = 0;
    for(int nn1 =0; nn1 < 2; nn1++)
      for(int nn2 =0; nn2 < 2; nn2++)
        for(int nn3=0; nn3 < 2; nn3++)
        {
          // check if (nn1,nn2,nn3) are lattice sites
          //  atom number, atom coordinates, species number of atom
          for(int i=0; i < numLattice; i++)
          {
            if(nn1 == 0 && nn2 == 0 && nn3 == 0)
            {
            // if(i != 4 || i != 5)
            // {
              fprintf(data_file_3,"%d   %f  %f   %f   %d\n",
                atom_count,
                lattice.a1[0]*nn1 + lattice.a2[0]*nn2 + lattice.a3[0]*nn3 + shiftVector[i][0],
                lattice.a1[1]*nn1 + lattice.a2[1]*nn2 + lattice.a3[1]*nn3 + shiftVector[i][1],
                lattice.a1[2]*nn1 + lattice.a2[2]*nn2 + lattice.a3[2]*nn3 + shiftVector[i][2],
                i);
              atom_count += 1;
            // }
            }
            else
            {
              fprintf(data_file_3,"%d   %f  %f   %f   %d\n",
                atom_count,
                lattice.a1[0]*nn1 + lattice.a2[0]*nn2 + lattice.a3[0]*nn3 + shiftVector[0][0],
                lattice.a1[1]*nn1 + lattice.a2[1]*nn2 + lattice.a3[1]*nn3 + shiftVector[0][1],
                lattice.a1[2]*nn1 + lattice.a2[2]*nn2 + lattice.a3[2]*nn3 + shiftVector[0][2],
                0);
              atom_count += 1;
            }            
          }
        }
    fclose(data_file_3); 
  }
  printf("lattice a1 = (%f, %f, %f)\n", lattice.a1[0],lattice.a1[1],lattice.a1[2]);
  printf("lattice a2 = (%f, %f, %f)\n", lattice.a2[0],lattice.a2[1],lattice.a2[2]);
  printf("lattice a3 = (%f, %f, %f)\n", lattice.a3[0],lattice.a3[1],lattice.a3[2]);
  printf("lattice l_start = (%d, %d, %d)\n", lattice.l_start[0],lattice.l_start[1],lattice.l_start[2]);
  printf("lattice l_end = (%d, %d, %d)\n", lattice.l_end[0],lattice.l_end[1],lattice.l_end[2]);  
  for(int i=0; i < numLattice; i++)
    printf("shift_%d = (%f, %f, %f)\n", i, shiftVector[i][0],shiftVector[i][1],shiftVector[i][2]);

  /**
   *
   */
  printf("Number sites %d number nodes %d\n", n_sites, n_nodes);

  /**
   * close
   */
  for(int i=0; i < MAX_LATTICE_NUMBER; i++)
    free(material_name[i]);
  
  exit(EXIT_SUCCESS);
  return(EXIT_SUCCESS);
}
