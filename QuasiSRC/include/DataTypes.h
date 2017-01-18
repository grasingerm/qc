#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdbool.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>

#if defined(__cplusplus)
namespace quasicontinuum {
extern "C" {
#endif /* __cplusplus */

/**
 * defines
 */

#define RETURN_SUCCESS EXIT_SUCCESS
#define RETURN_FAILURE EXIT_FAILURE

/* boundary condition mask  */
#define FREE_MASK 0x000
#define FIX_X_MASK 0x001
#define FIX_Y_MASK 0x002
#define FIX_Z_MASK 0x004

/* surface mask  */
#define SURFACE_MASK 0x020
#define SURFACE_X_MASK 0x040
#define SURFACE_Y_MASK 0x080
#define SURFACE_Z_MASK 0x100

#define LOAD_MASK 0x008
#define NEW_MASK 0x010

/* crack masks */
#define CRACK_TOP_MASK 0x1000
#define CRACK_BOT_MASK 0x2000

/* void masks  */
#define VOID_SURFACE_MASK 0x4000
#define VOID_MASK 0x8000

/* new mask for frequency */
#define FIX_W_MASK 0x001

/**
 * util
 */

enum loc_t { INSIDE, OUTSIDE };
enum flag_t { DOWN, UP };
enum found_t { FOUND, NOT_FOUND };

struct energy_t {
  double kinetic;
  double potential;
  double total;
#define ENERGY_INITIALIZER                                                     \
  { 0, 0, 0 }
};

// new data types added from other c files

// fcc_neighbor_shells.
struct shell_t {
  int (*site)[3];
  int number_sites;
#define SHELL_INITIALIZER                                                      \
  { NULL, 0 }
};

// indentor.h
struct indentor_t {
  double radius;
  double constant;
  double displacement_incr[3];
  double initial_position[3];
  double position[3];
  double force[3];
  struct energy_t energy;
  pthread_mutex_t lock;
#define INDENTOR_INITIALIZER                                                   \
  {                                                                            \
    0, 0, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},               \
        PTHREAD_MUTEX_INITIALIZER                                              \
  }
};

// monitor.h
enum mt_version_t { SINGLE_THREADED, MULTI_THREADED };

// internal_void.h
enum void_t { CIRCLE, OVAL, VPIT };

// used in CreateMesh.cc and Element.cc
struct mesh_data_t {
  int (*tets)[4];
  int number_tets;
  double eta_min;
  double rho_min;
  double sigma_min;
#define MESHDATA_INITIALIZER                                                   \
  { NULL, 0, 0.0, 0.0, 0.0 }
};

// boundingbox.h
struct boundingbox_t {
  double vertex_1[3];
  double vertex_2[3];
};

/**
 * basic types
 */

struct quad_point_t {
  union {
    struct atom_t *P_atom;
    struct node_t *P_node;
  } P_point;
  double weight;
};

enum switch_t { OFF, ON }; /* ON/OFF switch */

struct qc_options_t {
  enum switch_t restart;
#define QC_OPTIONS_INITIALIZER                                                 \
  { OFF }
};

struct site_list_t {
  int (*sites)[3];
  int number_sites;
  int alloc;
  pthread_mutex_t lock;
#define SITE_LIST_INITIALIZER                                                  \
  { NULL, 0, 0, PTHREAD_MUTEX_INITIALIZER }
};

/**
 * element types
 */

#define ELEM_CLEAN_MASK 0x00
#define ELEM_REMESH_MASK 0x01

struct element_t {
  unsigned int mask;
  int number;             /* element number                        */
  struct node_t *node[4]; /* list of all nodes                     */
  struct quad_point_t quad_points[6]; /* pointers to quad atoms             */
  double center[3];    /* center of a tetrahedron               */
  double cir_radius;   /* radius of a circumsphere              */
  double radii_ratio;  /* ratio of insphere to circumsphere     */
  double strain[6];    /* xx, yy, zz, xy, xz, yz                */
  double F[3][3];      /* deformation gradient                  */
  double number_sites; /* estimated number of sites inside      */
  enum sum_t { EXPLICIT, QUAD } sum_type; /*                                  */
  pthread_mutex_t lock;
};

/**
 * Node list
 */

struct node_list_t {
  struct node_t **nodes;
  struct node_t **inactive_nodes;
  int number_nodes;
  int number_inactive_nodes;
  int alloc_space;
  int alloc_inactive_space;
  pthread_mutex_t lock;
#define NODE_LIST_INITIALIZER                                                  \
  { NULL, NULL, 0, 0, 0, 0, PTHREAD_MUTEX_INITIALIZER }
};

/**
 * Element list
 */

struct element_list_t {
  struct element_t **elements;
  struct element_t **inactive_elements;
  int number_elements;
  int number_inactive_elements;
  int alloc_space;
  int alloc_inactive_space;
  pthread_mutex_t lock;
#define ELEMENT_LIST_INITIALIZER                                               \
  { NULL, NULL, 0, 0, 0, 0, PTHREAD_MUTEX_INITIALIZER }
};

/**
 * Node types
 */

struct node_t {
  int number;                     /* node number                      */
  unsigned int fix_mask;          /*3 bit field (z,y,x) 0-free 1-fixed*/
  unsigned int void_mask;         /* 0-outside 1-inside               */
  unsigned int void_surface_mask; /* 0-not on surface 1-on surface */
  unsigned int fix_w_mask;        /*1 bit mask for w 0-free 1-fixed*/
  double weight;                  /* nodal weight                     */
  double position[3];             /* nodal position                   */
  double temperature;             /* temperature dof of node          */
  double frequency;               /* frequency of node                */
  double tau; /* std dev from position of node                */
  double initial_position[3];
  double initial_temperature;
  double initial_frequency;
  double initial_tau;
  double velocity[3];     /* nodal velocity                   */
  double acceleration[3]; /* nodal acceleration               */
  double force_frequency; /* force due to frequency dof       */
  double body_force[3];   /* body force                       */
  double external_force[3];
  double external_freq_force;
  struct energy_t energy; /* energy                           */
  double density;
  int l[3];
  double *mass;
  struct site_list_t site_cluster;
  // struct site_list_t   *neighbor_list;  /* not using it */
  struct node_list_t node_list;
  struct element_list_t element_list; /* list of adjacent elements        */
  int bin_number;
  pthread_mutex_t lock;
};

/*
 * structure for all nodes (new and old)
 */

struct all_node_list_t {
  struct node_list_t node_list;
  struct node_list_t new_node_list;
  struct energy_t energy;
  struct energy_t initial_energy;
#define ALL_NODE_LIST_INITIALIZER                                              \
  {                                                                            \
    NODE_LIST_INITIALIZER, NODE_LIST_INITIALIZER, ENERGY_INITIALIZER,          \
        ENERGY_INITIALIZER                                                     \
  }
};

/**
 * lattice type
 */
enum lattice_type_t { FCC, BCC };
struct lattice_t {
#define LATTICE_NAME_MAX 128
#define LATTICE_INITIALIZER                                                    \
  {                                                                            \
    "", FCC, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, { 0, 0, 0 }           \
  }
  char name[LATTICE_NAME_MAX + 1];
  enum lattice_type_t type;
  double a1[3];
  double a2[3];
  double a3[3];
  int l_start[3];
  int l_end[3];
};

/**
 * face types
 */

struct face_t {
  struct node_t *node[3];
};

struct face_list_t {
  struct face_t **faces;
  int number_faces;
  int number_alloc_faces;
  pthread_mutex_t lock;
};

/**
 * generic object list
 */

struct object_list_t {
#define OBJECT_LIST_INITIALIZER                                                \
  { NULL, 0, 0 }
  void **objects;
  int n_objects;
  int alloc_space;
};

/**
 * periodic data type
 */

struct periodic_t {
#define NO_PERIODIC_MASK 0x0
#define PERIODIC_X_MASK 0x1
#define PERIODIC_Y_MASK 0x2
#define PERIODIC_Z_MASK 0x4
  unsigned mask;
  double box_size[3];
  double half_box_size[3];
};

#define PERIODIC_INITIALIZER                                                   \
  {                                                                            \
    NO_PERIODIC_MASK, {0, 0, 0}, { 0, 0, 0 }                                   \
  }

/**
 * input data
 */

struct options_t {
  char init_file_name[PATH_MAX];
  char data_file_name[PATH_MAX];
  int number_steps;
  enum { MD, CG } minimizer;
};

#define OPTIONS_INITIALIZER                                                    \
  { "", "", -1, MD }

/**
 *
 */

#if defined(__cplusplus)
}
}
#endif /* __cplusplus */

#endif /* _TYPES_H_ */
