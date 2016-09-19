//
// Indent.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#if !defined(_REENTRANT)
#define _REENTRANT
#endif
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#ifdef STDC_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* Added for strstr */
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#ifdef HAVE_CTYPE_H
#include <ctype.h>
#else
#error ctype.h not found
#endif /* HAVE_CTYPE_H */

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

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#else
#error fcntl.h not found
#endif /* HAVE_FCNTL_H */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found
#endif /* HAVE_SYS_TYPES_H */

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#else
#error sys/stat.h not found.
#endif /* HAVE_SYS_STAT_H */

#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#else
#error sys/mman.h not found
#endif /* HAVE_SYS_MMAN_H */

#include "DataTypes.h"
#include "Error.h"
#include "Indent.h"
#include "Node.h"
#include "Quasicontinua.h"

#include "monitor.h"
#include "threads.h"

// global variables
#define INDENTOR_DISPL -0.1
#define DEFAULT_INDENTOR_INIT_FILE_NAME "indentor.ini"

static const char *indentor_init_keywords[] = {

#define RADIUS 0
#define RADIUS_FLAG 0x1
    "radius",

#define CONSTANT 1
#define CONSTANT_FLAG 0x2
    "constant",

#define DISPLACEMENT_STEP 2
#define DISPLACEMENT_STEP_FLAG 0x4
    "displacement",

#define POSITION 3
#define POSITION_FLAG 0x8
    "position"

#define EMPTY_FLAG 0x0
};

//
//
//

namespace quasicontinuum {

//
//  namespace for indentorNidesInteractionWorker()
//
namespace {
struct indentorNodesInteractionData_t {
  struct node_list_t *P_node_list;
  struct indentor_t *P_indentor;
  int iQuasi;
};

//
//  indentorNodePotentialLocal()
//
void indentorNodePotentialLocal(double *force, double *energy, const double r,
                                const struct indentor_t *P_indentor) {
  double dr;

  dr = P_indentor->radius - r;

  *energy = P_indentor->constant * dr * dr * dr;
  *force = -3.0 * P_indentor->constant * dr * dr;

  return;
} // end of indentorNodePotentialLocal()

//
//  indentorNodesInteractionWorker()
//
void *indentorNodesInteractionWorker(void *arg) {

  struct node_list_t *P_node_list =
      ((struct indentorNodesInteractionData_t *)arg)->P_node_list;
  struct indentor_t *P_indentor =
      ((struct indentorNodesInteractionData_t *)arg)->P_indentor;
  const int iQuasi = ((struct indentorNodesInteractionData_t *)arg)->iQuasi;

  double indentor_force[3] = {0, 0, 0};
  double r_ij[3];
  double r_cut_sqr = P_indentor->radius * P_indentor->radius;
  double r;
  double energy;
  double indenter_energy = 0.0;
  double force;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_node_thread;
  int number_node_start;
  int number_node_end;

  int i_node;

  /**
    * get my share of data
    */

  get_share(my_id, number_threads, P_node_list->number_nodes,
            &number_node_thread, &number_node_start, &number_node_end);

  /**
    * vector to hold relative shift
    */
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  /**
    * loop over my share of nodes
    */

  for (i_node = number_node_start; i_node <= number_node_end; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];
    // get difference of node and indentor position

    r_ij[0] = P_node->position[0] - P_indentor->position[0];
    r_ij[1] = P_node->position[1] - P_indentor->position[1];
    r_ij[2] = P_node->position[2] - P_indentor->position[2];

    /**
      * update r_ij for global shift
      */
    r_ij[0] += shift[0];
    r_ij[1] += shift[1];
    r_ij[2] += shift[2];

    /**
      * calculate r^2
      */
    r = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

    /**
      * compute force if node is inside
      */

    if (r < r_cut_sqr) {
      r = sqrt(r);

      indentorNodePotentialLocal(&force, &energy, r, P_indentor);

      force = -force / r;

      force *= P_node->weight;
      /**
        * update force on node
        */

      P_node->acceleration[0] += r_ij[0] * force;
      P_node->acceleration[1] += r_ij[1] * force;
      P_node->acceleration[2] += r_ij[2] * force;

      /**
        * update force on indentor
        */

      indentor_force[0] -= r_ij[0] * force;
      indentor_force[1] -= r_ij[1] * force;
      indentor_force[2] -= r_ij[2] * force;

      /**
        * save energy
        */

      indenter_energy += energy;
    }
  }

  /**
    * update total indentor force
    */

  pthread_mutex_lock(&P_indentor->lock);

  P_indentor->force[0] += indentor_force[0];
  P_indentor->force[1] += indentor_force[1];
  P_indentor->force[2] += indentor_force[2];

  P_indentor->energy.potential += indenter_energy;

  pthread_mutex_unlock(&P_indentor->lock);

  return ((void *)NULL);
} // end of indentorNodesInteractionWorker()
} // end of namespace for indentorNodesInteractionWorker()

//
//
//

Indent *Indent::_instance = NULL;

//
// constructor
//

Indent::Indent() {

  //
  //
  //
  return;
}

//
// constructor
//

Indent::Indent(const int indent_enable) {
  d_enableNanoIndentation = indent_enable;

  //
  //
  //
  return;
}

//
// destructor
//

Indent::~Indent() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Indent *Indent::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Indent(0);
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Indent::destroyInstance() {

  //
  // delete instance
  //
  delete _instance;

  //
  //
  //
  return;
}

//
//  enableNanoIndentation()
//
//  0 : no indentation
//  1 : indentation is active
//
void Indent::setNanoIndentation(const int enable_nano_indentation_flag) {
  d_enableNanoIndentation = enable_nano_indentation_flag;

  return;
}

//
//  isIndentEnable()
//
//  0 : not enebaled
//  1 : enabled
//
int Indent::isIndentEnable() { return d_enableNanoIndentation; }

//
//  moveIndentor()
//
void Indent::moveIndentor(struct indentor_t *P_indentor) {
  P_indentor->position[0] += P_indentor->displacement_incr[0];
  P_indentor->position[1] += P_indentor->displacement_incr[1];
  P_indentor->position[2] += P_indentor->displacement_incr[2];

  return;
}

//
//  indentorNodesInteraction()
//
void Indent::indentorNodesInteraction(int iQuasi) {
  //
  //  get data of iQuasi
  //
  Quasicontinuum iQuasicontinuum =
      Quasicontinua::getInstance()->getQuasi(iQuasi);

  struct node_list_t node_list = iQuasicontinuum.getNodeList().node_list;

  struct indentor_t indentor = iQuasicontinuum.getIndentor();

  enum mt_version_t mt_version = SINGLE_THREADED;

  if (get_max_number_threads() > 1)
    mt_version = MULTI_THREADED;

  // data for thread worker
  struct indentorNodesInteractionData_t indentorNodesInteractionData;

  /**
    * initialize indenter data
    */

  indentor.force[0] = 0.0;
  indentor.force[1] = 0.0;
  indentor.force[2] = 0.0;

  indentor.energy.potential = 0.0;

  /**
    *
    */

  switch (MULTI_THREADED) {
  case MULTI_THREADED:

    indentorNodesInteractionData.P_node_list = &node_list;
    indentorNodesInteractionData.P_indentor = &indentor;
    indentorNodesInteractionData.iQuasi = iQuasi;

    // call thread_monitor() of monitor.c
    thread_monitor(indentorNodesInteractionWorker,
                   (void *)&indentorNodesInteractionData,
                   get_max_number_threads());
    break;

  case SINGLE_THREADED:

    fprintf(stderr, "[%s, %d] Feature not yet implemented\n", __FILE__,
            __LINE__);
    exit(EXIT_FAILURE);
    break;
  }

  return;
}

//
//  findIndentorContactRadius()
//
//  find contact zone between the indenter and the sample
//  arguments:
//  *P_R - pointer to a location where the radius of the contact zone
//         in the initial configuration is stored
//  *P_r - pointre to a location where the radius of the contact zone
//         in the current configuration is stored
//
void Indent::findIndentorContactRadius(double *P_R, double *P_r,
                                       const struct node_list_t *P_node_list,
                                       const struct indentor_t *P_indentor,
                                       const int iQuasi,
                                       enum mt_version_t mt_version) {
  /**
    * this interface is not MT, as the routine is called infrequently
    */

  double r_ij[3];
  double r;
  double r_cut_sqr = P_indentor->radius * P_indentor->radius;

  int i_node;

  /**
    * if both P_R and P_r ar NULL there is nothing to do
    */

  if ((P_R == NULL) && (P_r == NULL))
    return;

  /**
    * initialize radii
    */

  if (P_R)
    *P_R = 0.0;
  if (P_r)
    *P_r = 0.0;

  /**
    * vector to hold relative shift
    */
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  /**
    * scan all nodes
    */

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

    r_ij[0] = P_node->position[0] - P_indentor->position[0];
    r_ij[1] = P_node->position[1] - P_indentor->position[1];
    r_ij[2] = P_node->position[2] - P_indentor->position[2];

    /**
      * update r_ij for global shift
      */
    r_ij[0] += shift[0];
    r_ij[1] += shift[1];
    r_ij[2] += shift[2];

    /**
      * calculate r^2
      */
    r = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];

    if (r < r_cut_sqr) {

      if (P_R) {
        r_ij[0] = P_node->initial_position[0] - P_indentor->initial_position[0];
        r_ij[1] = P_node->initial_position[1] - P_indentor->initial_position[1];

        r = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1];

        if (r > (*P_R))
          *P_R = r;
      }

      if (P_r) {
        r_ij[0] = P_node->position[0] - P_indentor->position[0];
        r_ij[1] = P_node->position[1] - P_indentor->position[1];

        r = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1];

        if (r > (*P_r))
          *P_r = r;
      }
    }
  }

  /**
    *
    */

  if (P_R)
    *P_R = sqrt(*P_R);
  if (P_r)
    *P_r = sqrt(*P_r);

  return;
}

//
//  fixIndentBoundary()
//
void Indent::fixIndentBoundary(struct node_list_t *P_node_list,
                               struct lattice_t *P_lattice) {
  int i_node;

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

    /**
      * fix bottom
      */

    if (P_node->l[2] == P_lattice->l_start[2])
      P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

    /**
     * fix planes with normals in x
     */

    if ((P_node->l[0] == P_lattice->l_start[0]) ||
        (P_node->l[0] == P_lattice->l_end[0]))
      P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK;

    /**
     * fix planes with normals in y
     */

    if ((P_node->l[1] == P_lattice->l_start[1]) ||
        (P_node->l[1] == P_lattice->l_end[1]))
      P_node->fix_mask |= FIX_X_MASK | FIX_Y_MASK;

    /**
     * set surface data
     */

    if ((P_node->l[0] == P_lattice->l_start[0]) ||
        (P_node->l[0] == P_lattice->l_end[0]))
      P_node->fix_mask |= SURFACE_MASK | SURFACE_X_MASK;

    if ((P_node->l[1] == P_lattice->l_start[1]) ||
        (P_node->l[1] == P_lattice->l_end[1]))
      P_node->fix_mask |= SURFACE_MASK | SURFACE_Y_MASK;

    if ((P_node->l[2] == P_lattice->l_start[2]) ||
        (P_node->l[2] == P_lattice->l_end[2]))
      P_node->fix_mask |= SURFACE_MASK | SURFACE_Z_MASK;
  }

  return;
}
}