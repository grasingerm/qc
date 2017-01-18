/**
 * Provide map site -> element
 *
 * $Log: site_element_map.c,v $
 * Revision 1.3  2002/03/07 23:53:09  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.2  2001/10/03 22:29:52  knap
 * Checked in the latest changes.
 *
 * Revision 1.1  2000/04/12 22:07:21  knap
 * Removed references to atoms.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "C_Interface.h"
#include "DataTypes.h"
#include "Error.h"

#include "range_search_generic.h"
#include "site_element_map.h"
#include "site_element_map_cache.h"

#if !defined(lint)
static const char rcsid[] =
    "$Id: site_element_map.c,v 1.3 2002/03/07 23:53:09 knap Exp $";
#endif /* !lint */

/**
 * local data
 */

static struct site_element_data_t {
  struct BSP_tree_t *BSP_tree; /* BSP tree used in mapping                   */
  double h_min;                /* min element size                           */
  double h_max;                /* max element size                           */
#define INITIAL_WINDOW_SIZE 4.0
#define MAX_WINDOW_SIZE 50.0
} * site_element_data;

static int site_element_allocated = 0;
static const size_t site_element_size = sizeof(struct site_element_data_t);

/**
 * given an element return its position == barycenter coordinates
 */

static struct point_t get_element_position(const void *object) {

  struct point_t point;
  struct element_t *P_element = *((struct element_t **)object);

  point.x = P_element->center[0];
  point.y = P_element->center[1];
  point.z = P_element->center[2];

  return (point);
}

/**
 * create BSP tree and place all elements
 */

static struct BSP_tree_t *
create_element_BSP_tree(struct element_list_t *P_element_list) {

  return (BSP_tree_create((void *)P_element_list->elements,
                          (size_t)P_element_list->number_elements,
                          sizeof(struct element_t *), get_element_position));
}

/**
 * initialize site -> element map data
 */

void initialize_site_element_map_data(const int iQuasi,
                                      struct element_list_t *P_element_list) {
  /**
   * check if currentId has been allocated yet
   */
  if (iQuasi >= site_element_allocated) {

    /**
     * check if this is first allocation
     */
    if (iQuasi == 0) {

      /**
       * malloc size of site_element_data
       */
      site_element_data = malloc(site_element_size);

    } else {

      /**
       * realloc size of site_element_data
       */
      site_element_data =
          realloc(site_element_data, (iQuasi + 1) * site_element_size);
    }

    /**
     * check to make sure site_element_data allocated
     */
    if (site_element_data == NULL) {
      PERROR_MALLOC();
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    /**
     * increment site_element_allocated
     */
    ++site_element_allocated;
  }

  int i_elem;

  /**
   * build element tree
   */

  site_element_data[iQuasi].BSP_tree = create_element_BSP_tree(P_element_list);

  /**
   * locate smallest and largest element size
   */

  site_element_data[iQuasi].h_min = P_element_list->elements[0]->cir_radius;
  site_element_data[iQuasi].h_max = P_element_list->elements[0]->cir_radius;

  for (i_elem = 0; i_elem < P_element_list->number_elements; i_elem++) {

    struct element_t *P_element = P_element_list->elements[i_elem];

    if (P_element->cir_radius < site_element_data[iQuasi].h_min)
      site_element_data[iQuasi].h_min = P_element->cir_radius;

    if (P_element->cir_radius > site_element_data[iQuasi].h_max)
      site_element_data[iQuasi].h_max = P_element->cir_radius;
  }

  site_element_data[iQuasi].h_min = sqrt(site_element_data[iQuasi].h_min);
  site_element_data[iQuasi].h_max = sqrt(site_element_data[iQuasi].h_max);

  return;
}

/**
 * clean site_element_data
 */

void clean_site_element_map_data(const int iQuasi) {

  BSP_tree_destroy(site_element_data[iQuasi].BSP_tree);

  return;
}

/**
 * process object list
 */

#define ELEMENT_RADIUS_SCALE_COEFF 1.1

static struct element_t *
process_object_list(struct object_list_t *P_object_list, const double site[3]) {

  double r[3];
  double dr;

  int i_object;

  /**
   * if object list is empty return here
   */

  if (P_object_list->n_objects == 0)
    return (NULL);

  /**
   * loop over all objects
   */

  for (i_object = 0; i_object < P_object_list->n_objects; i_object++) {

    struct element_t *P_element =
        *((struct element_t **)P_object_list->objects[i_object]);

    /**
     * check first if the site is within sphere centered at the
     * barycentric coordinates of the element
     */

    r[0] = P_element->center[0] - site[0];
    r[1] = P_element->center[1] - site[1];
    r[2] = P_element->center[2] - site[2];

    dr = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

    if (dr <= P_element->cir_radius * ELEMENT_RADIUS_SCALE_COEFF) {

      if (checkPointInTetra(P_element->node[0]->initial_position,
                            P_element->node[1]->initial_position,
                            P_element->node[2]->initial_position,
                            P_element->node[3]->initial_position, site) == 1)
        return (P_element);
    }
  }

  return (NULL);
}

/**
 * given a site locate the element it is in
 * this does not depend on current id
 */

struct element_t *locate_site_element(double *shape, const int l[3],
                                      struct lattice_t *P_lattice,
                                      const int iQuasi) {

  /* struct object_list_t object_list = OBJECT_LIST_INITIALIZER; */
  struct object_list_t *P_object_list;
  struct element_t *P_element = NULL;
  struct rectangle_t rectangle;

  // iQuasi is currentId
  double r[3];
  double h = site_element_data[iQuasi].h_min * INITIAL_WINDOW_SIZE;

  int i_iter = 1;

  /**
   * check cache first
   */

  if ((P_element = hash_table_search_element(shape, l, iQuasi)) != NULL)
    return (P_element);

  /**
   * get coordinates of the site
   */

  getSiteInitialPosition(r, l, P_lattice, iQuasi);

  /**
   * loop until found
   */

  while (P_element == NULL) {

    /**
     * initialize window
     */

    rectangle.point1.x = r[0] - h;
    rectangle.point1.y = r[1] - h;
    rectangle.point1.z = r[2] - h;

    rectangle.point2.x = r[0] + h;
    rectangle.point2.y = r[1] + h;
    rectangle.point2.z = r[2] + h;

    /**
     * reinitialize number of objects on the list
     */

    /* object_list.n_objects = 0; */

    /**
     * locate all elements inside
     */

    P_object_list = range_search(site_element_data[iQuasi].BSP_tree, &rectangle,
                                 get_element_position);
    /* range_search(&object_list, site_element_data.BSP_tree, &rectangle,
       get_element_position); */

    /**
     * check all elements on the list
     */

    if ((P_element = process_object_list(P_object_list, r)) != NULL)
      break;

    /**
     * if element has not been found increase window size
     */

    h *= exp(i_iter);

    /**
     * some sanity checks; if h > alpha*h_max sth is probably wrong
     */

    if (h > MAX_WINDOW_SIZE * site_element_data[iQuasi].h_max) {
      fprintf(stderr, "No element found but h has grown to %8.6e\n", h);
      abort();
      exit(EXIT_FAILURE);
    }

    /**
     * increase iteration count
     */

    i_iter++;
  }

  /**
   * cleanup
   */

  /*  free(object_list.objects); */

  /**
   * add element to the cache
   */

  hash_table_add_site(shape, l, P_element, P_lattice, iQuasi);

  /**
   *
   */

  return (P_element);
}

/**
 * given site and element list find element containing it
 */

struct element_t *find_site_element(const int l[3], struct lattice_t *P_lattice,
                                    const struct element_list_t *P_element_list,
                                    const int iQuasi) {

  int i_element;

  for (i_element = 0; i_element < P_element_list->number_elements;
       ++i_element) {

    struct element_t *P_element = P_element_list->elements[i_element];
    double r[3];

    /**
     * compute position of the site
     */

    getSiteInitialPosition(r, l, P_lattice, iQuasi);

    /**
     * check the element
     */

    if (checkPointInTetra(P_element->node[0]->initial_position,
                          P_element->node[1]->initial_position,
                          P_element->node[2]->initial_position,
                          P_element->node[3]->initial_position, r) == 1)
      return P_element;
  }

  return NULL;
}

#if defined(_QC_SITE_ELEMENT_MAP_DEBUG_)

/**
 * take all lattice sites and make sure that they really belong to elements
 * returned by site_element_map
 */

void site_element_map_debug(struct lattice_t *P_lattice, const int iQuasi) {

  struct element_t *P_element;
  double r[3];
  int l[3];

  /**
   * generate lattice
   */

  for (l[0] = P_lattice->l_start[0]; l[0] <= P_lattice->l_end[0]; l[0]++)
    for (l[1] = P_lattice->l_start[1]; l[1] <= P_lattice->l_end[1]; l[1]++)
      for (l[2] = P_lattice->l_start[2]; l[2] <= P_lattice->l_end[2]; l[2]++) {

        /**
         * check if the (i,j,k) result in a site
         */

        if (isLatticeSite(l, P_lattice, iQuasi) == RETURN_FAILURE)
          continue;

        /**
         * find element the site is in
         */

        P_element = locate_site_element(NULL, l, P_lattice, iQuasi);

        /**
         * check if the site really is inside P_element
         */

        getSiteInitialPosition(r, l, P_lattice, iQuasi);

        if (checkPointInTetra(P_element->node[0]->initial_position,
                              P_element->node[1]->initial_position,
                              P_element->node[2]->initial_position,
                              P_element->node[3]->initial_position, r) == 0) {
          printf("locate_site_element() shows that site (%d,%d,%d) is "
                 "inside %d, but check_tetra() fails...\n",
                 l[0], l[1], l[2], P_element->number);
        }
      }

  return;
}

#endif /* _QC_SITE_ELEMENT_MAP_DEBUG_ */
