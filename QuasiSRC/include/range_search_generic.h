/**
 * $Id: range_search_generic.h,v 1.2 2003/08/14 20:59:29 knap Exp $
 *
 * $Log: range_search_generic.h,v $
 * Revision 1.2  2003/08/14 20:59:29  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.1  2000/04/12 22:06:18  knap
 * Removed references to atoms.
 *
 */
#ifndef _QC_RANGE_SEARCH_GENERIC_H_
#define _QC_RANGE_SEARCH_GENERIC_H_

#include "DataTypes.h"

#if defined(__cplusplus)
namespace quasicontinuum {
extern "C" {
#endif /* __cplusplus */

struct BSP_tree_t {
  struct tnode_t *P_head;
  int n_nodes;
  int alloc_space;
};

struct tnode_t {
  void *object;
  enum { X, Y, Z } plane;
  struct tnode_t *P_tnode_left;
  struct tnode_t *P_tnode_right;
};

struct point_t {
  double x;
  double y;
  double z;
};

struct rectangle_t {
  struct point_t point1;
  struct point_t point2;
};

/**
 * prototypes
 */

extern struct BSP_tree_t *
BSP_tree_create(void *base, size_t n_elem, size_t width,
                struct point_t (*get_object_pos)(const void *));
extern void BSP_tree_destroy(struct BSP_tree_t *BSP_tree);

/*
extern void
range_search(struct object_list_t     *P_object_list,
             struct BSP_tree_t        *BSP_tree,
             const struct rectangle_t *P_rectangle,
             struct point_t          (*get_object_pos)(const void *));
*/

extern struct object_list_t *
range_search(struct BSP_tree_t *BSP_tree, const struct rectangle_t *P_rectangle,
             struct point_t (*get_object_pos)(const void *));

#if defined(__cplusplus)
}
}
#endif /* __cplusplus */

#endif /* _QC_RANGE_SEARCH_GENERIC_H_ */
