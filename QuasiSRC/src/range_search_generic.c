/**
 * Generic range search routines in 3D
 *
 * $Log: range_search_generic.c,v $
 * Revision 1.4  2002/03/07 23:53:08  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.3  2001/10/03 22:29:52  knap
 * Checked in the latest changes.
 *
 * Revision 1.2  2000/04/20 01:09:09  knap
 * Ported to HP-UX.
 *
 * Revision 1.1  2000/04/12 22:07:15  knap
 * Removed references to atoms.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
/* No _REENTRANT */
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#ifdef STDC_HEADERS
#include <stdlib.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#else
#error sys/types.h not found.
#endif /* HAVE_SYS_TYPES_H */

#if defined(__QC_SGI)
#include "numa.h"
#endif /* sgi */

#include "monitor.h"
#include "range_search_generic.h"
#include "threads.h"

#include "DataTypes.h"
#include "Error.h"


#if !defined(lint)
static const char rcsid[] = "$Id: range_search_generic.c,v 1.4 2002/03/07 23:53:08 knap Exp $";
#endif /* !lint */

/**
 * local defs
 */

enum side_t { LEFT, RIGHT };


/**
 * add object to the list of objects
 */

#define OBJECT_ALLOC_BUCKET 10

static void
range_search_add_object(struct object_list_t *P_object_list,
			void                 *P_object)
{
  
  const size_t P_object_len=sizeof(P_object);
  int          i_object;

  /* i_object = P_object_list->n_objects + 1; */
  
  i_object = ++(P_object_list->n_objects);

  /**
   * check for allocated space
   */

  if( P_object_list->alloc_space == 0 ){
    if( (P_object_list->objects = malloc(OBJECT_ALLOC_BUCKET*P_object_len)) ==
	NULL ){
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }
    
    P_object_list->alloc_space += OBJECT_ALLOC_BUCKET;

  }
  else if( i_object > P_object_list->alloc_space ){
    
    P_object_list->alloc_space += OBJECT_ALLOC_BUCKET;
    
    P_object_list->objects = realloc(P_object_list->objects, P_object_len*
				     P_object_list->alloc_space);
    if( P_object_list->objects == NULL ){
      ERROR("realloc()");
      exit(EXIT_FAILURE);
    }

  }

  /*
   * insert atom
   */

  /*
  P_object_list->n_objects++;
  */

  i_object--;

  P_object_list->objects[i_object] = P_object;

  return;


}


/**
 * create node in a tree
 */

static struct tnode_t *
node_create(struct BSP_tree_t *BSP_tree)
{

  struct tnode_t *P_tnode;

  P_tnode=&(BSP_tree->P_head[BSP_tree->n_nodes]);

  /**
   * check for overflow
   */

  if (++(BSP_tree->n_nodes) > BSP_tree->alloc_space){
    fprintf(stderr, "Not enough space allocated for tree\n");
    exit(EXIT_FAILURE);
  }

  return(P_tnode);

}

/**
 * insert object
 */

static void
insert_object(struct BSP_tree_t *BSP_tree,
	      struct tnode_t    *P_tnode,
	      void              *object,
	      struct point_t   (*get_object_pos)(const void *))
{

  enum side_t     side;
  struct tnode_t *P_child_tnode;
  struct point_t  point = get_object_pos(P_tnode->object);
  struct point_t  new_point = get_object_pos(object);

  switch(P_tnode->plane){

  case X:
    
    if (new_point.x < point.x )
      side=LEFT;
    else
      side=RIGHT;
    
    break;

  case Y:

    if (new_point.y < point.y)
      side=LEFT;
    else
      side=RIGHT;

    break;

  case Z:

    if (new_point.z < point.z)
      side=LEFT;
    else                                       
      side=RIGHT;

    break;

  }

  /**
   * determine P_child_tnode
   */

  switch(side){

  case LEFT:
    
    P_child_tnode=P_tnode->P_tnode_left;
    break;

  case RIGHT:

    P_child_tnode=P_tnode->P_tnode_right;
    break;

  }

  /**
   * if P_child_tnode is the external node insert point here if not
   * recursively call itself with the child node
   */

  if( P_child_tnode == NULL ){

    switch(side){

    case LEFT:
      
      P_tnode->P_tnode_left=node_create(BSP_tree);
      P_child_tnode=P_tnode->P_tnode_left;
      break;
      
    case RIGHT:

      P_tnode->P_tnode_right=node_create(BSP_tree);
      P_child_tnode=P_tnode->P_tnode_right;
      break;

    }
    
    /**
     * fill in P_child_tnode
     */

    P_child_tnode->object = object;
    
    if(P_tnode->plane == X) P_child_tnode->plane = Y;
    if(P_tnode->plane == Y) P_child_tnode->plane = Z;
    if(P_tnode->plane == Z) P_child_tnode->plane = X;

  }
  else
    insert_object(BSP_tree, P_child_tnode, object, get_object_pos);

  return;

}




/**
 * create BSP tree 
 */

struct BSP_tree_t *
BSP_tree_create(void    *base,
		size_t   n_elem,
		size_t   width,
		struct point_t (*get_object_pos)(const void *))
{

  struct BSP_tree_t *BSP_tree;
  caddr_t address;
  int     i_object;

  /**
   * initialize BSP tree
   */

  if ((BSP_tree=malloc(sizeof(struct BSP_tree_t))) == NULL) {
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  if ((BSP_tree->P_head=calloc(n_elem, sizeof(struct tnode_t))) == NULL) {
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

#if defined(__QC_SGI)
  change_memory_placement(BSP_tree, sizeof(struct BSP_tree_t), 0);
  change_memory_placement(BSP_tree->P_head, n_elem*sizeof(struct tnode_t), 0);
#endif /* sgi */

  BSP_tree->n_nodes = 1;
  BSP_tree->alloc_space = n_elem;

  /**
   * initialize the head of the tree
   */

  BSP_tree->P_head->object = base;
  BSP_tree->P_head->plane = X;
  BSP_tree->P_head->P_tnode_left = NULL;
  BSP_tree->P_head->P_tnode_right = NULL;

  /**
   * build the tree
   */

  for (i_object = 1; i_object < n_elem; i_object++) {
    
    address = (caddr_t) base +i_object*width;
    insert_object(BSP_tree, BSP_tree->P_head, (void *) address, 
		  get_object_pos);

  }

  return(BSP_tree);

}

/**
 * destroy BSP tree
 */

void
BSP_tree_destroy(struct BSP_tree_t *BSP_tree)
{

  /**
   * free nodes
   */

  free(BSP_tree->P_head);

  /**
   * free tree
   */

  free(BSP_tree);

  return;

}

/**
 * check if object is inside a rectangle
 */

static int 
is_inside_rectangle(const void               *object,
		    const struct rectangle_t *P_rectangle,
		    struct point_t          (*get_object_pos)(const void *))
{

  struct point_t P_point = get_object_pos(object);

  if( P_rectangle->point1.x <= P_point.x && 
      P_rectangle->point1.y <= P_point.y &&
      P_rectangle->point1.z <= P_point.z &&
      P_rectangle->point2.x >= P_point.x &&
      P_rectangle->point2.y >= P_point.y &&
      P_rectangle->point2.z >= P_point.z) return(1);

  return(0);

}

/**
 * range search; generate list of objects inside a rectangle
 */

static void 
_range_search(struct object_list_t     *P_object_list,
	      const struct tnode_t     *P_tnode,
	      const struct rectangle_t *P_rectangle,
	      struct point_t          (*get_object_pos)(const void *))
{

  struct point_t P_point;

  /**
   * check if P_tnode is NOT an external node
   */
  
  if( P_tnode == NULL ) return;

  /**
   * get position of object at the node
   */

  P_point = get_object_pos(P_tnode->object);

  /**
   * check if node contains the point that is inside rectangle
   */

  if( is_inside_rectangle(P_tnode->object, P_rectangle, get_object_pos) )
    range_search_add_object(P_object_list, P_tnode->object);


  switch(P_tnode->plane){

  case X:

    if( P_rectangle->point1.x < P_point.x )
      _range_search(P_object_list, P_tnode->P_tnode_left, P_rectangle,
		    get_object_pos);

    if( P_rectangle->point2.x >= P_point.x )
      _range_search(P_object_list, P_tnode->P_tnode_right, P_rectangle,
		    get_object_pos);
    break;

  case Y:

    if( P_rectangle->point1.y < P_point.y )
      _range_search(P_object_list, P_tnode->P_tnode_left, P_rectangle,
		    get_object_pos);

    if( P_rectangle->point2.y >= P_point.y )
      _range_search(P_object_list, P_tnode->P_tnode_right, P_rectangle,
		    get_object_pos);

    break;

  case Z:

    if( P_rectangle->point1.z < P_point.z )
      _range_search(P_object_list, P_tnode->P_tnode_left, P_rectangle,
		    get_object_pos);

    if( P_rectangle->point2.z >= P_point.z )
      _range_search(P_object_list, P_tnode->P_tnode_right, P_rectangle,
		    get_object_pos);

    break;

  }

  return;

}

/**
 * initialize range search routine
 */

static pthread_once_t range_search_once = PTHREAD_ONCE_INIT;

static struct object_list_t *P_object_list;

void 
range_search_init(void)
{

  int number_threads = get_max_number_threads();

  if ((P_object_list = calloc(number_threads, sizeof(struct object_list_t)))
      == NULL) {
    ERROR("realloc()");
    exit(EXIT_FAILURE);
  }
  
#if defined(__QC_SGI)
  change_memory_placement(P_object_list, 
			  number_threads*sizeof(struct object_list_t), 0);
#endif /* sgi */  

  return;

}



/**
 * begin range search
 */

/*
void
range_search(struct object_list_t     *P_object_list,
	     struct BSP_tree_t        *BSP_tree,
	     const struct rectangle_t *P_rectangle,
	     struct point_t          (*get_object_pos)(const void *))
*/

struct object_list_t *
range_search(struct BSP_tree_t        *BSP_tree,
	     const struct rectangle_t *P_rectangle,
	     struct point_t          (*get_object_pos)(const void *))
{

  int tid = get_my_tid();

  /**
   * initialize thread specific data
   */

  pthread_once(&range_search_once, range_search_init);

  /**
   * 
   */

  P_object_list[tid].n_objects = 0;

  /**
   * start range search from the trunk
   */

  _range_search(&(P_object_list[tid]), BSP_tree->P_head, P_rectangle, 
		get_object_pos);

  /**
   * 
   */

  return(&(P_object_list[tid]));

}


