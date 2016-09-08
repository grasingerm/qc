/**
 * $Id: site_element_map_cache.h,v 1.2 2003/08/14 20:59:29 knap Exp $
 *
 * $Log: site_element_map_cache.h,v $
 * Revision 1.2  2003/08/14 20:59:29  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.2  2001/10/03 22:29:49  knap
 * Checked in the latest changes.
 *
 * Revision 1.1  2000/04/12 22:06:24  knap
 * Removed references to atoms.
 *
 */

#ifndef _QC_SITE_ELEMENT_MAP_CACHE_H_
#define _QC_SITE_ELEMENT_MAP_CACHE_H_

#include "DataTypes.h"

#if defined(__cplusplus)
namespace quasicontinuum {
  extern "C" {
#endif /* __cplusplus */

extern void check_hash_table_allocation(const int iQuasi);

extern void
hash_table_add_site(double                 *shape,
		    const int               l[3],
		    struct element_t       *P_element,
		    struct lattice_t *P_lattice,
		    const int 						  iQuasi);

extern struct element_t *
hash_table_search_element(double   *shape,
			  const int l[3],
			  const int  iQuasi);

extern void hash_table_clean(const int iQuasi);
/* extern void check_hash_table(struct node_list_t *P_node_list); */

#if defined(__cplusplus)
  }
}
#endif /* __cplusplus */

#endif /* _QC_SITE_ELEMENT_MAP_CACHE_H_ */

