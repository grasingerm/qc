/**
 * $Id: site_element_map.h,v 1.2 2003/08/14 20:59:29 knap Exp $
 *
 * $Log: site_element_map.h,v $
 * Revision 1.2  2003/08/14 20:59:29  knap
 * Added C++ support.
 *
 * Revision 1.1.1.1  2003/06/03 23:27:36  jaime
 * Imported sources
 *
 * Revision 1.2  2001/10/03 22:29:48  knap
 * Checked in the latest changes.
 *
 * Revision 1.1  2000/04/12 22:06:24  knap
 * Removed references to atoms.
 *
 */

#ifndef _QC_SITE_ELEMENT_MAP_H_
#define _QC_SITE_ELEMENT_MAP_H_

#include "DataTypes.h"
#include "range_search_generic.h"

#define _QC_SITE_ELEMENT_MAP_DEBUG_
// #if defined(_QC_SITE_ELEMENT_MAP_DEBUG_)
// #include "lattice.h"
// #endif /* _QC_SITE_ELEMENT_MAP_DEBUG_ */

#if defined(__cplusplus)
namespace quasicontinuum {
extern "C" {
#endif /* __cplusplus */

extern void
initialize_site_element_map_data(const int iQuasi,
                                 struct element_list_t *P_element_list);
extern void clean_site_element_map_data(const int iQuasi);

extern struct element_t *locate_site_element(double *shape, const int l[3],
                                             struct lattice_t *P_lattice,
                                             const int iQuasi);

extern struct element_t *
find_site_element(const int l[3], struct lattice_t *P_lattice,
                  const struct element_list_t *P_element_list,
                  const int iQuasi);

#if defined(_QC_SITE_ELEMENT_MAP_DEBUG_)
extern void site_element_map_debug(struct lattice_t *P_lattice,
                                   const int iQuasi);
#endif /* _QC_SITE_ELEMENT_MAP_DEBUG_ */

#if defined(__cplusplus)
}
}
#endif /* __cplusplus */

#endif /* _QC_SITE_ELEMENT_MAP_H_ */
