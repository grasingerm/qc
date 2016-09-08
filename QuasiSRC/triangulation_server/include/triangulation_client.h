/**
 * $Id: triangulation_client.h,v 1.1 1999/12/18 17:50:28 knap Exp $
 *
 * $Log: triangulation_client.h,v $
 * Revision 1.1  1999/12/18 17:50:28  knap
 * Initial rev.
 *
 */

#ifndef _TRIANGULATION_CLIENT_H_
#define _TRIANGULATION_CLIENT_H_

#include "DataTypes.h"

#if defined(__cplusplus)
namespace quasicontinuum {
  extern "C" {
#endif /* __cplusplus */

extern void create_elements_external(struct node_list_t    *P_node_list,
				      struct mesh_data_t    *P_mesh_data);

#if defined(__cplusplus)
  }
}
#endif /* __cplusplus */
#endif /* _TRIANGULATION_CLIENT_H_ */
