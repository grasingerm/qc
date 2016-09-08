/**
 *------------------------------------------------------------------------------
 *                                       
 *               Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
 *                     California Institute of Technology
 *                        (C) 2001 All Rights Reserved
 *                                       
 *------------------------------------------------------------------------------
 * $Id: f77.h,v 1.1 2002/03/08 00:02:47 knap Exp $
 *
 * $Log: f77.h,v $
 * Revision 1.1  2002/03/08 00:02:47  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.2  2001/10/02 22:43:57  knap
 * Imported sources.
 *
 * Revision 1.1  2001/07/17 05:43:26  knap
 * Added autoconf macros.
 *
 * Revision 1.1  2001/05/14 04:17:39  knap
 * Initial sources.
 *
 */

#if !defined(_system_F77_f77_h)
#define _system_F77_f77_h

#if defined(HAVE_F77_SYMBOL_UNDERSCORE)
#define FORTRAN_FUNCTION(f) f##_
#else /* ! HAVE_F77_SYMBOL_UNDERSCORE */
#define FORTRAN_FUNCTION(f) f
#endif /* HAVE_F77_SYMBOL_UNDERSCORE */

#endif /* _system_F77_f77_h */ 
