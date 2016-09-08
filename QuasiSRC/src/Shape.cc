//
// Shape.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "Shape.h"
#include "Lattice.h"
#include "DataTypes.h"


//
//
//

namespace quasicontinuum {
  //
  //
  //

  Shape* Shape::_instance = NULL;

  //
  // constructor
  //

  Shape::Shape()
  {

    //
    //
    //
    return;

  }

  //
  // destructor
  //

  Shape::~Shape()
  {

    //
    //
    //
    return;

  }

  //
  // getInstance method
  //

  Shape*
  Shape::getInstance()
  {

    //
    // if not created, create
    //
    if(_instance == NULL){
      _instance = new Shape();
    }

    //
    // return instance
    //
    return _instance;

  }

  //
  // destroy instance method
  //

  void
  Shape::destroyInstance()
  {

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
  //  computeShapeFunction()
  //
  double Shape::computeShapeFunction(const double p1[3],
            const double p2[3],
            const double p3[3],
            const double p4[3],
            const double p[3])  
  {
    double shape;

    shape = ( p2[0]*p3[2]*p4[1] - p2[0]*p3[1]*p4[2] - p3[2]*p4[1]*p[0] 
            + p3[1]*p4[2]*p[0]  - p2[0]*p3[2]*p[1]  + p3[2]*p4[0]*p[1]
            + p2[0]*p4[2]*p[1]  - p3[0]*p4[2]*p[1]  
            + p2[2]*( -p3[0]*p4[1]  + p3[1]*( p4[0] - p[0] ) 
                      + p4[1]*p[0] + p3[0]*p[1] - p4[0]*p[1] )
            + p2[0]*p3[1]*p[2] - p3[1]*p4[0]*p[2] - p2[0]*p4[1]*p[2]
            + p3[0]*p4[1]*p[2] 
            + p2[1]*( - p3[2]*p4[0] + p3[0]*p4[2] + p3[2]*p[0] 
                      - p4[2]*p[0]  - p3[0]*p[2]  + p4[0]*p[2]) 
            )/
            (-p1[0]*p2[2]*p3[1] + p1[0]*p2[1]*p3[2] + p2[2]*p3[1]*p4[0]
            - p2[1]*p3[2]*p4[0] + p1[0]*p2[2]*p4[1] - p2[2]*p3[0]*p4[1]
            - p1[0]*p3[2]*p4[1] + p2[0]*p3[2]*p4[1] 
            + p1[2]*( p2[0]*p3[1] - p3[1]*p4[0] + p2[1]*( -p3[0] + p4[0])
                     -p2[0]*p4[1] + p3[0]*p4[1] )
            - p1[0]*p2[1]*p4[2] + p2[1]*p3[0]*p4[2] + p1[0]*p3[1]*p4[2]
            - p2[0]*p3[1]*p4[2] + p1[1]*( p2[2]*p3[0] - p2[0]*p3[2]
                                         -p2[2]*p4[0] + p3[2]*p4[0]
                                         +p2[0]*p4[2] - p3[0]*p4[2]) 
            );

    return(shape);
  }

  //
  //  computeSiteShapeFunction()
  //
  double Shape::computeSiteShapeFunction(const struct node_t    *P_node_0,
                     const struct node_t    *P_node_1,
                     const struct node_t    *P_node_2,
                     const struct node_t    *P_node_3,
                     const int               l[3],
                     struct lattice_t *P_lattice,
                     const int               iQuasi)
  {
      double X[3];
      double shape;

      Lattice::getInstance()->getSiteInitialPosition(X, l, P_lattice, iQuasi);

      shape = computeShapeFunction(P_node_0->initial_position, 
                P_node_1->initial_position,
                P_node_2->initial_position,
                P_node_3->initial_position,
                X );

      return(shape);
  }  

  //
  //  computeShapeFunctionGradient()
  //
  void Shape::computeShapeFunctionGradient(const double p1[3],
                   const double p2[3],
                   const double p3[3],
                   const double p4[3],
                   double       gradient[3])
  {
    gradient[0]=(- p3[2]*p4[1] + p2[2]*(-p3[1] + p4[1]) 
                 + p2[1]*(p3[2] - p4[2]) + p3[1]*p4[2]
                )/
                (- p1[0]*p2[2]*p3[1] + p1[0]*p2[1]*p3[2] + p2[2]*p3[1]*p4[0] 
                 - p2[1]*p3[2]*p4[0] + p1[0]*p2[2]*p4[1] - p2[2]*p3[0]*p4[1] 
                 - p1[0]*p3[2]*p4[1] + p2[0]*p3[2]*p4[1] + p1[2]*(p2[0]*p3[1]
                 - p3[1]*p4[0] + p2[1]*(-p3[0] + p4[0]) 
                 - p2[0]*p4[1] + p3[0]*p4[1]) - p1[0]*p2[1]*p4[2] 
                 + p2[1]*p3[0]*p4[2] + p1[0]*p3[1]*p4[2] - p2[0]*p3[1]*p4[2] 
                 + p1[1]*(  p2[2]*p3[0] - p2[0]*p3[2] - p2[2]*p4[0] 
                          + p3[2]*p4[0] + p2[0]*p4[2] - p3[0]*p4[2] )
                 );

    gradient[1]=(- p2[0]*p3[2] + p2[2]*(p3[0] - p4[0]) + p3[2]*p4[0] 
                 + p2[0]*p4[2] - p3[0]*p4[2]
                )/
                (- p1[0]*p2[2]*p3[1] + p1[0]*p2[1]*p3[2] + p2[2]*p3[1]*p4[0] 
                 - p2[1]*p3[2]*p4[0] + p1[0]*p2[2]*p4[1] - p2[2]*p3[0]*p4[1] 
                 - p1[0]*p3[2]*p4[1] + p2[0]*p3[2]*p4[1] 
                 + p1[2]*( p2[0]*p3[1] - p3[1]*p4[0] 
                          +p2[1]*(-p3[0] + p4[0]) - p2[0]*p4[1] 
                          + p3[0]*p4[1] ) 
                 - p1[0]*p2[1]*p4[2] + p2[1]*p3[0]*p4[2] + p1[0]*p3[1]*p4[2] 
                 - p2[0]*p3[1]*p4[2] 
                 + p1[1]*( p2[2]*p3[0] - p2[0]*p3[2] - p2[2]*p4[0] 
                          +p3[2]*p4[0] + p2[0]*p4[2] - p3[0]*p4[2] )
                 );

    gradient[2]=(  p2[0]*p3[1] - p3[1]*p4[0] + p2[1]*(-p3[0] + p4[0]) 
                 - p2[0]*p4[1] + p3[0]*p4[1]
                )/
                (- p1[0]*p2[2]*p3[1] + p1[0]*p2[1]*p3[2] + p2[2]*p3[1]*p4[0] 
                 - p2[1]*p3[2]*p4[0] + p1[0]*p2[2]*p4[1] - p2[2]*p3[0]*p4[1] 
                 - p1[0]*p3[2]*p4[1] + p2[0]*p3[2]*p4[1] 
                 + p1[2]*( p2[0]*p3[1] - p3[1]*p4[0] + p2[1]*(-p3[0] + p4[0])
                          -p2[0]*p4[1] + p3[0]*p4[1] ) 
                 - p1[0]*p2[1]*p4[2] + p2[1]*p3[0]*p4[2] + p1[0]*p3[1]*p4[2] 
                 - p2[0]*p3[1]*p4[2] + p1[1]*(p2[2]*p3[0] - p2[0]*p3[2] 
                 - p2[2]*p4[0] + p3[2]*p4[0] + p2[0]*p4[2] - p3[0]*p4[2])
                 );

    return;
  }

}