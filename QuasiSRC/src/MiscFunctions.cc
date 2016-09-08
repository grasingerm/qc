//
// MiscFunctions.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifndef _REENTRANT
#define _REENTRANT
#endif /* _REENTRANT */

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found.
#endif /* HAVE_UNISTD_H */

// vector
#ifdef HAVE_VECTOR
#include <vector>
#else
#ifdef HAVE_VECTOR_H
#include <vector.h>
#else
#error No vector or vector.h available
#endif // HAVE_VECTOR_H
#endif // HAVE_VECTOR

#include <math.h>
#include <cassert>

#include "MiscFunctions.h"
#include "Error.h"
#include "DataTypes.h"
#include "Lattice.h"

#define ATOM_NODE_ERROR 0.01
#define CHECK_ERROR     1.0e-3

//
//
//

namespace quasicontinuum {

  namespace
  {
    //
    //  getMidplaneLocal()
    //
    void getMidplaneLocal(const double x1[3], const double x2[3], double a[4])
    {
      double xm=(x1[0]+x2[0])/2.0;
      double ym=(x1[1]+x2[1])/2.0;
      double zm=(x1[2]+x2[2])/2.0;
      
      a[0]=x2[0]-x1[0];
      a[1]=x2[1]-x1[1];
      a[2]=x2[2]-x1[2];
      a[3]=a[0]*xm+a[1]*ym+a[2]*zm;

      return;
    }
  } // end of namespace
  //
  //
  //

  MiscFunctions* MiscFunctions::_instance = NULL;

  //
  // constructor
  //

  MiscFunctions::MiscFunctions()
  {

    //
    //
    //
    return;

  }

  //
  // destructor
  //

  MiscFunctions::~MiscFunctions()
  {

    //
    //
    //
    return;

  }

  //
  // getInstance method
  //

  MiscFunctions*
  MiscFunctions::getInstance()
  {

    //
    // if not created, create
    //
    if(_instance == NULL){
      _instance = new MiscFunctions();
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
  MiscFunctions::destroyInstance()
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
  //  freeE()
  //
  void * MiscFunctions::freeE(void * p)
  {
    free(p);
    return ((void *) NULL);
  } // end of freeE()

  //
  //  tetraVolume()
  //
  double MiscFunctions::tetraVolume(const double   a[3],
        const double   b[3],
        const double   c[3],
        const double   d[3])
  {
    double vol;

    vol =   ((a[0]-b[0]) * 
            ((a[1]-c[1]) * (a[2]-d[2]) - (a[2]-c[2]) * (a[1]-d[1]))
           - (a[1]-b[1]) * 
            ((a[0]-c[0]) * (a[2]-d[2]) - (a[0]-d[0]) * (a[2]-c[2]))
           + (a[2]-b[2]) *
            ((a[0]-c[0]) * (a[1]-d[1]) - (a[1]-c[1]) * (a[0]-d[0]))
            )/6.0;

    return vol;
  }

  //
  //  getUnitVector()
  //
  int MiscFunctions::getUnitVector( const double  v[3],
    double u[3] )  
  {
    double len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

    if( len <= 0.000000001 )
      return(RETURN_FAILURE);
  
    u[0]=v[0]/len;
    u[1]=v[1]/len;
    u[2]=v[2]/len;
  
    return(RETURN_SUCCESS);
  }

  //
  //  getDistanceBetPoints()
  //
  double MiscFunctions::getDistanceBetPoints(const double a[3],
          const double     b[3])
  {
    double ab_local[3];
    
    ab_local[0]=b[0]-a[0];
    ab_local[1]=b[1]-a[1];
    ab_local[2]=b[2]-a[2];

    return(sqrt(ab_local[0]*ab_local[0]+ab_local[1]*ab_local[1]+
    ab_local[2]*ab_local[2]));    
  }

  //
  //  getDistanceSqrBetweenPoints()
  //
  double MiscFunctions::getDistanceSqrBetweenPoints(const double  a[3],
          const double      b[3])
  {
    double ab_local[3];
    
    ab_local[0]=b[0]-a[0];
    ab_local[1]=b[1]-a[1];
    ab_local[2]=b[2]-a[2];

    return(ab_local[0]*ab_local[0]+ab_local[1]*ab_local[1]+
     ab_local[2]*ab_local[2]);    
  } 

  //
  //  computeDifferenceOfVectorAndAbsValue()
  //
  //  writes the difference of vector to third argument
  //  returns sqrt of abs difference
  //
  double MiscFunctions::computeDifferenceAndAbsValue(const double a[3],
            const double           b[3],
            double                 ab[3])
  {
    ab[0]=b[0]-a[0];
    ab[1]=b[1]-a[1];
    ab[2]=b[2]-a[2];

    return(sqrt(ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2]));    
  }

  //
  //  checkPointInTetra()
  //
  enum loc_t  MiscFunctions::checkPointInTetra(const double v1[3], 
           const double v2[3], 
           const double v3[3],
           const double v4[3], 
           const double pos[3])
  {
    double v21[3];
    double v31[3];
    double v41[3];
    double p[3];

    double det;

    double r;
    double s;
    double t;
    double u;

    v21[0]=v2[0]-v1[0];
    v21[1]=v2[1]-v1[1];
    v21[2]=v2[2]-v1[2];

    v31[0]=v3[0]-v1[0];
    v31[1]=v3[1]-v1[1];
    v31[2]=v3[2]-v1[2]; 

    v41[0]=v4[0]-v1[0];
    v41[1]=v4[1]-v1[1];
    v41[2]=v4[2]-v1[2];

    p[0]=pos[0]-v1[0];
    p[1]=pos[1]-v1[1];
    p[2]=pos[2]-v1[2];

    /**
     *
     */

    det=-v21[2]*v31[1]*v41[0] + v21[1]*v31[2]*v41[0] + v21[2]*v31[0]*v41[1] -
         v21[0]*v31[2]*v41[1] - v21[1]*v31[0]*v41[2] + v21[0]*v31[1]*v41[2];

    /*
     * compute r
     */

    r=(-p[2]*v31[1]*v41[0] + p[1]*v31[2]*v41[0] + p[2]*v31[0]*v41[1] - 
        p[0]*v31[2]*v41[1] - p[1]*v31[0]*v41[2] + p[0]*v31[1]*v41[2])/det;

    /**
     * compute s
     */

    s=(p[2]*v21[1]*v41[0] - p[1]*v21[2]*v41[0] - p[2]*v21[0]*v41[1] + 
       p[0]*v21[2]*v41[1] + p[1]*v21[0]*v41[2] - p[0]*v21[1]*v41[2])/det;

    /*
     * compute t
     */

    t=(-p[2]*v21[1]*v31[0] + p[1]*v21[2]*v31[0] + p[2]*v21[0]*v31[1] - 
        p[0]*v21[2]*v31[1] - p[1]*v21[0]*v31[2] + p[0]*v21[1]*v31[2])/det;

    /**
     * u
     */

    u=1.0-r-s-t;
    
    /**
     * check if all coordinates are in [0,1]
     */
    /* printf("r= %f s= %f t= %f u= %f\n", r, s, t, u); */
    if( r < -CHECK_ERROR || r > 1.0+CHECK_ERROR ) return(OUTSIDE);
    if( s < -CHECK_ERROR || s > 1.0+CHECK_ERROR ) return(OUTSIDE);
    if( t < -CHECK_ERROR || t > 1.0+CHECK_ERROR ) return(OUTSIDE);
    if( u < -CHECK_ERROR || u > 1.0+CHECK_ERROR ) return(OUTSIDE);
    
    return(INSIDE);
  }

  //
  //  checkSiteInTetra()
  //
  enum loc_t MiscFunctions::checkSiteInTetra(const double            v1[3], 
         const double            v2[3], 
         const double            v3[3],
         const double            v4[3],
         const int               l[3],
         struct lattice_t *P_lattice,
         const int                iQuasi)
  {
    double r[3];

    Lattice::getInstance()->getSiteInitialPosition(r, l, P_lattice, iQuasi);

    return(checkPointInTetra(v1, v2, v3, v4, r));
  }

  //
  //  findTetraCenter()
  //
  double MiscFunctions::findTetraCenter(const double vertice1[3], 
       const double vertice2[3],
       const double vertice3[3], 
       const double vertice4[3],
       double       center[3])
  {
    double a[3][4];
    double b[4][3];
    double radius=-1.0;

    int i;
    int j;

    getMidplaneLocal(vertice1, vertice2, a[0]);
    getMidplaneLocal(vertice1, vertice3, a[1]);
    getMidplaneLocal(vertice1, vertice4, a[2]);

    /** 
     * transpose a into b
     */

    for(i=0; i < 4; i++)
      for(j=0; j < 3; j++)
        b[i][j]=a[j][i];
    
    if( solveLinearThreeD(b[0], b[1], b[2], b[3], center) == RETURN_SUCCESS )
      radius= (vertice1[0]-center[0])*(vertice1[0]-center[0])+
              (vertice1[1]-center[1])*(vertice1[1]-center[1])+
              (vertice1[2]-center[2])*(vertice1[2]-center[2]);

    return(radius);    
  }

  //
  //  computeDeterminant()
  //
  double MiscFunctions::computeDeterminant(const double      a[3],
        const double      b[3],
        const double      c[3])
  {
    return( a[0]*(b[1]*c[2]-b[2]*c[1]) - 
            b[0]*(a[1]*c[2]-a[2]*c[1]) +
            c[0]*(a[1]*b[2]-a[2]*b[1]) );    
  }

  //
  //  solveLinearThreeD()
  //
  int MiscFunctions::solveLinearThreeD(const double    a[3],
            const double    b[3],
            const double    c[3],
            const double    d[3],
            double          x[3]) 
  {
    double det=computeDeterminant(a, b, c);

    if( det == 0.0 ) return(RETURN_FAILURE);

    x[0]=computeDeterminant(d, b, c)/det;
    x[1]=computeDeterminant(a, d, c)/det;
    x[2]=computeDeterminant(a, b, d)/det;

    return(RETURN_SUCCESS);    
  }

  //
  // calculate cross product of 2 vectors
  //
  std::vector<double>
  MiscFunctions::crossProduct3x3(std::vector<double> a,
         std::vector<double> b)
  {

    //
    // check size of vectors
    //
    assert(a.size() == b.size() && int(a.size()) == 3);

    //
    // result
    //
    std::vector<double> result(3,0.0);

    //
    // calculate cross product
    //
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];

    //
    //
    //
    return result;

  }

  //
  // calculate determinant of 3x3
  //
  double
  MiscFunctions::determinant3x3(std::vector<std::vector<double> > a)
  {

    //
    // return variable
    //
    double result;

    //
    // calculate determinant
    //
    result = 
      a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
      a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
      a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    //
    //
    //
    return result;

  }

  //
  // calculate L2 Norm of vector
  //
  double
  MiscFunctions::L2Norm(std::vector<double> a)
  {

    //
    // check size of vector
    //
    assert(a.size() > 0);

    //
    // result
    //
    double result = 0.0;

    //
    // calculate squares
    //
    for (unsigned int i = 0; i < a.size(); ++i)
      result += a[i] * a[i];

    //
    // calculate square root
    //
    result = sqrt(result);

    //
    //
    //
    return result;

  }

  //
  // calculate dot product of two vectors
  //
  double
  MiscFunctions::dotProduct(const std::vector<double> &a,
          const std::vector<double> &b)
  {

    //
    // check size of vector
    //
    assert((a.size() > 0) && (a.size() == b.size()));

    //
    // result
    //
    double result = 0.0;

    //
    // calculate squares
    //
    for (unsigned int i = 0; i < a.size(); ++i)
      result += a[i] * b[i];

    //
    //
    //
    return result;

  }

  //
  // calculate area of triangle
  //
  double 
  MiscFunctions::TriangleArea(std::vector<double> vertexA,
            std::vector<double> vertexB,
            std::vector<double> vertexC)
  {

    //
    // vectors u and v
    //
    std::vector<double> u;
    std::vector<double> v;

    //
    // calculate vector u and v
    //
    for(int iDof = 0; iDof < 3; ++iDof){

      //
      // get vectors u and v
      //
      u.push_back(vertexB[iDof] - vertexA[iDof]);
      v.push_back(vertexC[iDof] - vertexA[iDof]);

    }

    //
    // calculate cross product of u and v
    //
    std::vector<double> crossVector = crossProduct3x3(u,v);

    //
    // calculate magnitude
    //
    double area = 0.5 * L2Norm(crossVector);

    //
    // 
    //
    return area;

  }

  //
  // calculate area of triangle
  //
  std::vector<double> 
  MiscFunctions::TriangleCenter(std::vector<double> vertexA,
        std::vector<double> vertexB,
        std::vector<double> vertexC)
  {

    //
    // center location 
    //
    std::vector<double> center;
 
    //
    // calculate center
    //
    for(int iDof = 0; iDof < 3; ++iDof){

      //
      // get center vector
      //
      center.push_back((vertexA[iDof] + vertexB[iDof] + vertexC[iDof]) / 3.0 );

    }

    //
    // 
    //
    return center;

  }

  //
  // get sphere ray intersection point
  //
  std::vector<double> 
  MiscFunctions::getSphereRayIntersection(std::vector<double> center,
            std::vector<double> nodeInSphere,
            std::vector<double> nodeOutsideSphere,
            double              radiusSquared)
  {

    //
    // location of intersection
    //
    std::vector<double> location;

    //
    // change coordinate system to center is 0
    //
    for(int iDof = 0; iDof < 3; ++iDof){
      nodeInSphere[iDof] -= center[iDof];
      nodeOutsideSphere[iDof] -= center[iDof];
    }

    //
    // get direction vector
    //
    std::vector<double> direction;
    for(int iDof = 0; iDof < 3; ++iDof)
      direction.push_back(nodeOutsideSphere[iDof] - nodeInSphere[iDof]);

    //
    // get A, B, C, and discrimant coefficients
    //
    double A = MiscFunctions::dotProduct(direction, direction);
    double B = 2.0 * MiscFunctions::dotProduct(direction, nodeInSphere);
    double C = MiscFunctions::dotProduct(nodeInSphere, nodeInSphere) - radiusSquared;
    double disc = B * B - 4.0 * A * C;
    double t;

    //
    // check discriminant sign
    //
    assert(disc > 0.0);

    //
    // compute q
    //
    double discSqrt = sqrt(disc);
    double q;
    if(B < 0.0)
      q = (-B - discSqrt)/ 2.0;
    else
      q = (-B + discSqrt)/ 2.0;

    //
    // compute t0 and t1
    //
    double t0 = q / A;
    double t1 = C / q;

    //
    // make sure t0 is smaller
    //
    if (t0 > t1){

      double temp = t0;
      t0 = t1;
      t1 = temp;

    }

    //
    // make sure t1 is greater 0
    //
    assert(t1 > 0.0);

    //
    // if t0 is less than 0, intersection is at t1
    //
    if(t0 < 0.0)
      t = t1;
    else
      t = t0;

    //
    // get location
    //
    for(int iDof = 0; iDof < 3; ++iDof)
      location.push_back(nodeInSphere[iDof] + t * direction[iDof] + center[iDof]);

    //
    //
    //
    return location;

  }

  //
  // difference between vectors
  //
  std::vector<double>
  MiscFunctions::vectorDifference(const std::vector<double> &vectorA,
          const std::vector<double> &vectorB)
  {

    //
    // make sure vector sizes are the same
    //
    assert(int(vectorA.size()) == int(vectorB.size()));

    //
    // calculate difference
    //
    std::vector<double> diff;
    for(int iDof = 0; iDof < int(vectorA.size()); ++iDof)
      diff.push_back(vectorA[iDof]-vectorB[iDof]);

    //
    //
    //
    return diff;
    
  }

  //
  // add 2 vectors
  //
  std::vector<double>
  MiscFunctions::vectorAddition(const std::vector<double> &vectorA,
        const std::vector<double> &vectorB)
  {

    //
    // make sure vector sizes are the same
    //
    assert(int(vectorA.size()) == int(vectorB.size()));

    //
    // add vectors
    //
    std::vector<double> addition;
    for(int iDof = 0; iDof < int(vectorA.size()); ++iDof)
      addition.push_back(vectorA[iDof]+vectorB[iDof]);

    //
    //
    //
    return addition;
    
  }

  //
  // add 3 vectors
  //
  std::vector<double>
  MiscFunctions::vectorAddition(const std::vector<double> &vectorA,
        const std::vector<double> &vectorB,
        const std::vector<double> &vectorC)
  {
    // make sure vector sizes are the same
    assert(int(vectorA.size()) == int(vectorB.size()) &&
     int(vectorA.size()) == int(vectorC.size()));

    // add vectors
    std::vector<double> addition;
    for(int i_dof = 0; i_dof < int(vectorA.size()); ++i_dof)
      addition.push_back(vectorA[i_dof] + vectorB[i_dof] + vectorC[i_dof]);

    //
    return addition;
    
  }

  //
  // scale vector
  //
  std::vector<double>
  MiscFunctions::vectorScale(std::vector<double> vectorA,
           const double scale)
  {

    //
    // scale vector
    //
    for(int iDof = 0; iDof < int(vectorA.size()); ++iDof)
      vectorA[iDof] *= scale;

    //
    //
    //
    return vectorA;
    
  }

  //
  // triangle sphere intersection
  //
  bool
  MiscFunctions::triangleSphereIntersection(const std::vector<double> &vertexA,
              const std::vector<double> &vertexB,
              const std::vector<double> &vertexC,
              const std::vector<double> &sphereCenter,
              const double &sphereRadiusSquared)
  {

    //
    // code from http://realtimecollisiondetection.net/blog/?p=103
    //
    std::vector<double> A = vectorDifference(vertexA,sphereCenter);
    std::vector<double> B = vectorDifference(vertexB,sphereCenter);
    std::vector<double> C = vectorDifference(vertexC,sphereCenter);

    std::vector<double> V = crossProduct3x3(vectorDifference(B,A), vectorDifference(C,A));
    double d = dotProduct(A, V);
    double e = dotProduct(V, V);
    bool sep1 = d * d > sphereRadiusSquared * e;
    double aa = dotProduct(A, A);
    double ab = dotProduct(A, B);
    double ac = dotProduct(A, C);
    double bb = dotProduct(B, B);
    double bc = dotProduct(B, C);
    double cc = dotProduct(C, C);
    bool sep2 = (aa > sphereRadiusSquared) & (ab > aa) & (ac > aa);
    bool sep3 = (bb > sphereRadiusSquared) & (ab > bb) & (bc > bb);
    bool sep4 = (cc > sphereRadiusSquared) & (ac > cc) & (bc > cc);
    std::vector<double> AB = vectorDifference(B, A);
    std::vector<double> BC = vectorDifference(C, B);
    std::vector<double> CA = vectorDifference(A, C);
    double d1 = ab - aa;
    double d2 = bc - bb;
    double d3 = ac - cc;
    double e1 = dotProduct(AB, AB);
    double e2 = dotProduct(BC, BC);
    double e3 = dotProduct(CA, CA);
    std::vector<double> Q1 = vectorDifference(vectorScale(A,e1),vectorScale(AB,d1));
    std::vector<double> Q2 = vectorDifference(vectorScale(B,e2),vectorScale(BC,d2));
    std::vector<double> Q3 = vectorDifference(vectorScale(C,e3),vectorScale(CA,d3));
    std::vector<double> QC = vectorDifference(vectorScale(C,e1),Q1);
    std::vector<double> QA = vectorDifference(vectorScale(A,e2),Q2);
    std::vector<double> QB = vectorDifference(vectorScale(B,e3),Q3);
    bool sep5 =
      (dotProduct(Q1, Q1) > sphereRadiusSquared * e1 * e1) & (dotProduct(Q1, QC) > 0);
    bool sep6 =
      (dotProduct(Q2, Q2) > sphereRadiusSquared * e2 * e2) & (dotProduct(Q2, QA) > 0);
    bool sep7 =
      (dotProduct(Q3, Q3) > sphereRadiusSquared * e3 * e3) & (dotProduct(Q3, QB) > 0);
    bool separated = sep1 | sep2 | sep3 | sep4 | sep5 | sep6 | sep7;

    //
    //
    //
    return separated;
    
  }

  //
  // triangle sphere intersection
  //
  bool
  MiscFunctions::edgeSphereIntersection(const std::vector<double> &vertexA,
          const std::vector<double> &vertexB,
          const std::vector<double> &sphereCenter,
          const double &sphereRadiusSquared,
          std::vector<double> &intersectPoint1,
          std::vector<double> &intersectPoint2)
  {

    //
    // clear intersect points
    //
    intersectPoint1.clear();
    intersectPoint2.clear();

    //
    // calculations
    //
    const std::vector<double> dp = vectorDifference(vertexB,vertexA);
    const double a = dotProduct(dp,dp);
    const double b =
      2.0 * (dp[0] * (vertexA[0] - sphereCenter[0]) + dp[1] * (vertexA[1] - sphereCenter[1]) + dp[2] * (vertexA[2] - sphereCenter[2]));
    const double c =
      dotProduct(sphereCenter,sphereCenter) + dotProduct(vertexA,vertexA) - 2.0 * dotProduct(sphereCenter,vertexA) - sphereRadiusSquared;
    const double bb4ac = b * b - 4.0 * a * c;
    if (bb4ac <= 0.0)
      return false;

    //
    // calculate points
    //
    const double mu1 = (-b + sqrt(bb4ac)) / (2.0 * a);
    const double mu2 = (-b - sqrt(bb4ac)) / (2.0 * a);

    //
    // check if intersection is between two points
    //
    if(mu1 >= 0.0 && mu1 <= 1.0 && mu2 >= 0.0 && mu2 <= 1.0){

      for(int iDof = 0; iDof < 3; ++iDof){
  intersectPoint1.push_back(vertexA[iDof] + mu1 * dp[iDof]);
  intersectPoint2.push_back(vertexA[iDof] + mu2 * dp[iDof]);
      }

      //
      //
      //
      return true;

    }

    //
    //
    //
    return false;
    
  }

  //
  // check if point is inside sphere
  //
  bool MiscFunctions::IsPointInSphere(const std::vector<double> test_point,
              const std::vector<double> sphere_center,
              const double sphere_radius_squared)
  {
    // get difference between two vectors
    const std::vector<double> vector_difference =
      vectorDifference(test_point, sphere_center);

    // sum of difference of two points squared
    double radius_squared = 0.0;
    for(int i_dof = 0; i_dof < test_point.size(); ++i_dof)
      radius_squared += vector_difference[i_dof] * vector_difference[i_dof];

    // check radius squared against cutoff
    if(radius_squared <= sphere_radius_squared)
      return true;
    
    // if here, point is not in sphere
    return false;
    
  }  

}