//
// MiscFunctions.h
//

#if !defined(MISCFUNCTIONS_H)
#define MISCFUNCTIONS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

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

#include "DataTypes.h"


#if !defined(MAX) && !defined(MIN)
#  define MIN(a, b) ((a) < (b) ? (a) : (b))
#  define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif /* !MAX && !MIN */

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class MiscFunctions {

    //
    // public methods
    //

  public:
 
    /**
     * @brief getInstance.
     */
    static MiscFunctions * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    //
    //  freeE()
    //
    void * freeE(void * p);

    //
    //  tetraVolume()
    //
    double tetraVolume(const double   a[3],
          const double   b[3],
          const double   c[3],
          const double   d[3]);

    //
    //  getUnitVector()
    //
    int getUnitVector( const double  v[3],
      double u[3] );

    //
    //  getDistanceBetPoints()
    //
    double getDistanceBetPoints(const double a[3],
            const double     b[3]);

    //
    //  getDistanceSqrBetweenPoints()
    //
    double getDistanceSqrBetweenPoints(const double  a[3],
            const double      b[3]);

    //
    //  computeDifferenceOfVectorAndAbsValue()
    //
    //  writes the difference of vector to third argument
    //  returns sqrt of abs difference
    //
    double computeDifferenceAndAbsValue(const double a[3],
              const double           b[3],
              double                 ab[3]);

    //
    //  checkPointInTetra()
    //
    enum loc_t  checkPointInTetra(const double v1[3], 
             const double v2[3], 
             const double v3[3],
             const double v4[3], 
             const double p[3]);

    //
    //  checkSiteInTetra()
    //
    enum loc_t checkSiteInTetra(const double            v1[3], 
           const double            v2[3], 
           const double            v3[3],
           const double            v4[3],
           const int               l[3],
           struct lattice_t *P_lattice,
           const int                iQuasi);

    //
    //  findTetraCenter()
    //
    double findTetraCenter(const double vertice1[3], 
         const double vertice2[3],
         const double vertice3[3], 
         const double vertice4[3],
         double       center[3]);

    //
    //  computeDeterminant()
    //
    double computeDeterminant(const double      a[3],
          const double      b[3],
          const double      c[3]);

    //
    //  solveLinearThreeD()
    //
    int solveLinearThreeD(const double    a[3],
              const double    b[3],
              const double    c[3],
              const double    d[3],
              double          x[3]);

    /**
     * @brief Compute cross product of two 3x3 vectors.
     *
     * @param a 1st vector.
     * @param b 2nd vector.
     *
     * @return cross product result
     */
    std::vector<double> crossProduct3x3(std::vector<double> a,
                 std::vector<double> b);

    /**
     * @brief Compute determinant of 3x3 matrix.
     *
     * @param a matrix
     *
     * @return determinant
     */
    double determinant3x3(std::vector<std::vector<double> > a);

    /**
     * @brief Compute L2Norm of vector.
     *
     * @param a vector.
     *
     * @return norm
     */
    double L2Norm(std::vector<double> a);

   /**
     * @brief Compute dot product of two vectors.
     *
     * @param a vector.
     * @param b vector.
     *
     * @return dotProduct
     */
    double dotProduct(const std::vector<double> &a,
           const std::vector<double> &b);

    /**
     * @brief Area of triangle.
     *
     * @param vertexA Vertex A of triangle.
     * @param vertexB Vertex B of triangle.
     * @param vertexC Vertex C of triangle.
     *
     * @return area
     */
    double TriangleArea(std::vector<double> vertexA,
             std::vector<double> vertexB,
             std::vector<double> vertexC);

    /**
     * @brief Center of triangle.
     *
     * @param vertexA Vertex A of triangle.
     * @param vertexB Vertex B of triangle.
     * @param vertexC Vertex C of triangle.
     *
     * @return center
     */
    std::vector<double> TriangleCenter(std::vector<double> vertexA,
                std::vector<double> vertexB,
                std::vector<double> vertexC);

    /**
     * @brief Sphere ray intersection.
     *
     * @param center Center of sphere.
     * @param nodeInSphere Node that is in sphere.
     * @param nodeOutsideSphere Node that is outside sphere.
     * @param radiusSquared Radius of sphere squared
     *
     * @return location
     */
    std::vector<double> getSphereRayIntersection(std::vector<double> center,
              std::vector<double> nodeInSphere,
              std::vector<double> nodeOutsideSphere,
              double              radiusSquared);

    /**
     * @brief Vector A - Vector B.
     *
     * @param vectorA
     * @param vectorB
     *
     * @return difference
     */
    std::vector<double> vectorDifference(const std::vector<double> & vectorA,
            const std::vector<double> & vectorB);

    /**
     * @brief Vector A + Vector B.
     *
     * @param vectorA
     * @param vectorB
     *
     * @return addition
     */
    std::vector<double> vectorAddition(const std::vector<double> & vectorA,
                const std::vector<double> & vectorB);

    /**
     * @brief Vector A + Vector B + Vector C.
     *
     * @param vectorA
     * @param vectorB
     * @param vectorC
     *
     * @return addition
     */
    std::vector<double> vectorAddition(const std::vector<double> & vectorA,
                const std::vector<double> & vectorB,
                const std::vector<double> & vectorC);

    /**
     * @brief Vector A * scale.
     *
     * @param vectorA
     * @param scale
     *
     * @return scaled vector
     */
    std::vector<double> vectorScale(std::vector<double> vectorA,
             const double scale);

    /**
     * @brief Triangle sphere intersection.
     *
     * @param vertexA
     * @param vertexB
     * @param vertexC
     * @param sphereCenter
     * @param sphereRadiusSquared
     *
     * @return bool
     */
    bool triangleSphereIntersection(const std::vector<double> &vertexA,
             const std::vector<double> &vertexB,
             const std::vector<double> &vertexC,
             const std::vector<double> &sphereCenter,
             const double &sphereRadiusSquared);

    /**
     * @brief Edge sphere intersection.
     *
     * @param vertexA
     * @param vertexB
     * @param sphereCenter
     * @param intersection point 1 return value
     * @param intersection point 2 return value
     *
     * @return bool
     */
     bool edgeSphereIntersection(const std::vector<double> &vertexA,
          const std::vector<double> &vertexB,
          const std::vector<double> &sphereCenter,
          const double &sphereRadiusSquared,
          std::vector<double> &intersectPoint1,
          std::vector<double> &intersectPoint2);

     /**
     * @brief Is point inside sphere
     *
     * @param test_point Point to test whether it is in sphere.
     * @param sphere_center Center of sphere.
     * @param sphere_radius_squared Radius squared of sphere.
     *
     * @return bool
     */
     bool IsPointInSphere(const std::vector<double> test_point,
         const std::vector<double> sphere_center,
         const double sphere_radius_squared);    


    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    MiscFunctions();

    /**
     * @brief Copy constructor.
     */
    MiscFunctions(MiscFunctions const&);

    /**
     * @brief Assignment operator.
     */
    const MiscFunctions & operator=(const MiscFunctions &);

    /**
     * @brief Destructor.
     */
    ~MiscFunctions();

    //
    // private data types
    //
  private:

    static MiscFunctions*                        _instance;
  };

}

#endif // MISCFUNCTIONS_H