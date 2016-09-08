//
// Lattice.h
//

#if !defined(LATTICE_H)
#define LATTICE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#include <vector>

#include "DataTypes.h"

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class Lattice {

    //
    // public methods
    //

  public:
 
    /**
     * @brief getInstance.
     */
    static Lattice * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    //
    //  setInitialConfigFlag()
    //
    void setInitialConfigFlag(int    flag);    

    //
    //  getShell()
    //
    struct shell_t * getShell(int            shell_number, 
            const struct lattice_t          *P_lattice);

    //
    //  getShellRadius()
    //
    double getShellRadius(struct lattice_t *P_lattice,
                   const int               l[3],
                   const int               iQuasi);

    //
    //  initializeShells()
    //
    int initializeShells(struct lattice_t lattice);

    //
    //  findClosestSite()
    //
    void findClosestSite(int                l[3],
                  double                    X[3],
                  const struct lattice_t    *P_lattice,
                  const int                  iQuasi);

    //
    //  isSiteInsideLattice()
    //
    int isSiteInsideLattice(const struct lattice_t *P_lattice,
                  const int               l[3],
                  const int               iQuasi);

    //
    //  isLatticeSite()
    //
    int isLatticeSite(const int               l[3],
               const struct lattice_t *P_lattice,
               const int               iQuasi);

    //
    //  compareSites()
    //
    int compareSites(const int l1[3], const int l2[3]);


    //
    //  getSiteInitialPosition()
    //
    void getSiteInitialPosition(double                  r[3],
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);

    //
    //  getSiteInitialFrequency()
    //
    void getSiteInitialFrequency(double                 &w,
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);    

    //
    //  getSiteInitialTemperature()
    //
    void getSiteInitialTemperature(double                  &T,
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);

    //
    //  getSiteInitialState() : return (x,y,z,w)
    //
    void getSiteInitialState(double                  s[4],
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);

    //
    //  getSiteInitialStateInterpolate()
    //
    void getSiteInitialStateInterpolate(double                  S[4],
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);    


    //
    //  getSiteCurrentPosition()
    //
    void getSiteCurrentPosition(double                  r[3],
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);

    //
    //  getSiteICurrentFrequency()
    //
    void getSiteCurrentFrequency(double                  &w,
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);    

    //
    //  getSiteCurrentTemperature()
    //
    void getSiteCurrentTemperature(double                  &T,
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);

    //
    //  getSiteCurrentState() : return (x,y,z,w)
    //
    void getSiteCurrentState(double                  s[4],
           const int               l[3],
           struct lattice_t *P_lattice,
           const int               iQuasi);


    //
    //  getSiteCurrentStateAndNodeInfo()
    //
    //  Note : this takes std::vector type data as input
    //  This function writes position+freq dof on input data. it also 
    //  writes NodeInfo on input data
    //
    void getSiteCurrentStateAndNodeInfo(
      std::pair< std::pair< std::vector<int>, std::vector<double> >, std::vector<int> > &iC_data, 
        struct lattice_t   *P_lattice,
        const int           iQuasi);

    //
    //  getSiteInitialState() : return (x,y,z,w,t)
    //
    void getSiteCurrentStateNew(std::vector<double> &state,
      std::vector<int>    site,
      const int iQuasi);    

    //
    //  getSiteInitialPositionNew()
    //
    void getSiteInitialPositionNew(std::vector<double> &X,
      std::vector<int> site,
      const int iQuasi);

    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    Lattice();

    /**
     * @brief Copy constructor.
     */
    Lattice(Lattice const&);

    /**
     * @brief Assignment operator.
     */
    const Lattice & operator=(const Lattice &);

    /**
     * @brief Destructor.
     */
    ~Lattice();

    //
    //  GenerateNewShell() : used in initializeShells()
    //
    int GenerateNewShell(int               shell_number,
       struct lattice_t  lattice);

    //
    //  FillShell()
    //
    int FillShell(struct shell_t*      data_shell,
        int                            shell_number,
        struct lattice_t               lattice);

    //
    //  AddNewSite()
    //
    int AddNewSite(struct shell_t*      data_shell,
        int                             l1,
        int                             l2,
        int                             l3);    

    //
    // GetLatticeCoordinates()
    //
    void GetLatticeCoordinates(double                  r[3],
            double                  X[3],
            const struct lattice_t *P_lattice);

    //
    // private data types
    //
  private:

    static Lattice*                        _instance;
    static int                             d_bcc_number_shells;
    static int                             d_fcc_number_shells;
    static int                             d_fcc_p;
    static int                             d_bcc_p;
    static int                             d_initialize_bcc;
    static int                             d_initialize_fcc;
    std::vector<struct shell_t>            d_fcc_shells;
    std::vector<struct shell_t>            d_bcc_shells;
    static int                                    d_initialConfigFlag;
  };

}

#endif // LATTICE_H
