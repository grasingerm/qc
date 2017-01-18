#if !defined(LATTICE_H)
#define LATTICE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <pthread.h>
#include <vector>

#include "DataTypes.h"

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Lattice {

public:
  /**
   * @brief getInstance.
   */
  static Lattice *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  void setInitialConfigFlag(int flag);

  struct shell_t *getShell(int shell_number, const struct lattice_t *P_lattice);

  double getShellRadius(struct lattice_t *P_lattice, const int l[3],
                        const int iQuasi);

  int initializeShells(struct lattice_t lattice);

  void findClosestSite(int l[3], double X[3], const struct lattice_t *P_lattice,
                       const int iQuasi);

  int isSiteInsideLattice(const struct lattice_t *P_lattice, const int l[3],
                          const int iQuasi);

  int isLatticeSite(const int l[3], const struct lattice_t *P_lattice,
                    const int iQuasi);

  int compareSites(const int l1[3], const int l2[3]);

  void getSiteInitialPosition(double r[3], const int l[3],
                              struct lattice_t *P_lattice, const int iQuasi);

  void getSiteInitialFrequency(double &w, const int l[3],
                               struct lattice_t *P_lattice, const int iQuasi);

  void getSiteInitialTemperature(double &T, const int l[3],
                                 struct lattice_t *P_lattice, const int iQuasi);

  void getSiteInitialState(double s[4], const int l[3],
                           struct lattice_t *P_lattice, const int iQuasi);

  void getSiteInitialStateInterpolate(double S[4], const int l[3],
                                      struct lattice_t *P_lattice,
                                      const int iQuasi);

  void getSiteCurrentPosition(double r[3], const int l[3],
                              struct lattice_t *P_lattice, const int iQuasi);

  void getSiteCurrentFrequency(double &w, const int l[3],
                               struct lattice_t *P_lattice, const int iQuasi);

  void getSiteCurrentTemperature(double &T, const int l[3],
                                 struct lattice_t *P_lattice, const int iQuasi);

  void getSiteCurrentState(double s[4], const int l[3],
                           struct lattice_t *P_lattice, const int iQuasi);

  void getSiteCurrentStateAndNodeInfo(
      std::pair<std::pair<std::vector<int>, std::vector<double>>,
                std::vector<int>> &iC_data,
      struct lattice_t *P_lattice, const int iQuasi);

  void getSiteCurrentStateNew(std::vector<double> &state, std::vector<int> site,
                              const int iQuasi);

  void getSiteInitialPositionNew(std::vector<double> &X, std::vector<int> site,
                                 const int iQuasi);

private:
  Lattice();
  Lattice(Lattice const &);
  const Lattice &operator=(const Lattice &);
  ~Lattice();
  
  int GenerateNewShell(int shell_number, struct lattice_t lattice);
  int FillShell(struct shell_t *data_shell, int shell_number,
                struct lattice_t lattice);
  int AddNewSite(struct shell_t *data_shell, int l1, int l2, int l3);
  void GetLatticeCoordinates(double r[3], double X[3],
                             const struct lattice_t *P_lattice);

private:
  static Lattice *_instance;
  static int d_bcc_number_shells;
  static int d_fcc_number_shells;
  static int d_fcc_p;
  static int d_bcc_p;
  static int d_initialize_bcc;
  static int d_initialize_fcc;
  std::vector<struct shell_t> d_fcc_shells;
  std::vector<struct shell_t> d_bcc_shells;
  static int d_initialConfigFlag;
};
}

#endif // LATTICE_H
