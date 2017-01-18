#if !defined(QUASICONTINUA_H)
#define QUASICONTINUA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "Quasicontinuum.h"
#include <vector>

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Quasicontinua {

public:
  typedef unsigned int id_type;
  typedef unsigned int size_type;

public:
  /**
   * @brief getInstance.
   */
  static Quasicontinua *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  /**
   * @brief Insert a Quasicontinuum.
   *
   * @param num_quasi_insert Number of Quasicontinuum instances to insert.
   * @param temperature : temperature of Quasicontinuum.
   * @param ac Command line options.
   * @param av Command line options.
   *
   * @return Size of quasicontinua.
   */
  void insert(int num_quasi, int ac, char *av[]);

  /**
   * @brief Get an instance of Quasicontinuum.
   *
   * @param A Quasicontinuum id.
   *
   * @return A Quasicontinuum instance corresponding to the id.
   */
  const Quasicontinuum &get(id_type id) const;
  Quasicontinuum &get(id_type id);

  /**
   * @brief Get an instance of Quasicontinuum/
   *        Current Id safe version of get()
   *
   * @param A Quasicontinuum id.
   *
   * @return A Quasicontinuum instance corresponding to the id.
   */
  const Quasicontinuum &getQuasi(id_type id) const;
  Quasicontinuum &getQuasi(id_type id);

  //    /**
  //     * @brief Get the temperature of system.
  //     *
  //     * @return The temperature of Quasicontinua.
  //     */
  //    double getTemperature() const;

  /**
   * @brief Get the number of Quasicontinuum instances.
   *
   * @return The number of Quasicontinuum instances.
   */
  size_type size() const;

  /**
   * @brief Get the current Quasicontinuum id.
   *
   * @return Id of the current Quasicontinuum.
   */
  id_type getCurrentId() const;

  /**
   * @brief Insert a shift vector.
   *
   * @param quasicontinuumId Id of Quasicontinuum
   * @param globalVector Shift of Quasicontinuum w.r.t. global coordinates.
   */
  void insertShift(int quasicontinuumId, std::vector<double> shiftVector,
                   int coreShell);

  /**
    * @brief Get relative shift between 2 Quasicontinuum.
    *
    * @param iQuasi 1st Quasicontinuum ID.
    * @param jQuasi 2nd Quasicontinuum ID.
    *
    * @return Relative vector between Quasicontinuum.
    */
  std::vector<double> getRelativeShift(unsigned int iQuasi,
                                       unsigned int jQuasi);

  /**
     * @brief Get global shift vector of Quasicontinuum.
     *
     * @param QuasicontinuumID Quasicontinuum ID.
     *
     * @return Global shift vector of Quasicontinuum.
     */
  std::vector<double> getShift(unsigned int QuasicontinuumID);

  /**
     * @brief Get core shell interaction between Quasicontinuum.
     *
     * @param QuasicontinuumID Quasicontinuum ID.
     *
     * @return Core shell info.
     */
  int getCoreShell(unsigned int QuasicontinuumID);

  //
  //  setTemperature()
  //
  void setTemperature(const double temperature);

  //
  //  getTemperature()
  //
  double getTemperature();

  //
  //  getSigmaVector()
  //
  std::vector<double> getSigmaVector(void);

  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  Quasicontinua();

  /**
   * @brief Copy constructor.
   */
  Quasicontinua(Quasicontinua const &);

  /**
   * @brief Assignment operator.
   */
  const Quasicontinua &operator=(const Quasicontinua &);

  /**
   * @brief Destructor.
   */
  ~Quasicontinua();

private:
  static Quasicontinua *_instance;
  std::vector<quasicontinuum::Quasicontinuum> d_quasicontinua;
  mutable id_type d_quasicontinuum_id;
  std::vector< std::pair< int, std::vector<double> > > d_shifts;
  std::vector<int> d_coreShell;
  double d_temperature;
};
}

#endif // QUASICONTINUA_H
