#if !defined(INDENT_H)
#define INDENT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <pthread.h>
#include "DataTypes.h"

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Indent {

public:
  /**
   * @brief getInstance.
   */
  static Indent *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  void setNanoIndentation(const int enable_nano_indentation_flag);

  int isIndentEnable();

  void moveIndentor(struct indentor_t *P_indentor);

  void indentorNodesInteraction(int iQuasi);

  void findIndentorContactRadius(double *P_R, double *P_r,
                                 const struct node_list_t *P_node_list,
                                 const struct indentor_t *P_indentor,
                                 const int iQuasi,
                                 enum mt_version_t mt_version);

  void fixIndentBoundary(struct node_list_t *P_node_list,
                         struct lattice_t *P_lattice);

private:
  Indent();
  Indent(const int indent_flag);
  Indent(Indent const &);
  const Indent &operator=(const Indent &);
  ~Indent();

  void IndentorNodePotential(double *force, double *energy, const double r,
                             const struct indentor_t *P_indentor);

private:
  static Indent *_instance;
  int d_enableNanoIndentation;
};
}

#endif // INDENT_H
