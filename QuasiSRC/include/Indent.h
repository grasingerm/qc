//
// Indent.h
//

#if !defined(INDENT_H)
#define INDENT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#include "DataTypes.h"

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class Indent {

    //
    // public methods
    //

  public:
 
    /**
     * @brief getInstance.
     */
    static Indent * getInstance();

    /**
     * @brief destroyInstance.
     */
    static void destroyInstance();

    //
    //  setNanoIndentation()
    //
    void setNanoIndentation(const int  enable_nano_indentation_flag);

    //
    //  isIndentEnable()
    //
    int isIndentEnable();

    //
    //  moveIndentor()
    //
    void moveIndentor( struct indentor_t*       P_indentor);

    //
    //  indentorNodesInteraction()
    //
    void indentorNodesInteraction(int iQuasi);

    //
    //  findIndentorContactRadius()
    //
    void findIndentorContactRadius(double   *P_R,
                    double                   *P_r,
                    const struct node_list_t *P_node_list,
                    const struct indentor_t  *P_indentor,
                    const int                iQuasi,
                    enum mt_version_t         mt_version);

    //
    //  fixIndentBoundary()
    //
    void fixIndentBoundary(struct node_list_t     *P_node_list,
        struct lattice_t      *P_lattice);
    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    Indent();

    //
    //  Indent()
    //
    Indent(const int indent_flag);

    /**
     * @brief Copy constructor.
     */
    Indent(Indent const&);

    /**
     * @brief Assignment operator.
     */
    const Indent & operator=(const Indent &);

    /**
     * @brief Destructor.
     */
    ~Indent();

    //
    //  IndentorNodePotential()
    //
    void IndentorNodePotential(double                  *force,
            double                  *energy,
            const double             r,
            const struct indentor_t *P_indentor);

    //
    // private data types
    //
  private:

    static Indent*                        _instance;
    //global fixed parameters
    int                                   d_enableNanoIndentation;
  };

}

#endif // INDENT_H
