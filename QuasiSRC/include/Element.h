//
// Element.h
//

#if !defined(ELEMENT_H)
#define ELEMENT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "DataTypes.h"

//
//
//

namespace quasicontinuum{

  /**
   * @brief Singleton container for Quasicontinua instance
   */
  class Element {

    //
    // public methods
    //

  public:
 
   /**
    * @brief getInstance.
    */
   static Element * getInstance();

   /**
    * @brief destroyInstance.
    */
   static void destroyInstance();

  // ***************************************//
  // I = binelements.c
   
   //
   //  binElements()
   //
   void binElements(struct element_list_t      *P_element_list,
             struct node_list_t                 *P_node_list,
             enum mt_version_t                  mt_version);

  // ***************************************//
  // I = elements.c

   //
   // makeActiveElement()
   // 
   // get active element; if there is an inactive element return it
   // if not allocate a new one
   //
   void makeActiveElement(struct element_list_t    *P_element_list,
                  const int                        i_elem);

   //
   // makeInactiveElement()
   //
   void makeInactiveElement(struct element_list_t    *P_element_list,
                  const int                          i_elem);

  // ***************************************//
  // I = free_elements.c

   //
   // cleanElements()
   //
   void cleanElements(struct element_list_t   *P_element_list,
            enum mt_version_t        mt_version);

  // ***************************************//
  // I = site_element_map_cache.h

   //
   // checkHashTableAllocation()
   //
   void checkHashTableAllocation(const int iQuasi);

   //
   // hashTableClean()
   //
   void hashTableClean(const int iQuasi);

  // ***************************************//
  // I = site_element_map.h

  //
  //  initializeSiteElementMapData()
  //
  void initializeSiteElementMapData(const int iQuasi,
    struct element_list_t *P_element_list);

  //
  //  cleanSiteElementMapData()
  //
  void cleanSiteElementMapData(const int iQuasi);

  //
  //  locateSiteElement()
  //
  struct element_t * locateSiteElement(double                 *shape,
        const int               l[3],
        struct lattice_t *P_lattice,
        const int               iQuasi);

  // ***************************************//
  // I = delaunay_3d.h  

  //
  //  computeElementNumberSites()
  //
  void computeElementNumberSites(struct element_list_t  *P_element_list, 
              const struct lattice_t *P_lattice);

  // ***************************************//
  // I = triangulation_client.h

  //
  //  createElementsFromMeshData()
  //
  void 
  createElementsFromMeshData(struct element_list_t *P_element_list,
            struct node_list_t      *P_node_list,
            struct lattice_t        *P_lattice,
            struct mesh_data_t      *P_mesh_data);

  // ***************************************//
  // I = strain.c  

  //
  //  computeDeformGradInfinityNorm()
  //
  double computeDeformGradInfinityNorm(double F[3][3]);

  //
  //  computeDeformGrad()
  //
  void computeDeformGrad(struct element_t  *P_element);

  //
  //  computeElementStrain()
  //
  double computeElementStrain(struct element_t *P_element,
             int         i_load);


    //
    // private methods
    //

  private:

    /**
     * @brief Constructor.
     */
    Element();

    /**
     * @brief Copy constructor.
     */
    Element(Element const&);

    /**
     * @brief Assignment operator.
     */
    const Element & operator=(const Element &);

    /**
     * @brief Destructor.
     */
    ~Element();

    //
    // private data types
    //
  private:

    static Element*                        _instance;
  };

}

#endif // ELEMENT_H
