//
// Void.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#ifdef HAVE_ASSERT_H
#include <assert.h>
#else 
#error assert.h not found
#endif /* HAVE_ASSERT_H */

#ifdef HAVE_VECTOR
#include <vector>
#else
#ifdef HAVE_VECTOR_H
#include <vector.h>
#else
#error No vector or vector.h available
#endif // HAVE_VECTOR_H
#endif // HAVE_VECTOR

#include <iostream>
#include <utility>

#include "Node.h"
#include "DataTypes.h"
#include "Element.h"
#include "Void.h"
#include "Quasicontinua.h"
#include "MiscFunctions.h"
#include "C_Interface.h"

//
//
//

namespace quasicontinuum {
  //
  //
  //

  Void* Void::_instance = NULL;

  //
  // constructor
  //

  Void::Void()
  {

    //
    //
    //
    return;

  }

  //
  // destructor
  //

  Void::~Void()
  {

    //
    //
    //
    return;

  }

  //
  // new constructor
  //
  Void::Void(const int                      enable_void_default)
  {
    d_enableVoid = enable_void_default;

    d_voidCenter.resize(3);
    for(int i = 0; i< 3; i++)
    {
      d_voidCenter[i] = 0.0;
    }

    d_voidParameters.resize(1);
    d_voidParameters[0] = -1;

    d_voidRadius = -1.0;

    d_voidType = CIRCLE;

    return;
  }

  //
  // getInstance method
  //

  Void*
  Void::getInstance()
  {

    //
    // if not created, create
    //
    if(_instance == NULL){
      _instance = new Void(0);
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
  Void::destroyInstance()
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
  //  createVoid()
  //
  void Void::createVoid(struct node_list_t       *P_node_list,
            const struct lattice_t   *P_lattice,
            const int                 iQuasi)
  {
    double dr[3];
    double void_radius_sqr;
    int i_node;
    int j_node;
    int k_node;

    /**
      * clear void site cache 
      */
    if(iQuasi == 0)
      clearVoidSiteCache();

    /**
      * create void from different types
      */
    switch (d_voidType)
    {
      case CIRCLE:
      {

        /**
          * set void radius and void radius squared
          */
        d_voidRadius = d_voidParameters[0];
        void_radius_sqr = d_voidRadius*d_voidRadius;

        /**
          * vector to hold relative shift
          */
        std::vector<double> shift =
          Quasicontinua::getInstance()->getShift(iQuasi);
      
        /**
          * assign (void_mask = 1) to nodes inside void
          */
        for (i_node=0; i_node < P_node_list->number_nodes; i_node++)
        {
          dr[0]=P_node_list->nodes[i_node]->initial_position[0]+shift[0]
                -d_voidCenter[0];
          dr[1]=P_node_list->nodes[i_node]->initial_position[1]+shift[1]
                -d_voidCenter[1];
          dr[2]=P_node_list->nodes[i_node]->initial_position[2]+shift[2]
                -d_voidCenter[2];

          if ( dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] < void_radius_sqr)
          {
            /**
              * set void mask
              */
            P_node_list->nodes[i_node]->void_mask = 1;  

            /**
              * insert node into void cache
              */
            std::vector<int> site;
            site.resize(3);

            for(int i =0; i< 3; i++)
            {
              site[i] = P_node_list->nodes[i_node]->l[i];
            }

            insertVoidSiteCache(site, iQuasi);
          }

        }

      }     

      break;

      case VPIT:
      {
        /**
          * set void radius and void radius squared
          */
        const double sideLength = d_voidParameters[0];
        const double pitDepth = d_voidParameters[1];

        /**
          * vectors
          */
        double currentPosition[3];
        
        const double vertex1[3] = {d_voidCenter[0] - sideLength,
                  d_voidCenter[1], d_voidCenter[2]};


        const double vertex2[3] = {d_voidCenter[0] - 0.500*sideLength,
                  d_voidCenter[1] + 0.866*sideLength, d_voidCenter[2]};

        const double vertex3[3] = {d_voidCenter[0] + 0.500*sideLength,
                  d_voidCenter[1] + 0.866*sideLength, d_voidCenter[2]};

        const double vertex4[3] = {d_voidCenter[0] + sideLength,
                  d_voidCenter[1], d_voidCenter[2]};

        const double vertex5[3] = {d_voidCenter[0] + 0.500*sideLength,
                  d_voidCenter[1] - 0.866*sideLength, d_voidCenter[2]};

        const double vertex6[3] = {d_voidCenter[0] - 0.500*sideLength,
                  d_voidCenter[1] - 0.866*sideLength, d_voidCenter[2]};

        const double vertexBottom[3] = {d_voidCenter[0],
                  d_voidCenter[1], d_voidCenter[2] - pitDepth};
      
        /**
          * get shift
          */
        std::vector<double> shift = 
          Quasicontinua::getInstance()->getShift(iQuasi);

        /**
          * assign (void_mask = 1) to nodes inside void
          */
        for (i_node=0; i_node < P_node_list->number_nodes; i_node++)
        {
          currentPosition[0]=P_node_list->nodes[i_node]->initial_position[0]+shift[0];
          currentPosition[1]=P_node_list->nodes[i_node]->initial_position[1]+shift[1];
          currentPosition[2]=P_node_list->nodes[i_node]->initial_position[2]+shift[2];

          // create double voidCenter[3] to pass to checkPointInTetra()
          double voidCenter[3];
          for(int i=0; i<3; i++)
          {
            voidCenter[i]=d_voidCenter[i];
          }

          MiscFunctions* miscC = MiscFunctions::getInstance();

          if ((miscC->checkPointInTetra(vertex1, vertex2, voidCenter, 
                  vertexBottom, currentPosition) == INSIDE) ||
              (miscC->checkPointInTetra(vertex2, vertex3, voidCenter,
                  vertexBottom, currentPosition) == INSIDE) ||
              (miscC->checkPointInTetra(vertex3, vertex4, voidCenter,
                  vertexBottom, currentPosition) == INSIDE) ||
              (miscC->checkPointInTetra(vertex4, vertex5, voidCenter,
                  vertexBottom, currentPosition) == INSIDE) ||
              (miscC->checkPointInTetra(vertex5, vertex6, voidCenter,
                  vertexBottom, currentPosition) == INSIDE) ||
              (miscC->checkPointInTetra(vertex6, vertex1, voidCenter,
                  vertexBottom, currentPosition) == INSIDE))
          {
            /**
              * set void mask
              */
            P_node_list->nodes[i_node]->void_mask = 1;  

            /**
              * insert node into void cache
              */              
            std::vector<int> site;
            site.resize(3);

            for(int i =0; i< 3; i++)
            {
              site[i] = P_node_list->nodes[i_node]->l[i];
            }

            insertVoidSiteCache(site, iQuasi);
          }
        }
      } 

      break;

      case OVAL:

        d_print("Void type not implemented\n");
        exit(EXIT_FAILURE);

      break;

      default:
        d_print("Void type not implemented\n");
        exit(EXIT_FAILURE);
    }

    /**
      *
      */

    return;
  }

  //
  //  removeVoidElements()
  //
  void Void::removeVoidElements(struct element_list_t   *P_element_list) 
  {
    int  number_elements = 0;
    int i_elem;

    /**
      * loop over all elements
      */

    for (i_elem=0; i_elem < P_element_list->number_elements; i_elem++)
    {
      struct element_t *P_element=P_element_list->elements[i_elem];

      if ((P_element->node[0]->void_mask == 1) ||
          (P_element->node[1]->void_mask == 1) ||
          (P_element->node[2]->void_mask == 1) ||
          (P_element->node[3]->void_mask == 1)   )

        Element::getInstance()->makeInactiveElement(P_element_list, i_elem);
    }

    /**
      * compress element list
      */
  
    for (i_elem=0; i_elem < P_element_list->number_elements; i_elem++)
    {
      struct element_t *P_element=P_element_list->elements[i_elem];

      if (P_element != NULL)
      {
        P_element->number = number_elements;
        P_element_list->elements[number_elements] = P_element;

        number_elements++;
      }
    }

    P_element_list->number_elements = number_elements;

    return;
  } 

  //
  //  removeVoidNodes()
  //
  void Void::removeVoidNodes(struct node_list_t  *P_node_list)  
  {
    int i_node;
    int new_number_nodes = 0;

    for (i_node=0; i_node < P_node_list->number_nodes; i_node++)
    {
      if (P_node_list->nodes[i_node]->void_mask)
        Node::getInstance()->makeInactiveNode(P_node_list, i_node);
    }

    /**
      * clean node list
      */

    for (i_node=0; i_node < P_node_list->number_nodes; i_node++)
      if (P_node_list->nodes[i_node] != NULL) 
      {
        P_node_list->nodes[i_node]->number = new_number_nodes;
        P_node_list->nodes[new_number_nodes++] = P_node_list->nodes[i_node];
      }

    P_node_list->number_nodes = new_number_nodes;

    return;    
  }

  //
  //  insertVoidNodes()
  //  
  //  insert back removed void nodes; we assume that ALL nodes on the
  //  inactive node list have been removed using removeVoidNodes() above.
  //
  void Void::insertVoidNodes(struct node_list_t  *P_node_list) 
  {
    int i_node;

    while (P_node_list->number_inactive_nodes != 0)
      Node::getInstance()->makeActiveNode(P_node_list, 
                              P_node_list->number_nodes++);

    return;
  }

  //
  //  setVoid()
  //
  void Void::setVoid(const int enable_void_flag)  
  {
    d_enableVoid = enable_void_flag;

    return;
  }

  //
  //  isVoidEnable()
  //
  int Void::isVoidEnable()  
  {
    return d_enableVoid;
  }

  //
  //  setVoidParameters()
  //
  void Void::setVoidParameters(const double      center[3],
                  const enum void_t         voidType,
                  const int                 number_parameters,
                  std::vector<double>      param_list)
  {
    // enter the void center
    d_voidCenter[0] = center[0];
    d_voidCenter[1] = center[1];
    d_voidCenter[2] = center[2];

    // void Type
    d_voidType = voidType;

    // parameters
    d_voidParameters.resize(number_parameters);
    for(int i =0; i< number_parameters; i++)
    {
      d_voidParameters[i] = param_list[i];
    }

    return;
  }


  //
  //  insertVoidSiteCache()
  //
  void Void::insertVoidSiteCache(const std::vector<int> site,
                const int       iQuasi)
  {
    // check size of d_voidSiteCache
    if(iQuasi >= d_voidSiteCache.size())
      d_voidSiteCache.resize(iQuasi+1);

    //insert site to cache
    d_voidSiteCache[iQuasi].push_back(site);

    return;
  }

  //
  //  clearVoidCache()
  //
  void Void::clearVoidSiteCache()
  {
    d_voidSiteCache.clear();

    return;
  }

  //
  //  findSiteInVoidCache()
  //  0 : site not a void
  //  1 : site is void
  //
  int Void::findSiteInVoidCache(const std::vector<int> site,
            const int       iQuasi)  
  {
    //
    // loop over all sites in iQuasi
    //
    if(iQuasi < d_voidSiteCache.size())
    {
      for(unsigned int iSite = 0; iSite < d_voidSiteCache[iQuasi].size(); ++iSite)
      {
        //
        // see if site matches
        //
        if(site[0] == d_voidSiteCache[iQuasi][iSite][0] &&
           site[1] == d_voidSiteCache[iQuasi][iSite][1] &&
           site[2] == d_voidSiteCache[iQuasi][iSite][2])
           return 1;
      }
    }

    // site is not void
    return 0;
  }

  //
  //  printVoidCache()
  //  
  void Void::printVoidCache()  
  {
    //
    // loop over all void sites
    //
    for(int iQuasi = 0; iQuasi < d_voidSiteCache.size(); ++iQuasi)
      for(int iSite = 0; iSite < d_voidSiteCache[iQuasi].size(); ++iSite)
        std::cout<<"Quasi: "<<iQuasi
          <<" Site: "<<d_voidSiteCache[iQuasi][iSite][0]
          <<" "<<d_voidSiteCache[iQuasi][iSite][1]
          <<" "<<d_voidSiteCache[iQuasi][iSite][2]<<std::endl;

    //
    //
    //
    return;    
  }

  //
  //  voidAtomsMissed()
  //
  std::vector<std::pair<int,std::vector<int> > > 
  Void::voidAtomsMissed()  
  {
    //
    // atoms to return
    //
    std::vector<std::pair<int,std::vector<int> > > checkSites;

    //
    // loop over base lattice void sites
    //
    for(int iSite = 0; iSite < d_voidSiteCache[0].size(); ++iSite)
    {
      //
      // loop over non base lattices
      //
      for(int jQuasi = 1; jQuasi < d_voidSiteCache.size(); ++jQuasi)
      {
        //
        // site flag
        //
        bool addSite = true;

        //
        // loop over jSites
        //
        for(int jSite = 0; jSite < d_voidSiteCache[jQuasi].size(); ++jSite)
        {
          //
          // check sites
          //
          if(  d_voidSiteCache[0][iSite][0]==d_voidSiteCache[jQuasi][jSite][0] 
            && d_voidSiteCache[0][iSite][1]==d_voidSiteCache[jQuasi][jSite][1]
            && d_voidSiteCache[0][iSite][2]==d_voidSiteCache[jQuasi][jSite][2])
          {
            addSite = false;
            break;
          }
        }

        //
        // add site, if not found
        //
        if(addSite == true)
        {
          std::pair<int,std::vector<int> > addInfo =
            std::make_pair(jQuasi,d_voidSiteCache[0][iSite]);

          checkSites.push_back(addInfo);
        }
      }
    }

    return checkSites;    
  }

  //
  //  fixVoidBoundary()
  //
  void Void::fixVoidBoundary(struct node_list_t   *P_node_list,
    struct element_list_t    *P_element_list)
  {
    // call removeVoidELements() and removeVoidNodes()
    removeVoidElements(P_element_list);

    removeVoidNodes(P_node_list);     

    return;
  }

}