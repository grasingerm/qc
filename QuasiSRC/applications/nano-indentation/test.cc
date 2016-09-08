//
// Main complex non-local quasicontinuum file
//
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
/* no  _REENTRANT */
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found
#endif /* HAVE UNISTD_H */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#error limits.h not found
#endif /* HAVE_LIMITS_H */

#include <cmath>
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <vector>

#include "Input.h"

#include "CGNonLinearSolver.h"
#include "CrossNeighborList.h"
#include "CreateMesh.h"
#include "ForceEnergyCalculation.h"
#include "PairPotentials.h"
#include "Quasicontinua.h"
#include "QuasicontinuaForceFunction.h"
#include "Quasicontinuum.h"
#include "QuadraturePoints.h"
#include "Output.h"
#include "Electrostatics.h"
#include "Lattice.h"
// #include "C_Interface.h"
#include "Error.h"
#include "MiscFunctions.h"
#include "Element.h"
#include "DataTypes.h"
#include "RunData.h"
#include "read_pipe.h"

#if !defined(lint)
static char rcsid[]="$Id: qc3d.c,v 1.10 2002/03/07 23:53:08 knap Exp $";
#endif // !lint

//
//  tetrahedronVolume()
//
double tetrahedronVolume(double a[3],
          double b[3],
          double c[3],
          double d[3])
{
  //
  // call MiscFunctions to compute volume
  //
  using namespace quasicontinuum;

  double volume = MiscFunctions::getInstance()->tetraVolume(a, b, c, d);

  //
  // make sure volume is positive
  //
  if (volume < 0.0)
    volume = volume * -1.0;

  //
  //
  //
  return volume;
}

//
//  computeDeformGrad()
//
void computeDeformGrad(void)
{
  using namespace quasicontinuum;

  //
  //  computing deformation graident for each element and corresponding
  //  unit cell volume
  //
  Quasicontinua * quasicontinua = Quasicontinua::getInstance();
  Quasicontinuum zQuasi = quasicontinua->getQuasi(0);
  struct element_list_t * P_element_list = &(zQuasi.getElementList());
  struct quasicontinuum::lattice_t * P_lattice = &(quasicontinua->getQuasi(0).getLattice());

  Element * elementC = Element::getInstance();

  //
  //  first compute unit cell volume, assuming crystal lattice coordinates
  //
  double unitCellVol = 0.0;
  unitCellVol =
     (P_lattice->a1[1] * P_lattice->a2[2] - P_lattice->a1[2] * P_lattice->a2[1]) * P_lattice->a3[0] +
  (P_lattice->a1[2] * P_lattice->a2[0] - P_lattice->a1[0] * P_lattice->a2[2]) * P_lattice->a3[1] +
  (P_lattice->a1[0] * P_lattice->a2[1] - P_lattice->a1[1] * P_lattice->a2[0]) * P_lattice->a3[2];

  if(unitCellVol < 0.0)
    unitCellVol = unitCellVol*(-1.0);

  double unitCellVol_1 = 0.0;
  double unitCellVol_2 = 0.0;

  for(int iElem=0; iElem < P_element_list->number_elements; iElem++)
  {
    struct element_t * P_element = P_element_list->elements[iElem];

    //
    //  compute deformation gradient
    //
    elementC->computeDeformGrad(P_element);

    //
    //  Method 1 : Following the methodology in Electrostatics.cc
    //
    double initialTetraVol =
    tetrahedronVolume(P_element->node[0]->initial_position,
      P_element->node[1]->initial_position,
      P_element->node[2]->initial_position,
      P_element->node[3]->initial_position);

    double currentTetraVol =
    tetrahedronVolume(P_element->node[0]->position,
      P_element->node[1]->position,
      P_element->node[2]->position,
      P_element->node[3]->position);

    double changeInVol = currentTetraVol / initialTetraVol;

    unitCellVol_1 = changeInVol * unitCellVol;

    //
    //  Method 2 : Multiplying initial volume by determinant of F
    //
    double a[3];
    double b[3];
    double c[3];

    for(int dof =0; dof < 3; dof++)
    {
      a[dof] = P_element->F[0][dof];
      b[dof] = P_element->F[1][dof];
      c[dof] = P_element->F[2][dof];
    }

    double detF = MiscFunctions::getInstance()->computeDeterminant(a, b, c);

    unitCellVol_2 = detF * unitCellVol;

    d_print("Element : %d   %lf     %lf     %lf\n",
      iElem, unitCellVol, unitCellVol_1, unitCellVol_2);

    d_print("     %lf       %lf       %lf\n", a[0], a[1], a[2]);
    d_print("     %lf       %lf       %lf\n", b[0], b[1], b[2]);
    d_print("     %lf       %lf       %lf\n", c[0], c[1], c[2]);
    d_print("\n");
  }

  return;
}

int 
main(int argc, char **argv)
{
  // timing data
  timeval t1,t2;
  double elapsedTime;

  // start timer
  gettimeofday(&t1,NULL);

  // quasi specific namespace
  using namespace quasicontinuum;

  // load complex inputs
  int    NUM_QUASI;
  int    NUMBER_LOAD;
  int    boundaryFlag;
  std::vector<std::vector<double> > loadDeformation;
  double cgTolerance;
  int    cgMaxIterations;
  int    cgDebugLevel;
  double cgLineSearchTolerance;
  int    cgLineSearchMaxIterations;
  int    getDefects;
  int    addNodeCutoff;
  int    addNodeFlag;
  int    addNodeRestartNum;
  std::vector<double> defectCenter;
  std::vector<double> defectBox;
  double defectLat;
  std::vector<std::pair<double, std::vector<double> > > charges;
  double temperature;
  int    statistics;
  int    quadMethod;
  int    output_dir_flag;
  char*  output_directory;
  int    test_flag;
  int    minMethod;
  int    minMethodMaxIterations;
  double minMethodTolerance;
  int    initialConfigFlag;

  // get path_max for allocating memory to the char * variables
  long int path_max;
  path_max = pathconf(".", _PC_PATH_MAX);   
  output_directory = (char *) malloc(path_max*sizeof(char));  

  printf("Main : Reading main input parameters\n");
  Input::getInstance()->mainInput(NUM_QUASI,
            NUMBER_LOAD,
            boundaryFlag,
            loadDeformation,
            cgTolerance,
            cgMaxIterations,
            cgDebugLevel,
            cgLineSearchTolerance,
            cgLineSearchMaxIterations,
            getDefects,
            addNodeCutoff,
            addNodeFlag,
            addNodeRestartNum,
            defectCenter,
            defectBox,
            defectLat,
            charges,
            temperature,
            statistics,
            quadMethod,
            output_dir_flag,
            output_directory,
            test_flag,
            minMethod,
            minMethodMaxIterations,
            minMethodTolerance,
            initialConfigFlag,
            argc,
            argv);
  printf("Main : Done\n");

  //
  //  check if it is quasi harmonic approximation
  //
  if(quadMethod == 1)
  {
    if(minMethod == 2)
      d_print("*** frequency minimization with quasi harmonic approximation ***\n");
    else
      d_print("*** quasi harmonic approximation ***\n");

    if(PairPotentials::getInstance()->d_EAMFlag == true)
    {
      d_print("EAM with quasi harmonic approximation is not implemented\n");
      d_print("exiting the code\n");
      exit(EXIT_SUCCESS);
    }
  }

  FILE * debugFile;
  debugFile = Output::getInstance()->d_debugFile;

  //
  // parameters used in remeshing
  //
  int   total_nodes_added = 0;
  std::vector<int>  nodes_added;
  //
  // initialize nodes_added
  //
  for(int i = 0; i < NUM_QUASI; i++)
    nodes_added.push_back(0);
  
  // get 0 Quasi -0 Quasi interaction number
  int potentialNumber = PairPotentials::getInstance()->doQuasicontinuumInteract(0,0); 

  // define rebuild_cutoff
  const double rebuild_cutoff =  PairPotentials::getInstance()->getCutoffRadiusNeighborList(potentialNumber);

  //
  // if getDefects = 1 ---> exit, as it is not implemented currently.
  //
  if(getDefects == 1)
  {
    d_print("getDefects == %d\n", getDefects);
    d_print("get defects is not implemented\n");
    d_print("change value of getDefects flag to 0\n");

    exit(EXIT_SUCCESS);
  }

  //
  // get Quasicontinua instance
  //
  Quasicontinua * quasicontinua = Quasicontinua::getInstance();

  // create Quasicontinuum objects in Quasicontinua class
  d_print("Main : Creating %i number of Quasicontinuum \n", NUM_QUASI);

  quasicontinua->insert(NUM_QUASI, argc, argv);

  d_print("Main : Done\n");

  //
  //  Dump the information of problem to file
  //
  Input::getInstance()->systemInfoOutput();


  //
  //  test flag
  //
  //  -1 - do nothing
  //  0 - get force and energy from P.E., Entropy, Electrostatics
  //      separately and output the data
  //  1 - test crossneighborlist
  //
  if(test_flag == -1)
  {
    d_print("Nothing to do in test debug of code\n");
    d_print("exiting the code\n");
    exit(EXIT_SUCCESS);

    return EXIT_SUCCESS;
  }

  switch(test_flag)
  {
    case 0:
    {
      //
      //  compute the forces and energy for each cluster site from
      //  P.E., Entropy, Electrostatics separately and output
      //  the data
      //
      d_print("Main : Computing force energy separately for interatomic potential, entropy and electrostatics\n");

      //
      //  call create elements : compute forces in createElementsCentral()
      //  only if also want output of total forces on nodes.
      //
      //  0 : compute forces in Create Mesh
      //  1 : do not compute forces in Create Mesh
      //
      //  DO NOT USE computeForce=1. IT GIVES SEG FAULT ERROR.
      //
      int computeForce = 0;
      d_print("Main : Meshing all Quasis\n");
      CreateMesh * createMesh = CreateMesh::getInstance();
      createMesh->createElementsCentral(computeForce);
      d_print("Main : Done\n");

      if(computeForce == 0)
      {
        //
        //  output
        //
        unsigned int flags;
        flags |= NODE_OUTPUT_FLAG;

        // 
        //  create output_flag_vector
        //
        std::vector<int> output_flag_vector(1,0);

        d_print("Main : Performing output for all quasis\n");

        Output::getInstance()->performOutputCentral(flags, output_flag_vector);
        d_print("Main : Done\n");
      }

      //
      //  Loop over each quasi, compute forces individually for pairwise,
      //  entropy, electrostatics, and store it in data handle
      //
      //  Note : reset the cluster force data in ForceEnergyCalculation
      //         before calling function for force calculation
      //
      ForceEnergyCalculation * forceC = ForceEnergyCalculation::getInstance();

      // force data
      // pair<vector<double>, vector<double> >
      // --> pair<vector of forces, vector of energy>
      //
      //  vector of forces
      //  --> <0 to 3 for  EAM (if enabled),
      //       4 to 7 for Pairwise,
      //       8 for entropy freq force,
      //       9 to 12 for electrostatics force (if enabled)>
      //
      //  vector of energy
      //  --> < energy for EAM(if enabled),
      //        energy for pairwise,
      //        energy for entropy,
      //        energy for electrostatics (if enabled)>
      //
      std::vector< std::vector< std::vector< std::pair< std::vector< double >, std::vector<double> > > > > test_cluster_force_energy;

      std::vector< std::vector< std::vector< std::pair< std::vector< double >, double > > > > test_total_cluster_force_energy;

      //  node force energy data
      std::vector< std::vector< std::pair< std::vector<double>, std::vector<double> > > > test_node_force_energy;
      std::vector< std::vector< std::pair< std::vector<double>, double > > > test_total_node_force_energy;

      //
      //  call ForceEnergyCalculation function which computes electric field
      //  and also initializes node dat
      //
      d_print("Main : Computing Forces\n");
      int rebuild_neighbor_flag = 0;
      int rebuild_cluster_flag = 1;
      forceC->testForceEnergyCalculation(rebuild_neighbor_flag,
        rebuild_cluster_flag,
        test_cluster_force_energy,
        test_node_force_energy,
        test_total_cluster_force_energy,
        test_total_node_force_energy);
      d_print("Main : Done\n");

      // perform output by calling testEnergyOutput()
      d_print("Main : Performing output\n");
      Output::getInstance()->testEnergyOutput(&test_cluster_force_energy,
        &test_node_force_energy,
        &test_total_cluster_force_energy,
        &test_total_node_force_energy);
      d_print("Main : Done\n");

      //
      //  use below to check if results are same when
      //  neighbor list is fresh built, rebuilt, and updates
      //
      int test_neighbor_list = 1;
      if(test_neighbor_list == 1)
      {
        // checking if crossneighborlist gives same result
        d_print("Main : Computing Forces\n");
        rebuild_neighbor_flag = 1;
        rebuild_cluster_flag = 0;

        // test_cluster_force_energy.clear();
        // test_total_cluster_force_energy.clear();
        // test_node_force_energy.clear();
        // test_total_node_force_energy.clear();

        forceC->testForceEnergyCalculation(rebuild_neighbor_flag,
          rebuild_cluster_flag,
          test_cluster_force_energy,
          test_node_force_energy,
          test_total_cluster_force_energy,
          test_total_node_force_energy);
        d_print("Main : Done\n");

        // perform output by calling testEnergyOutput()
        d_print("Main : Performing output\n");
        Output::getInstance()->testEnergyOutput(&test_cluster_force_energy,
          &test_node_force_energy,
          &test_total_cluster_force_energy,
          &test_total_node_force_energy);
        d_print("Main : Done\n");

        // checking if crossneighborlist gives same result
        d_print("Main : Computing Forces\n");
        rebuild_neighbor_flag = 0;
        rebuild_cluster_flag = 0;

        // test_cluster_force_energy.clear();
        // test_total_cluster_force_energy.clear();
        // test_node_force_energy.clear();
        // test_total_node_force_energy.clear();

        forceC->testForceEnergyCalculation(rebuild_neighbor_flag,
          rebuild_cluster_flag,
          test_cluster_force_energy,
          test_node_force_energy,
          test_total_cluster_force_energy,
          test_total_node_force_energy);
        d_print("Main : Done\n");

        // perform output by calling testEnergyOutput()
        d_print("Main : Performing output\n");
        Output::getInstance()->testEnergyOutput(&test_cluster_force_energy,
          &test_node_force_energy,
          &test_total_cluster_force_energy,
          &test_total_node_force_energy);
        d_print("Main : Done\n");
      }

      d_print("Main : Done\n");

      exit(EXIT_SUCCESS);
      return EXIT_SUCCESS;
    }
    break;

    case 1:
    {
      //
      //  test of crossneighborlist data
      //

      // not implemented yet
      d_print("test of crossneighborlist is not implemented yet\n");

      exit(EXIT_SUCCESS);
      return EXIT_SUCCESS;
    }
    break;

    case 2:
    {
      //  call create elements : don't compute forces in
      //  createElementsCentral()
      int computeForce = 1;
      d_print("Main : Meshing all Quasis\n");
      CreateMesh * createMesh = CreateMesh::getInstance();
      createMesh->createElementsCentral(computeForce);
      d_print("Main : Done\n");

      //
      //
      //
      d_print("computing unit cell volume and also volume of unit cell after deformation\n");
      for(int iQuasi =0; iQuasi < NUM_QUASI; iQuasi++)
      {
        const std::vector<double> iShift =
          quasicontinua->getShift(iQuasi);

        quasicontinua->getQuasi(iQuasi).applyDeformation(boundaryFlag,
                  loadDeformation,
                  iShift);
      }

      //
      //  compute deformation gradient and volume of unit cell for quasi 0
      //
      computeDeformGrad();

      d_print("Main : Done\n");

      exit(EXIT_SUCCESS);
      return EXIT_SUCCESS;
    }
    break;

    case 3:
    {
      //
      //  output node data of all quasi's
      //
      d_print("Main : Outputting node data of quasi's\n");

      char *filename;
      filename = (char *) malloc(path_max*sizeof(char));
      for(int iQuasi=0; iQuasi < NUM_QUASI; iQuasi++)
      {
        struct node_list_t *P_node_list =&(quasicontinua->getQuasi(iQuasi).getNodeList().node_list);

        //
        // call function to write the data
        //
        sprintf(filename, "input_data_node_quasi_%i", iQuasi);

        Output::getInstance()->performNodeListOutput(filename,P_node_list);

        d_print("Main : Done\n");

        exit(EXIT_SUCCESS);
        return EXIT_SUCCESS;
      }
    }
    break;

    case 4:
    {
      //
      //  meshes the quasi and outputs the node and mesh data
      //
      //
      //  call create elements : don't compute forces in
      //  createElementsCentral()
      //
      int computeForce = 1;
      d_print("Main : Meshing all Quasis\n");
      CreateMesh * createMesh = CreateMesh::getInstance();
      createMesh->createElementsCentral(computeForce);
      d_print("Main : Done\n");

      //
      //  output
      //
      unsigned int flags;
      flags |= NODE_OUTPUT_FLAG;

      // 
      //  create output_flag_vector
      //
      std::vector<int> output_flag_vector(1,0);

      d_print("Main : Performing output for all quasis\n");

      Output::getInstance()->performOutputCentral(flags, output_flag_vector);
      d_print("Main : Done\n");

      exit(EXIT_SUCCESS);
      return EXIT_SUCCESS;
    }
    break;

    case 5:
    {
      //
      //  outputs the initial forces
      //  Note
      //  1: The forces are after removing the residual forces
      //  2: It removes the atomistic loads
      //  3: The residual forces at initial configuration
      //     depends on remove_residual_flag
      //     If flag asks to remove initial forces in quasi.ini
      //     it removes the forces.
      //
      d_print("Main : Initial force and energy\n");

      // setting d_EAMFlag to false for now
      // PairPotentials::getInstance()->d_EAMFlag = false;

      //
      //  call create elements : compute forces in createElementsCentral()
      //  only if also want output of total forces on nodes.
      //
      //  0 : compute forces in Create Mesh
      //  1 : do not compute forces in Create Mesh
      //
      int computeForce = 0;
      d_print("Main : Meshing all Quasis\n");
      CreateMesh * createMesh = CreateMesh::getInstance();
      createMesh->createElementsCentral(computeForce);
      d_print("Main : Done\n");

      //
      //  output
      //
      unsigned int flags;
      flags |= NODE_OUTPUT_FLAG;

      // 
      //  create output_flag_vector
      //
      std::vector<int> output_flag_vector(1,0);

      d_print("Main : Performing output for all quasis\n");

      Output::getInstance()->performOutputCentral(flags, output_flag_vector);
      d_print("Main : Done\n");

      d_print("Main : Done\n");

      exit(EXIT_SUCCESS);
      return EXIT_SUCCESS;
    }
    break;

    case 6:
    {
      //
      //  Outputs the time required for each component of interaction
      //  Note
      //  1: The forces are after removing the residual forces
      //  2: It removes the atomistic loads
      //  3: The residual forces at initial configuration
      //     depends on remove_residual_flag
      //     If flag asks to remove initial forces in quasi.ini
      //     it removes the forces.
      //
      d_print("Main : Collecting time data for the QC code\n");

      //  set RunData flag to true and initialize its data
      RunData::getInstance()->d_timeOut = true;
      RunData::getInstance()->initializeData();

      // setting d_EAMFlag to false for now
      // PairPotentials::getInstance()->d_EAMFlag = false;

      //
      //  call create elements : compute forces in createElementsCentral()
      //  only if also want output of total forces on nodes.
      //
      //  0 : compute forces in Create Mesh
      //  1 : do not compute forces in Create Mesh
      //
      int computeForce = 0;
      d_print("Main : Meshing all Quasis\n");
      CreateMesh * createMesh = CreateMesh::getInstance();
      createMesh->createElementsCentral(computeForce);
      d_print("Main : Done\n");

      // collect time
      gettimeofday(&t2,NULL);
      elapsedTime = (t2.tv_sec - t1.tv_sec);
      elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;

      if(RunData::getInstance()->d_timeOut == true)
        RunData::getInstance()->addInitialToFirstForceCalculationTime(0,
          elapsedTime);
      //
      //  output
      //
      unsigned int flags;
      flags |= NODE_OUTPUT_FLAG;

      // 
      //  create output_flag_vector
      //
      std::vector<int> output_flag_vector(1,0);


      d_print("Main : Performing output for all quasis\n");

      Output::getInstance()->performOutputCentral(flags, output_flag_vector);
      d_print("Main : Done\n");

      d_print("Main : Done\n");

      exit(EXIT_SUCCESS);
      return EXIT_SUCCESS;
    }
    break;

    default:
    {
      d_print("check test_flag in quasi.ini\n");
      d_print("exiting the code\n");
    }
    break;
  }

  d_print("Main : Done\n");

  exit(EXIT_SUCCESS);
  return EXIT_SUCCESS;
} // end of main




