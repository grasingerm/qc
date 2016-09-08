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
#include "C_Interface.h"
#include "Error.h"

#if !defined(lint)
static char rcsid[]="$Id: qc3d.c,v 1.10 2002/03/07 23:53:08 knap Exp $";
#endif // !lint

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

  // get 0 Quasi -0 Quasi interaction number
  int potentialNumber = PairPotentials::getInstance()->doQuasicontinuumInteract(0,0); 

  // define rebuild_cutoff
  const double rebuild_cutoff =  PairPotentials::getInstance()->getCutoffRadiusNeighborList(potentialNumber);

  //
  // initialize nodes_added
  //
  for(int i = 0; i < NUM_QUASI; i++)
    nodes_added.push_back(0);

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
  //  meshing the quasis
  //
  d_print("Main : Meshing all Quasis\n");
  CreateMesh * createMesh = CreateMesh::getInstance();
  createMesh->createElementsCentral(0);
  d_print("Main : Done\n");

  //
  //  Dump the information of problem to file
  //
  Input::getInstance()->systemInfoOutput();

  //
  //  initialize CGNonLinearSolver
  //
  CGNonLinearSolver CGSolver( cgTolerance,
                              cgMaxIterations,
                              cgDebugLevel,
                              cgLineSearchTolerance,
                              cgLineSearchMaxIterations,
                              rebuild_cutoff,
                              minMethod,
                              minMethodMaxIterations,
                              minMethodTolerance); 

  //
  // initialConfigFlag 
  // 0 - use initial configuration provided by input file.
  // 1 - first equilibrate the configuration provided by input 
  //     file and then set the equilibrated configuration
  //     as initial configuration
  //
  if(initialConfigFlag == 1)
  {
    // create QuasicontinuaForceFunction
    QuasicontinuaForceFunction forceFunction;

    // equilibrate
    d_print("CreateMesh : Equilibrating system to define initial configuration\n");
    
    int output_flag = 0;
    int loadNumber = 0;

    CGSolver.solve(forceFunction,
      cgTolerance,
      cgMaxIterations,
      cgDebugLevel,
      output_flag,
      loadNumber,
      minMethod,
      minMethodMaxIterations,
      minMethodTolerance);
    d_print("CreateMesh : Done\n");

    d_print("CreateMesh : Redefining Initial Configuration for all quasi\n");

    for(int iQuasi = 0; iQuasi < NUM_QUASI; iQuasi++)
      quasicontinua->getQuasi(iQuasi).modifyInitialConfiguration();

    d_print("CreateMesh : Done\n");

    d_print("CreateMesh : performOutput on converged configuration\n");

    unsigned int flags;
    flags |= NODE_OUTPUT_FLAG;
    // flags |= RESTART_FLAG;
    
    std::vector<int> output_flag_vector;
    output_flag_vector.push_back(0);

    Output::getInstance()->performOutputCentral(flags, output_flag_vector);     

    d_print("CreateMesh : Done\n");        

    d_print("CreateMesh : Setting the initial configuration flag in Lattice Class\n");

    Lattice::getInstance()->setInitialConfigFlag(initialConfigFlag);

    d_print("CreateMesh : Done\n");      
  }

  //
  // get Instance of Electrostatics
  //
  Electrostatics * electrostaticsC = Electrostatics::getInstance();

  // initialize usage timer
  // initialize_usage_timer();

  // 
  // loop over loading
  //
  int   iLoad;
  // for( iLoad = 0; iLoad < 0; iLoad++)
  for( iLoad = 0; iLoad <= NUMBER_LOAD; iLoad++)
  {
    //
    //  set output flag for solve()
    //
    int output_flag = 0;
    if(iLoad == 0 || iLoad == 1)
      output_flag = 1;

    d_print("****************** Loading *********************\n");
    d_print("Main : Apply load = %i\n ", iLoad);

    // apply loading to Quasicontinuun, Electrostatics
    if (iLoad != 0)
    {
      // 
      // apply deformation to Quasicontinuum
      //      
      if(boundaryFlag > -1)
      {
        d_print("Main : Apply deformation gradient to all quasis\n");
        for(int iQuasi =0; iQuasi < NUM_QUASI; iQuasi++)
        {
          const std::vector<double> iShift =
            quasicontinua->getShift(iQuasi); 
          
          quasicontinua->getQuasi(iQuasi).applyDeformation(boundaryFlag,
                    loadDeformation,
                    iShift);
        }
      }

      //
      // apply external charge and electric field
      //
      std::vector<std::pair<double, std::vector<double> > > temp_charges;

      for(int iCharge = 0; iCharge < charges.size(); iCharge++)
      {
        temp_charges.push_back(std::make_pair(charges[iCharge].first*iLoad,
          charges[iCharge].second));
      } 

      if(electrostaticsC->isElectrostaticEnable() == 1)
      {
        electrostaticsC->clearExternalPointCharges();
        d_print("Main : Insert external point charges\n");
        electrostaticsC->insertExternalPointCharges(temp_charges);
      }
    } // apply deformation

    d_print("Main : Done\n");

    //
    //  Create QuasicontinuaForceFunction object
    //
    d_print("Main : Creating QuasicontinuaForceFunction Object\n");
    QuasicontinuaForceFunction forceFunction_1;

    //
    //  equilibrate the system
    //
    d_print("Main : Calling CGNonLinearSolver.solve()\n");
    CGSolver.solve(forceFunction_1,
        cgTolerance,
        cgMaxIterations,
        cgDebugLevel,
        output_flag,
        iLoad,
        minMethod,
        minMethodMaxIterations,
        minMethodTolerance);
    d_print("Main : Done\n");

    //
    //  restart
    //
    bool restart = false;
    int restartNum = 0;

    d_print("Main : Entering meshing loop for load step = %i\n", iLoad);
    while(restart == true)
    {
      d_print("Main : Remesh iteration = %i\n", restartNum);
      restartNum++;

      // reset total nodes added from remeshing
      total_nodes_added = 0;

      // 
      //  remesh
      //
      d_print("Main : Remesh all the quasis\n");
      total_nodes_added = createMesh->remeshCentral(iLoad);

      //
      //  create new QuasicontinuaForceFunction object
      //
      d_print("Main : Creating QuasicontinuaForceFunction Object\n");
      QuasicontinuaForceFunction forceFunction_2;

      //
      //  quilibrate system
      //        
      d_print("Main : Calling CGNonLinearSolver.solve()\n");
      CGSolver.solve(forceFunction_2,
            cgTolerance,
            cgMaxIterations,
            cgDebugLevel,        
            output_flag,
            iLoad,
            minMethod,
            minMethodMaxIterations,
            minMethodTolerance);
      d_print("Main : Done\n");

      //
      //  remesh information check
      //
      d_print("Main : Checking remesh criteria\n");
      if(addNodeFlag == 0)
      {
        if(total_nodes_added > addNodeCutoff)
          restart = true;
        else
          restart = false;
      }
      else if(addNodeFlag == 1)
      {
        if(restartNum < addNodeRestartNum)
          restart = true;
        else
          restart = false;
      }
    } // end of while loop
    d_print("Main : Done with remesh iterations\n");

    //
    //  output
    //
    unsigned int flags;
    flags |= NODE_OUTPUT_FLAG;
    flags |= RESTART_FLAG;
    // flags |= FULL_LATTICE_OUTPUT_FLAG;

    // 
    //  create output_flag_vector
    //
    std::vector<int> output_flag_vector(6,0);
    output_flag_vector[0] = minMethod; // pass the minimization method
    output_flag_vector[5] = 3; // pass the cg method status

    d_print("Main : Performing output for all quasis\n");

    Output::getInstance()->performOutputCentral(flags, output_flag_vector);
    d_print("Main : Done\n");      
  } // end of loop over load

  // stop timer
  gettimeofday(&t2,NULL);

  // compute and print elapsed time in seconds
  elapsedTime = (t2.tv_sec - t1.tv_sec);
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
  d_print("Main : Program executed in = %f s. \n", elapsedTime);
  //std::cout<<"Program executed in: "<<elapsedTime<<" s."<<std::endl;  

  exit(EXIT_SUCCESS);
  return EXIT_SUCCESS;
} // end of main




