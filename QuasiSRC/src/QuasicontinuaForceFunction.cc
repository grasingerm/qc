//
// File:      QuadraticFunction.cc
//
// Solvers tests.
//
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
/* No _REENTRANT */
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#include <cassert>
#include <cstdlib>
#include <iostream>

// c++ files
#include "CrossNeighborList.h"
#include "DataTypes.h"
#include "Electrostatics.h"
#include "Error.h"
#include "ForceEnergyCalculation.h"
#include "PairPotentials.h"
#include "QuadraturePoints.h"
#include "Quasicontinua.h"
#include "QuasicontinuaForceFunction.h"
#include "C_Interface.h"
#include "Error.h"

// c files
#include "threads.h"
#include "monitor.h"


//
//
//
namespace quasicontinuum {

	//
	// Constructor.
	//
	QuasicontinuaForceFunction::QuasicontinuaForceFunction(): 
     d_lock(PTHREAD_MUTEX_INITIALIZER),
     d_version(quasicontinuum::SINGLE_THREADED) // Note: Multithreaded not implemented
	{
	  #if 0
  		//
  		// determine whether multithreaded
  		//
  		if(quasicontinuum::get_max_number_threads() == 1)
  	  		d_version = quasicontinuum::SINGLE_THREADED;
	  #endif
  	//
  	// switch between thread versions
 		//
  	switch(d_version)
  	{
  		case quasicontinuum::MULTI_THREADED:
 		 	{
  			//
  			// Not implemented
  			//
  			d_print("MULTI_THREADED not implemented\n");
 		 		exit(EXIT_FAILURE);
  		}

  	  break;

  		case quasicontinuum::SINGLE_THREADED:
  	  {

      		//
      		// get quasicontinua solution
      		//
      		GetQuasicontinuaSolution();
		  }

    	break;
  	}

  	//
  	//
  	//
  	return;
	}
  
  //
  // Destructor.
  //
  QuasicontinuaForceFunction::~QuasicontinuaForceFunction()
  {

  	//
   	//
   	//
   	return;
  }

  //
  // Obtain number of unknowns.
  //
  int
	QuasicontinuaForceFunction::getNumberUnknowns() const
  {
  	//
   	//
   	//
   	return d_solution.size();
  }

	//
	// Compute function value (aka energy).
	//
	double
	QuasicontinuaForceFunction::value() const
	{
    //
  	// total value or energy
  	//
		double totalEnergy = 0.0;

		Quasicontinua * quasicontinua = Quasicontinua::getInstance();

		int numQuasi = quasicontinua->size();

		// loop over all quasis
		for(int iQuasi = 0; iQuasi < numQuasi; iQuasi++)
		{
			// get all_node_list of iQuasi
			struct all_node_list_t * P_all_node_list = &(quasicontinua->getQuasi(iQuasi).getNodeList());

			totalEnergy += P_all_node_list->energy.total;
		}

  	return totalEnergy;
	}

	//
	// Compute function gradient (aka forces).
	//
	void
	QuasicontinuaForceFunction::gradient(std::vector<double> & gradient,
			       const int             rebuild_neighbor_flag) const
	{
    // Quasicontinua 
    Quasicontinua * quasicontinua = Quasicontinua::getInstance();

    // get size of Quasicontinua
    int numQuasi = quasicontinua->size();

		// get force from ForceEnergyCalculation.cc
    bool compute_in_reference = false;
    d_print("Force...");
		ForceEnergyCalculation::getInstance()->computeForceEnergy(rebuild_neighbor_flag, compute_in_reference);
		d_print("complete\n");

    // loop over all quasicontinuum and update
    for(int iQuasi = 0; iQuasi < numQuasi; ++iQuasi)
    {
      // vector to hold Forces values
      const std::vector<double> iForces = 
        quasicontinua->getQuasi(iQuasi).getForces();

      // start point of each Quasicontinuum Forces
      int startForces;
      if (iQuasi == 0)
        startForces = 0;
      else
        startForces += d_numQuasicontinuumUnknowns[iQuasi - 1];

      // populate gradient
      for(unsigned int i = 0; i < iForces.size(); ++i)
      {
        gradient[startForces + i] = iForces[i];
      }
    }

		return;
	}

	//
	// Apply preconditioner to function gradient.
	//
	void
	QuasicontinuaForceFunction::precondition(std::vector<double>       & s,
				   const std::vector<double> & gradient) const
	{
		//
		// switch between versions
		//

		switch(d_version)
		{
			case quasicontinuum::MULTI_THREADED:
			{
				//
				// Not implemented
				//

				d_print("MULTI_THREADED not implemented\n");
				exit(EXIT_FAILURE);
			}

			break;

			case quasicontinuum::SINGLE_THREADED:
			{
				//
				// check sizes
				//

				assert(s.size() == gradient.size());

				//
				// get quasicontinua
				//
				Quasicontinua * quasicontinua = Quasicontinua::getInstance();

				//
				// get number of Quasicontinuums
				//
				const int numQuasi = quasicontinua->size();

				//
				// loop over all quasicontinuum and update
				//
				for(int iQuasi = 0; iQuasi < numQuasi; ++iQuasi)
				{
					//
					// vector to hold Precondition values
					//
					const std::vector<double> iPrecondition = 
					    quasicontinua->getQuasi(iQuasi).getPreconditionValues();

					//
					// start point of each Quasicontinuum Precondition
					//

					int startPrecondition;
					if (iQuasi == 0)
						startPrecondition = 0;
					else
						startPrecondition += 
							d_numQuasicontinuumUnknowns[iQuasi - 1];

					//
					// populate preconditionedvalue
					//
					for(unsigned int i = 0; i < iPrecondition.size(); ++i)
					{
						//
						// check precondition value
						//
						assert(iPrecondition[i] != 0.0);

						//
						// update s
						//
						s[startPrecondition + i] = 
							-gradient[startPrecondition + i] 
							/ iPrecondition[i];
					}
				}
			}

			break;

		}

		//
		//
		//
		return;
	}

	//
	// Update solution.
	//
	void
	QuasicontinuaForceFunction::update(const std::vector<double> & solution)
	{
		//
		// switch between versions
		//
		switch(d_version)
		{
			case quasicontinuum::MULTI_THREADED:
			{
				//
				// Not implemented
				//
				d_print("MULTI_THREADED not implemented\n");
				exit(EXIT_FAILURE);
			}

			break;

			case quasicontinuum::SINGLE_THREADED:
			{
				//
				// check sizes
				// 
				assert(d_solution.size() == solution.size());

				//
				// clear solution
				//
				d_solution.clear();

				//
				// store solution
				//
				d_solution.insert(d_solution.end(), solution.begin(),
					  solution.end());

				//
				// get quasicontinua
				//
				Quasicontinua * quasicontinua = Quasicontinua::getInstance();

				//
				// get number of Quasicontinuums
				//

				const int numQuasi = quasicontinua->size();

				//
				// loop over all quasicontinuum and update
				//
				for(int iQuasi = 0; iQuasi < numQuasi; ++iQuasi)
				{
					//
					// vector to hold update
					//
					std::vector<double> iUpdate;

					//
					// start point of each Quasicontinuum solution
					//

					int startUpdate;
					if (iQuasi == 0)
						startUpdate = 0;
					else
						startUpdate += d_numQuasicontinuumUnknowns[iQuasi - 1];

					//std::cout<<startUpdate<<std::endl;
					//
					// get update for individual Quasicontinuum
					//

					for(int iCount = 0 ; 
						iCount < d_numQuasicontinuumUnknowns[iQuasi]; ++iCount)
						iUpdate.push_back(solution[startUpdate + iCount]);

					//
					// update Quasicontinuum
					//
					quasicontinua->getQuasi(iQuasi).updateState(iUpdate);
				}
			}

			break;

		}

		//
		//
		//
		return;
	}

	//
	// Obtain current solution.
	//
	void
	QuasicontinuaForceFunction::solution(std::vector<double> & solution) const
	{
    // update the data d_solution
    GetQuasicontinuaSolution();

		//
		// clear solution
		//
		solution.clear();

		//
		// insert d_solution in solution
		//
		solution.insert(solution.end(), d_solution.begin(),
			d_solution.end());
	
	  //
   	//
   	//
    return;
	}

	//
	// For each unknown indicate whether that unknown should be
	// accumulated. This functionality is needed in the case of parallel
	// execution when domain decomposition is employed. Unknowns
	// residing on processor boundary should only be accumulated once
	// when dot products of vertex fields are computed (e.g. residual).
	//
	std::vector<int>
	QuasicontinuaForceFunction::getUnknownCountFlag() const
	{	
		//
	  //
	  //
	  return std::vector<int>(d_solution.size(), 1); // FIXME: not implemented for MPI
	}

	//
	// get solution size
	//
	int
	QuasicontinuaForceFunction::getSolutionSize() const
	{
 		//
 		//
  	//
  	return d_solution.size(); 
  }

  //
  // get Quasicontinua solutions
  //
  void
  QuasicontinuaForceFunction::GetQuasicontinuaSolution() const
  {
		//
		// clear solution
		//
		d_solution.clear();

		//
		// get quasicontinua
		//
		Quasicontinua * quasicontinua = Quasicontinua::getInstance();

		//
		// get number of Quasicontinuums
		//
		const unsigned int numQuasi = quasicontinua->size();

		//
		// check size of d_numQuasicontinuumUnknowns
		//
		if (d_numQuasicontinuumUnknowns.size() != numQuasi)
	  		d_numQuasicontinuumUnknowns.resize(numQuasi);

		//
		// loop over all quasicontinuum and append solution
		//
		for(unsigned int iQuasi = 0; iQuasi < numQuasi; ++iQuasi)
		{
			//
			// get Quasicontinuum solution
			//
			const std::vector<double> iSolution = 
			    quasicontinua->getQuasi(iQuasi).getState();

			//
			// append Quasicontinuum solution to d_solution
			//
			d_solution.insert(d_solution.end(), iSolution.begin(),
				   iSolution.end());

			//
			// set number of unknowns for each Quasicontinuum
			//
			d_numQuasicontinuumUnknowns[iQuasi] = iSolution.size();
		}

  	//
  	//
  	//
  	return; 
  }

  //
  //	setPositionFixity()
  //
  void QuasicontinuaForceFunction::setPositionFixity(const int  fixity_flag) const
  {
		//
		// get quasicontinua
		//
		Quasicontinua * quasicontinua = Quasicontinua::getInstance();

		//
		// get number of Quasicontinuums
		//
		const unsigned int numQuasi = quasicontinua->size();

		//
		//	loop over Quasicontinuum and call the function which sets the
		// 	fixity
		//
		for(unsigned int iQuasi = 0; iQuasi < numQuasi; ++iQuasi)
		{
			quasicontinua->getQuasi(iQuasi).setPositionFixity(fixity_flag);
		}

		return;
  }

  //
  //	setFrequencyFixity()
  //
  void QuasicontinuaForceFunction::setFrequencyFixity(const int  fixity_flag) const
  {
		//
		// get quasicontinua
		//
		Quasicontinua * quasicontinua = Quasicontinua::getInstance();

		//
		// get number of Quasicontinuums
		//
		const unsigned int numQuasi = quasicontinua->size();

		//
		//	loop over Quasicontinuum and call the function which sets the
		// 	fixity
		//
		for(unsigned int iQuasi = 0; iQuasi < numQuasi; ++iQuasi)
		{
			quasicontinua->getQuasi(iQuasi).setFrequencyFixity(fixity_flag);
		}

		return;
  }


  //
  // recomputes the data of SolverFunction instance
  //
  void QuasicontinuaForceFunction::recomputeData() const
  {
  	//
  	//	recompute all the data of this instance. 
  	//
  	//	Note: This is done because we may change the fixity mask of 
  	//	freq or position dof for minimization purpose. This
  	//  changes the number of unknown in the system.
  	//
		GetQuasicontinuaSolution();

		return;
  }

}
