//
// File:      CGNonLinearSolver.cc
// Package:   solvers
//
// Various generic solvers.
//
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H


#if 0 //defined(HAVE_MPI)
#include <LengthError.h>
#include <Accumulator.h>
#undefine HAVE_MPI
#endif // HAVE_MPI

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#include <cmath>
#include <iostream>

// C++ source
#include "CGNonLinearSolver.h"
#include "SolverFunction.h"
#include "DataTypes.h"
#include "Output.h"
#include "C_Interface.h"

// C source
#include "monitor.h"
#include "threads.h"

// max cg iterations

static int alternateIter = 0;
static int minimizationMethod = 0;
static int freq_or_position_min = 0;

//
//
//
namespace quasicontinuum {

    //
    //
    //
    namespace{

      //
      // local types, data and functions
      //
      struct initializeDirection_argument_t {
       double                d_deltaNew;
       std::vector<double> * d_sOld;
       std::vector<double> * d_direction;
       std::vector<int>    * d_unknownCountFlag;
       std::vector<double> * d_s;
       std::vector<double> * d_gradient;
       size_t                size;
       pthread_mutex_t     * lock;
      };

     struct computeDeltaD_argument_t {
       double                deltaD;
       double                etaP;
       std::vector<double> * d_direction;
       std::vector<int>    * d_unknownCountFlag;
       std::vector<double> * d_gradient;
       size_t                size;
       pthread_mutex_t     * lock;
      };

     struct computeEta_argument_t {
       double                eta;
       std::vector<double> * d_direction;
       std::vector<int>    * d_unknownCountFlag;
       std::vector<double> * d_gradient;
       size_t                size;
       pthread_mutex_t     * lock;
      };

     struct computeDeltas_argument_t {
       double                d_deltaNew;
       double                d_deltaMid;
       std::vector<double> * d_s;
       std::vector<int>    * d_unknownCountFlag;
       std::vector<double> * d_sOld;
       std::vector<double> * d_gradient;
       size_t                size;
       pthread_mutex_t     * lock;
      };

     struct updateDirection_argument_t {
       double                d_beta;
       std::vector<double> * d_s;
       std::vector<double> * d_direction;
       size_t                size;
      };

     struct computeResidualL2Norm_argument_t {
       double                norm;
       std::vector<int>    * d_unknownCountFlag;
       std::vector<double> * d_gradient;
       size_t                size;
       pthread_mutex_t     * lock;
      };

     struct updateSolution_argument_t {
       double                      alpha;
       const std::vector<int>    * unknownCountFlag;
       const std::vector<double> * direction;
       std::vector<double>       * d_solution;
       size_t                      size;
      };

      //
      // distribute initialize direction work among threads 
      //
      void * 
      threadInitializeDirectionWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct initializeDirection_argument_t * initializeDirection_argument = 
          static_cast<struct initializeDirection_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                initializeDirection_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

         //
         // get handles
         //
         std::vector<double> & d_sOld = 
           *(initializeDirection_argument->d_sOld);
         std::vector<double> & d_direction = 
           *(initializeDirection_argument->d_direction);
         std::vector<double> & d_gradient = 
           *(initializeDirection_argument->d_gradient);
         std::vector<double> & d_s = 
           *(initializeDirection_argument->d_s);
          #if 0 //defined(HAVE_MPI)
           std::vector<int> & d_unknownCountFlag = 
             *(initializeDirection_argument->d_unknownCountFlag);
          #endif
          //
          // local variable
          //
          double localDeltaNew = 0.0;

          //
         // loop over all elements of solution
          //
          for (int i = number_start; i <= number_end; ++i){


          #if 0 //defined(HAVE_MPI)
              const double factor = d_unknownCountFlag[i];
          #else
          const double factor = 1.0;
          #endif // HAVE_MPI
    
          const double r = -d_gradient[i];
          const double s = d_s[i];

          d_sOld[i]        = s;
          d_direction[i]   = s;
          localDeltaNew   += factor * r * d_direction[i];
          }
    
        //
        // accumulate the value
        //
        pthread_mutex_lock(initializeDirection_argument->lock);
      
        initializeDirection_argument->d_deltaNew += localDeltaNew;

        pthread_mutex_unlock(initializeDirection_argument->lock);

        //
        //
        //
        return NULL;

      }

      //
      // distribute computeDeltaD work among threads 
      //
      void * 
      threadComputeDeltaDWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct computeDeltaD_argument_t * computeDeltaD_argument = 
          static_cast<struct computeDeltaD_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                computeDeltaD_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

       //
       // get handle
       //
       std::vector<double> & d_gradient = 
         *(computeDeltaD_argument->d_gradient);
       std::vector<double> & d_direction = 
         *(computeDeltaD_argument->d_direction);
        #if 0 //defined(HAVE_MPI)
         std::vector<int> & d_unknownCountFlag = 
         *(computeDeltaD_argument->d_unknownCountFlag);
        #endif
        //
        // local variable
        //
        double localdeltaD = 0.0;
        double localetaP = 0.0;

        //
        // loop over all elements of solution
        //
        for (int i = number_start; i <= number_end; ++i){
   
        #if 0 //defined(HAVE_MPI)
         const double factor = d_unknownCountFlag[i];
        #else
         const double factor = 1.0;
        #endif // HAVE_MPI
      
         const double direction = d_direction[i];
    
         localdeltaD += factor * direction * direction;
         localetaP   += factor * d_gradient[i] * direction;
        }
    
        //
        // accumulate the values
        //
        pthread_mutex_lock(computeDeltaD_argument->lock);
      
        computeDeltaD_argument->deltaD += localdeltaD;
        computeDeltaD_argument->etaP += localetaP;

        pthread_mutex_unlock(computeDeltaD_argument->lock);

        //
        //
        //
        return NULL;

      }

      //
      // distribute computeEta work among threads 
      //
      void * 
      threadComputeEtaWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct computeEta_argument_t * computeEta_argument = 
          static_cast<struct computeEta_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                computeEta_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

         //
         // get handle
         //
         std::vector<double> & d_gradient = 
          *(computeEta_argument->d_gradient);
         std::vector<double> & d_direction = 
          *(computeEta_argument->d_direction);
         #if 0 //defined(HAVE_MPI)
         std::vector<int> & d_unknownCountFlag = 
         *(computeEta_argument->d_unknownCountFlag);
         #endif
        //
        // local variable
        //
        double localEta = 0.0;

        //
        // loop over all elements of solution
        //
        for (int i = number_start; i <= number_end; ++i){

         #if 0 //defined(HAVE_MPI)
            const double factor = d_unknownCountFlag[i];
         #else
         const double factor = 1.0;
         #endif // HAVE_MPI
  
         const double direction = d_direction[i];

         localEta += factor * d_gradient[i] * direction;
         }   
    
        //
        // accumulate the values
        //
        pthread_mutex_lock(computeEta_argument->lock);
      
        computeEta_argument->eta += localEta;

        pthread_mutex_unlock(computeEta_argument->lock);

        //
        //
        //
        return NULL;

      }

      //
      // distribute computeDeltas work among threads 
      //
      void * 
      threadComputeDeltasWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct computeDeltas_argument_t * computeDeltas_argument = 
          static_cast<struct computeDeltas_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                computeDeltas_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

         //
         // get handle
         //
         std::vector<double> & d_sOld = 
          *(computeDeltas_argument->d_sOld);
         std::vector<double> & d_gradient = 
          *(computeDeltas_argument->d_gradient);
         std::vector<double> & d_s = 
          *(computeDeltas_argument->d_s);
         #if 0 //defined(HAVE_MPI)
         std::vector<int> & d_unknownCountFlag = 
         *(computeDeltas_argument->d_unknownCountFlag);
         #endif
        //
        // local variable
        //
        double localDeltaMid = 0.0;
        double localDeltaNew = 0.0;

        //
        // loop over all elements of solution
        //
        for (int i = number_start; i <= number_end; ++i){

         #if 0 //defined(HAVE_MPI)
            const double factor = d_unknownCountFlag[i];
         #else
            const double factor = 1.0;
         #endif // HAVE_MPI

         const double r    = -d_gradient[i];
         const double s    = d_s[i];
         const double sOld = d_sOld[i];

         //
         // compute delta mid
         //
         localDeltaMid += factor * r * sOld;

         //
         // save gradient old
         //
         d_sOld[i] = s;

         //
         // compute delta new
          //
         localDeltaNew += factor * r * s;

         }

    
        //
        // accumulate the values
        //
        pthread_mutex_lock(computeDeltas_argument->lock);
      
        computeDeltas_argument->d_deltaMid += localDeltaMid;
        computeDeltas_argument->d_deltaNew += localDeltaNew;

        pthread_mutex_unlock(computeDeltas_argument->lock);

        //
        //
        //
        return NULL;

      }

      //
      // distribute updateDirection work among threads 
      //
      void * 
      threadUpdateDirectionWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct updateDirection_argument_t * updateDirection_argument = 
          static_cast<struct updateDirection_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                updateDirection_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

         //
         // get handle
         //
         std::vector<double> & d_direction = 
          *(updateDirection_argument->d_direction);
          std::vector<double> & d_s = 
          *(updateDirection_argument->d_s);

        //
        // loop over all elements of direction
        //
        for (int i = number_start; i <= number_end; ++i){

         d_direction[i] *= updateDirection_argument->d_beta;
        d_direction[i] += d_s[i];

        }

        //
        //
        //
        return NULL;

      }

     //
      // distribute computeResidualL2Norm work among threads 
      //
      void * 
      threadComputeResidualL2NormWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct computeResidualL2Norm_argument_t * computeResidualL2Norm_argument = 
          static_cast<struct computeResidualL2Norm_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                computeResidualL2Norm_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

        //
        // get handles
        //
        std::vector<double> & d_gradient  = 
        *(computeResidualL2Norm_argument->d_gradient);
        #if 0 //defined(HAVE_MPI)
          std::vector<int> & d_unknownCountFlag = 
          *(computeResidualL2Norm_argument->d_unknownCountFlag);
        #endif
        //
        // local variable
        //
        double localNorm = 0.0;

        //
        // loop over all elements of gradient
        //
        for (int i = number_start; i <= number_end; ++i){
   
        #if 0 //defined(HAVE_MPI)
            const double factor = d_unknownCountFlag[i];
        #else
            const double factor = 1.0;
        #endif // HAVE_MPI

            const double gradient = d_gradient[i];

            localNorm += factor * gradient * gradient;

        }

    
        //
        // accumulate the values
        //
        pthread_mutex_lock(computeResidualL2Norm_argument->lock);
      
        computeResidualL2Norm_argument->norm += localNorm;

        pthread_mutex_unlock(computeResidualL2Norm_argument->lock);

        //
        //
        //
        return NULL;

      }

      //
      // distribute updateSolution work among threads 
      //
      void * 
      threadUpdateSolutionWorker(void * arg)
      {

        //
        // cast back to get initializeDirection_argument_t
        //
        struct updateSolution_argument_t * updateSolution_argument = 
          static_cast<struct updateSolution_argument_t *>(arg);
      
        //
        // get share of data
        //
        const int number_threads = quasicontinuum::get_number_threads();
        const int my_id = quasicontinuum::get_my_tid();
        int number_data_thread;
        int number_start;
        int number_end;
  
        //
        // get share of work
        //
        quasicontinuum::get_share(my_id,
                number_threads,
                updateSolution_argument->size,
                &number_data_thread,
                &number_start,
                &number_end);

        //
        // get handle
        //
        std::vector<double> & d_solution = 
        *(updateSolution_argument->d_solution);
        const std::vector<double> & direction  = 
        *(updateSolution_argument->direction);
        const std::vector<int> & unknownCountFlag = 
        *(updateSolution_argument->unknownCountFlag);

        //
        // loop over all elements of solution
        //
        for (int i = number_start; i <= number_end; ++i){
    
        //
        // update solution
        //
        d_solution[i] = 
        (d_solution[i] + updateSolution_argument->alpha * direction[i]) * unknownCountFlag[i];

        }
        //
        //
        //
        return NULL;

      }

    }

    //
    // Constructor.
    //
    CGNonLinearSolver::CGNonLinearSolver(double tolerance,
             int    maxNumberIterations,
             int    debugLevel,
             double lineSearchTolerance,
             int    lineSearchMaxIterations,
             double rebuild_cutoff,
             int    minMethod,
             int    minMethodMaxIterations,
             double minMethodTolerance) :
      NonLinearSolver(tolerance,
                      maxNumberIterations,
                      debugLevel),
      d_lineSearchTolerance(lineSearchTolerance),
      d_lineSearchMaxIterations(lineSearchMaxIterations),
      d_version(quasicontinuum::MULTI_THREADED),
      REBUILD(1),
      NO_REBUILD(0),
      d_accumulated_r(0.0),
      d_rebuild_cutoff(rebuild_cutoff),
      d_minMethod(minMethod),
      d_minMethodMaxIterations(minMethodMaxIterations),
      d_minMethodTolerance(minMethodTolerance)
    {

      #ifdef HAVE_PTHREAD_H
      //
      // initialize lock
      //
      d_lock = PTHREAD_MUTEX_INITIALIZER;
      #endif /* HAVE_PTHREAD_H */

      //
      // check if threaded or not
      //
      if(quasicontinuum::get_max_number_threads() == 1)
      d_version = quasicontinuum::SINGLE_THREADED;

      //
      //
      //
      return;

    }

    //
    // Destructor.
    //
    CGNonLinearSolver::~CGNonLinearSolver()
    {

      //
      //
      //
      return;

    }

    //
    // initialize direction
    //
    void
    CGNonLinearSolver::initializeDirection()
    {

      //
      // firewalls
      //
      #if 0 //defined(HAVE_MPI)
      if(d_unknownCountFlag.empty() == true) {

        const std::string message("For MPI the dot product fuctor must be "
                                  "initialized.");
        throw LengthError(message);

      }
      #endif // HAVE_MPI

      //
      // initialize d_deltaNew
      //
      d_deltaNew = 0.0;

      //
      // switch between thread versions
      //
      switch(d_version){

      case quasicontinuum::MULTI_THREADED:
      {
        //
        // instantiate initializeDirection_argument_t
        //
        initializeDirection_argument_t initializeDirection_argument;

        //
        // initialize initializeDirection_argument
        //
        initializeDirection_argument.d_deltaNew = 0.0;
        initializeDirection_argument.d_sOld = &(d_sOld);
        initializeDirection_argument.d_direction = &(d_direction);
        initializeDirection_argument.d_unknownCountFlag = &(
          d_unknownCountFlag);
        initializeDirection_argument.d_s = &(d_s);
        initializeDirection_argument.d_gradient = &(d_gradient);
        initializeDirection_argument.size = d_numberUnknowns;
        initializeDirection_argument.lock = &(d_lock);

        //
        // distribute work among workers
        //
        quasicontinuum::thread_monitor(threadInitializeDirectionWorker,
             static_cast<void *>(&initializeDirection_argument),
             quasicontinuum::get_max_number_threads());

        //
        // assign value to d_deltaNew
        //
        d_deltaNew = initializeDirection_argument.d_deltaNew;
      }

      break;

      case quasicontinuum::SINGLE_THREADED:
      {
        //
        // iterate over unknowns
        //
        for (int i = 0; i < d_numberUnknowns; ++i) 
        {
          #if 0 //defined(HAVE_MPI)
          const double factor = d_unknownCountFlag[i];
          #else
          const double factor = 1.0;
          #endif // HAVE_MPI

          const double r = -d_gradient[i];
          const double s = d_s[i];

          d_sOld[i]        = s;
          d_direction[i]   = s;
          d_deltaNew      += factor*r*d_direction[i];

        }

        #if 0 //defined(HAVE_MPI)
        //
        // accumulate d_deltaNew
        //
        d_deltaNew = Accumulator().accumulate(d_deltaNew);
        #endif // HAVE_MPI
      }

      break;
      }

      //
      //
      //
      return;
    }

    //
    // Compute delta_d and eta_p.
    //
    std::pair<double, double>
    CGNonLinearSolver::computeDeltaD()
    {

      //
      // firewalls
      //
      #if 0 //defined(HAVE_MPI)
        if(d_unknownCountFlag.empty() == true) {

        const std::string message("For MPI the dot product fuctor must be "
                                  "initialized.");
        throw LengthError(message);

      }
      #endif // HAVE_MPI

      //
      // initialize delta_d and eta_p
      //
      double deltaD = 0.0;
      double etaP   = 0.0;

      //
      // switch between thread versions
      //
      switch(d_version)
      {
        case quasicontinuum::MULTI_THREADED:
        {
          //
          // instantiate computeDeltaD_argument_t
          //
          computeDeltaD_argument_t computeDeltaD_argument;

          //
          // initialize computeDeltaD_argument
          //
          computeDeltaD_argument.deltaD = 0.0;
          computeDeltaD_argument.etaP = 0.0;
          computeDeltaD_argument.d_direction = &(d_direction);
          computeDeltaD_argument.d_unknownCountFlag = &(d_unknownCountFlag);
          computeDeltaD_argument.d_gradient = &(d_gradient);
          computeDeltaD_argument.size = d_numberUnknowns;
          computeDeltaD_argument.lock = &(d_lock);

          //
          // distribute work among workers
          //
          quasicontinuum::thread_monitor(threadComputeDeltaDWorker,
             static_cast<void *>(&computeDeltaD_argument),
             quasicontinuum::get_max_number_threads());

          //
          // assign value to deltaD
          //
          deltaD = computeDeltaD_argument.deltaD;

          //
          // assign value to etaP
          //
          etaP = computeDeltaD_argument.etaP;
        }
        break;

        case quasicontinuum::SINGLE_THREADED:
        {
          //
          // iterate over unknowns
          //
          for (int i = 0; i < d_numberUnknowns; ++i) 
          {
            #if 0 //defined(HAVE_MPI)
              const double factor = d_unknownCountFlag[i];
            #else
              const double factor = 1.0;
            #endif // HAVE_MPI

            const double direction = d_direction[i];

            deltaD += factor*direction*direction;
            etaP   += factor*d_gradient[i]*direction;
          }
        }

        break;
      }

      #if 0 //defined(HAVE_MPI)
      //
      // instantiate Accumulator
      //
      Accumulator accumulator;

      //
      // accumulate deltaD
      //
      deltaD = accumulator.accumulate(deltaD);

      //
      // accumulate etaP
      //
      etaP = accumulator.accumulate(etaP);

      #endif // HAVE_MPI

      //
      //
      //
      return std::make_pair(deltaD, etaP);

    }

    //
    // Compute eta.
    //
    double
    CGNonLinearSolver::computeEta()
    {

      //
      // firewalls
      //
      #if 0 //defined(HAVE_MPI)
      if(d_unknownCountFlag.empty() == true) {

        const std::string message("For MPI the dot product fuctor must be "
                                  "initialized.");
        throw LengthError(message);

      }
      #endif // HAVE_MPI
      
      //
      // initialize eta
      //
      double eta = 0.0;

      //
      // switch between thread versions
      //
      switch(d_version)
      {
        case quasicontinuum::MULTI_THREADED:
        {
          //
          // instantiate computeEta_argument_t
          //
          computeEta_argument_t computeEta_argument;

          //
          // initialize computeEta_argument
          //
          computeEta_argument.eta = 0.0;
          computeEta_argument.d_direction = &(d_direction);
          computeEta_argument.d_unknownCountFlag = &(d_unknownCountFlag);
          computeEta_argument.d_gradient = &(d_gradient);
          computeEta_argument.size = d_numberUnknowns;
          computeEta_argument.lock = &(d_lock);

          //
          // distribute work among workers
          //
          quasicontinuum::thread_monitor(threadComputeEtaWorker,
             static_cast<void *>(&computeEta_argument),
             quasicontinuum::get_max_number_threads());

          //
          // assign value to etaP
          //

          eta = computeEta_argument.eta;
        }

        break;

        case quasicontinuum::SINGLE_THREADED:
        {
          //
          // iterate over unknowns
          //
          for (int i = 0; i < d_numberUnknowns; ++i) 
          {
            #if 0 //defined(HAVE_MPI)
              const double factor = d_unknownCountFlag[i];
            #else
              const double factor = 1.0;
            #endif // HAVE_MPI

            const double direction = d_direction[i];
            eta += factor*d_gradient[i]*direction;
          }
        }

        break;
      }

      #if 0 //defined(HAVE_MPI)
      //
      // accumulate eta
      //
      eta = Accumulator().accumulate(eta);

      #endif // HAVE_MPI

      //
      //
      //
      return eta;

    }

    //
    // Compute delta new and delta mid.
    //
    void
    CGNonLinearSolver::computeDeltas()
    {

      //
      // firewalls
      //
      #if 0 //defined(HAVE_MPI)
      if(d_unknownCountFlag.empty() == true) {

        const std::string message("For MPI the dot product fuctor must be "
                                  "initialized.");
        throw LengthError(message);

      }
      #endif // HAVE_MPI

      //
      // initialize delta new and delta mid.
      //
      d_deltaMid = 0.0;
      d_deltaNew = 0.0;

      //
      // switch between thread versions
      //
      switch(d_version)
      {
        case quasicontinuum::MULTI_THREADED:
        {
          //
          // instantiate computeDeltas_argument_t
          //
          computeDeltas_argument_t computeDeltas_argument;

          //
          // initialize computeDeltas_argument
          //
          computeDeltas_argument.d_deltaNew = 0.0;
          computeDeltas_argument.d_deltaMid = 0.0;
          computeDeltas_argument.d_s = &(d_s);
          computeDeltas_argument.d_unknownCountFlag = &(d_unknownCountFlag);
          computeDeltas_argument.d_sOld = &(d_sOld);
          computeDeltas_argument.d_gradient = &(d_gradient);
          computeDeltas_argument.size = d_numberUnknowns;
          computeDeltas_argument.lock = &(d_lock);

          //
          // distribute work among workers
          //
          quasicontinuum::thread_monitor(threadComputeDeltasWorker,
             static_cast<void *>(&computeDeltas_argument),
             quasicontinuum::get_max_number_threads());

          //
          // assign values
          //

          d_deltaMid = computeDeltas_argument.d_deltaMid;
          d_deltaNew = computeDeltas_argument.d_deltaNew;
        }

        break;

        case quasicontinuum::SINGLE_THREADED:
        {
          //
          // iterate over unknowns
          //
          for (int i = 0; i < d_numberUnknowns; ++i) 
          {
            #if 0 //defined(HAVE_MPI)
            const double factor = d_unknownCountFlag[i];
            #else
            const double factor = 1.0;
            #endif // HAVE_MPI

            const double r    = -d_gradient[i];
            const double s    = d_s[i];
            const double sOld = d_sOld[i];

            //
            // compute delta mid
            //
            d_deltaMid += factor*r*sOld;

            //
            // save gradient old
            //
            d_sOld[i] = s;

            //
            // compute delta new
            //
            d_deltaNew += factor*r*s;

          }
        }
        break;
      }

      #if 0 //defined(HAVE_MPI)
      //
      // instantiate Accumulator
      //
      Accumulator accumulator;

      //
      // accumulate d_deltaMid
      //
      d_deltaMid = accumulator.accumulate(d_deltaMid);

      //
      // accumulate d_deltaNew
      //
      d_deltaNew = accumulator.accumulate(d_deltaNew);

      #endif // HAVE_MPI

      //
      //
      //
      return;

    }

    //
    // Update direction.
    //
    void
    CGNonLinearSolver::updateDirection()
    {

      //
      // switch between thread versions
      //
      switch(d_version){

      case quasicontinuum::MULTI_THREADED:
      {
        //
        // instantiate updateDirection_argument_t
        //
        updateDirection_argument_t updateDirection_argument;

        //
        // initialize updateDirection_argument
        //
        updateDirection_argument.d_beta = d_beta;
        updateDirection_argument.d_s = &(d_s);
        updateDirection_argument.d_direction = &(d_direction);
        updateDirection_argument.size = d_numberUnknowns;

        //
        // distribute work among workers
        //
        quasicontinuum::thread_monitor(threadUpdateDirectionWorker,
             static_cast<void *>(&updateDirection_argument),
             quasicontinuum::get_max_number_threads());
      }

      break;

      case quasicontinuum::SINGLE_THREADED:
      {
        //
        // iterate over unknowns
        //
        for (int i = 0; i < d_numberUnknowns; ++i) 
        {
          d_direction[i] *= d_beta;
          d_direction[i] += d_s[i];
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
    // Compute residual L2-norm.
    //
    double
    CGNonLinearSolver::computeResidualL2Norm() const
    {

      //
      // firewalls
      //
      #if 0 //defined(HAVE_MPI)
      if(d_unknownCountFlag.empty() == true) {

        const std::string message("For MPI the dot product fuctor must be "
                                  "initialized.");
        throw LengthError(message);

      }
      #endif // HAVE_MPI

      //
      // initialize norm
      //
      double norm = 0.0;

      //
      // switch between thread versions
      //
      switch(d_version){

      case quasicontinuum::MULTI_THREADED:
      {
        //
        // instantiate computeResidualL2Norm_argument_t
        //
        computeResidualL2Norm_argument_t computeResidualL2Norm_argument;

        //
        // initialize computeL2ResidualNorm_argument
        //
        computeResidualL2Norm_argument.norm = 0.0;
        computeResidualL2Norm_argument.d_unknownCountFlag = &(d_unknownCountFlag);
        computeResidualL2Norm_argument.d_gradient = &(d_gradient);
        computeResidualL2Norm_argument.size = d_numberUnknowns;
        computeResidualL2Norm_argument.lock = &(d_lock);

        //
        // distribute work among workers
        //
        quasicontinuum::thread_monitor(threadComputeResidualL2NormWorker,
             static_cast<void *>(&computeResidualL2Norm_argument),
             quasicontinuum::get_max_number_threads());

        //
        // assign values
        //

        norm = computeResidualL2Norm_argument.norm;
        }

        break;

        case quasicontinuum::SINGLE_THREADED:
        {
          //
          // iterate over unknowns
          //
          for (int i = 0; i < d_numberUnknowns; ++i) 
          {
            #if 0 //defined(HAVE_MPI)
            const double factor = d_unknownCountFlag[i];
            #else
            const double factor = 1.0;
            #endif // HAVE_MPI

            const double gradient = d_gradient[i];
            norm += factor*gradient*gradient;
          }
        }
        break;
        }

        #if 0 //defined(HAVE_MPI)
        //
        // accumulate eta
        //
        norm = Accumulator().accumulate(norm);

        #endif // HAVE_MPI

        //
        // take square root
        // 
        norm = std::sqrt(norm);

        //
        //
        //
        return norm;

    }

    //
    // Compute the total number of unknowns in all processors.
    //
    int
    CGNonLinearSolver::computeTotalNumberUnknowns() const
    {

      #if 0 //defined(HAVE_MPI)
      //
      // initialize total number of unknowns
      //
      int totalNumberUnknowns = 0;

      //
      // switch between thread versions
      //
      switch(d_version){

      case quasicontinuum::MULTI_THREADED:
       {

       //
       // FIXME: Implement
       //
       std::cout 
         << "Need to implement threaded CGNonLinearSolver::computeTotalNumberUnknowns" 
          << std::endl;
       exit (1);

         }

       break;

      case quasicontinuum::SINGLE_THREADED:
       {

       //
       // iterate over unknowns
       //
       for (int i = 0; i < d_numberUnknowns; ++i) {

       const int factor = d_unknownCountFlag[i];

       totalNumberUnknowns += factor;

     }

     }

     break;

      }

      //
      // accumulate totalnumberUnknowns
      //
      totalNumberUnknowns = static_cast<int>(Accumulator().accumulate(totalNumberUnknowns));

      //
      //
      //
      return totalNumberUnknowns;

      #else
      return d_numberUnknowns;
      #endif // HAVE_MPI

    }

    //
    // Update solution x -> x + \alpha direction.
    //
    void
    CGNonLinearSolver::updateSolution(double                      alpha,
                                      const std::vector<double> & direction,
                                      SolverFunction            & function)
    {

      //
      // get the solution from solver function
      //
      function.solution(d_solution);

      //
      // get unknownCountFlag
      //
      const std::vector<int> unknownCountFlag = function.getUnknownCountFlag();

      //
      // get the size of solution
      //
      const std::vector<double>::size_type solutionSize = d_solution.size();

      //
      // switch between thread versions
      //
      switch(d_version)
      { 
        case quasicontinuum::MULTI_THREADED:
        {
          //
          // instantiate updateSolution_argument_t
          //
          updateSolution_argument_t updateSolution_argument;

          //
          // initialize updateSolution_argument
          //
          updateSolution_argument.alpha = alpha;
          updateSolution_argument.unknownCountFlag = &(unknownCountFlag);
          updateSolution_argument.direction = &(direction);
          updateSolution_argument.d_solution = &(d_solution);
          updateSolution_argument.size = solutionSize;

          //
          // distribute work among workers
          //
          quasicontinuum::thread_monitor(threadUpdateSolutionWorker,
             static_cast<void *>(&updateSolution_argument),
             quasicontinuum::get_max_number_threads());
        }
        break;

        case quasicontinuum::SINGLE_THREADED:
        {
          //
          // update solution
          //
          for (std::vector<double>::size_type i = 0; i < solutionSize; ++i)
            d_solution[i] = (d_solution[i] + alpha*direction[i])*unknownCountFlag[i];
        }

        break;
      }

      //
      // store solution
      //
      function.update(d_solution);

      //
      //
      //
      return;

    }

    //
    // Perform line search.
    //
    CGNonLinearSolver::ReturnValueType
    CGNonLinearSolver::lineSearch(SolverFunction & function,
                                  double           tolerance,
                                  int              maxNumberIterations,
                                  int              debugLevel)
    {
      //
      // local data
      //
      const double toleranceSqr = tolerance*tolerance;

      //
      // value of sigma0
      //
      const double sigma0 = 0.1;

      //
      // set the initial value of alpha
      //
      double alpha = -sigma0;

      //
      // update unknowns
      //
      CGNonLinearSolver::updateSolution(sigma0,
                                        d_direction,
                                        function);

      //
      // evaluate function gradient
      //
      function.gradient(d_gradient, NO_REBUILD);

      //
      // compute delta_d and eta_p
      //
      std::pair<double, double> deltaDReturnValue = 
        CGNonLinearSolver::computeDeltaD();
      double deltaD = deltaDReturnValue.first;
      double etaP   = deltaDReturnValue.second;
     
      //
      // update unknowns removing earlier update
      //
      CGNonLinearSolver::updateSolution(-sigma0,
                                        d_direction,
                                        function);

      //
      // begin iteration
      //
      for (int iter = 0; iter < maxNumberIterations; ++iter) 
      {
        //
        // evaluate function gradient
        //
        if(d_accumulated_r > d_rebuild_cutoff)
        {
          //
          // reset accumulated r
          //
          d_accumulated_r = 0.0;

          //
          // output
          //
          // std::cout<<"Rebuilding neighbor list..."<<std::endl;

          //
          // calculate gradient with neighbor rebuild
          //

          function.gradient(d_gradient, REBUILD);

        }

        else
        {
          //
          // calculate gradient without neighbor rebuild
          //
          function.gradient(d_gradient, NO_REBUILD);

        }

        //
        // compute eta
        //
        const double eta = CGNonLinearSolver::computeEta();

        //
        // update alpha
        //
        alpha *= eta/(etaP - eta);

        //
        // update accumulated r
        //
        double direction_max = 0.0;

        for(unsigned int iDirection = 0; iDirection < d_direction.size(); ++iDirection)
          if(std::abs(d_direction[iDirection]) > direction_max)
            direction_max = std::abs(d_direction[iDirection]);

        double direction_length = std::sqrt(3.0 * direction_max * direction_max);

        d_accumulated_r += 2.0*alpha*direction_length;

        //
        // update unknowns
        //
        CGNonLinearSolver::updateSolution(alpha,
                                          d_direction,
                                          function);

        //
        // update etaP
        //
        etaP = eta;

        //
        // output 
        //
        if (debugLevel >= 2)
          d_print("iteration : %d    alpha : %lf \n", iter, alpha);

        //
        // check for convergence
        //
        if (alpha*alpha*deltaD < toleranceSqr)
          return SUCCESS;

      }

      //
      //
      //
      return MAX_ITER_REACHED;

    }



    //
    // Perform function minimization.
    //
    NonLinearSolver::ReturnValueType
    CGNonLinearSolver::solve(SolverFunction & function,
                             double           tolerance,
                             int              maxNumberIterations,
                             int              debugLevel,
                             int              output_flag,
                             int              loadNumber,
                             int              minMethod,
                             int              minMethodMaxIterations,
                             double           minMethodTolerance)
    {
      #if 0 //defined(HAVE_MPI)
      //
      // get MPIController
      //
      const MPIController & mpiController = 
        MPIControllerSingleton::getInstance();

      //
      // get root task id
      //
      const MPIController::mpi_task_id_type rootTaskId = 
        mpiController.getRootId();

      //
      // get task id
      //
      const MPIController::mpi_task_id_type taskId = 
        mpiController.getId();

      //
      // get dot product factor
      //
      d_unknownCountFlag = function.getUnknownCountFlag();

      #endif // HAVE_MPI

      switch(minMethod)
      {
        case 0:
        {
          //
          //  for debug output in between iterations of conjugate
          //  gradient method, let the code know what kind of
          //  minimization we are running.
          //
          //  0 - simultaneously minimizing energy wrt u and w
          //  1 - wrt w
          //  2 - wrt u
          //
          minimizationMethod = 0;

          // minimization of free energy wrt u and w simultaneously
          d_print("CG Minimization: Minimizing free energy wrt u and w simultaneously\n");
          ReturnValueType returnValue = 
            CGMinimization(function,
              tolerance,
              maxNumberIterations,
              debugLevel,
              output_flag,
              loadNumber);

          d_print("CG Minimization: Done\n");

          if(returnValue == SUCCESS)
          {
            d_print("CG Minimization: Minimization converged\n");
            return returnValue;          
          }
          else if(returnValue == MAX_ITER_REACHED)
          {
            d_print("CG Minimization: Max iteration reached for CG Minimization\n");
            d_print("CG Minimization: Exiting the code\n");
            exit(EXIT_FAILURE);            
          }      
        }
        break;

        case 1:
        {
          // minimization of free energy wrt u and w alternatively
          d_print("CG Minimization: Minimizing free energy wrt u and w in an alternate way\n");

          ReturnValueType returnValue = 
            CGAlternateMinimization(function,
              tolerance,
              maxNumberIterations,
              debugLevel,
              output_flag,
              loadNumber,
              minMethodMaxIterations,
              minMethodTolerance);
          
          d_print("CG Minimization: Done\n");            

          return returnValue;          
        }
        break;

        case 2:
        {
          // used for debug output
          minimizationMethod = 2;
          
          // minimization of free energy wrt w keeping u fixed

          //
          //  first fix position degree of freedom of all nodes
          // 
          //  fixity_flag = 0 - reset the fixity mask
          //                1 - set fixity mask so that all 
          //                    position dof are fixed
          //                2 - set fixity mask so that all 
          //                    position dof are free
          //
          int fixity_flag = 1;
          function.setPositionFixity(fixity_flag);

          //
          //  call CGMinimization()
          //  Note: Since position dof is fixed, it will only
          //  minimize the frequency dof
          //
          d_print("CG Minimization: Minimizing free energy wrt w keeping u fixed\n");

          ReturnValueType returnValue = 
            CGMinimization(function,
              tolerance,
              maxNumberIterations,
              debugLevel,
              output_flag,
              loadNumber);

          d_print("CG Minimization: Done\n");            

          //
          //  set the fixity mask of position dof all nodes
          //  back to original value
          //
          fixity_flag = 0;
          function.setPositionFixity(fixity_flag);

          if(returnValue == SUCCESS)
          {
            d_print("CG Minimization: Frequency minimization converged\n");
            return returnValue;          
          }
          else if(returnValue == MAX_ITER_REACHED)
          {
            d_print("CG Minimization: Max iteration reached for frequency CG Minimization\n");
            d_print("CG Minimization: Exiting the code\n");
            exit(EXIT_FAILURE);            
          }      
        }
        break;

        case 3:
        {
          // used for debug output
          minimizationMethod = 3;

          // minimization of free energy wrt u keeping w fixed

          //
          //  first fix frequency degree of freedom of all nodes
          // 
          //  fixity_flag = 0 - reset the fixity mask
          //                1 - set fixity mask so that all 
          //                    frequency dof are fixed
          //                2 - set fixity mask so that all 
          //                    frequency dof are free
          //
          int fixity_flag = 1;
          function.setFrequencyFixity(fixity_flag);

          //
          //  call CGMinimization()
          //  Note: Since frequency dof is fixed, it will only
          //  minimize the position dof
          //
          d_print("CG Minimization: Minimizing free energy wrt u keeping w fixed\n");

          ReturnValueType returnValue = 
            CGMinimization(function,
              tolerance,
              maxNumberIterations,
              debugLevel,
              output_flag,
              loadNumber);

          //
          //  set the fixity mask of position dof all nodes
          //  back to original value
          //
          fixity_flag = 0;
          function.setFrequencyFixity(fixity_flag);

          if(returnValue == SUCCESS)
          {
            d_print("CG Minimization: Position minimization converged\n");
            return returnValue;          
          }
          else if(returnValue == MAX_ITER_REACHED)
          {
            d_print("CG Minimization: Max iteration reached for position CG Minimization\n");
            d_print("CG Minimization: Exiting the code\n");
            exit(EXIT_FAILURE);            
          }

        }
        break;

        default:
        {
          // invalid value of minMethod
          d_print("Check minMethod = %d in quasi init file\n", minMethod);
          exit(EXIT_FAILURE);
        }
        break;        
      } 

      //
      // should not reach here
      //
      d_print("Check CGNonLinearSolver::solve()\n");
      exit(EXIT_FAILURE);
    }


    //
    //  CGMinimization()  
    //
    NonLinearSolver::ReturnValueType
    CGNonLinearSolver::CGMinimization(SolverFunction & function,
                                 double           tolerance,
                                 int              maxNumberIterations,
                                 int              debugLevel,
                                 int              output_flag,
                                 int              loadNumber)
    {
      //
      // method const data
      //
      const double toleranceSqr = tolerance*tolerance;

      //
      //  recompute the data of SolverFunction instance, as
      //  we may have changed the fixity mask for minimization
      //
      function.recomputeData();

      //
      // get total number of unknowns in the problem.
      //
      d_numberUnknowns = function.getNumberUnknowns();

      //
      // get total number of unknowns
      //
      const int totalnumberUnknowns = 
        CGNonLinearSolver::computeTotalNumberUnknowns();

      //
      // allocate space for direction, gradient and old gradient
      // values.
      //
      d_direction.resize(d_numberUnknowns);
      d_gradient.resize(d_numberUnknowns);
      d_sOld.resize(d_numberUnknowns);
      d_s.resize(d_numberUnknowns);

      //
      // compute initial values of function and function gradient
      //
      function.gradient(d_gradient, NO_REBUILD);

      double functionValue = function.value();
      
      //
      // apply preconditioner
      //
      function.precondition(d_s, d_gradient);

      //
      // initialize delta new and direction
      //
      CGNonLinearSolver::initializeDirection();
      
      //
      // check for convergence
      //
      #if 0
      if ( d_deltaNew < toleranceSqr*totalnumberUnknowns*totalnumberUnknowns)
        return SUCCESS;
      #endif

      if (computeResidualL2Norm() < tolerance)
      {
        d_print("Conjugate Gradient converged at the beginning.\n");
        return SUCCESS;
      }

      //
      // main CG iteration loop
      //
      #if 0 //defined(HAVE_MPI)
      if (taskId == rootTaskId)
      #endif // HAVE_MPI 
      if (debugLevel >= 1)
        d_print("Iteration no. |  delta new  |  resdiual norm   |   residual norm avg \n");
        
      for (d_iter = 0; d_iter < maxNumberIterations; ++d_iter) 
      {
        //
        // compute L2-norm of the residual (gradient)
        //
        const double residualNorm = computeResidualL2Norm();

        //
        // output at the begining of the iteration
        //
        #if 0 //defined(HAVE_MPI)
        if (taskId == rootTaskId)
        #endif // HAVE_MPI
          
        if (debugLevel >= 1) 
          d_print("CG Iteration : %d               %lf         %lf           %lf \n",
            d_iter,
            d_deltaNew,
            residualNorm,
            residualNorm/totalnumberUnknowns);

        // if(residualNorm > 1000.0)
        // {
        //   d_print("Residual norm is above acceptable value\n");
        //   d_print("exiting code\n");
        //   exit(EXIT_SUCCESS);
        // }

        //
        // perform line search along direction
        //
        ReturnValueType lineSearchReturnValue = 
          CGNonLinearSolver::lineSearch(function,
                                        d_lineSearchTolerance,
                                        d_lineSearchMaxIterations,
                                        debugLevel);
          
        //
        // evaluate gradient
        //
        function.gradient(d_gradient, NO_REBUILD);

        //
        //  check if we need to output data
        //
        if(output_flag == 1)
        {
          if(d_iter >= 4 && d_iter % 9 == 0)
          {
            // since we are making an restart data output
            // also, set the fixity to original value and reset it back
            if(minimizationMethod == 1)
            {
              if(freq_or_position_min == 0)
              {
                // for freq minimization, we may have set the
                // position fixity to fixed. Reset it to original
                function.setPositionFixity(0);
              }
              else
              {
                // for position minimization, we may have set the
                // frequency fixity to fixed. Reset it to original
                function.setFrequencyFixity(0);
              }
            }
            else if(minimizationMethod == 2)
            {
              // this is freq minimization hence position fixity
              // are fixed. Set it to original
              function.setPositionFixity(0);
            }
            else if(minimizationMethod == 3)
            {
              // this is position minimization hence frequency fixity
              // are fixed. Set it to original
              function.setFrequencyFixity(0);
            }

            // set flag            
            std::vector<int> output_flag_vector(6,0);
            output_flag_vector[0] = minimizationMethod;
            output_flag_vector[1] = freq_or_position_min;
            output_flag_vector[2] = d_iter;
            output_flag_vector[3] = alternateIter;
            output_flag_vector[4] = 0; // inside iteration loop
            output_flag_vector[5] = loadNumber;

            unsigned int flags;
            flags |= NODE_OUTPUT_FLAG;
            flags |= RESTART_FLAG;

            d_print("CG Minimization : Performing output for all quasis\n");
            Output::getInstance()->performOutputCentral(flags, output_flag_vector);
            d_print("CG Minimization : Done\n");

            // since we are making an restart data output
            // also, set the fixity to the value it was before
            // making an output
            if(minimizationMethod == 1)
            {
              if(freq_or_position_min == 0)
              {
                // for freq minimization, we may have set the
                // position fixity to fixed. Reset it to original
                function.setPositionFixity(1);
              }
              else
              {
                // for position minimization, we may have set the
                // frequency fixity to fixed. Reset it to original
                function.setFrequencyFixity(1);
              }
            }
            else if(minimizationMethod == 2)
            {
              // this is freq minimization hence position fixity
              // are fixed. Set it to original
              function.setPositionFixity(1);
            }
            else if(minimizationMethod == 3)
            {
              // this is position minimization hence frequency fixity
              // are fixed. Set it to original
              function.setFrequencyFixity(1);
            }
          }
        }
        else
        {
          // make restart data output every 15 iteration
          if(d_iter >= 4 && d_iter % 15 == 0)
          {
            // since we are making an restart data output
            // also, set the fixity to original value and reset it back
            if(minimizationMethod == 1)
            {
              if(freq_or_position_min == 0)
              {
                // for freq minimization, we may have set the
                // position fixity to fixed. Reset it to original
                function.setPositionFixity(0);
              }
              else
              {
                // for position minimization, we may have set the
                // frequency fixity to fixed. Reset it to original
                function.setFrequencyFixity(0);
              }
            }
            else if(minimizationMethod == 2)
            {
              // this is freq minimization hence position fixity
              // are fixed. Set it to original
              function.setPositionFixity(0);
            }
            else if(minimizationMethod == 3)
            {
              // this is position minimization hence frequency fixity
              // are fixed. Set it to original
              function.setFrequencyFixity(0);
            }

            // set flag            
            std::vector<int> output_flag_vector(6,0);
            output_flag_vector[0] = minimizationMethod;
            output_flag_vector[1] = freq_or_position_min;
            output_flag_vector[2] = d_iter;
            output_flag_vector[3] = alternateIter;
            output_flag_vector[4] = 0; // inside iteration loop
            output_flag_vector[5] = loadNumber;

            unsigned int flags;
            flags |= NODE_OUTPUT_FLAG;
            flags |= RESTART_FLAG;

            d_print("CG Minimization : Performing output for all quasis\n");
            Output::getInstance()->performOutputCentral(flags, output_flag_vector);
            d_print("CG Minimization : Done\n");

            // since we are making an restart data output
            // also, set the fixity to the value it was before
            // making an output
            if(minimizationMethod == 1)
            {
              if(freq_or_position_min == 0)
              {
                // for freq minimization, we may have set the
                // position fixity to fixed. Reset it to original
                function.setPositionFixity(1);
              }
              else
              {
                // for position minimization, we may have set the
                // frequency fixity to fixed. Reset it to original
                function.setFrequencyFixity(1);
              }
            }
            else if(minimizationMethod == 2)
            {
              // this is freq minimization hence position fixity
              // are fixed. Set it to original
              function.setPositionFixity(1);
            }
            else if(minimizationMethod == 3)
            {
              // this is position minimization hence frequency fixity
              // are fixed. Set it to original
              function.setFrequencyFixity(1);
            }
          }          
        }          
        
        //
        // apply preconditioner
        //
        function.precondition(d_s, d_gradient);

        //
        // update values of delta_new and delta_mid
        //
        d_deltaOld = d_deltaNew;
        CGNonLinearSolver::computeDeltas();

        //
        // compute beta
        //
        d_beta = (lineSearchReturnValue == SUCCESS) ? 
           (d_deltaNew - d_deltaMid)/d_deltaOld : 0.0;

        //
        // update direction
        //
        CGNonLinearSolver::updateDirection();

        //
        // check for convergence
        //
        #if 0
        if (d_deltaNew < toleranceSqr*totalnumberUnknowns*totalnumberUnknowns)
          break;
        #endif 

        if (residualNorm < tolerance)
          break;

        // // max number of iterations
        // if(d_iter > 100)
        //   break;
      }

      // if(output_flag == 1)
      // {
        // since we are making an restart data output
        // also, set the fixity to original value and reset it back
        if(minimizationMethod == 1)
        {
          if(freq_or_position_min == 0)
          {
            // for freq minimization, we may have set the
            // position fixity to fixed. Reset it to original
            function.setPositionFixity(0);
          }
          else
          {
            // for position minimization, we may have set the
            // frequency fixity to fixed. Reset it to original
            function.setFrequencyFixity(0);
          }
        }
        else if(minimizationMethod == 2)
        {
          // this is freq minimization hence position fixity
          // are fixed. Set it to original
          function.setPositionFixity(0);
        }
        else if(minimizationMethod == 3)
        {
          // this is position minimization hence frequency fixity
          // are fixed. Set it to original
          function.setFrequencyFixity(0);
        }

        unsigned int flags;
        flags |= NODE_OUTPUT_FLAG;
        flags |= RESTART_FLAG;

        std::vector<int> output_flag_vector(6,0);
        output_flag_vector[0] = minimizationMethod;
        output_flag_vector[1] = freq_or_position_min;
        output_flag_vector[2] = d_iter;
        output_flag_vector[3] = alternateIter;
        output_flag_vector[4] = 1; // converged cg
        output_flag_vector[5] = loadNumber;

        d_print("CG Minimization : Performing output for all quasis\n");
        Output::getInstance()->performOutputCentral(flags, output_flag_vector);
        d_print("CG Minimization : Done\n"); 
      // }

      //
      // set error condition
      //
      ReturnValueType returnValue = SUCCESS;
      
      if(d_iter == maxNumberIterations) 
        returnValue = MAX_ITER_REACHED;

      //
      // compute function value
      //
      functionValue = function.value();
      
      //
      // final output
      //
      #if 0 //defined(HAVE_MPI)
      if (taskId == rootTaskId)
      #endif // HAVE_MPI 
      
      if (debugLevel >= 1) 
      {
        if (returnValue == SUCCESS) 
        {
          d_print("Conjugate Gradient converged after %d iterations.\n", d_iter);
        } 
        else 
        {
          d_print("Conjugate Gradient failed to converge after d_iter iterations.\n", d_iter);
        }

        d_print("Final function value: %lf\n", functionValue);
      }

      //
      //
      //
      return SUCCESS;  
    }  

    //
    //  Minimizing free energy wrt u and w in an alternate way
    //
    NonLinearSolver::ReturnValueType
    CGNonLinearSolver::CGAlternateMinimization(SolverFunction & function,
                     double           tolerance,
                     int              maxNumberIterations,
                     int              debugLevel,
                     int              output_flag,
                     int              loadNumber,
                     int              minMethodMaxIterations,
                     double           minMethodTolerance)
    {
      //
      //  compute residual norm in the beginning
      //

      //
      // recompute the data of SolverFunctio instance
      // 
      function.recomputeData();

      //
      // get total number of unknowns in the problem.
      //
      d_numberUnknowns = function.getNumberUnknowns();

      //
      // get total number of unknowns
      //
      const int totalnumberUnknowns = 
        CGNonLinearSolver::computeTotalNumberUnknowns();

      //
      // allocate space for direction, gradient and old gradient
      // values.
      //
      d_direction.resize(d_numberUnknowns);
      d_gradient.resize(d_numberUnknowns);
      d_sOld.resize(d_numberUnknowns);
      d_s.resize(d_numberUnknowns);

      //
      // compute initial values of function and function gradient
      //
      function.gradient(d_gradient, NO_REBUILD);

      double functionValue = function.value();
      
      //
      // apply preconditioner
      //
      function.precondition(d_s, d_gradient);

      //
      // initialize delta new and direction
      //
      CGNonLinearSolver::initializeDirection();
      
      //
      // check for convergence
      //
      if (computeResidualL2Norm() < minMethodTolerance)
      {
        d_print("Conjugate Gradient converged at the beginning.\n");
        d_print("Residual Norm = %lf \n", computeResidualL2Norm());
        return SUCCESS;
      }

      //
      //  loop over minimization wrt frequency and position
      //  in an alternate way
      //
      for(int d_minIter = 0; d_minIter < minMethodMaxIterations; d_minIter++)
      {
        alternateIter = d_minIter;
        minimizationMethod = 1;
        freq_or_position_min = 0;
        
        d_print("CG Minimization: Iteration = %d\n", d_minIter);

        //
        //  first minimize wrt to frequency
        //
        // fix position dof for freq minimization
        int fixity_flag = 1;
        function.setPositionFixity(fixity_flag);    

        // minimize wrt frequency
        d_print("Iteration = %d - Frequency minimization\n", d_minIter);

        ReturnValueType returnValue = 
          CGMinimization(function,
            tolerance,
            maxNumberIterations,
            debugLevel,
            output_flag,
            loadNumber);
        
        // reset the position dof
        fixity_flag = 0;
        function.setPositionFixity(fixity_flag);

        if(returnValue == SUCCESS)
        {  
          d_print("Iteration = %d - Frequency minimization converged\n", d_minIter);
        }
        else if(returnValue == MAX_ITER_REACHED)
        {
          d_print("Iteration = %d - Frequency minimization max iteration reached\n", d_minIter); 
          d_print("Iteration = %d - Exiting the minimization\n", d_minIter); 
          exit(EXIT_FAILURE);
        }

        //
        //  minimize wrt to position
        //
        // fix frequency dof for freq minimization
        fixity_flag = 1;
        function.setFrequencyFixity(fixity_flag);    

        // minimize wrt frequency
        d_print("Iteration = %d - Position minimization\n", d_minIter);

        // set the parameter for output
        freq_or_position_min = 1;

        returnValue = 
          CGMinimization(function,
            tolerance,
            maxNumberIterations,
            debugLevel,
            output_flag,
            loadNumber);
        
        // reset the frequency dof
        fixity_flag = 0;
        function.setFrequencyFixity(fixity_flag);  

        if(returnValue == SUCCESS)
        {  
          d_print("Iteration = %d - Position minimization converged\n", d_minIter);
        }
        else if(returnValue == MAX_ITER_REACHED)
        {
          d_print("Iteration = %d - Position minimization max iteration reached\n", d_minIter); 
          d_print("Iteration = %d - Exiting the minimization\n", d_minIter); 
          exit(EXIT_FAILURE);
        }

        //
        //  defer the checking of convergence till the end of 
        //  second loop
        //
        if(d_minIter > 0)
        {
          //
          //  check if the minimized state has forces below 
          //  tolerance
          //

          //
          // reset the SolverFunctio instance, as fixity mask is modified
          //
          function.recomputeData();

          //
          // get total number of unknowns in the problem.
          //
          d_numberUnknowns = function.getNumberUnknowns();

          //
          // get total number of unknowns
          //
          const int totalnumberUnknowns = 
            CGNonLinearSolver::computeTotalNumberUnknowns();

          //
          // allocate space for direction, gradient and old gradient
          // values.
          //
          d_direction.resize(d_numberUnknowns);
          d_gradient.resize(d_numberUnknowns);
          d_sOld.resize(d_numberUnknowns);
          d_s.resize(d_numberUnknowns);

          //
          // compute initial values of function and function gradient
          //
          function.gradient(d_gradient, NO_REBUILD);

          double functionValue = function.value();
          
          //
          // apply preconditioner
          //
          function.precondition(d_s, d_gradient);
        
          //
          // check for convergence
          //
          d_print("Iteration = %d - Checking convergence for CG Alternate Minimization method\n", d_minIter);
          if (computeResidualL2Norm() < minMethodTolerance)
          {
            d_print("Iteration = %d - Converged\n", d_minIter);
            d_print("Residual Norm = %lf \n", computeResidualL2Norm());
            d_print("CG Minimization: Final function value: %lf\n", functionValue);

            // output
            unsigned int flags;
            flags |= NODE_OUTPUT_FLAG;
            flags |= RESTART_FLAG;

            std::vector<int> output_flag_vector(6,0);
            output_flag_vector[0] = minimizationMethod;
            output_flag_vector[1] = freq_or_position_min;
            output_flag_vector[2] = d_iter;
            output_flag_vector[3] = alternateIter;
            output_flag_vector[4] = 2; // converged cg alternate
            output_flag_vector[5] = loadNumber;

            d_print("CG Minimization : Performing output for all quasis\n");
            Output::getInstance()->performOutputCentral(flags, output_flag_vector);
            d_print("CG Minimization : Done\n");
            return SUCCESS;            
          }
        } // if d_minIter > 0
      }

      //
      //  reached here means we have reached max number of iterations
      //  and code did not converge
      //
      ReturnValueType returnValue = MAX_ITER_REACHED;

      d_print("CG Minimization: Max iteration reached for Alternate CG Minimization\n");
      d_print("CG Minimization: Exiting the code\n");
      exit(EXIT_FAILURE);
    }        

}
