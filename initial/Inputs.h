/*****************************************************************************************
CTSat -- Copyright (c) 2020, Marc Hartung
                        Zuse Institute Berlin, Germany

Maple_LCM_Dist_Chrono -- Copyright (c) 2018, Vadim Ryvchin, Alexander Nadel

GlucoseNbSAT -- Copyright (c) 2016,Chu Min LI,Mao Luo and Fan Xiao
                           Huazhong University of science and technology, China
                           MIS, Univ. Picardie Jules Verne, France

MapleSAT -- Copyright (c) 2016, Jia Hui Liang, Vijay Ganesh

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/
#ifndef SOURCES_INITIAL_INPUTS_H_
#define SOURCES_INITIAL_INPUTS_H_

#include "utils/Options.h"

namespace ctsat
{
class Inputs
{
 public:


   static StringOption database;
   static StringOption branch;
   static StringOption reduce;
   static StringOption restart;
   static StringOption exchange;
   static StringOption propagate;
   static StringOption analyze;

   static BoolOption LAA_alwaySwap;
   static IntOption LAA_levelDiffEnforce;
   static IntOption LAA_numInitialConflicts;
   static IntOption LAA_levelQueueSz;


   static DoubleOption step_size;
   static DoubleOption step_size_dec;
   static DoubleOption min_step_size;
   static DoubleOption var_decay;
   static DoubleOption clause_decay;
   static DoubleOption random_seed;
   static IntOption ccmin_mode;
   static IntOption maxEntendedBinaryResolutionSz;
   static IntOption maxFullImplicationMinLbd;
   static IntOption phase_saving;
   static BoolOption rnd_init_act;
   static BoolOption rnd_polarity;
   static BoolOption initVarPolZero;
   static IntOption restart_first;
   static DoubleOption restart_inc;
   static DoubleOption garbage_frac;
   static IntOption chrono;
   static IntOption conf_to_chrono;
   static IntOption conflict_budget;
   static IntOption propagation_budget;

   static IntOption first_reduce_db;
   static IntOption inc_reduce_db;
   static IntOption spec_inc_reduce_db;
   static IntOption maxProtectableLbd;

   static BoolOption useVivification;
   static BoolOption use_elim;
   static IntOption grow;
   static IntOption clause_lim;
   static IntOption subsumption_lim;
   static DoubleOption simp_garbage_frac;

   static IntOption verb;
   static DoubleOption print_interval;
   static IntOption cpu_lim;
   static IntOption mem_lim;

   static BoolOption verifySat;
   static BoolOption model;
   static BoolOption drup;
   static StringOption drup_file;

   static BoolOption onlyExportWhenMin;
   static BoolOption minimize_import_cl;
   static IntOption nThreads;
   static IntOption max_export_lbd;
   static IntOption max_export_sz;
   static IntOption numConflictsToDelete;

   static DoubleOption mbExchangeBufferPerThread;

   static DoubleOption mpiMbBufferSize;
   static DoubleOption mpi_send_interval;
   static IntOption nMpiPartitions;
   static BoolOption mpiAutoThreads;
   static BoolOption pinSolver;
   static BoolOption log;
   static BoolOption mpiHashClauseFilter;
   static StringOption logDirectory;

   static int * argc;
   static char *** argv;

   static void setArgs(int * argc, char *** argv)
   {
      Inputs::argc = argc;
      Inputs::argv = argv;
   }

};

}

#endif /* SOURCES_INITIAL_INPUTS_H_ */
