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
#include "initial/Inputs.h"

using namespace ctsat;

const char * _logger = "LOGGER";
BoolOption Inputs::log(_logger, "logging", "Logs events of each solver", false);
StringOption Inputs::logDirectory(_logger, "log-dir", "Directory in which log file will be created",
                                  "");

const char * _analyze = "ANALYZE";
StringOption Inputs::analyze(_analyze, "analyze", "possible options are 'fuip', 'muip' or 'laa'",
                             "laa");

BoolOption Inputs::LAA_alwaySwap(
      _analyze,
      "laa-always-swap",
      "If an additional learnt clause is asserting and smaller than the fuip clause, they will be used",
      true);

IntOption Inputs::LAA_levelDiffEnforce(
      _analyze,
      "laa-level-diff",
      "when the average conflict level is this higher than the current, additional clauses will be learnt, -1 disables additional clause learning",
      -1, IntRange(-1, INT32_MAX));

IntOption Inputs::LAA_numInitialConflicts(
      _analyze, "laa-initial-confl",
      "Number of conflicts to learn additional clauses after solver start", 10000,
      IntRange(0, INT32_MAX));

IntOption Inputs::LAA_levelQueueSz(
      _analyze, "laa-queue-sz",
      "Length of the queue to determine the current average conflict level", 20,
      IntRange(1, INT32_MAX));

const char * _mpi = "MPI - has only effect when mpi build is used";
BoolOption Inputs::mpiAutoThreads(
      _mpi, "mpi-auto-threads",
      "Automatically chooses different thread numbers for different mpi ranks", false);

BoolOption Inputs::mpiHashClauseFilter(
      _mpi, "mpi-hash-filter",
      "Uses a hash filter to exclude already seen clauses from being imported", true);

DoubleOption Inputs::mpiMbBufferSize(_mpi, "mpi-buffer-sz",
                                     "Number of MB allocated for send buffers", 0.5,
                                     DoubleRange(0.1, true, HUGE_VAL, false));

DoubleOption Inputs::mpi_send_interval(_mpi, "mpi-interval", "Seconds between mpi send operations",
                                       0.05, DoubleRange(0.0, true, HUGE_VAL, false));

IntOption Inputs::nMpiPartitions(
      _mpi,
      "mpi-partition",
      "Number of partitions to distribute mpi processes over. Nodes in the same partition exchange more clauses than nodes from different partitions",
      0, IntRange(0, INT32_MAX));

const char * _parallel = "PARALLEL - no effect when serial solver is used -";
StringOption Inputs::exchange(_parallel, "exchange", "possible options are 'none' or 'importbuff'",
                              "importbuff");

IntOption Inputs::nThreads(_parallel, "nthreads", "number of threads to run in parallel", 2,
                           IntRange(1, 512));
DoubleOption Inputs::mbExchangeBufferPerThread(
      _parallel, "mb-exchange", "number of mega bytes per thread to use for clause exchange buffer",
      1.0, DoubleRange(0.01, true, HUGE_VAL, false));
IntOption Inputs::max_export_lbd(_parallel, "exp-lbd",
                                 "Maximal allowed lbd of a clause to be exported", 5,
                                 IntRange(1, 10));
IntOption Inputs::max_import_lbd(_parallel, "exp-lbd",
                                 "Maximal allowed lbd of a clause to be imported", 5,
                                 IntRange(1, 10));

IntOption Inputs::max_export_sz(_parallel, "exp-sz",
                                "Maximal allowed size of a clause to be exported", 30,
                                IntRange(1, 1000));
IntOption Inputs::numConflictsToDelete(
      _parallel, "nconfl-to-delete",
      "Number of conflicts after import the clause is deleted when unused", 15000,
      IntRange(1, INT32_MAX));
BoolOption Inputs::minimize_import_cl(_parallel, "min-import-cl",
                                      "Allows minimization of imported clauses", false);

BoolOption Inputs::onlyExportWhenMin(_parallel, "min-export-cl",
                                     "Only export clauses when tried to vifify them", false);

BoolOption Inputs::pinSolver(_parallel, "pin-solvers", "Pins solvers to cores", false);

const char * _bra = "BRANCH";
StringOption Inputs::branch(_bra, "branch",
                            "possible options are 'dist_mixed', 'dist', 'vsids' or 'lrb",
                            "dist_mixed");
DoubleOption Inputs::step_size(_bra, "step-size", "Initial step size", 0.40,
                               DoubleRange(0, false, 1, false));
DoubleOption Inputs::step_size_dec(_bra, "step-size-dec", "Step size decrement", 0.000001,
                                   DoubleRange(0, false, 1, false));
DoubleOption Inputs::min_step_size(_bra, "min-step-size", "Minimal step size", 0.06,
                                   DoubleRange(0, false, 1, false));
DoubleOption Inputs::vsids_var_decay(_bra, "var-decay", "The variable activity decay factor", 0.80,
                                     DoubleRange(0, false, 1, false));
DoubleOption Inputs::vsids_max_var_decay(_bra, ",max-var-decay",
                                         "The maximum variable activity decay factor", 0.95,
                                         DoubleRange(0, false, 1, false));
IntOption Inputs::phase_saving(_bra, "phase-saving",
                               "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2,
                               IntRange(0, 2));
IntOption Inputs::vsids_var_decay_timer(
      _bra, "phase-saving", "Defines after how many conflicts the variable decay factor will be increased", 5000,
      IntRange(1000, INT32_MAX));

BoolOption Inputs::rnd_init_act(_bra, "rnd-activ", "Randomize the initial activity", false);
BoolOption Inputs::rnd_polarity(_bra, "rnd-pol", "Randomize the initial var polarity", false);
BoolOption Inputs::initVarPolZero(
      _bra, "var-init-false",
      "If not random all vars polarity will be set to false. When disabled to true", true);

DoubleOption Inputs::random_seed(_bra, "rnd-seed", "Used by the random variable selection",
                                 91648253, DoubleRange(0, false, HUGE_VAL, false));

const char * _min = "MINIMIZE";
IntOption Inputs::ccmin_mode(_min, "ccmin-mode",
                             "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2,
                             IntRange(0, 2));
IntOption Inputs::maxEntendedBinaryResolutionSz(
      _min, "max-ext-bin-res",
      "Maximal size to use extended binary resolution minimization during conflict analyzes", 5,
      IntRange(0, 30));
IntOption Inputs::maxFullImplicationMinLbd(_min, "max-full-impl",
                                           "Maximal lbd to use full implication graph minimization",
                                           3, IntRange(0, 30));

BoolOption Inputs::useVivification(_min, "vivi",
                                   "Uses vivification during restart to minimize clauses", true);

const char * _rest = "RESTART";
StringOption Inputs::restart(_rest, "restart", "possible options are 'mixed', 'luby' or 'glucose'",
                             "mixed");
IntOption Inputs::luby_base_factor(_rest, "luby-base", "The base restart interval for luby restarts", 100,
                                IntRange(1, INT32_MAX));
DoubleOption Inputs::luby_inc_factor(_rest, "rinc", "Restart interval increase factor", 2,
                                 DoubleRange(1, false, HUGE_VAL, false));

const char* _back = "BACKTRACK";
IntOption Inputs::chrono(_back, "chrono", "Controls if to perform chrono backtrack", 100,
                         IntRange(-1, INT32_MAX));
IntOption Inputs::conf_to_chrono(_back, "confl-to-chrono",
                                 "Controls number of conflicts to perform chrono backtrack", 4000,
                                 IntRange(-1, INT32_MAX));

const char* _red = "REDUCE";
StringOption Inputs::reduce(_red, "reduce", "possible options are 'chanseok' or 'glucose'",
                            "chanseok");
DoubleOption Inputs::clause_decay(_red, "cla-decay", "The clause activity decay factor", 0.999,
                                  DoubleRange(0, false, 1, false));

IntOption Inputs::first_reduce_db(_red, "firstReduceDB",
                                  "The number of conflicts before the first reduce DB", 2000,
                                  IntRange(0, INT32_MAX));
IntOption Inputs::inc_reduce_db(_red, "incReduceDB", "Increment for reduce DB", 300,
                                IntRange(0, INT32_MAX));
IntOption Inputs::spec_inc_reduce_db(_red, "specialIncReduceDB", "Special increment for reduce DB",
                                     1000, IntRange(0, INT32_MAX));
IntOption Inputs::maxProtectableLbd(
      _red,
      "max-prot-lbd",
      "Clauses with lbd an updated LBD equal or less than this will be protected during the next reduce",
      30, IntRange(0, INT32_MAX));

const char* _simp = "SIMP";
BoolOption Inputs::use_elim(_simp, "elim", "Perform variable elimination.", true);
IntOption Inputs::grow(_simp, "grow",
                       "Allow a variable elimination step to grow by a number of clauses.", 0);
IntOption Inputs::clause_lim(
      _simp,
      "cl-lim",
      "Variables are not eliminated if it produces a resolvent with a length above this limit. -1 means no limit",
      20, IntRange(-1, INT32_MAX));
IntOption Inputs::subsumption_lim(
      _simp, "sub-lim",
      "Do not check if subsumption against a clause larger than this. -1 means no limit.", 1000,
      IntRange(-1, INT32_MAX));
DoubleOption Inputs::simp_garbage_frac(
      _simp,
      "simp-gc-frac",
      "The fraction of wasted memory allowed before a garbage collection is triggered during simplification.",
      0.4, DoubleRange(0, false, HUGE_VAL, false));

const char* _main = "MAIN";
IntOption Inputs::verb(_main, "verb", "Verbosity level (0=silent, 1=some, 2=more).", 1,
                       IntRange(0, 2));
IntOption Inputs::secToSwitchHeuristic(_main, "sec-heu-switch", "Seconds until the solver switches heuristics", 2500,
                       IntRange(10, INT32_MAX));

DoubleOption Inputs::print_interval(_main, "print-interval",
                                    "Elapsed time between solver state prints", 5.0,
                                    DoubleRange(1.0, false, HUGE_VAL, false));
const char* _budg = "RESOURCES";
IntOption Inputs::conflict_budget(_budg, "budget-confl",
                                  "Maximal number of conflicts before abort solving", -1,
                                  IntRange(-1, INT32_MAX));
IntOption Inputs::propagation_budget(_budg, "budget-prop",
                                     "Maximal number of propagations before abort solving", -1,
                                     IntRange(-1, INT32_MAX));

IntOption Inputs::cpu_lim(_budg, "cpu-lim", "Limit on CPU time allowed in seconds.\n", INT32_MAX,
                          IntRange(0, INT32_MAX));
IntOption Inputs::mem_lim(_budg, "mem-lim", "Limit on memory usage in megabytes.\n", INT32_MAX,
                          IntRange(0, INT32_MAX));

StringOption Inputs::database(_main, "database", "possible options are 'minisat' and 'sticky'",
                              "minisat");
DoubleOption Inputs::garbage_frac(
      _main, "gc-frac",
      "The fraction of wasted memory allowed before a garbage collection is triggered", 0.20,
      DoubleRange(0, false, HUGE_VAL, false));
BoolOption Inputs::verifySat(_main, "verify-sat",
                             "On sat answere, the solution is checked against the dimacs file.",
                             true);
BoolOption Inputs::model(_main, "model", "Print the SAT model", false);
BoolOption Inputs::drup(_main, "drup", "Generate DRUP UNSAT proof.", false);
StringOption Inputs::drup_file(_main, "drup-file", "DRUP UNSAT proof ouput file.", "");

int * Inputs::argc = nullptr;
char *** Inputs::argv = nullptr;
