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

#include "mtl/Sort.h"
#include "utils/System.h"
#include "utils/Timer.h"
#include "Preprocessor.h"

#include <cstdio>
#include <zlib.h>

namespace ctsat {

//=================================================================================================
// Options:

//=================================================================================================
// Constructor/Destructor:

Preprocessor::Preprocessor(SolverConfig const & config) :
		elim(config.elim), ok(true), verb(config.verbosity), grow(config.grow), clause_lim(
				config.clause_lim), subsumption_lim(config.subsumption_lim), nRedundantAssignedLits(
				0), bwdsub_assigns(0), n_touched(0), eliminated_vars(0), bwdsub_tmpunit(
				ClauseAllocator::npos()), simp_garbage_frac(
				config.simp_garbage_frac), occurs(ClauseDeleted(ca)), elim_heap(
				ElimLt(n_occ)), randEngine(config.rnd_seed), drat(
				config.drat_file), ca(), ig(ca), branch(
				Branch<ClauseAllocator>::BranchInputArgs(config, smode,
						randEngine, stat, ca, ig)), propEngine(stat, ca, ig) {
	vec<Lit> dummy(1, Lit::Undef());
	ca.extra_clause_field = true;
	bwdsub_tmpunit = ca.alloc(dummy, false);
}

Preprocessor::~Preprocessor() {
}

bool Preprocessor::readInstance_main(StreamBuffer & in) {
	vec<Lit> lits;
	int vars = 0;
	int clauses = 0;
	int cnt = 0;
	for (;;) {
		skipWhitespace(in);
		if (*in == EOF)
			break;
		else if (*in == 'p') {
			if (eagerMatch(in, "p cnf")) {
				vars = parseInt(in);
				clauses = parseInt(in);
			} else
				printf("PARSE ERROR! Unexpected char: %c\n", *in), throw InputException();
		} else if (*in == 'c' || *in == 'p')
			skipLine(in);
		else {
			cnt++;
			readInstance_clause(in, lits);
			if (!addClause(lits, true))
				return false;
		}
	}
	if (vars != ig.nVars())
		fprintf(stderr,
				"c WARNING! DIMACS header mismatch: wrong number of variables.\n");
	if (cnt != clauses)
		fprintf(stderr,
				"c WARNING! DIMACS header mismatch: wrong number of clauses.\n");
	return true;
}
void Preprocessor::readInstance_clause(StreamBuffer & in, vec<Lit> & lits) {
	int parsed_lit, var;
	lits.clear();
	for (;;) {
		parsed_lit = parseInt(in);
		if (parsed_lit == 0)
			break;
		var = abs(parsed_lit) - 1;
		while (var >= ig.nVars())
			newVar();
		lits.push(Lit(var, parsed_lit < 0));
	}
}

bool Preprocessor::readInstance(std::string const & filename) {
	double initial_time = cpuTime();
	gzFile in = gzopen(filename.c_str(), "rb");
	if (in == NULL)
		printf("c ERROR! Could not open file: %s\n", filename.c_str()), throw InputException();

	StreamBuffer strb(in);
	bool res = readInstance_main(strb);

	gzclose(in);
	if (verb > 0) {
		printf(
				"c ##############################   Initial  ############################\n");
		size_t pos = filename.find_last_of('/') + 1;
		if (pos == std::string::npos)
			pos = 0;
		std::string satName = filename.substr(pos, filename.size());
		printf("c name: %s \n", satName.c_str());
		printf("c nVars: %12d nCls:%12d IO time:%3.2fs\n", ig.nVars(),
				clauses.size(), cpuTime() - initial_time);
	}
	return res;
}

SatInstance Preprocessor::getInstance(std::string const & filename) {
	Timer preprocessTime;
	if (!readInstance(filename) || (elim && !eliminate())
			|| !removeRedundant(true))
		return SatInstance();

	// Force full cleanup (this is safe and desirable since it only happens once):
	garbageCollect(true);
	vec<lbool> model(ig.nVars(), lbool::Undef());
	vec<bool> isDecisionVar(ig.nVars(), false);
	for (int i = 0; i < ig.nVars(); ++i) {
		Var const v = i;
		model[v] = ig.value(v);
		isDecisionVar[v] = !eliminated[v] && ig.value(v).isUndef();
		eliminated_vars -= eliminated[v] && !ig.value(v).isUndef();
	}

	SatInstance res(std::move(model), std::move(isDecisionVar),
			std::move(clauses), std::move(ca), std::move(elimDb),
			std::move(drat));
	assert(res.isClean());
	printf(
			"c #########################  after Preprocessor  #######################\n");
	printf("c nVars: %12d nCls:%12d    time:%3.2fs\n",
			ig.nVars() - ig.nAssigns() - eliminated_vars, res.clauses.size(),
			preprocessTime.getPassedTime());
	printf(
			"c ######################################################################\n\n");
	return res;
}

typename Preprocessor::Var Preprocessor::newVar() {
	Var v = ig.nVars();
	branch.newVar(false);
	ig.newVar();
	propEngine.newVar();

	eliminated.push((char) false);

	if (elim) {
		n_occ.push(0);
		n_occ.push(0);
		occurs.init(v);
		touched.push(0);
		elim_heap.insert(v);
	}
	return v;
}

bool Preprocessor::addClause(vec<Lit>& ps, bool const initial) {
#ifndef NDEBUG
	for (int i = 0; i < ps.size(); i++)
		assert(!isEliminated(ps[i].var()));
#endif
	assert(ig.decisionLevel() == 0);
	if (!isOk())
		return false;

	if (initial && drat.isActive())
		ps.copyTo(add_oc);

	int const initSize = ps.size();
	if (ig.removeRedundant(ps))
		return true;

	if (initSize != ps.size() || !initial) {
		drat.addClause(ps);
		if (initial)
			drat.removeClause(add_oc);
	}

	if (ps.size() == 0) {
		return setOk(false);
	} else if (ps.size() == 1) {
		propEngine.uncheckedEnqueue(branch, ps[0]);
		setOk(propEngine.propagate(branch) == CRef_Undef);
		return true;
	}
	CRef cr = ca.alloc(ps, false);
	clauses.push(cr);
	propEngine.attachClause(cr);

	if (elim && cr != CRef_Undef) {
		const Clause& c = ca[cr];

		// NOTE: the clause is added to the queue immediately and then
		// again during 'gatherTouchedClauses()'. If nothing happens
		// in between, it will only be checked once. Otherwise, it may
		// be checked twice unnecessarily. This is an unfortunate
		// consequence of how backward subsumption is used to mimic
		// forward subsumption.
		subsumption_queue.insert(cr);
		for (int i = 0; i < c.size(); i++) {
			Var const v = c[i].var();
			occurs[v].push(cr);
			n_occ[toInt(c[i])]++;
			touched[v] = 1;
			n_touched++;
			if (elim_heap.inHeap(v))
				elim_heap.increase(v);
		}
	}

	return true;
}

void Preprocessor::removeClause(CRef const cr) {
	Clause& c = ca[cr];
	assert(c.mark() != 1);
	if (elim)
		for (int i = 0; i < c.size(); i++) {
			n_occ[toInt(c[i])]--;
			updateElimHeap(c[i].var());
			occurs.smudge(c[i].var());
		}

	propEngine.detachClause(cr);
	if (ig.locked(c)) {
		Lit implied =
				c.size() != 2 ? c[0] : (ig.value(c[0]).isTrue() ? c[0] : c[1]);
		ig.reason((implied.var())) = CRef_Undef;
	}
	ca.remove(cr);
}

bool Preprocessor::strengthenClause(CRef cr, Lit l) {
	Clause& c = ca[cr];
	assert(ig.decisionLevel() == 0);
	assert(elim);

	// FIX: this is too inefficient but would be nice to have (properly implemented)
	// if (!find(subsumption_queue, &c))
	subsumption_queue.insert(cr);

	drat.addClauseExcludeLit(c, l);

	if (c.size() == 2) {
		removeClause(cr);
		c.strengthen(l);
		return propEngine.enqueue(branch, c[0])
				&& propEngine.propagate(branch) == CRef_Undef;
	} else {
		drat.removeClause(c);

		propEngine.detachClause(cr, true);
		c.strengthen(l);
		propEngine.attachClause(cr);
		remove(occurs[l.var()], cr);
		n_occ[toInt(l)]--;
		updateElimHeap(l.var());
		return true;
	}

}

// Returns FALSE if clause is always satisfied ('out_clause' should not be used).

bool Preprocessor::merge(const Clause& _ps, const Clause& _qs, Var v,
		vec<Lit>& out_clause) {
	out_clause.clear();

	bool ps_smallest = _ps.size() < _qs.size();
	const Clause& ps = ps_smallest ? _qs : _ps;
	const Clause& qs = ps_smallest ? _ps : _qs;

	for (int i = 0; i < qs.size(); i++) {
		if (qs[i].var() != v) {
			for (int j = 0; j < ps.size(); j++)
				if (ps[j].var() == qs[i].var())
					if (ps[j] == ~qs[i])
						return false;
					else
						goto next;
			out_clause.push(qs[i]);
		}
		next: ;
	}

	for (int i = 0; i < ps.size(); i++)
		if (ps[i].var() != v)
			out_clause.push(ps[i]);

	return true;
}

// Returns FALSE if clause is always satisfied.

bool Preprocessor::merge(const Clause& _ps, const Clause& _qs, Var v,
		int& size) {
	bool ps_smallest = _ps.size() < _qs.size();
	const Clause& ps = ps_smallest ? _qs : _ps;
	const Clause& qs = ps_smallest ? _ps : _qs;
	const Lit* __ps = (const Lit*) ps;
	const Lit* __qs = (const Lit*) qs;

	size = ps.size() - 1;

	for (int i = 0; i < qs.size(); i++) {
		if (__qs[i].var() != v) {
			for (int j = 0; j < ps.size(); j++)
				if (__ps[j].var() == __qs[i].var())
					if (__ps[j] == ~__qs[i])
						return false;
					else
						goto next;
			size++;
		}
		next: ;
	}

	return true;
}

void Preprocessor::gatherTouchedClauses() {
	if (n_touched == 0)
		return;

	int i, j;
	for (i = j = 0; i < subsumption_queue.size(); i++) {
		Clause & c = ca[subsumption_queue[i]];
		if (c.mark() == 0)
			c.mark(2);
	}

	for (i = 0; i < touched.size(); i++)
		if (touched[i]) {
			const vec<CRef>& cs = occurs.lookup(i);
			for (j = 0; j < cs.size(); j++) {
				Clause & c = ca[cs[j]];
				if (c.mark() == 0) {
					subsumption_queue.insert(cs[j]);
					c.mark(2);
				}
			}
			touched[i] = 0;
		}

	for (i = 0; i < subsumption_queue.size(); i++) {
		Clause & c = ca[subsumption_queue[i]];
		if (c.mark() == 2)
			c.mark(0);
	}

	n_touched = 0;
}

// Backward subsumption + backward subsumption resolution

bool Preprocessor::backwardSubsumptionCheck(bool verbose) {
	int cnt = 0;
	int subsumed = 0;
	int deleted_literals = 0;
	assert(ig.decisionLevel() == 0);

	while (subsumption_queue.size() > 0 || bwdsub_assigns < ig.nAssigns()) {

		// Empty subsumption queue and return immediately on user-interrupt:
		if (hasInterrupt()) {
			subsumption_queue.clear();
			bwdsub_assigns = ig.nAssigns();
			break;
		}

		// Check top-level assignments by creating a dummy clause and placing it in the queue:
		if (subsumption_queue.size() == 0 && bwdsub_assigns < ig.nAssigns()) {
			Lit l = ig.getTrailLit(bwdsub_assigns++);
			ca[bwdsub_tmpunit][0] = l;
			ca[bwdsub_tmpunit].calcAbstraction();
			subsumption_queue.insert(bwdsub_tmpunit);
		}

		CRef cr = subsumption_queue.peek();
		subsumption_queue.pop();
		Clause& c = ca[cr];

		if (c.mark())
			continue;

		if (verbose && verb >= 2 && cnt++ % 1000 == 0)
			printf(
					"c subsumption left: %10d (%10d subsumed, %10d deleted literals)\r",
					subsumption_queue.size(), subsumed, deleted_literals);

		assert(c.size() > 1 || ig.value(c[0]).isTrue()); // Unit-clauses should have been propagated before this point.

				// Find best variable to scan:
		Var best = c[0].var();
		for (int i = 1; i < c.size(); i++)
			if (occurs[c[i].var()].size() < occurs[best].size())
				best = c[i].var();

		// Search all candidates:
		vec<CRef>& _cs = occurs.lookup(best);
		CRef* cs = (CRef*) _cs;

		for (int j = 0; j < _cs.size(); j++)
			if (c.mark())
				break;
			else if (!ca[cs[j]].mark() && cs[j] != cr
					&& (subsumption_lim == -1
							|| ca[cs[j]].size() < subsumption_lim)) {
				Lit l = c.subsumes(ca[cs[j]]);

				if (l == Lit::Undef())
					subsumed++, removeClause(cs[j]);
				else if (l != Lit::Error()) {
					deleted_literals++;

					if (!strengthenClause(cs[j], ~l))
						return false;

					// Did current candidate get deleted from cs? Then check candidate at index j again:
					if (l.var() == best)
						j--;
				}
			}
	}

	return true;
}

bool Preprocessor::eliminateVar(Var const v) {
	assert(!isEliminated(v));
	assert(ig.value(v).isUndef());

// Split the occurrences into positive and negative:
//
	const vec<CRef>& cls = occurs.lookup(v);
	vec<CRef> pos, neg;
	for (int i = 0; i < cls.size(); i++)
		(find(ca[cls[i]], Lit(v, false)) ? pos : neg).push(cls[i]);

// Check wether the increase in number of clauses stays within the allowed ('grow'). Moreover, no
// clause must exceed the limit on the maximal clause size (if it is set):
//
	int cnt = 0;
	int clause_size = 0;

	for (int i = 0; i < pos.size(); i++)
		for (int j = 0; j < neg.size(); j++)
			if (merge(ca[pos[i]], ca[neg[j]], v, clause_size)
					&& (++cnt > cls.size() + grow
							|| (clause_lim != -1 && clause_size > clause_lim)))
				return true;

// Delete and store old clauses:
	eliminated[v] = true;
	branch.setDecisionVar(v, false);
	eliminated_vars++;

	if (pos.size() > neg.size()) {
		for (int i = 0; i < neg.size(); i++)
			elimDb.addElimClause(v, ca[neg[i]]);
		elimDb.addElimUnit(Lit(v, false));
	} else {
		for (int i = 0; i < pos.size(); i++)
			elimDb.addElimClause(v, ca[pos[i]]);
		elimDb.addElimUnit(Lit(v, true));
	}

// Produce clauses in cross product:
	vec<Lit>& resolvent = add_tmp;
	for (int i = 0; i < pos.size(); i++)
		for (int j = 0; j < neg.size(); j++)
			if (merge(ca[pos[i]], ca[neg[j]], v, resolvent)
					&& !addClause(resolvent))
				return false;

	for (int i = 0; i < cls.size(); i++)
		removeClause(cls[i]);

// Free occurs list for this variable:
	occurs[v].clear(true);

// Free watchers lists for this variable, if possible:
	propEngine.removeVar(v);

	return backwardSubsumptionCheck();
}

bool Preprocessor::substitute(Var v, Lit x) {
	assert(!isEliminated(v));
	assert(ig.value(v).isUndef());

	if (!isOk())
		return false;

	eliminated[v] = true;
	branch.setDecisionVar(v, false);
	const vec<CRef>& cls = occurs.lookup(v);

	vec<Lit>& subst_clause = add_tmp;
	for (int i = 0; i < cls.size(); i++) {
		Clause& c = ca[cls[i]];

		subst_clause.clear();
		for (int j = 0; j < c.size(); j++) {
			Lit p = c[j];
			subst_clause.push(p.var() == v ? x ^ p.sign() : p);
		}

		if (!addClause(subst_clause))
			return setOk(false);

		removeClause(cls[i]);
	}

	return true;
}

int Preprocessor::nFreeVars() const {
	return (int) branch.nDecVars()
			- (ig.decisionLevel() == 0 ? ig.nAssigns() : ig.levelEnd(0));
}
bool Preprocessor::removeRedundant(bool const removeFalseLits) {
	if (!isOk() || propEngine.propagate(branch) != ClauseAllocator::npos())
		return setOk(false);
	if (removeFalseLits) {
		int i = 0, j = 0;
		if (!drat.isActive()) {
			for (; i < clauses.size(); ++i) {
				Clause & c = ca[clauses[i]];
				if (c.mark() == 0) {
					if (!ig.removeSetLits(c)) {
						assert(c.size() != 0);
						if (c.size() > 1)
							clauses[j++] = clauses[i];
					} else
						removeClause(clauses[i]);
				}
			}
		} else {
			vec<Lit> oldC;
			for (; i < clauses.size(); ++i) {
				Clause & c = ca[clauses[i]];
				if (c.mark() == 0) {
					lbool state = ig.checkLits(c);
					if (!state.isTrue()) {
						if (state.isFalse()) {
							oldC.copyFrom(c);
							ig.removeSetLits(c);
							assert(c.size() < oldC.size());
							drat.addClause(c);
							drat.removeClause(oldC);
							assert(c.size() > 0);
							if(c.size() < 2)
								continue;
						}
						clauses[j++] = clauses[i];
					}
				}
			}
		}
		clauses.shrink(i - j);

	} else if (nRedundantAssignedLits < ig.nAssigns()) {

		int i = 0, j = 0;
		for (; i < clauses.size(); ++i) {
			Clause const & c = ca[clauses[i]];
			if (c.mark() == 0) {
				if (c.mark() == 0 && !ig.satisfied(c))
					clauses[j++] = clauses[i];
				else
					removeClause(clauses[i]);
			}
		}
		clauses.shrink(i - j);
	}
	nRedundantAssignedLits = ig.nAssigns();
	return true;
}

bool Preprocessor::eliminateVarMinimal(Var const v) {
	assert(!isEliminated(v));
	assert(ig.value(v).isUndef());

// Split the occurrences into positive and negative:
//
	const vec<CRef>& cls = occurs.lookup(v);
	vec<CRef> pos, neg;
	for (int i = 0; i < cls.size(); i++)
		(find(ca[cls[i]], Lit(v, false)) ? pos : neg).push(cls[i]);

// Delete and store old clauses:
	eliminated[v] = true;
	branch.setDecisionVar(v, false);
	eliminated_vars++;

	if (pos.size() > neg.size()) {
		for (int i = 0; i < neg.size(); i++)
			elimDb.addElimClause(v, ca[neg[i]]);
		elimDb.addElimUnit(Lit(v, false));
	} else {
		for (int i = 0; i < pos.size(); i++)
			elimDb.addElimClause(v, ca[pos[i]]);
		elimDb.addElimUnit(Lit(v, true));
	}

// Produce clauses in cross product:
	vec<Lit>& resolvent = add_tmp;
	for (int i = 0; i < pos.size(); i++)
		for (int j = 0; j < neg.size(); j++)
			if (merge(ca[pos[i]], ca[neg[j]], v, resolvent)
					&& !addClause(resolvent))
				return false;

	for (int i = 0; i < cls.size(); i++)
		removeClause(cls[i]);

// Free occurs list for this variable:
	occurs[v].clear(true);

// Free watchers lists for this variable, if possible:
	propEngine.removeVar(v);

	return true;
}

bool Preprocessor::eliminateMinimalNiver() {
	if (!removeRedundant())
		return setOk(false);

	while (!elim_heap.empty()) {
		Var const elim = elim_heap.removeMin();

		uint32_t const npos = n_occ[Lit(elim, false).toInt()], nneg = n_occ[Lit(
				elim, true).toInt()];

		if (npos > 1 && nneg > 1)
			break;
		if (hasInterrupt())
			break;

		if (isEliminated(elim) || !ig.value(elim).isUndef())
			continue;

		// At this point, the variable may have been set by assymetric branching, so check it
		// again. Also, don't eliminate frozen variables:
		if (ig.value(elim).isUndef() && !eliminateVarMinimal(elim)) {
			setOk(false);
			break;
		}

		garbageCollect();
	}
	return true;
}

// The technique and code are by the courtesy of the GlueMiniSat team. Thank you!
// It helps solving certain types of huge problems tremendously.

bool Preprocessor::eliminate() {
	assert(elim);
	int iter = 0;
	int n_cls, n_cls_init, n_vars;

	if (ig.nVars() == 0)
		goto cleanup;
// User disabling preprocessing.

	n_cls_init = nClauses();
	if (nClauses() > 480000) {
		if (verb > 0)
			std::cout
					<< "c using Niver on variables with less occurances\nc too many clauses in instance for SatElite!\n";
		eliminateMinimalNiver();
		goto cleanup;
	}
	eliminate_();  // The first, usual variable elimination of MiniSat.
	if (!isOk())
		goto cleanup;

	n_cls = nClauses();
	n_vars = nFreeVars();

	if ((double) n_cls / n_vars >= 10 || n_vars < 10000) {
		goto cleanup;
	}

	grow = grow ? grow * 2 : 8;
	for (; grow < 10000; grow *= 2) {
// Rebuild elimination variable heap.
		for (int i = 0; i < clauses.size(); i++) {
			const Clause& c = ca[clauses[i]];
			for (int j = 0; j < c.size(); j++)
				if (!elim_heap.inHeap(c[j].var()))
					elim_heap.insert(c[j].var());
				else
					elim_heap.update(c[j].var());
		}

		int n_cls_last = nClauses();
		int n_vars_last = nFreeVars();

		eliminate_();
		if (!isOk() || n_vars_last == nFreeVars())
			break;
		iter++;

		int n_cls_now = nClauses();
		int n_vars_now = nFreeVars();

		double cl_inc_rate = (double) n_cls_now / n_cls_last;
		double var_dec_rate = (double) n_vars_last / n_vars_now;

		if (n_cls_now > n_cls_init || cl_inc_rate > var_dec_rate)
			break;
	}

	cleanup:

	if (!isOk()) {
		drat.addEmptyClause();
		drat.flush();
	}
	return isOk();
}

void Preprocessor::clear() {
	touched.clear(true);
	occurs.clear(true);
	n_occ.clear(true);
	elim_heap.clear(true);
	subsumption_queue.clear(true);
	// FXIME: implement and use clear functions
//   ca.~ClauseAllocator();
//   branch.~Branch<ClauseAllocator>();
//   propEngine.~MinisatPropagate<ClauseAllocator, Branch<ClauseAllocator>>();
//   ig.~ImplicationGraph<ClauseAllocator>();
}

bool Preprocessor::eliminate_() {
	assert(elim);
	if (!removeRedundant())
		return setOk(false);

	int trail_size_last = ig.nAssigns();

// Main simplification loop:
//
	while (n_touched > 0 || bwdsub_assigns < ig.nAssigns()
			|| elim_heap.size() > 0) {

		gatherTouchedClauses();
// printf("  ## (time = %6.2f s) BWD-SUB: queue = %d, trail = %d\n", cpuTime(), subsumption_queue.size(), trail.size() - bwdsub_assigns);
		if ((subsumption_queue.size() > 0 || bwdsub_assigns < ig.nAssigns())
				&& !backwardSubsumptionCheck(true)) {
			setOk(false);
			goto cleanup;
		}

// Empty elim_heap and return immediately on user-interrupt:
		if (hasInterrupt()) {
			assert(bwdsub_assigns == ig.nAssigns());
			assert(subsumption_queue.size() == 0);
			assert(n_touched == 0);
			elim_heap.clear();
			goto cleanup;
		}

// printf("  ## (time = %6.2f s) ELIM: vars = %d\n", cpuTime(), elim_heap.size());
		for (int cnt = 0; !elim_heap.empty(); cnt++) {
			Var elim = elim_heap.removeMin();

			if (hasInterrupt())
				break;

			if (isEliminated(elim) || !ig.value(elim).isUndef())
				continue;

			if (verb >= 2 && cnt % 100 == 0)
				printf("c elimination left: %10d\r", elim_heap.size());
			// At this point, the variable may have been set by assymetric branching, so check it
			// again. Also, don't eliminate frozen variables:
			if (ig.value(elim).isUndef() && !eliminateVar(elim)) {
				setOk(false);
				goto cleanup;
			}

			garbageCollect();
		}

		assert(subsumption_queue.size() == 0);
	}
	cleanup:
// To get an accurate number of clauses.
	if (trail_size_last != ig.nAssigns())
		removeRedundant();
	garbageCollect();

	return isOk();
}

//=================================================================================================
// Garbage Collection methods:

void Preprocessor::relocAll(ClauseAllocator& to, bool const finalGarbage) {
// All occurs lists:
//
	for (int i = 0; i < clauses.size(); ++i)
		ca.reloc(clauses[i], to);

	if (!finalGarbage) {
		propEngine.relocAll(to);

		occurs.cleanAll();
		for (int i = 0; i < nVars(); i++) {
			vec<CRef>& cs = occurs[i];
			for (int j = 0; j < cs.size(); j++)
				ca.reloc(cs[j], to);
		}

// Subsumption queue:
//
		for (int i = 0; i < subsumption_queue.size(); i++)
			ca.reloc(subsumption_queue[i], to);

// Temporary clause:
//
		ca.reloc(bwdsub_tmpunit, to);
	}
}

void Preprocessor::garbageCollect(bool const finalGarbage) {
	if (!finalGarbage && ca.wasted() < ca.size() * simp_garbage_frac)
		return;
// Initialize the next region to a size corresponding to the estimated utilization degree. This
// is not precise but should avoid some unnecessary reallocations for the new region:
	ClauseAllocator to(ca.size() - ca.wasted());
	to.extra_clause_field = !finalGarbage && ca.extra_clause_field;	// NOTE: this is important to keep (or lose) the extra fields.
	relocAll(to, finalGarbage);
	to.moveTo(ca);
}
}
