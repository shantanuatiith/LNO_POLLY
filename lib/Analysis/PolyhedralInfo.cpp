//===--------- PolyhedralInfo.cpp  - Create Scops from LLVM IR-------------===//
///
///                     The LLVM Compiler Infrastructure
///
/// This file is distributed under the University of Illinois Open Source
/// License. See LICENSE.TXT for details.
///
//===----------------------------------------------------------------------===//
///
/// An interface to the Polyhedral analysis engine(Polly) of LLVM.
///
/// This pass provides an interface to the polyhedral analysis performed by
/// Polly.
///
/// This interface provides basic interface like isParallel, isVectorizable
/// that can be used in LLVM transformation passes.
///
/// Work in progress, this file is subject to change.
//===----------------------------------------------------------------------===//

#include "polly/PolyhedralInfo.h"
#include "polly/DependenceInfo.h"
#include "polly/LinkAllPasses.h"
#include "polly/Options.h"
#include "polly/ScopInfo.h"
#include "polly/Support/GICHelper.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Support/Debug.h"
#include <isl/map.h>
#include <isl/union_map.h>

using namespace llvm;
using namespace polly;

#define DEBUG_TYPE "polyhedral-info"

STATISTIC(LoopsAnalyzedP, "Number of Loops Analyzed for Parallelism");
STATISTIC(LoopsAnalyzedV, "Number of Loops Analyzed for Vectorization");
STATISTIC(ParallelLoops, "Number of Parallel Loops");
STATISTIC(VectorizableLoops, "Number of Vectorizable Loops");

static cl::opt<bool> CheckParallel("polly-check-parallel",
                                   cl::desc("Check for parallel loops"),
                                   cl::Hidden, cl::init(false), cl::ZeroOrMore,
                                   cl::cat(PollyCategory));

static cl::opt<bool> CheckVectorizable("polly-check-vectorizable",
                                       cl::desc("Check for vectorizable loops"),
                                       cl::Hidden, cl::init(false),
                                       cl::ZeroOrMore, cl::cat(PollyCategory));

static cl::opt<bool>
    CheckDependences("polly-check-dependences",
                     cl::desc("Polly instruction wise dependences"), cl::Hidden,
                     cl::init(false), cl::ZeroOrMore, cl::cat(PollyCategory));

static cl::opt<bool> PrintDependenceGraph("print-dependence-graph",
			cl::desc("Prints the dependence graph constructed from polly dependence analysis"),
			cl::Hidden, cl::init(false), cl::ZeroOrMore, cl::cat(PollyCategory));

// DependenceDirectionVector methods

DependenceDirectionVector::DependenceDirectionVector(Instruction *Source,
                               Instruction *Destination,
                               bool PossiblyLoopIndependent,
                               unsigned CommonLevels)
    : Dependence(Source, Destination), Levels(CommonLevels),
      LoopIndependent(PossiblyLoopIndependent) {
  Consistent = true;
  if (CommonLevels)
    DV = make_unique<DVEntry[]>(CommonLevels);
}

// The rest are simple getters that hide the implementation.

// getDirection - Returns the direction associated with a particular level.
unsigned DependenceDirectionVector::getDirection(unsigned Level) const {
  assert(Level <= Levels && "Level out of range");
  return DV[Level].Direction;
}


// Returns the distance (or NULL) associated with a particular level.
const llvm::SCEV *DependenceDirectionVector::getDistance(unsigned Level) const {
  assert(0 < Level && Level <= Levels && "Level out of range");
  return DV[Level - 1].Distance;
}

// Temporary fix untill SCEV version is fixed.
int DependenceDirectionVector::getDepDistance(unsigned Level) const {
  assert(Level < Levels && "Level out of range");
  return DV[Level].DepDistance;
}


// Returns true if a particular level is scalar; that is,
// if no subscript in the source or destination mention the induction
// variable associated with the loop at this level.
bool DependenceDirectionVector::isScalar(unsigned Level) const {
  assert(0 < Level && Level <= Levels && "Level out of range");
  return DV[Level - 1].Scalar;
}


// Returns true if peeling the first iteration from this loop
// will break this dependence.
bool DependenceDirectionVector::isPeelFirst(unsigned Level) const {
  assert(0 < Level && Level <= Levels && "Level out of range");
  return DV[Level - 1].PeelFirst;
}


// Returns true if peeling the last iteration from this loop
// will break this dependence.
bool DependenceDirectionVector::isPeelLast(unsigned Level) const {
  assert(0 < Level && Level <= Levels && "Level out of range");
  return DV[Level - 1].PeelLast;
}


// Returns true if splitting this loop will break the dependence.
bool DependenceDirectionVector::isSplitable(unsigned Level) const {
  assert(0 < Level && Level <= Levels && "Level out of range");
  return DV[Level - 1].Splitable;
}

void PolyhedralInfo::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.addRequiredTransitive<DependenceInfoWrapperPass>();
  AU.addRequired<LoopInfoWrapperPass>();
  AU.addRequiredTransitive<ScopInfoWrapperPass>();
  AU.setPreservesAll();
}

bool PolyhedralInfo::runOnFunction(Function &F) {
  DI = &getAnalysis<DependenceInfoWrapperPass>();
  SI = getAnalysis<ScopInfoWrapperPass>().getSI();
  return false;
}

void PolyhedralInfo::print(raw_ostream &OS, const Module *) const {
  auto &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
  int VF = 2;
   for (auto *TopLevelLoop : LI) {
     for (auto *L : depth_first(TopLevelLoop)) {
       OS.indent(2) << L->getHeader()->getName() << ":\t";
       if (CheckParallel && isParallel(L))
         OS << "Loop is parallel.\n";
       else if (CheckParallel)
         OS << "Loop is not parallel.\n";
       if (CheckVectorizable && isVectorizable(L, &VF)) {
         OS << "Loop is vectorizable with max VF = ";
         if (VF == INT_MAX)
           OS << "infinite";
         else
           OS << VF;
         OS << "\n";
       } else if (CheckVectorizable)
         OS << "Loop is not vectorizable.\n";
     }
   }

	 // Printing dependence graph
	 if(PrintDependenceGraph) {
		 std::vector<DependenceDirectionVector *> *DDV;
		 OS << "\n\n";
	   for (auto *TopLevelLoop : LI) {
		   DDV = getDDVs(getScopContainingLoop(TopLevelLoop));
			 for(unsigned i = 0; i < DDV->size(); i++) {
			   for(unsigned level = 0; level < DDV->at(i)->getLevels(); level++) {
				   if(DDV->at(i)->getDepDistance(level) != 0)
					    OS.indent(4) << "\n" << *DDV->at(i)->getSrc() << "  ----(distance " << DDV->at(i)->getDepDistance(level)
									 << ")----(level " << level << ")----(direction " << DDV->at(i)->getDirection(level) <<")----"
									 << *DDV->at(i)->getDst() << "\n";
					 }
			 }
		   OS << "\n\n";
		 }

	 }
	 
}
/// @brief Get the minimum dependence distance from the quasi-affine expression
///        @p Dist and update it to @p User
static isl_stat getMinDependenceDistance(__isl_take isl_set *Set,
                                         __isl_take isl_aff *Dist, void *User) {
  int *MinDepDistance = static_cast<int *>(User);

  isl_val *Distance = isl_aff_get_constant_val(Dist);
  int Flag = isl_val_cmp_si(Distance, *MinDepDistance);
  if (Flag == -1)
    *MinDepDistance = isl_val_get_num_si(Distance);

  isl_val_free(Distance);
  isl_set_free(Set);
  isl_aff_free(Dist);
  return isl_stat_ok;
}

static isl_stat getDependenceDistance(__isl_take isl_set *Set,
                                      __isl_take isl_aff *Dist, void *User) {
  int *DepDistance = static_cast<int *>(User);
  isl_bool isConst = isl_aff_is_cst(Dist);
  if(isConst == isl_bool_true) {
    isl_val *Distance = isl_aff_get_constant_val(Dist);
    *DepDistance = isl_val_get_num_si(Distance);
    isl_val_free(Distance);

  }else if(isConst == isl_bool_false)
    *DepDistance = INT_MAX;    

  isl_set_free(Set);
  isl_aff_free(Dist);
  return isl_stat_ok;
}

bool PolyhedralInfo::checkParallel(Loop *L, isl_pw_aff **MinDepDistPtr) const {
  bool IsParallel;
  const Scop *S = getScopContainingLoop(L);
  if (!S)
    return false;
  const Dependences &D =
      DI->getDependences(const_cast<Scop *>(S), Dependences::AL_Access);
  if (!D.hasValidDependences())
    return false;
  DEBUG(dbgs() << "Loop :\t" << L->getHeader()->getName() << ":\n");

  isl_union_map *Deps =
      D.getDependences(Dependences::TYPE_RAW | Dependences::TYPE_WAW |
                       Dependences::TYPE_WAR | Dependences::TYPE_RED);
  DEBUG(dbgs() << "Dependences :\t" << stringFromIslObj(Deps) << "\n");

  isl_union_map *Schedule = getScheduleForLoop(S, L);
  DEBUG(dbgs() << "Schedule: \t" << stringFromIslObj(Schedule) << "\n");

  IsParallel = D.isParallel(Schedule, Deps, MinDepDistPtr);
  isl_union_map_free(Schedule);
  return IsParallel;
}

bool PolyhedralInfo::isParallel(Loop *L) const { return checkParallel(L); }

const Scop *PolyhedralInfo::getScopContainingLoop(Loop *L) const {
  assert((SI) && "ScopInfoWrapperPass is required by PolyhedralInfo pass!\n");
  for (auto &It : *SI) {
    Region *R = It.first;
    if (R->contains(L))
      return It.second.get();
  }
  return nullptr;
}

bool PolyhedralInfo::isVectorizable(Loop *L, int *VF) const {
  ++LoopsAnalyzedV;
  isl_pw_aff *MinDistancePtr = nullptr;
  bool IsVectorizable = checkParallel(L, &MinDistancePtr);

  if (IsVectorizable || !MinDistancePtr) {
    // Trivially Vectorizable. Set VF to infinity
    *VF = IsVectorizable ? INT_MAX : *VF;
    isl_pw_aff_free(MinDistancePtr);
    if (IsVectorizable)
      ++VectorizableLoops;
    return IsVectorizable;
  }

  // Currently handle only constant distances
  if (!isl_pw_aff_is_empty(MinDistancePtr) &&
      !isl_pw_aff_involves_nan(MinDistancePtr) &&
      isl_pw_aff_is_cst(MinDistancePtr)) {
    int MinDepDistance = INT_MAX;

    isl_pw_aff_foreach_piece(MinDistancePtr, getMinDependenceDistance,
                             &MinDepDistance);
    IsVectorizable = (*VF <= MinDepDistance) ? true : false;

    // If minimum dependence distance is greater than vectorization factor @p VF
    // provided by user, update VF to maximum permissible value, which is the
    // minimum dependence distance.
    *VF = IsVectorizable ? MinDepDistance : *VF;
    DEBUG(dbgs() << "Min Dep Distance(constant):\t" << MinDepDistance << "\n");
  }

  isl_pw_aff_free(MinDistancePtr);
  if (IsVectorizable)
    ++VectorizableLoops;
  return IsVectorizable;
}


//  Given a Loop and the containing SCoP, we compute the partial schedule
//  by taking union of individual schedules of each ScopStmt within the loop
//  and projecting out the inner dimensions from the range of the schedule.
//   for (i = 0; i < n; i++)
//      for (j = 0; j < n; j++)
//        A[j] = 1;  //Stmt
//
//  The original schedule will be
//    Stmt[i0, i1] -> [i0, i1]
//  The schedule for the outer loop will be
//    Stmt[i0, i1] -> [i0]
//  The schedule for the inner loop will be
//    Stmt[i0, i1] -> [i0, i1]
__isl_give isl_union_map *PolyhedralInfo::getScheduleForLoop(const Scop *S,
                                                             Loop *L) const {
  isl_union_map *Schedule = isl_union_map_empty(S->getParamSpace());
  int CurrDim = S->getRelativeLoopDepth(L);
  DEBUG(dbgs() << "Relative loop depth:\t" << CurrDim << "\n");
  assert(CurrDim >= 0 && "Loop in region should have at least depth one");

  for (auto *BB : L->blocks()) {
    auto *SS = S->getStmtFor(BB);
    if (!SS)
      continue;

    unsigned int MaxDim = SS->getNumIterators();
    DEBUG(dbgs() << "Maximum depth of Stmt:\t" << MaxDim << "\n");
    auto *ScheduleMap = SS->getSchedule();
    assert(ScheduleMap &&
           "Schedules that contain extension nodes require special handling.");

    ScheduleMap = isl_map_project_out(ScheduleMap, isl_dim_out, CurrDim + 1,
                                      MaxDim - CurrDim - 1);
    ScheduleMap =
        isl_map_set_tuple_id(ScheduleMap, isl_dim_in, SS->getDomainId());
    Schedule =
        isl_union_map_union(Schedule, isl_union_map_from_map(ScheduleMap));
  }

  Schedule = isl_union_map_coalesce(Schedule);
  return Schedule;
}

static isl_stat getMapVector(__isl_take isl_map *map, void *user) {
  std::vector<isl_map* > *depMapVector = (std::vector<isl_map *> *)user;
	depMapVector->push_back(map);
  return isl_stat_ok;
}

std::vector<DependenceDirectionVector *> *PolyhedralInfo::getDDVs(const Scop *S) const{

	const Dependences &D =
      DI->getDependences(const_cast<Scop*>(S), Dependences::AL_Access);

  isl_union_map *Deps = D.getWrappedDependences(Dependences::TYPE_RAW
		  | Dependences::TYPE_WAR | Dependences::TYPE_WAW);

  DEBUG(dbgs() << isl_union_map_to_str(Deps) );
  
  std::vector<DependenceDirectionVector*> *Dependence =
				 new std::vector<DependenceDirectionVector*>();

  //isl_union_map_foreach_map(Deps, getMapDDV, (void*)Dependence);
	std::vector<isl_map *> *depMapVector = new std::vector<isl_map *>();
  isl_union_map_foreach_map(Deps, getMapVector, (void*)depMapVector);
	for(auto *map : *depMapVector) {
	  getMapDDV(S, map, Dependence);
		map = nullptr; 
	}

	// clean up of manual union map
	free(depMapVector);
	

  unsigned numMaps = isl_union_map_n_map(Deps);
  //assert(numMaps == Dependences.size() && 
  //          "All dependences are not captured to DDV object");
  isl_union_map_free(Deps);
  DEBUG(dbgs() << "No of Dependences : " << Dependence->size() << "\n");
/*
	for(unsigned i=0; i<Dependences.size(); i++) {
     if(i>0)
      Dependences[i]->setNextPredecessor(Dependences[i-1]);
    if(i<numMaps-1)
      Dependences[i]->setNextSuccessor(Dependences[i+1]);
  }
*/
  if(Dependence->size() != 0)
    return Dependence;
  else
    return nullptr;
}


void PolyhedralInfo::getMapDDV(const Scop *Scp, isl_map *map,
				std::vector<llvm::DependenceDirectionVector *> *DDV) const {
	const Dependences &D =
      DI->getDependences(const_cast<Scop*>(Scp), Dependences::AL_Access);
  //get the pointer to the dependences vector
  //std::vector<DependenceDirectionVector *> *Dependences =
  //                 (std::vector<DependenceDirectionVector *> *)user;
  isl_space *S = nullptr;
  // Source Instruction
  S = isl_space_domain(isl_map_get_space(map));
  if (!isl_space_is_wrapping(S)) {
    DEBUG(dbgs() << "Domain space is not wrapping\n");
    dbgs() << "Domain space is not wrapping\n";
    isl_space_free(S);
    isl_map_free(map);
    return;
  }
  S = isl_space_unwrap(S);

  isl_id *id_Source = isl_space_get_tuple_id(S, isl_dim_out);
  MemoryAccess *MA1 = (MemoryAccess *)(isl_id_get_user(id_Source));
  isl_id_free(id_Source);
  isl_space_free(S);

  // Sink Instruction
  S = isl_space_range(isl_map_get_space(map));
  if (!isl_space_is_wrapping(S)) {
    DEBUG(dbgs() << "Range space is not wrapping\n");
    dbgs() << "Range space is not wrapping\n";
    isl_space_free(S);
    isl_map_free(map);
    return;
  }

  S = isl_space_unwrap(S);

  isl_id *id_Sink = isl_space_get_tuple_id(S, isl_dim_out);

  MemoryAccess *MA2 = (MemoryAccess *)(isl_id_get_user(id_Sink));
  isl_id_free(id_Sink);
  isl_space_free(S);

  // insert dependence information into map
  Instruction *Src = MA1->getAccessInstruction();
  Instruction *Dst = MA2->getAccessInstruction();

  // get number of levels in the DDV
  unsigned domainDims = isl_map_dim(map, isl_dim_in);
  unsigned rangeDims = isl_map_dim(map, isl_dim_out);
  // first get number of possible levels
  unsigned Levels = (domainDims>rangeDims)?rangeDims:domainDims;
  for(unsigned i = 0; i < Levels; i++) {
    // check that the loops are actually the same
    isl_id *id_domain = isl_map_get_dim_id(map, isl_dim_in, i);
    isl_id *id_range = isl_map_get_dim_id(map, isl_dim_out, i);
    bool equal = (isl_id_get_user(id_domain) == isl_id_get_user(id_range));
    isl_id_free(id_domain);
    isl_id_free(id_range);
    if(!equal) {
      // different loops - recalculate and exit
      Levels = i+1;
      break;
    }
  }

  //remove extra dimensions in the map
  map = isl_map_remove_dims(map, isl_dim_in, Levels, domainDims-Levels);
  map = isl_map_remove_dims(map, isl_dim_out, Levels, rangeDims-Levels);

  map = isl_map_zip(map);
  isl_set *set = isl_map_domain(map);
  map = isl_set_unwrap(set);

  //FIXME:
  if(isl_space_tuple_is_equal(isl_map_get_space(map), isl_dim_in,
             isl_map_get_space(map), isl_dim_out) == isl_bool_false) {
    isl_map_free(map);
    return;
  }

  //get the final dependence polyhedron
  isl_set *DP = isl_map_deltas(isl_map_copy(map));

  DependenceDirectionVector *ddv =
            new DependenceDirectionVector(Src, Dst, false, Levels);


	bool IsParallel = false;
	unsigned depth = 0;
	int distance = 0;
  auto &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();

  for (auto *TopLevelLoop : LI) {
	  depth = -1;
  	if(Scp->contains(TopLevelLoop) &&
						TopLevelLoop->contains(const_cast<Instruction *> (Src)) &&
						TopLevelLoop->contains(const_cast<Instruction *> (Dst))) {
      for (auto *L : depth_first(TopLevelLoop)) {
				depth++;
				distance = 0;
				isl_pw_aff *islDepDist = nullptr;
        isl_union_map *Schedule = getScheduleForLoop(Scp, L);
			  IsParallel = D.isParallel(Schedule, isl_union_map_from_map(isl_map_copy(map)),
												&islDepDist);
				if(IsParallel)
				  ddv->DV[depth].DepDistance = 0;
				else {
          isl_pw_aff_foreach_piece(islDepDist, getDependenceDistance, &distance);
				  ddv->DV[depth].DepDistance = distance;
				}
			}
		
		}
	}

  for(unsigned i = 0; i < Levels; i++) {
    // Direction calculation //
    // which regions is the polyhedron in?
    bool inGreater = false, inEqual = false, inLesser = false;
    isl_set *subset;
    isl_constraint *cnst;
    isl_space *tempSpace = isl_set_get_space(DP);
    isl_local_space *localSpace = isl_local_space_from_space(tempSpace);
    // check if DP is present in lesser than region
    cnst = isl_constraint_alloc_inequality(isl_local_space_copy(localSpace));
    cnst = isl_constraint_set_coefficient_si(cnst, isl_dim_set, i, 1);
    cnst = isl_constraint_set_constant_si(cnst, -1);
    subset = isl_set_add_constraint(isl_set_copy(DP), cnst);
    inLesser = (isl_set_is_empty(subset) == isl_bool_true);
    isl_set_free(subset);
    // check if DP is present in greater than region
    cnst = isl_constraint_alloc_inequality(isl_local_space_copy(localSpace));
    cnst = isl_constraint_set_coefficient_si(cnst, isl_dim_set, i, -1);
    cnst = isl_constraint_set_constant_si(cnst, -1);
    subset = isl_set_add_constraint(isl_set_copy(DP), cnst);
    inGreater = (isl_set_is_empty(subset) == isl_bool_true);
    isl_set_free(subset);
    // check if DP is present in equal to region
    cnst = isl_constraint_alloc_equality(localSpace);
    cnst = isl_constraint_set_coefficient_si(cnst, isl_dim_set, i, 1);
    subset = isl_set_add_constraint(isl_set_copy(DP), cnst);
    inEqual = (isl_set_is_empty(subset) == isl_bool_true);
    isl_set_free(subset);
    ddv->DV[i].Direction = DependenceDirectionVector::DVEntry::NONE;
    if(inLesser)
      ddv->DV[i].Direction |= DependenceDirectionVector::DVEntry::LT;
    if(inGreater)
      ddv->DV[i].Direction |= DependenceDirectionVector::DVEntry::GT;
    if(inEqual)
      ddv->DV[i].Direction |= DependenceDirectionVector::DVEntry::EQ;
  }

  // add this DDV to the vector
  DDV->push_back(ddv);

  // cleanup
  isl_set_free(DP);
  isl_map_free(map);

  // if we are here, everything went fine
  return;// isl_stat_ok;
}


char PolyhedralInfo::ID = 0;

Pass *polly::createPolyhedralInfoPass() { return new PolyhedralInfo(); }

INITIALIZE_PASS_BEGIN(PolyhedralInfo, "polyhedral-info",
                      "Polly - Interface to polyhedral analysis engine", false,
                      false);
INITIALIZE_PASS_DEPENDENCY(DependenceInfoWrapperPass);
INITIALIZE_PASS_DEPENDENCY(LoopInfoWrapperPass);
INITIALIZE_PASS_DEPENDENCY(ScopInfoWrapperPass);
INITIALIZE_PASS_END(PolyhedralInfo, "polyhedral-info",
                    "Polly - Interface to polyhedral analysis engine", false,
                    false)
/*

struct DDVPrinter
    : public DOTGraphTraitsPrinter<PolyhedralInfo, false> {
  static char ID;
  ScopPrinter()
      : DOTGraphTraitsPrinter<PolyhedralInfo, false>("DDV", ID) {}
};
char DDVPrinter::ID = 0;

static RegisterPass<DDVPrinter> X("dot-ddv",
                                   "Polly - Print DDVs of function");
*/
