//===- polly/PolyhedralInfo.h - PolyhedralInfo class definition -*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
///
/// This file contains the declaration of the PolyhedralInfo class, which will
/// provide an interface to expose polyhedral analysis information of Polly.
///
/// This is work in progress. We will add more API's as and when deemed
/// required.
//===----------------------------------------------------------------------===///

#ifndef POLLY_POLYHEDRAL_INFO_H
#define POLLY_POLYHEDRAL_INFO_H

#include "llvm/Pass.h"
#include "isl/ctx.h"
#include "isl/union_map.h"
#include "llvm/Analysis/DependenceAnalysis.h"

namespace llvm {
class Loop;
class Instruction;
/// based on the FullDependence class in llvm/Analysis/DependenceAnalysis.h
class DependenceDirectionVector : public Dependence {
public:
  DependenceDirectionVector(Instruction *Src, Instruction *Dst,
                 bool LoopIndependent, unsigned Levels);

  DependenceDirectionVector(DependenceDirectionVector &&RHS)
      : Dependence(std::move(RHS)), Levels(RHS.Levels),
        LoopIndependent(RHS.LoopIndependent), Consistent(RHS.Consistent),
        DV(std::move(RHS.DV)) {}

  /// isLoopIndependent - Returns true if this is a loop-independent
  /// dependence.
  bool isLoopIndependent() const override { return LoopIndependent; }

  /// isConfused - Returns true if this dependence is confused
  /// (the compiler understands nothing and makes worst-case
  /// assumptions).
  bool isConfused() const override { return false; }

  /// isConsistent - Returns true if this dependence is consistent
  /// (occurs every time the source and destination are executed).
  bool isConsistent() const override { return Consistent; }

  /// getLevels - Returns the number of common loops surrounding the
  /// source and destination of the dependence.
  unsigned getLevels() const override { return Levels; }

  /// getDirection - Returns the direction associated with a particular
  /// level.
  unsigned getDirection(unsigned Level) const override;

  /// getDistance - Returns the distance (or NULL) associated with a
  /// particular level.
  const llvm::SCEV *getDistance(unsigned Level) const override;

  /// isPeelFirst - Returns true if peeling the first iteration from
  /// this loop will break this dependence.
  bool isPeelFirst(unsigned Level) const override;

  /// isPeelLast - Returns true if peeling the last iteration from
  /// this loop will break this dependence.
  bool isPeelLast(unsigned Level) const override;

  /// isSplitable - Returns true if splitting the loop will break
  /// the dependence.
  bool isSplitable(unsigned Level) const override;

  /// isScalar - Returns true if a particular level is scalar; that is,
  /// if no subscript in the source or destination mention the induction
  /// variable associated with the loop at this level.
  bool isScalar(unsigned Level) const override;

	/// getDistance - Returs the integer value of distance.
	int getDepDistance(unsigned Level) const;

private:
  unsigned short Levels;
  bool LoopIndependent;
  bool Consistent; // Init to true, then refine.
	int DepDistance; //Temporary integer version
  friend class llvm::DependenceInfo;

public:
  std::unique_ptr<DVEntry[]> DV; 
};
} // namespace llvm

namespace polly {

class Scop;
class ScopInfo;
class DependenceInfoWrapperPass;

class PolyhedralInfo : public llvm::FunctionPass {
public:
  static char ID; // Pass identification, replacement for typeid

  /// Construct a new PolyhedralInfo pass.
  PolyhedralInfo() : FunctionPass(ID) {}
  ~PolyhedralInfo() {}

  /// Check if a given loop is parallel.
  ///
  /// @param L The loop.
  ///
  /// @return  Returns true, if loop is parallel false otherwise.
  bool isParallel(llvm::Loop *L) const;

 /// @brief   Check if a given loop is vectorizable and compute the minimum
  ///          dependence distance.
  /// 
  /// @param L  The loop.
  /// @param VF The vectorization factor. Update if bigger VF possible.
  ///           If loop is trivially vectorizable, set VF to INT_MAX.
  ///           If loop is vectorizable for the given VF, update the VF to the
  ///           minimum dependence distance, which is the maximum permissible
  ///           VF.
  /// 
  /// @return  Returns true if loop is vectorizable for the given @p VF, false
  ///          otherwise.
  bool isVectorizable(llvm::Loop *L, int *VF) const;

  /// Return the SCoP containing the @p L loop.
  ///
  /// @param L The loop.
  ///
  /// @return  Returns the SCoP containing the given loop.
  ///          Returns null if the loop is not contained in any SCoP.
  const Scop *getScopContainingLoop(llvm::Loop *L) const;

  /// Computes the partial schedule for the given @p L loop.
  ///
  /// @param S The SCoP containing the given loop
  /// @param L The loop.
  ///
  /// @return  Returns the partial schedule for the given loop
  __isl_give isl_union_map *getScheduleForLoop(const Scop *S,
                                               llvm::Loop *L) const;

  /// Get the SCoP and dependence analysis information for @p F.
  bool runOnFunction(llvm::Function &F) override;

  /// Release the internal memory.
  void releaseMemory() override {}

  /// Print to @p OS if each dimension of a loop nest is parallel or not.
  void print(llvm::raw_ostream &OS,
             const llvm::Module *M = nullptr) const override;

  /// Register all analyses and transformation required.
  void getAnalysisUsage(llvm::AnalysisUsage &AU) const override;

  /// Checks memory dependence between two LLVM instructions
  /// @param Inst1    1st LLVM memory instruction
  /// @param Inst2    2nd LLVM memory instruction
  /// @return         Returns true if dependenent, false otherwise
  bool isDependent(llvm::Instruction *Inst1, llvm::Instruction *Inst2);

  /// Get the dependence direction vectors
	std::vector<llvm::DependenceDirectionVector*> *getDDVs(const Scop *S) const;

private:

	/// Check if a given loop is parallel or vectorizable.
  ///
  /// @param L             The loop.
  /// @param MinDepDistPtr If not nullptr, the minimal dependence distance will
  ///                      be returned at the address of that pointer
  ///
  /// @return  Returns true if loop is parallel or vectorizable, false
  ///          otherwise.
  bool checkParallel(llvm::Loop *L,
                     __isl_give isl_pw_aff **MinDepDistPtr = nullptr) const;

/// Calculate the dependence direction vector for a map
  ///
  /// @param map  The map to extract a DDV from
  /// @param user The user pointer, in this case, a pointer to a vector of
  ///             DependenceDirectionVector*
  ///
  /// @return isl_stat_ok is returned if no error occurs, else isl_stat_error
  ///         is returned
	//FIXME: The above comment doesn't match.
  void getMapDDV(const Scop *S, __isl_take isl_map *map,
									std::vector<llvm::DependenceDirectionVector *> *DDV) const;

  ScopInfo *SI;
  DependenceInfoWrapperPass *DI;
};

} // end namespace polly

namespace llvm {
class PassRegistry;
void initializePolyhedralInfoPass(llvm::PassRegistry &);
} // namespace llvm

#endif
