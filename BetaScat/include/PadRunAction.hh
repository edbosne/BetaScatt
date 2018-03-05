#ifndef PadRunAction_h
#define PadRunAction_h 1

#include "G4UserRunAction.hh"

class PadAnalysisManager;

class PadRunAction : public G4UserRunAction {
public:
  PadRunAction(PadAnalysisManager* = 0);
  ~PadRunAction();
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
private:
  PadAnalysisManager* fAnalysisManager;
};

#endif
