#ifndef PadSteppingAction_h
#define PadSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class PadAnalysisManager;

class PadSteppingAction : public G4UserSteppingAction {
public:
  PadSteppingAction(PadAnalysisManager* = 0);
  virtual ~PadSteppingAction();
  virtual void UserSteppingAction(const G4Step*);
private:
  PadAnalysisManager* fAnalysisManager;
};

#endif
