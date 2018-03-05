#ifndef PadEventAction_h
#define PadEventAction_h 1

#include "G4UserEventAction.hh"

class PadAnalysisManager;

class PadEventAction : public G4UserEventAction {
public:
  PadEventAction(PadAnalysisManager* = 0);
  virtual ~PadEventAction();
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
private:
  PadAnalysisManager* fAnalysisManager;
};

#endif
