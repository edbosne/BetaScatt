#ifdef G4ANALYSIS_USE
#include "PadAnalysisManager.hh"
#endif

#include "PadRunAction.hh"

PadRunAction::PadRunAction(PadAnalysisManager* aAnalysisManager)
:fAnalysisManager(aAnalysisManager){
}

PadRunAction::~PadRunAction(){}

void PadRunAction::BeginOfRunAction(const G4Run* aRun){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->BeginOfRun(aRun);
#endif
}

void PadRunAction::EndOfRunAction(const G4Run* aRun){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->EndOfRun(aRun);
#endif
}
