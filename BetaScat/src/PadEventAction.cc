#ifdef G4ANALYSIS_USE
#include "PadAnalysisManager.hh"
#endif

#include "PadEventAction.hh"

PadEventAction::PadEventAction(
 PadAnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

PadEventAction::~PadEventAction(){}

void PadEventAction::BeginOfEventAction(const G4Event* aEvent){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->BeginOfEvent(aEvent);
#endif
}

void PadEventAction::EndOfEventAction(const G4Event* aEvent) {
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->EndOfEvent(aEvent);
#endif
}
