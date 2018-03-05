#ifdef G4ANALYSIS_USE
#include "PadAnalysisManager.hh"
#endif

#include "PadSteppingAction.hh"

PadSteppingAction::PadSteppingAction(
 PadAnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

PadSteppingAction::~PadSteppingAction(){}
void PadSteppingAction::UserSteppingAction(const G4Step* aStep){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->Step(aStep);
#endif
}

