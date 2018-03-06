#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ios.hh"


#ifdef G4ANALYSIS_USE
//#include <AIDA/IAnalysisFactory.h>
#include <TFile.h>
#include <TDirectory.h>
#endif

#include "PadAnalysisManager.hh"
#include "PadDetectorConstruction.hh"
//#include "PadPhysicsList.hh"
#include "PadPrimaryGeneratorAction.hh"
#include "PadRunAction.hh"
#include "PadEventAction.hh"
#include "PadSteppingAction.hh"
#include "PadCentralData.hh"
#include "globals.hh"

#include "Randomize.hh"

#ifdef GEOMETRY_DEBUG
#include "PadVisManager.hh"
#endif

void WriteHelp() {
  G4cout << "Usage: BetaScat <input-file> <output-path>" << G4endl << G4endl;
  G4cout << "This program simulates the scattering of electrons in" << G4endl;
  G4cout << "the emission channeling setups of IKS and ITN." << G4endl << G4endl;
  
  G4cout << "If command-line arguments are passed to the program (not" << G4endl;
  G4cout << "mandatory), the FIRST argument should be the INPUT FILE" << G4endl;
  G4cout << "(give full path).  The SECOND optional argument should" << G4endl;
  G4cout << "contain the PATH where the other in- and output files" << G4endl;
  G4cout << "can be found/written to.  Caution: if set, this path" << G4endl;
  G4cout << "argument will override the FilePath entry in the input" << G4endl;
  G4cout << "file." << G4endl << G4endl;
  G4cout << "More information on how to run the program can be found" << G4endl;
  G4cout << "in manual.pdf." << G4endl << G4endl;

  G4cout << "(c) Bart De Vries - 2006" << G4endl; 
}


int main(int argc, char* argv[])
{

#ifdef GDEBUG
std::cout << " begin program! " << std::endl;
#endif

   // First test if the short help text should be printed
  if ( (argc == 2) && (argv[1][0] == '-') ) {
    WriteHelp();
    return 0;
  }

#ifdef GDEBUG
std::cout << " help text. " << std::endl;
#endif

  // Construct default run manager
  G4RunManager* runManager = new G4RunManager;

  // Construct object to hold data to be shared among classes
  // The constructor will read the default input file, or the
  // one given on the command line
  PadCentralData *dataObject = new PadCentralData(argc,argv);

  // Before reading the input file, check if there are command
  // line arguments, which contain alternative files and paths
  // If arguments are given, pass them on to *dataObject.
  // First argument should be the input filename (full path!!), 
  // the second one the path of the output files.

  // Set the Random number generator engine and
  // read previously saved seed file
  G4String SeedFileName = dataObject->GetSeedFileName();
  CLHEP::HepRandom::createInstance();
  CLHEP::RanluxEngine* theRanluxEngine = new CLHEP::RanluxEngine();
  theRanluxEngine->restoreStatus(SeedFileName);
  CLHEP::HepRandom::setTheEngine(theRanluxEngine);

#ifdef GDEBUG
std::cout << " CLHEP intitialization " << std::endl;
#endif

  // Construct analysis object if variable is set
#ifdef G4ANALYSIS_USE
//   AIDA::IAnalysisFactory* analysisFactory = AIDA_createAnalysisFactory();
//   PadAnalysisManager* analysisManager = new PadAnalysisManager(dataObject,analysisFactory);
   TFile * file = new TFile (dataObject->GetHistogramFileName(),"RECREATE");
   PadAnalysisManager* analysisManager = new PadAnalysisManager(dataObject,file);
#endif

#ifdef GDEBUG
std::cout << " Created file and analysisManager " << std::endl;
#endif

  // Mandatory initialization classes
  // Detector Construction
  runManager->SetUserInitialization(new PadDetectorConstruction);

  // Physics List
  G4int verbose = 1;
  G4PhysListFactory factory;
  G4VModularPhysicsList* physlist = factory.GetReferencePhysList("QBBC_EMY"); //QBBC_EMZ
  physlist->SetVerboseLevel(verbose);
  runManager->SetUserInitialization(physlist);
  //runManager->SetUserInitialization(new PadPhysicsList);

#ifdef GDEBUG
std::cout << " Set user initialization " << std::endl;
#endif

#ifdef G4ANALYSIS_USE
  // Mandatory user action class
  runManager->SetUserAction(new PadPrimaryGeneratorAction(dataObject, analysisManager));
#endif
  // Analysis routines to make histograms

#ifdef G4ANALYSIS_USE
  runManager->SetUserAction(new PadRunAction(analysisManager));
  runManager->SetUserAction(new PadEventAction(analysisManager));
  runManager->SetUserAction(new PadSteppingAction(analysisManager));
#else
  runManager->SetUserAction(new PadRunAction());
  runManager->SetUserAction(new PadEventAction());
  runManager->SetUserAction(new PadSteppingAction());
#endif

#ifdef GDEBUG
std::cout << " runManager setting user action " << std::endl;
#endif

  // Initialize G4 kernel
  runManager->Initialize();
#ifdef GDEBUG
std::cout << " runManager initialize " << std::endl;
#endif

#ifdef GEOMETRY_DEBUG
  G4UIsession* session = new G4UIterminal(new G4UItcsh);
#else
  // get the pointer to UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose "+dataObject->GetVerboseRun());
  UI->ApplyCommand("/event/verbose "+dataObject->GetVerboseEvent());
  UI->ApplyCommand("/tracking/verbose "+dataObject->GetVerboseTrack());
std::cout << " G4UI verbosity "  << std::endl;
#endif



#ifdef GEOMETRY_DEBUG
  G4VisManager* visManager = new PadVisManager;
  visManager->Initialize();
#endif

  // start a run

#ifdef GEOMETRY_DEBUG
  session->SessionStart();
  delete session;
#else
  runManager->BeamOn(dataObject->GetNumberOfEvents());
  std::cout << " runManager Beam on. "  << std::endl;
#endif
  
  //job termination
  
  if(dataObject->GetWriteSeed() == 1) {
    theRanluxEngine->saveStatus(SeedFileName);
  }
  delete theRanluxEngine;
  
//#ifdef GDEBUG
std::cout << " Ending job. " << std::endl;
//#endif

#ifdef G4ANALYSIS_USE
  delete analysisManager;


//  delete analysisFactory;
  file->Write(0, TObject::kOverwrite);
  file->Close();
  delete file;
#endif

std::cout << " Ending job. " << std::endl;
#ifdef GEOMETRY_DEBUG
  delete visManager;
#endif
  delete dataObject;
  delete runManager;

  return 0;
}



