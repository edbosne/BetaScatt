#ifdef G4ANALYSIS_USE

#include <vector>
#include <iomanip>
#include "PadAnalysisManager.hh"
#include "PadCentralData.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Track.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TDirectory.h>

/*LA
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/ITreeFactory.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IAxis.h>
*/

/*PadAnalysisManager::PadAnalysisManager(PadCentralData *aCentralData, 
				       AIDA::IAnalysisFactory* aAIDA) 
  : dataPointer(aCentralData), fAIDA(aAIDA) {
*/

PadAnalysisManager::PadAnalysisManager(PadCentralData *aCentralData, TFile* rfile)
  : dataPointer(aCentralData), file(rfile) {

  fstream& logFile = dataPointer->GetOutputStream();
  try{

  // Initialize histograms (new form=AIDA)!!!!
  // Start with creating the necessary factories

  // Could fail if no AIDA implementation found :
/*  if(!fAIDA) {
    G4cout << "AIDA analysis factory not found." << G4endl;
    return;
  }
*/

   if (file->GetDirectory("histograms") == 0){
    dHistograms = file->mkdir ("histograms");
    file->cd("histograms");
    std::cout<<"creating histograms!" << std::endl;
   } else dHistograms = file->GetDirectory("histograms");
  } catch (...) {}

/*
  AIDA::ITreeFactory* treeFactory = fAIDA->createTreeFactory();
  if(!treeFactory) return;

  // Create a tree-like container to handle histograms.
  // This tree is associated to the output file (=first parameter).
  std::string opts = "compress=yes";
  //std::string opts = "compress=no";
  fTree = treeFactory->create(dataPointer->GetHistogramFileName(),
			      "xml",false,true,opts);
  //std::string opts = "export=root";
  //fTree = treeFactory->create(dataPointer->GetHistogramFileName(),
  //                            "ROOT",false,true,opts);

  // Factories are not "managed" by an AIDA analysis system.
  // They must be deleted by the AIDA user code.
  delete treeFactory;

  if(!fTree) return;

  fTree->mkdir("histograms");
  fTree->cd("histograms");
*/  


  G4String max_energy_reason = "Maximum e- energy set from file (keV) = ";


  // If CE should be simulated, readjust the maximum energy so that the
  // highest energy line still shows up in the energy spectra
  
  if (dataPointer->GetElecType() != 1) {  // CE or mixed beta/CE
    // If there is a conversion line with higher energy than MaximumHEnergy,
    // take this value as the new setting
    G4double CE_line_energy;
    for (G4int line=0; line < dataPointer->GetNumberOfConvElecLines(); line++) {
      CE_line_energy = dataPointer->GetConvElecEnergies()[line]/keV;
      if (CE_line_energy + 1 > dataPointer->GetMaximumHEnergy()) {
	dataPointer->SetMaximumHEnergy(CE_line_energy+1);
	max_energy_reason = "Maximum energy set from conversion e- lines (keV) = ";
      }
    }
  }
  
  
  // Check if we have to load a beta spectrum (ElecType=0 or 1)
  // This will automatically set the maximum energy for the energy
  // spectra
  
  if (dataPointer->GetElecType() != 2) {  //beta or mixed beta/CE
    G4double CE_max_energy = dataPointer->GetMaximumHEnergy();
    InitializeBetaSpectrum();
    if (CE_max_energy > dataPointer->GetMaximumHEnergy())
      max_energy_reason = "Maximum beta energy from beta-spectrum file (keV) = ";
  }
  
  MaximumEnergy = (int)(dataPointer->GetMaximumHEnergy());
  logFile << max_energy_reason << MaximumEnergy << G4endl << G4endl;


  // Get energy windows from input file

  NumberOfWindows = dataPointer->GetNumberOfEnergyWindows();
  WindowBoundaries = dataPointer->GetEnergyWindow();


  // Before registering histograms, first set the boundary and name 
  // information needed for creating the histograms
  // All histograms use distances in cm and energies in keV internally
  // This means that the histogram file contains everything in cm and keV
  // Do not mix with different units or you'll get in trouble !!!

  G4String DetectorFragment[3], ElectronKind[4];
  
  DetectorFragment[0] = "Full detector: ";
  DetectorFragment[1] = "Central 3x3 pixels: ";
  DetectorFragment[2] = "Central pixel: ";
  
  ElectronKind[0] = " (all e-)";
  ElectronKind[1] = " (non-scattered e-)";
  ElectronKind[2] = " (e- scattered in tubes)";
  ElectronKind[3] = " (backscattered e-)";
  


  G4double DetectorLowerBound, DetectorUpperBound;
  
  switch (dataPointer->GetPadSetup()) {
  case 1:
  default:
    DetectorLowerBound = 28.6;    DetectorUpperBound = 28.65;
    break;
  case 2:
    DetectorLowerBound = 30.175;  DetectorUpperBound = 30.225;
    break;
  case 3:
    DetectorLowerBound = 30.87;  DetectorUpperBound = 30.92;
    break;
  case 4:
    DetectorLowerBound = 30.87;  DetectorUpperBound = 30.92;
    break;
  case 5:
    DetectorLowerBound = 35.785;  DetectorUpperBound = 35.835;
    break;
  case 6:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 7:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 8:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 9:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 10:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 11:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 12:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 13:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 14:
    DetectorLowerBound = 31.575;  DetectorUpperBound = 31.625;
    break;
  case 256:
    DetectorLowerBound = 31.6 - 0.5 * 0.03;  DetectorUpperBound = 31.6 + 0.5 * 0.03;
    break;
  case 512:
	DetectorLowerBound = 31.6 - 0.5 * 0.03;  DetectorUpperBound = 31.6 + 0.5 * 0.03;
    break;
  };
  
  G4double DetectorXYSize[3] = {1.5, 0.2, 0.065};

/*LA*/
/*
  // Create an histo factory that will create histo in the tree :
  AIDA::IHistogramFactory* histoFactory =
    fAIDA->createHistogramFactory(*fTree);
  if(histoFactory) {

    // Histogram which contains starting energy of electrons
    HStartingEnergy = 
      histoFactory->createHistogram1D("Start Energy", MaximumEnergy,
				      0., MaximumEnergy);

//    AIDA::IHistogram1D* test = histoFactory->createHistogram1D("Test", MaximumEnergy, 0., MaximumEnergy);

the end of this if is after the histoFactory delete*/
  HStartingEnergy = new TH1D ("Start Energy", "Start Energy",MaximumEnergy, 0., MaximumEnergy);

  //std::cout << HStartingEnergy->FindBin(MaximumEnergy-0.2) << "=binMax-1 " <<HStartingEnergy->FindBin(MaximumEnergy) << "=binMax binMin=" <<HStartingEnergy->FindBin(0.) << " nbins=" << HStartingEnergy->GetXaxis()->GetNbins() << std::endl;
    // Register the three kinds of histograms (depth profile, energy spectrum 
    // and xy-profile) for the three detector "sizes", i.e. the full detector (i=0),
    // the central 3x3 pixels (i=1), and only the central pixel (i=2).

    // Loop over the three kinds of electrons, i.e. all e- (j=0), 
    // non-scattered electrons (j=1), e- scattered in the tubes (j=2) and
    // e- backscattered in the sample (j=3)

    G4String HistogramName;
    for (G4int i=0; i<3; i++) {
      for (G4int j=0; j<4; j++) {

	HistogramName = DetectorFragment[i] + "Depth" + ElectronKind[j];
    HDetectorDepth[i][j] = new TH1D (HistogramName, HistogramName, 100, DetectorLowerBound, DetectorUpperBound);
//	 histoFactory->createHistogram1D(HistogramName, 100, DetectorLowerBound, DetectorUpperBound);
	
	HistogramName = DetectorFragment[i] + "Spectrum" + ElectronKind[j];
	HDetectorSpectrum[i][j] = new TH1D (HistogramName, HistogramName,MaximumEnergy, 0., MaximumEnergy);
//	  histoFactory->createHistogram1D(HistogramName, MaximumEnergy, 0., MaximumEnergy);

	HistogramName = DetectorFragment[i] + "XY plot" + ElectronKind[j];
	HDetectorXY[i][j] = new TH2D (HistogramName, HistogramName,100, -DetectorXYSize[i],DetectorXYSize[i], 100, -DetectorXYSize[i],DetectorXYSize[i]);
//	  histoFactory->createHistogram2D(HistogramName, 100, -DetectorXYSize[i],DetectorXYSize[i], 100, -DetectorXYSize[i],DetectorXYSize[i]);

	HistogramName = DetectorFragment[i] + "X plot" + ElectronKind[j];
	HDetectorX[i][j] = new TH1D (HistogramName, HistogramName,500, -DetectorXYSize[i],DetectorXYSize[i]);
//	  histoFactory->createHistogram1D(HistogramName, 500, -DetectorXYSize[i],DetectorXYSize[i]);

      }
    }


    // Register the histograms for the thin film

    HFilmDepth = new TH1D ("Sample: Depth profile","Sample: Depth profile",100, -0.0002, 0.);
//      histoFactory->createHistogram1D("Sample: Depth profile", 100, -0.0002, 0.);

    HFilmSpectrum = new TH1D ("Sample: Energy Spectrum","Sample: Energy Spectrum",MaximumEnergy, 0., MaximumEnergy);
//      histoFactory->createHistogram1D("Sample: Energy Spectrum", MaximumEnergy, 0., MaximumEnergy);

    HFilmXY = new TH2D ("Sample: XY plot","Sample: XY plot",100, -0.25, 0.25, 100, -0.25, 0.25);
//      histoFactory->createHistogram2D("Sample: XY plot", 100, -0.25, 0.25, 100, -0.25, 0.25);


    // Register the histograms for the substrate

    HSubstrateDepth = new TH1D ("Substrate: Depth profile","Substrate: Depth profile", 100, -0.04, 0.);
//      histoFactory->createHistogram1D("Substrate: Depth profile", 100, -0.04, 0.);
    HSubstrateSpectrum = new TH1D ("Substrate: Energy Spectrum", "Substrate: Energy Spectrum", MaximumEnergy, 0., MaximumEnergy);
//      histoFactory->createHistogram1D("Substrate: Energy Spectrum", MaximumEnergy, 0., MaximumEnergy);
    HSubstrateXY = new TH2D ("Substrate: XY plot", "Substrate: XY plot",100, -0.25, 0.25, 100, -0.25, 0.25);
//      histoFactory->createHistogram2D("Substrate: XY plot", 100, -0.25, 0.25, 100, -0.25, 0.25);


    // Register the histograms for the sample holder

    HHolderDepth = new TH1D ("Holder: Depth profile", "Holder: Depth profile", 100, -0.14, -0.04);
//      histoFactory->createHistogram1D("Holder: Depth profile", 100, -0.14, -0.04);
    HHolderSpectrum = new TH1D ("Holder: Energy Spectrum", "Holder: Energy Spectrum", MaximumEnergy, 0., MaximumEnergy);
//      histoFactory->createHistogram1D("Holder: Energy Spectrum", MaximumEnergy, 0., MaximumEnergy);
    HHolderXY = new TH2D ("Holder: XY plot", "Holder: XY plot", 100, -0.5, 0.5, 100, -0.5, 0.5);
//      histoFactory->createHistogram2D("Holder: XY plot", 100, -0.5, 0.5, 100, -0.5, 0.5);

  //  delete histoFactory;
 // }  // end of histogram creation

}

PadAnalysisManager::~PadAnalysisManager() {
} 

void PadAnalysisManager::BeginEnergy(G4double energy) {
  HStartingEnergy->Fill(energy,1.0);
}

void PadAnalysisManager::BeginOfRun(const G4Run*){
}

void PadAnalysisManager::EndOfRun(const G4Run*){
  fstream& logFile = dataPointer->GetOutputStream();
  logFile << G4endl << "At end of Run" << G4endl;
  WriteOut();
}


void PadAnalysisManager::WriteOut(){
  using namespace std;
  
// Save histogram file without versions
  //fTree->commit();
  file->Write(0, TObject::kOverwrite);

  // End write out the summary for the energy windows to file
  fstream& logFile = dataPointer->GetOutputStream();

  //apply energy windows to get the count rate and fractions
  std::vector<vector<vector<G4double> > > DetectorCounts;
  DetectorCounts.resize(NumberOfWindows);

  // loop over energy windows
  for (G4int window=0; window<NumberOfWindows; window++) {
    DetectorCounts[window].resize(3);

    logFile << "Count in energy window " << window << " = " 
	    << WindowBoundaries[window][0] << " keV to " 
	    << WindowBoundaries[window][1] << " keV" << G4endl;
    
    logFile << "DetSegm      Total e-    Dir e-   Tube e-   Back e-  BackFrac  CorrFact" << G4endl;
    
    // loop over detector (i)
    for (G4int i=0; i<3; i++) {
      DetectorCounts[window][i].resize(4,0);
      switch (i) {
      case 0:
	logFile << setw(10) << left << "Full det" << right;
	break;
      case 1:
	logFile << setw(10) << left << "3x3 pix" << right;
	break;
      case 2:
	logFile << setw(10) << left << "1x1 pix" << right;
	break;
      }      
      // loop over electron kind (k)
      for (G4int k=0; k<4; k++) {
	//Integrate counts within energy window
	DetectorCounts[window][i][k] = HistoIntegrate(HDetectorSpectrum[i][k],
						      WindowBoundaries[window][0],
						      WindowBoundaries[window][1]);
	
	logFile << " " << setw(9) << DetectorCounts[window][i][k];
      }

      G4double BackFrac, CorrFact;

      if (DetectorCounts[window][i][0] != 0) {
	BackFrac = 1 - DetectorCounts[window][i][1]/DetectorCounts[window][i][0];
      } else BackFrac = 0;
      if (DetectorCounts[window][i][1] != 0) {
	CorrFact = DetectorCounts[window][i][0]/DetectorCounts[window][i][1];
      } else CorrFact = 0;

      logFile << "  " << fixed << setw(9) << BackFrac 
	      << " " << setw(9) << CorrFact << G4endl;
      logFile.unsetf(ios::fixed);
    }
  }
  dataPointer->FlushOutputFile();
}

void PadAnalysisManager::BeginOfEvent(const G4Event*){
  G4int i,j;

  for (i=0;i<3;++i) DetectorFlag[i] = 0; // loop over detector segment
  for (i=0;i<4;++i) ScatteringFlag[i] = 0; // loop over possible scattering events
  for (i=0;i<3;++i)
    for (j=0;j<4;++j) DepositedEnergyDetector[i][j] = 0.;

  DepositedEnergyFilm = 0.;
  DepositedEnergySubstrate = 0.;
  DepositedEnergyHolder = 0.;
}

void PadAnalysisManager::EndOfEvent(const G4Event* aEvent){

  fstream& logFile = dataPointer->GetOutputStream();

  G4int i, k;

  if (DepositedEnergyFilm > 1*eV) {
    HFilmSpectrum->Fill(DepositedEnergyFilm,1.0);
  }
  if (DepositedEnergySubstrate > 1*eV) {
    HSubstrateSpectrum->Fill(DepositedEnergySubstrate,1.0);
  }
  if (DepositedEnergyHolder > 1*eV) {
    HHolderSpectrum->Fill(DepositedEnergyHolder,1.0);
  }

  for (i=0; i<3; i++) { // loop over detector segment
    for (k=0; k<4; k++) { // loop over possible scattering events
      if (DepositedEnergyDetector[i][k] > 1*eV){
        HDetectorSpectrum[i][k]->Fill(DepositedEnergyDetector[i][k],1.0);
      }
    }
  }

  G4int eventNr = aEvent->GetEventID();
  if (eventNr != 0) {
    if ((eventNr % dataPointer->GetEventWriteOut()) == 0) {
      logFile << G4endl << "Event Number: " << eventNr << G4endl;
      WriteOut();
    }
  }
}

void PadAnalysisManager::Step(const G4Step* aStep){
  G4double EnergyStep = ((aStep->GetTotalEnergyDeposit())/keV);
  G4double EnergyStepWeight = EnergyStep/(dataPointer->GetNumberOfEvents());
  
  G4ThreeVector deltaPosition = aStep->GetDeltaPosition();

  G4Track* aTrack = aStep->GetTrack();
  G4VPhysicalVolume* aVolume = aTrack->GetVolume();
  G4String NameVolume = aVolume->GetName();
  G4ThreeVector position = aTrack->GetPosition();

  // get position of electron and transform to virtual box and sample coordinates
  G4double univX = position.x()/cm;
  G4double univY = position.y()/cm;
  G4double univZ = position.z()/cm;

  G4double phi = dataPointer->GetPhiRotation();
  G4double theta = dataPointer->GetThetaRotation();

  G4double boxtempX =   univX * cos(theta) + univZ * sin(theta);
  G4double boxtempY =   univY;
  G4double boxtempZ = - univX * sin(theta) + univZ * cos(theta);

  G4double boxX =   boxtempX * cos(phi) + boxtempY * sin(phi);
  G4double boxY = - boxtempX * sin(phi) + boxtempY * cos(phi);
  G4double boxZ =   boxtempZ;

  //std::cout << aStep->GetTrack()->GetGlobalTime() << " " << EnergyStep << "-" << dataPointer->GetNumberOfEvents() <<" in (" <<univX<<","<<univY<<","<< univZ << ")"<< std::endl;

  // delta position in coordinates of box
  G4double deltaZ = deltaPosition.z() * cos(theta) - deltaPosition.x() * sin(theta);


  if ((NameVolume == "Film" ) && (DetectorFlag[0] == 0)) {
    HFilmDepth->Fill(boxZ,EnergyStepWeight);
    HFilmXY->Fill(boxX,boxY,EnergyStepWeight);
    DepositedEnergyFilm += EnergyStep;
    if (deltaZ < 0) {
      ScatteringFlag[0] = 1;
      ScatteringFlag[3] = 1;
    }
  }

  if ((NameVolume == "Substrate")  && (DetectorFlag[0] == 0)) {
    HSubstrateDepth->Fill(boxZ,EnergyStepWeight);
    HSubstrateXY->Fill(boxX,boxY,EnergyStepWeight);
    DepositedEnergySubstrate += EnergyStep;
    ScatteringFlag[0] = 1;
    ScatteringFlag[3] = 1;
  }

  if ((NameVolume == "SampleHolder") && (DetectorFlag[0] == 0)) {
    HHolderDepth->Fill(boxZ,EnergyStepWeight);
    HHolderXY->Fill(boxX,boxY,EnergyStepWeight);
    DepositedEnergyHolder += EnergyStep;
    ScatteringFlag[0] = 1;
    ScatteringFlag[3] = 1;
  }

  if ((NameVolume == "TaHolderPart1" || NameVolume == "TaHolderPart2" ||
       NameVolume == "Shield" || NameVolume == "Tube1" || 
       NameVolume == "Tube2" || NameVolume == "Tube3" ||
       NameVolume == "Tube4" || NameVolume == "Tube5" ||
       NameVolume == "Tube5" || NameVolume == "Tube6" ||
       NameVolume == "Tube7" || NameVolume == "Tube8" ||
       NameVolume == "Tube9" || NameVolume == "Valve1" || 
       NameVolume == "Valve2" || NameVolume == "Valve3" || 
       NameVolume == "Shield" || NameVolume == "Shield2" ||
       NameVolume == "Collimator" || NameVolume == "Cone" || NameVolume == "Cone2" ||
	NameVolume == "Shieldplate" || NameVolume == "ShieldplateA" ||
	NameVolume == "ShieldplateB" || NameVolume == "ShieldplateC" ||
	NameVolume == "ShieldplateD") 
      && DetectorFlag[0] == 0) {
    ScatteringFlag[0] = 1;
    ScatteringFlag[2] = 1;
  }


  if (NameVolume == "Detector1" || NameVolume == "Detector2" || 
      NameVolume == "Detector3") {
    DetectorFlag[0] = 1;
    HDetectorDepth[0][0]->Fill(univZ,EnergyStepWeight);
    HDetectorXY[0][0]->Fill(univX,univY,EnergyStepWeight);
    HDetectorX[0][0]->Fill(univX,EnergyStepWeight);
    DepositedEnergyDetector[0][0] += EnergyStep;

    
    if (ScatteringFlag[2] == 1) {
      HDetectorDepth[0][2]->Fill(univZ,EnergyStepWeight);
      HDetectorXY[0][2]->Fill(univX,univY,EnergyStepWeight);
      HDetectorX[0][2]->Fill(univX,EnergyStepWeight);
      DepositedEnergyDetector[0][2] += EnergyStep;
    };
    if (ScatteringFlag[3] == 1) {
      HDetectorDepth[0][3]->Fill(univZ,EnergyStepWeight);
      HDetectorXY[0][3]->Fill(univX,univY,EnergyStepWeight);
      HDetectorX[0][3]->Fill(univX,EnergyStepWeight);
      DepositedEnergyDetector[0][3] += EnergyStep;
    };
    if (ScatteringFlag[0] == 0) {
      HDetectorDepth[0][1]->Fill(univZ,EnergyStepWeight);
      HDetectorXY[0][1]->Fill(univX,univY,EnergyStepWeight);
      HDetectorX[0][1]->Fill(univX,EnergyStepWeight);
      DepositedEnergyDetector[0][1] += EnergyStep;
    };
    
    if (NameVolume == "Detector2" || NameVolume == "Detector3") {
      DetectorFlag[1] = 1;
      HDetectorDepth[1][0]->Fill(univZ,EnergyStepWeight);
      HDetectorXY[1][0]->Fill(univX,univY,EnergyStepWeight);
      HDetectorX[1][0]->Fill(univX,EnergyStepWeight);
      DepositedEnergyDetector[1][0] += EnergyStep;
      
      if (ScatteringFlag[2] == 1) {
        HDetectorDepth[1][2]->Fill(univZ,EnergyStepWeight);
        HDetectorXY[1][2]->Fill(univX,univY,EnergyStepWeight);
        HDetectorX[1][2]->Fill(univX,EnergyStepWeight);
	DepositedEnergyDetector[1][2] += EnergyStep;
      };
      if (ScatteringFlag[3] == 1) {
        HDetectorDepth[1][3]->Fill(univZ,EnergyStepWeight);
        HDetectorXY[1][3]->Fill(univX,univY,EnergyStepWeight);
        HDetectorX[1][3]->Fill(univX,EnergyStepWeight);
	DepositedEnergyDetector[1][3] += EnergyStep;
      };
      if (ScatteringFlag[0] == 0) {
        HDetectorDepth[1][1]->Fill(univZ,EnergyStepWeight);
        HDetectorXY[1][1]->Fill(univX,univY,EnergyStepWeight);
        HDetectorX[1][1]->Fill(univX,EnergyStepWeight);
	DepositedEnergyDetector[1][1] += EnergyStep;
      };
      
      if (NameVolume == "Detector3") {
	DetectorFlag[2] = 1;
        HDetectorDepth[2][0]->Fill(univZ,EnergyStepWeight);
        HDetectorXY[2][0]->Fill(univX,univY,EnergyStepWeight);
        HDetectorX[2][0]->Fill(univX,EnergyStepWeight);
	DepositedEnergyDetector[2][0] += EnergyStep;
	
	if (ScatteringFlag[2] == 1) {
          HDetectorDepth[2][2]->Fill(univZ,EnergyStepWeight);
          HDetectorXY[2][2]->Fill(univX,univY,EnergyStepWeight);
          HDetectorX[2][0]->Fill(univX,EnergyStepWeight);
	  DepositedEnergyDetector[2][2] += EnergyStep;
	};
	if (ScatteringFlag[3] == 1) {
          HDetectorDepth[2][3]->Fill(univZ,EnergyStepWeight);
          HDetectorXY[2][3]->Fill(univX,univY,EnergyStepWeight);
          HDetectorX[2][3]->Fill(univX,EnergyStepWeight);
	  DepositedEnergyDetector[2][3] += EnergyStep;
	};
	if (ScatteringFlag[0] == 0) {
          HDetectorDepth[2][1]->Fill(univZ,EnergyStepWeight);
          HDetectorXY[2][1]->Fill(univX,univY,EnergyStepWeight);
          HDetectorX[2][1]->Fill(univX,EnergyStepWeight);
	  DepositedEnergyDetector[2][1] += EnergyStep;
	};
		
      }
      
    }
    
  }  
  
}

void PadAnalysisManager::InitializeBetaSpectrum() {

  vector<double> betaList, betaIntList;

  ifstream betaFile;
  betaFile.open(dataPointer->GetBetaSpectrumFileName());

  G4double temp=0.;
  G4double sum=0.;
  while(!betaFile.eof()) {
    betaFile >> temp;
    sum += temp;
    betaList.push_back(temp);
    betaIntList.push_back(sum);
  };
  betaFile.close();

  // If a larger maximum energy is given, then fill up the vectors until
  // those values
  G4int size = betaList.size();
  G4int max_energy = (int)dataPointer->GetMaximumHEnergy();
  if (size < max_energy) { 
    for (G4int i = size; i < max_energy+1; i++) {
      betaList.push_back(0);
      betaIntList.push_back(sum);
    }
    size = max_energy;
  } 
  else dataPointer->SetMaximumHEnergy(size);

/*LA*/
/*  AIDA::IHistogramFactory* histoFactory = 
    fAIDA->createHistogramFactory(*fTree);
  
  HBetaSpectrum = histoFactory->createHistogram1D("BetaSpectrum",size,0,size);
  HBetaIntSpectrum = histoFactory->createHistogram1D("BetaIntSpectrum",size,0,size);

  delete histoFactory;
*/

  HBetaSpectrum = new TH1D ("BetaSpectrum", "BetaSpectrum", size, 0, size);
  HBetaIntSpectrum = new TH1D ("BetaIntSpectrum", "BetaIntSpectrum", size, 0, size);

  //std::cout << "Beta Spectrum" << std::endl;
  for (G4int i=0; i<size; i++){
      HBetaSpectrum->Fill(i,betaList[i]);
      HBetaIntSpectrum->Fill(i,betaIntList[i]/sum);
      //std::cout << HBetaIntSpectrum->FindBin(i)-1 << "-" << betaIntList[i]/sum << std::endl;
  };
}

G4double PadAnalysisManager::HistoIntegrate(TH1D* aHistogram, G4double lowerCoord, G4double upperCoord) {
//G4double PadAnalysisManager::HistoIntegrate(AIDA::IHistogram1D* aHistogram, 
//					    G4double lowerCoord, G4double upperCoord) {
  G4int lower, upper;
  if (lowerCoord < 0.) {
    lower = 0;
  } else {
//    lower = aHistogram->coordToIndex(lowerCoord);
    lower = aHistogram->FindBin(lowerCoord)-1;
  };

  if (upperCoord >= MaximumEnergy) {
//    upper = aHistogram->coordToIndex(aHistogram->axis().upperEdge()-0.5);
    upper = aHistogram->FindBin(aHistogram->GetXaxis()->GetXmax()-0.5)-1; //Last bin is nbins in root and nbins-1 in aida
  } else {
//    upper = aHistogram->coordToIndex(upperCoord);
    upper = aHistogram->FindBin(upperCoord)-1;
  };
  
/*  G4double sum = 0;
  for (G4int binnr=lower; binnr<=upper; binnr++) {
    sum += aHistogram->binHeight(binnr);
  };
  return sum;*/

  return aHistogram->Integral (lower+1, upper+1);
  //return aHistogram->Integral (lower+1, upper+1, "width");
}



G4double PadAnalysisManager::GenerateRandomBeta() {
  G4double rand;

  if (dataPointer->GetElecType() != 2) {
//    return GetRandom(HBetaIntSpectrum);
    rand = GetRandom(HBetaIntSpectrum);
//    std::cout<<rand << "=rand"<<std::endl;
    return rand;
  } else {
    return 0;
  }
}

//G4double PadAnalysisManager::GetRandom(AIDA::IHistogram1D* aHistogram) {
G4double PadAnalysisManager::GetRandom(TH1D* aHistogram) {

  G4double r = CLHEP::RandFlat::shoot(); // r is in [0,1]
  G4int ibin = BinarySearch(aHistogram,r);

//  std::cout<<ibin << "=ibin r=" << r<<std::endl;
  if(ibin==(-1)) { //r in [0,aHistogram.binHeight(0)]
//    G4double dr = aHistogram->binHeight(0);
// in GetBinContent - bin=0 underflow bin
      G4double dr = aHistogram->GetBinContent(1);
      if(dr==0) {
//      return (G4double)(aHistogram->axis().binLowerEdge(0) + aHistogram->axis().binWidth(0)*r);
        return (G4double)(aHistogram->GetXaxis()->GetBinLowEdge(1) 
                        + aHistogram->GetXaxis()->GetBinWidth(1)*r);
      } else {
//      return (G4double)(aHistogram->axis().binLowerEdge(0) + aHistogram->axis().binWidth(0)* r/dr);
        return (G4double)(aHistogram->GetXaxis()->GetBinLowEdge(1) 
                        + aHistogram->GetXaxis()->GetBinWidth(1)* r/dr);
      }
//  } else if(ibin >= (G4int)(aHistogram->axis().bins()-1)) {
  } else if(ibin >= (G4int)(aHistogram->GetXaxis()->GetNbins()-1)) {
    // We pass here when r is stricly equal to 1.
//    return (G4double)(aHistogram->axis().upperEdge());
      return (G4double)(aHistogram->GetXaxis()->GetXmax());
    // If histogrammed with same binning than aHistogram, 
    // it should go in the overflow.
  } else {
    // (ibin+1) < aHistogram.bins()
//    double dr = aHistogram->binHeight(ibin+1) - aHistogram->binHeight(ibin);
      double dr = aHistogram->GetBinContent(ibin+2) 
                  - aHistogram->GetBinContent(ibin+1);
    if(dr==0) {
//      return (G4double)(aHistogram->axis().binLowerEdge(ibin+1) + aHistogram->axis().binWidth(ibin+1)	* r);
      return (G4double)(aHistogram->GetXaxis()->GetBinLowEdge(ibin+2) 
                        + aHistogram->GetXaxis()->GetBinWidth(ibin+2) * r);
    } else {
//      return (G4double)(aHistogram->axis().binLowerEdge(ibin+1) +  aHistogram->axis().binWidth(ibin+1) * (r-aHistogram->binHeight(ibin))/dr);
        return (G4double)(aHistogram->GetXaxis()->GetBinLowEdge(ibin+2) 
                          + aHistogram->GetXaxis()->GetBinWidth(ibin+2) 
                          * (r - aHistogram->GetBinContent(ibin+1))/dr);
    }
  }
}


//G4int PadAnalysisManager::BinarySearch(AIDA::IHistogram1D* aHistogram, G4double aValue)
G4int PadAnalysisManager::BinarySearch(TH1D* aHistogram, G4double aValue)
  // Binary search in an array of n values to locate value.
  //
  // Array is supposed  to be sorted prior to this call.
  // If match is found, function returns position of element.
  // If no match found, function gives nearest element smaller than value.
  // If aArray is empty (-1) is returned.
  // If aValue is stricly below first element (-1) is returned.
  //
  // From ROOT/TMath code.
{
//  G4int nabove = aHistogram->axis().bins() + 1;
  G4int nabove = aHistogram->GetXaxis()->GetNbins() + 1;
//  G4int nbelow = 0;
  G4int nbelow = 0;
  while(nabove-nbelow > 1) {
    G4int middle = (nabove+nbelow)/2;
//    G4double value = aHistogram->binHeight(middle-1);
    G4double value = aHistogram->GetBinContent(middle);
    if (aValue == value) return middle-1;
    if (aValue  < value) nabove = middle;
    else nbelow = middle;
  }
   return nbelow-1;
}


#endif
