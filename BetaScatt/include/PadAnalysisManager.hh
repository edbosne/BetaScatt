#ifndef PadAnalysisManager_h
#define PadAnalysisManager_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"
#include <fstream>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TDirectory.h>

#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"


class G4Run;
class G4Event;
class G4Step;
class PadCentralData;
/*namespace AIDA {
    class IHistogram1D;
    class IHistogram2D;
    class IAnalysisFactory;
    class ITree;
}*/

class PadAnalysisManager {
public:
//  PadAnalysisManager(PadCentralData*, AIDA::IAnalysisFactory*);
  PadAnalysisManager(PadCentralData*, TFile*);
//  PadAnalysisManager(PadCentralData*, std::string);
  virtual ~PadAnalysisManager();

public:
  virtual void BeginOfRun(const G4Run*);
  virtual void EndOfRun(const G4Run*);
  virtual void BeginOfEvent(const G4Event*);
  virtual void EndOfEvent(const G4Event*);
  virtual void Step(const G4Step*);
  void WriteOut();
  void BeginEnergy(G4double);
  G4double GenerateRandomBeta();

private:
  void InitializeBetaSpectrum();
//  G4double GetRandom(G4double *);
//  G4double GetRandom(AIDA::IHistogram1D*);
  G4double GetRandom(TH1D*);
//  G4int BinarySearch(G4double *, G4double);
//  G4int BinarySearch(AIDA::IHistogram1D*, G4double);
  G4int BinarySearch(TH1D*, G4double);
//  G4double HistoIntegrate(G4double *, G4double, G4double);
//  G4double HistoIntegrate(AIDA::IHistogram1D*, G4double, G4double);
  G4double HistoIntegrate(TH1D*, G4double, G4double);
  std::pair<G4double,G4double>  GetEngCentroid(G4THitsMap<G4double>*) const;

private:
  PadCentralData *dataPointer;
//  AIDA::IAnalysisFactory *fAIDA;
//  AIDA::ITree *fTree;
  TFile * file;
  TDirectory * dAnalysis;
  TDirectory * dHistograms;

private:
  std::string fileName;
  G4int MaximumEnergy;
  G4int DetectorFlag[3], ScatteringFlag[4], BeforeDetectorFlag;
  G4double DepositedEnergyDetector[3][4];
  G4double DepositedEnergyFilm;
  G4double DepositedEnergySubstrate;
  G4double DepositedEnergyHolder;
  G4int NumberOfWindows;
  G4double KineticEnergy;
  std::vector<std::vector<double> > WindowBoundaries;

  // Histograms containing the data from the simulation
//  G4double * HBetaSpectrum, *HBetaIntSpectrum, *HStartingEnergy;
//  G4double * HStartingEnergy, *HBetaIntSpectrumD, *HBetaSpectrumD;
//  AIDA::IHistogram1D *HBetaSpectrum, *HBetaIntSpectrum;
//  AIDA::IHistogram1D* HStartingEnergy;
  TH1D *HBetaSpectrum, *HBetaIntSpectrum;
  TH1D *HStartingEnergy;

  // First index specifies the detector segment, second index
  // specifies the scattering of electrons
/*  AIDA::IHistogram1D* HDetectorDepth[3][4];
  AIDA::IHistogram1D* HDetectorSpectrum[3][4];
  AIDA::IHistogram2D* HDetectorXY[3][4];
  AIDA::IHistogram1D* HDetectorX[3][4];*/
  TH1D * HDetectorDepth[3][4];
  TH1D * HDetectorSpectrum[3][4];
  TH2D * HDetectorXY[3][4];
  TH1D * HDetectorX[3][4];
  TH2D * PixelDetectorPatt[4];
//  G4double * HDetectorDepth[3][4];
//  G4double * HDetectorSpectrum[3][4];
//  G4double * HDetectorX[3][4];
  

/*  AIDA::IHistogram1D* HFilmDepth;
  AIDA::IHistogram1D* HFilmSpectrum;
  AIDA::IHistogram2D* HFilmXY;*/
  TH1D * HFilmDepth;
  TH1D * HFilmSpectrum;
  TH2D * HFilmXY;
//  G4double * HFilmSpectrum;

/*  AIDA::IHistogram1D* HSubstrateDepth;
  AIDA::IHistogram1D* HSubstrateSpectrum;
  AIDA::IHistogram2D* HSubstrateXY;*/
  TH1D * HSubstrateDepth;
  TH1D * HSubstrateSpectrum;
  TH2D * HSubstrateXY;
//  G4double * HSubstrateSpectrum;

/*  AIDA::IHistogram1D* HHolderDepth;
  AIDA::IHistogram1D* HHolderSpectrum;
  AIDA::IHistogram2D* HHolderXY;*/
  TH1D * HHolderDepth;
  TH1D * HHolderSpectrum;
  TH2D * HHolderXY;
//  G4double * HHolderSpectrum;

//  AIDA::IHistogram1D* HTotalEnergyBeforeDetector;
};

#endif

#endif
