#ifndef PadCentralData_h
#define PadCentralData_h 1

#include "globals.hh"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class PadCentralData;

class PadCentralData
{
public:
  PadCentralData(int, char**);
  ~PadCentralData();

  void ReadInputFile();

  void SetInputFileName(G4String);
  G4String GetInputFileName();
  void SetOutputFileName(G4String fileName);
  G4String GetOutputFileName();
  fstream& GetOutputStream();
  void FlushOutputFile();
  void SetSeedFileName(G4String fileName);
  G4String GetSeedFileName();
  void SetHistogramFileName(G4String fileName);
  G4String GetHistogramFileName();
  void SetBetaSpectrumFileName(G4String fileName);
  G4String GetBetaSpectrumFileName();
  void SetNumberOfEvents(G4int intvalue);
  G4int GetNumberOfEvents();
  void SetEventWriteOut(G4int intvalue);
  G4int GetEventWriteOut();
  void SetMaximumHEnergy(G4double value);
  G4double GetMaximumHEnergy();
  void SetPadSetup(G4int);
  G4int GetPadSetup();
  void SetFilmThickness(G4double value);
  G4double GetFilmThickness();
  void SetFilmMaterial(G4String value);
  G4String GetFilmMaterial();
  void SetSampleThickness(G4double value);
  G4double GetSampleThickness();
  void SetSubstrateMaterial(G4String value);
  G4String GetSubstrateMaterial();
  void SetMaxRadius(G4double value);
  G4double GetMaxRadius();
  void SetSigmaDepth(G4double value);
  G4double GetSigmaDepth();  
  void SetMeanDepth(G4double value);
  G4double GetMeanDepth();
  void SetThetaRotation(G4double value);
  G4double GetThetaRotation();
  void SetPhiRotation(G4double value);
  G4double GetPhiRotation();
  void SetElecType(G4int);
  G4int GetElecType();
  void SetNumberOfConvElecLines(G4int);
  G4int GetNumberOfConvElecLines();
  void SetConvElecEnergies(vector<double>);
  vector<double> GetConvElecEnergies();
  void SetConvElecWeights(vector<double>);
  vector<double> GetConvElecWeights();
  void SetConvElecAbsoluteWeight(G4double);
  G4double GetConvElecAbsoluteWeight();
  void SetNumberOfEnergyWindows(G4int);
  G4int GetNumberOfEnergyWindows();
  void SetEnergyWindow(vector< vector<double> >);
  vector< vector<double> > GetEnergyWindow();
  void SetDefaultCutValue(G4double);
  G4double GetDefaultCutValue();
  void SetWriteSeed(G4int);
  G4int GetWriteSeed();
  void SetVerboseRun(G4String);
  G4String GetVerboseRun();
  void SetVerboseEvent(G4String);
  G4String GetVerboseEvent();
  void SetVerboseTrack(G4String);
  G4String GetVerboseTrack();
  void SetVerbosePhysics(G4int);
  G4int GetVerbosePhysics();
  void SetStartingPosition(vector<double>);
  vector<double> GetStartingPosition();
  void SetStartingMomentum(vector<double>);
  vector<double> GetStartingMomentum();
  void SetMaxFilmStep(G4double);
  G4double GetMaxFilmStep();
  void SetFilePath(G4String);
  G4String GetFilePath();

private:
  G4bool OnlySpaces(G4String);

public:
  static PadCentralData* GetDataPointer();
  
private:
  static PadCentralData* fPointer;
  
private:
  G4String InputFileName;
  G4String OutputFileName;
  G4String FilePath;
  fstream OutputStream;
  G4String SeedFileName;
  G4String HistogramFileName;
  G4String BetaSpectrumFileName;
  G4int NumberOfEvents;
  G4int EventWriteOut;
  G4double MaximumHEnergy;
  G4int PadSetup;
  G4double FilmThickness;
  G4String FilmMaterial;
  G4double SampleThickness;
  G4String SubstrateMaterial;
  G4double MaxRadius;
  G4double MeanDepth;
  G4double SigmaDepth;
  G4double ThetaRotation;
  G4double PhiRotation;
  G4int ElecType;
  G4int NumberOfConvElecLines;
  vector<double> ConvElecEnergies;
  vector<double> ConvElecWeights;
  G4double ConvElecAbsoluteWeight;
  G4int NumberOfEnergyWindows;
  vector< vector<double> > EnergyWindow;
  G4double DefaultCutValue;
  G4int WriteSeed;
  G4String VerboseRun;
  G4String VerboseEvent;
  G4String VerboseTrack;
  G4int VerbosePhysics;
  vector<double> StartingPosition;
  vector<double> StartingMomentum;
  G4double MaxFilmStep;
  ofstream OutFile;

};

#endif
