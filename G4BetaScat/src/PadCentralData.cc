#include "PadCentralData.hh"
#include "globals.hh"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

PadCentralData* PadCentralData::fPointer = 0;

PadCentralData* PadCentralData::GetDataPointer() {
  return fPointer;
}

PadCentralData::PadCentralData(int argc, char* argv[]) {
  fPointer = this;

  StartingPosition.resize(4,0.);
  StartingMomentum.resize(4,0.);

  // First set the input file and path.  Use default values, unless
  // command line arguments are given: the first argument defines the
  // input file that should be used, the second the directory of 
  // in- and output

  FilePath = ".";
  InputFileName = "input.data";

  if (argc > 1) {
    InputFileName = argv[1];
    if (argc > 2) FilePath = argv[2];
  }
  

  // This is the list of variables that STRICTLY NEED TO BE initialized

  MaximumHEnergy = 1;
  NumberOfConvElecLines = 0;
  NumberOfEnergyWindows = 1;

  vector<double> window;
  window.push_back(0); window.push_back(1);
  EnergyWindow.push_back(window);

  DefaultCutValue = -1;
  MaxFilmStep = -1;
  WriteSeed = 0;
  VerboseRun = "0";
  VerboseEvent = "0";
  VerboseTrack = "0";
  VerbosePhysics = 0;
  StartingPosition[0]=0;
  StartingMomentum[0]=0;

  // Default values for variables from input file

  OutputFileName = "output.data";
  SeedFileName = "Ranlux.conf";
  HistogramFileName = "pad.aida";
  NumberOfEvents = 1; 
  EventWriteOut = 1000000;
  PadSetup = 1;
  FilmThickness = 1.6 * um;
  FilmMaterial = "GaN";
  SampleThickness = 0.04 * cm;
  SubstrateMaterial = "Sapphire";
  MaxRadius = 0.05*cm; 
  MeanDepth = 156.*angstrom;
  SigmaDepth = 53.2*angstrom;
  ThetaRotation = 0*deg; 
  PhiRotation = 0*deg; 
  ElecType = 2;

  // Now the defaults are set, we can read the input file
  ReadInputFile();

}

PadCentralData::~PadCentralData() {
  OutputStream.close();
}

void PadCentralData::ReadInputFile() {

  using namespace std;

  int const NrOfInputs = 32;
  string* InputFields = new string[NrOfInputs];

  InputFields[0] = "OutputFileName       :";
  InputFields[1] = "SeedFileName         :";
  InputFields[2] = "HistogramFileName    :";
  InputFields[3] = "BetaSpectrumFileName :";
  InputFields[4] = "NumberOfEvents       :";
  InputFields[5] = "MaximumHEnergy       :";
  InputFields[6] = "FilmThickness        :";
  InputFields[7] = "FilmMaterial         :";
  InputFields[8] = "SampleThickness      :";
  InputFields[9] = "SubstrateMaterial    :";
  InputFields[10] = "MaxRadius            :";
  InputFields[11] = "MeanDepth            :";
  InputFields[12] = "SigmaDepth           :";
  InputFields[13] = "ThetaRotation        :";
  InputFields[14] = "PhiRotation          :";
  InputFields[15] = "ElecType             :";
  InputFields[16] = "NumberOfConvElecLines:";
  InputFields[17] = "ConvElecEnergies     :";
  InputFields[18] = "ConvElecWeights      :";
  InputFields[19] = "EnergyWindow         :";
  InputFields[20] = "DefaultCutValue      :";
  InputFields[21] = "WriteSeed            :";
  InputFields[22] = "VerboseRun           :";
  InputFields[23] = "VerboseEvent         :";
  InputFields[24] = "VerboseTrack         :";
  InputFields[25] = "VerbosePhysics       :";
  InputFields[26] = "PadSetup             :";
  InputFields[27] = "EventWriteOut        :";
  InputFields[28] = "StartingPosition     :";
  InputFields[29] = "StartingMomentum     :";
  InputFields[30] = "MaxFilmStep          :";
  InputFields[31] = "FilePath             :";

  ifstream InFile (InputFileName,ios::in);

  G4cout << "InputFileName        : " << InputFileName << G4endl;

  string ReadLine, dummy_str;
  double dummy_dbl, dummy_sum;
  int FieldNr;
  double lower, upper;
  vector<double> limits;

  // main while loop to read data records from file sequentially
  // data is transferred line per line to string ReadLine

  while(getline(InFile,ReadLine,'\n')) {
    
    // Default to illegal output value
    FieldNr = -1;
    
    // Check if the first part of ReadLine contains a known data field
    for (G4int i=0; i<NrOfInputs; i++) {
      if (ReadLine.substr(0,22) == InputFields[i]) {
	FieldNr = i;
      }
    };

    if (FieldNr == -1)
     // If no known data field is found:

      if (ReadLine[0] == '#' || OnlySpaces(ReadLine)) {
	// This case will happen when the line is blank or 
	// contains a comment -> ignore 
      } else {
	// If line contains just unrecognizable garbage
	G4cout << "Illegal Input Field: " << ReadLine << G4endl;
      }
    
    else {
      // If a known data field is found
      
      // First check if the string contains more characters than just
      // the data field itself.  Also check if the extra characters
      // are not just spaces.
      if (ReadLine.size() > 22)
	if (OnlySpaces(ReadLine.substr(22,ReadLine.size()-22))){
	  G4cout << InputFields[FieldNr] << "no value given" << G4endl;
	}
	else {
	  istringstream InputValue(ReadLine.substr(22,ReadLine.size()-22));
	  
	  G4int j = 0;
	  
	  switch(FieldNr) {
	  case 0:
	    InputValue >> OutputFileName;
	    G4cout << InputFields[FieldNr] << " " << OutputFileName << G4endl;
	    OutputFileName = FilePath + "/" + OutputFileName;
	    break;
	  case 1:
	    InputValue >> SeedFileName;
	    G4cout << InputFields[FieldNr] << " " << SeedFileName << G4endl;
	    SeedFileName = FilePath + "/" + SeedFileName;
	    break;
	  case 2:
	    InputValue >> HistogramFileName;
	    G4cout << InputFields[FieldNr] << " " << HistogramFileName << G4endl;
	    HistogramFileName = FilePath + "/" + HistogramFileName;
	    break;
	  case 3:
	    InputValue >> BetaSpectrumFileName;
	    G4cout << InputFields[FieldNr] << " " << BetaSpectrumFileName << G4endl;
	    BetaSpectrumFileName = FilePath + "/" + BetaSpectrumFileName;
	    break;
	  case 4:
	    InputValue >> NumberOfEvents;
	    G4cout << InputFields[FieldNr] << " " << NumberOfEvents << G4endl;
	    break;
	  case 5:
	    InputValue >> MaximumHEnergy;
	    // Also set the upper bound of the first energy window
	    EnergyWindow[0][1] = MaximumHEnergy+1;
	    G4cout << InputFields[FieldNr] << " " << MaximumHEnergy << G4endl;
	    break;
	  case 6:
	    InputValue >> dummy_dbl;
	    FilmThickness = dummy_dbl * cm;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 7:
	    InputValue >> FilmMaterial;
	    G4cout << InputFields[FieldNr] << " " << FilmMaterial << G4endl;
	    break;
	  case 8:
	    InputValue >> dummy_dbl;
	    SampleThickness = dummy_dbl * cm;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 9:
	    InputValue >> SubstrateMaterial;
	    G4cout << InputFields[FieldNr] << " " << SubstrateMaterial << G4endl;
	    break;
	  case 10:
	    InputValue >> dummy_dbl;
	    MaxRadius = dummy_dbl * cm;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 11:
	    InputValue >> dummy_dbl;
	    MeanDepth = dummy_dbl * cm;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 12:
	    InputValue >> dummy_dbl;
	    SigmaDepth = dummy_dbl * cm;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 13:
	    InputValue >> dummy_dbl;
	    ThetaRotation = dummy_dbl * deg;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 14:
	    InputValue >> dummy_dbl;
	    PhiRotation = dummy_dbl * deg;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 15:
	    InputValue >> ElecType;
	    G4cout << InputFields[FieldNr] << " " << ElecType << G4endl;
	    break;
	  case 16:
	    InputValue >> NumberOfConvElecLines;
	    G4cout << InputFields[FieldNr] << " " << NumberOfConvElecLines << G4endl;
	    break;
	  case 17:
	    G4cout << "Conversion Electron Energies in keV          : ";
	    for(j=0; j<NumberOfConvElecLines; j++) {
	      InputValue >> dummy_dbl;
	      ConvElecEnergies.push_back(dummy_dbl*keV);
	      G4cout << setw(9) << dummy_dbl << "  ";
	    }
	    G4cout << G4endl;
	    break;
	  case 18:
	    dummy_dbl = 0;
	    dummy_sum = 0;
	    G4cout << "Conversion Electron Weights (NOT integral)   : ";
	    for (j=0; j<NumberOfConvElecLines; j++) {
	      InputValue >> dummy_dbl; 
	      dummy_sum += dummy_dbl;
	      G4cout << setw(9) << dummy_dbl << "  ";
	      ConvElecWeights.push_back(dummy_sum);
	    }
	    G4cout << G4endl << "Conversion Electron Weights (integral)       : ";
	    for (j=0; j<NumberOfConvElecLines; j++) {
	      ConvElecWeights[j] = ConvElecWeights[j] / dummy_sum;
	      G4cout << setw(9) << ConvElecWeights[j] << "  ";
	    }
	    G4cout << G4endl;
	    ConvElecAbsoluteWeight = dummy_sum;
	    G4cout << "Total absolute weight of conversion electrons: "
		   << setw(9) << ConvElecAbsoluteWeight << G4endl;
	    break;
	  case 19:
	    limits.clear();
	    InputValue >> lower >> upper;
	    limits.push_back(lower);
	    limits.push_back(upper);
	    EnergyWindow.push_back(limits);
	    NumberOfEnergyWindows++;
	    G4cout << "EnergyWindow " << NumberOfEnergyWindows-1 << "       : "
		   << setw(4) << EnergyWindow[NumberOfEnergyWindows-1][0] 
		   << setw(4) << EnergyWindow[NumberOfEnergyWindows-1][1] 
		   << G4endl;
	    break;
	  case 20:
	    InputValue >> dummy_dbl;
	    DefaultCutValue = dummy_dbl * angstrom;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 21:
	    InputValue >> WriteSeed;
	    G4cout << InputFields[FieldNr] << " " << WriteSeed << G4endl;
	    break;
	  case 22:
	    InputValue >> VerboseRun;
	    G4cout << InputFields[FieldNr] << " " << VerboseRun << G4endl;
	    break;
	  case 23:
	    InputValue >> VerboseEvent;
	    G4cout << InputFields[FieldNr] << " " << VerboseEvent << G4endl;
	    break;
	  case 24:
	    InputValue >> VerboseTrack;
	    G4cout << InputFields[FieldNr] << " " << VerboseTrack << G4endl;
	    break;
	  case 25:
	    InputValue >> VerbosePhysics;
	    G4cout << InputFields[FieldNr] << " " << VerbosePhysics << G4endl;
	    break;
	  case 26:
	    InputValue >> PadSetup;
	    G4cout << InputFields[FieldNr] << " " << PadSetup << G4endl;
	    break;
	  case 27:
	    InputValue >> EventWriteOut;
	    G4cout << InputFields[FieldNr] << " " << EventWriteOut << G4endl;
	    break;
	  case 28:
	    G4cout << InputFields[FieldNr] << " ";
	    StartingPosition[0] = 1;
	    for (j=1; j<4; j++) {
	      InputValue >> dummy_dbl;
	      StartingPosition[j] = dummy_dbl * cm;
	      G4cout << dummy_dbl << "  ";
	    };
	    G4cout << G4endl;
	    break;
	  case 29:
	    G4cout << InputFields[FieldNr] << " ";
	    StartingMomentum[0] = 1;
	    for (j=1; j<4; j++) {
	      InputValue >> StartingMomentum[j];
	      G4cout << StartingMomentum[j] << "  ";
	    };
	    G4cout << G4endl;
	    break;
	  case 30:
	    InputValue >> dummy_dbl;
	    MaxFilmStep = dummy_dbl * angstrom;
	    G4cout << InputFields[FieldNr] << " " << dummy_dbl << G4endl;
	    break;
	  case 31:
	    // Only read the new file path if it was not already set as a
	    // command line parameter (i.e. if it is now != ".")
	    if (FilePath == ".") {
	      InputValue >> FilePath;
	      G4cout << InputFields[FieldNr] << " " << FilePath << G4endl;
	      OutputFileName = FilePath + "/" + OutputFileName;
	      HistogramFileName = FilePath + "/" + HistogramFileName;
	      SeedFileName = FilePath + "/" + SeedFileName;
	      BetaSpectrumFileName = FilePath + "/" + BetaSpectrumFileName;
	    } else {
	      InputValue >> dummy_str; // ignore the value
	      G4cout << "FilePath (cmd line)  : " << FilePath << G4endl;
	    }
	    break;
	    
	  } // end switch-statement

	} // end else (if a good value is present on input-file)

    } // end else (if a know data field is found)

  } // end while
  
  G4cout << "Defined starting position: ";
  if (StartingPosition[0] == 0) G4cout << "no"; else G4cout << "yes";
  G4cout << G4endl;

  G4cout << "Defined starting momentum: ";
  if (StartingMomentum[0] == 0) G4cout << "no"; else G4cout << "yes"; 
  G4cout << G4endl;

  InFile.close();
  OutputStream.open(OutputFileName,ios::out);

  delete[] InputFields;
}

G4bool PadCentralData::OnlySpaces(G4String str) {
  G4bool only_spaces = true;
  for (unsigned int i=0; i<str.size(); i++) {
    if (str[i] != ' ') only_spaces = false;
  }
  return only_spaces;
}

void PadCentralData::SetInputFileName(G4String fileName) {
  InputFileName = fileName;
}

G4String PadCentralData::GetInputFileName() {
  return InputFileName;
}

void PadCentralData::SetHistogramFileName(G4String fileName) {
  HistogramFileName = fileName;
}

G4String PadCentralData::GetHistogramFileName() {
  return HistogramFileName;
}

void PadCentralData::SetBetaSpectrumFileName(G4String fileName) {
  BetaSpectrumFileName = fileName;
}

G4String PadCentralData::GetBetaSpectrumFileName() {
  return BetaSpectrumFileName;
}

void PadCentralData::SetFilePath(G4String path) {
  FilePath = path;
}

G4String PadCentralData::GetFilePath() {
  return FilePath;
}

void PadCentralData::SetOutputFileName(G4String fileName) {
  OutputFileName = fileName;
}

G4String PadCentralData::GetOutputFileName() {
  return OutputFileName;
}

void PadCentralData::FlushOutputFile()
{
  OutputStream.close();
  OutputStream.open(OutputFileName,ios::out | ios::app);
}

fstream& PadCentralData::GetOutputStream()
{
  return OutputStream;
}

void PadCentralData::SetSeedFileName(G4String fileName) {
  SeedFileName = fileName;
}

G4String PadCentralData::GetSeedFileName() {
  return SeedFileName;
}

void PadCentralData::SetNumberOfEvents(G4int intvalue) {
  NumberOfEvents = intvalue;
}

G4int PadCentralData::GetNumberOfEvents() {
  return NumberOfEvents;
}

void PadCentralData::SetMaximumHEnergy(G4double value) {
  MaximumHEnergy = value;
  
  // don't forget to set the upper bound of the first energy
  // window appropriately
  EnergyWindow[0][1] = MaximumHEnergy + 1;
}

G4double PadCentralData::GetMaximumHEnergy() {
  return MaximumHEnergy;
}

void PadCentralData::SetPadSetup(G4int intvalue) {
  PadSetup = intvalue;
}

G4int PadCentralData::GetPadSetup() {
  return PadSetup;
}

void PadCentralData::SetFilmThickness(G4double value) {
  FilmThickness = value;
}

G4double PadCentralData::GetFilmThickness() {
  return FilmThickness;
}

void PadCentralData::SetFilmMaterial(G4String value) {
  FilmMaterial = value;
}

G4String PadCentralData::GetFilmMaterial() {
  return FilmMaterial;
}

void PadCentralData::SetSampleThickness(G4double value) {
  SampleThickness = value;
}

G4double PadCentralData::GetSampleThickness() {
  return SampleThickness;
}

void PadCentralData::SetSubstrateMaterial(G4String value) {
  SubstrateMaterial = value;
}

G4String PadCentralData::GetSubstrateMaterial() {
  return SubstrateMaterial;
}

void PadCentralData::SetMaxRadius(G4double value) {
  MaxRadius = value;
}

G4double PadCentralData::GetMaxRadius() {
  return MaxRadius;
}

void PadCentralData::SetMeanDepth(G4double value) {
  MeanDepth = value;
}

G4double PadCentralData::GetMeanDepth() {
  return MeanDepth;
}

void PadCentralData::SetSigmaDepth(G4double value) {
  SigmaDepth = value;
}

G4double PadCentralData::GetSigmaDepth() {
  return SigmaDepth;
}

void PadCentralData::SetThetaRotation(G4double value) {
  ThetaRotation = value;
}

G4double PadCentralData::GetThetaRotation() {
  return ThetaRotation;
}

void PadCentralData::SetPhiRotation(G4double value) {
  PhiRotation = value;
}

G4double PadCentralData::GetPhiRotation() {
  return PhiRotation;
}

void PadCentralData::SetElecType(G4int intvalue) {
  ElecType = intvalue;
}

G4int PadCentralData::GetElecType() {
  return ElecType;
}

void PadCentralData::SetNumberOfConvElecLines(G4int intvalue) {
  NumberOfConvElecLines = intvalue;
}

G4int PadCentralData::GetNumberOfConvElecLines() {
  return NumberOfConvElecLines;
}

void PadCentralData::SetConvElecEnergies(vector<double> energies) {
  ConvElecEnergies = energies;
}

vector<double> PadCentralData::GetConvElecEnergies() {
  return ConvElecEnergies;
}

void PadCentralData::SetConvElecWeights(vector<double> weights) {
  ConvElecWeights = weights;
}

vector<double> PadCentralData::GetConvElecWeights() {
  return ConvElecWeights;
}

void PadCentralData::SetConvElecAbsoluteWeight(G4double value) {
  ConvElecAbsoluteWeight = value;
}

G4double PadCentralData::GetConvElecAbsoluteWeight() {
  return ConvElecAbsoluteWeight;
}

void PadCentralData::SetNumberOfEnergyWindows(G4int intvalue) {
  NumberOfEnergyWindows = intvalue;
}

G4int PadCentralData::GetNumberOfEnergyWindows() {
  return NumberOfEnergyWindows;
}

void PadCentralData::SetEnergyWindow(vector< vector<double> > window) {
  EnergyWindow = window;
}

vector< vector<double> > PadCentralData::GetEnergyWindow() {
  return EnergyWindow;
}

void PadCentralData::SetDefaultCutValue(G4double value) {
  DefaultCutValue = value;
}

G4double PadCentralData::GetDefaultCutValue() {
  return DefaultCutValue;
}

void PadCentralData::SetWriteSeed(G4int intvalue) {
  WriteSeed = intvalue;
}

G4int PadCentralData::GetWriteSeed() {
  return WriteSeed;
}

void PadCentralData::SetVerboseRun(G4String intvalue) {
  VerboseRun = intvalue;
}

G4String PadCentralData::GetVerboseRun() {
  return VerboseRun;
}

void PadCentralData::SetVerboseEvent(G4String intvalue) {
  VerboseEvent = intvalue;
}

G4String PadCentralData::GetVerboseEvent() {
  return VerboseEvent;
}

void PadCentralData::SetVerboseTrack(G4String intvalue) {
  VerboseTrack = intvalue;
}

G4String PadCentralData::GetVerboseTrack() {
  return VerboseTrack;
}

void PadCentralData::SetVerbosePhysics(G4int intvalue) {
  VerbosePhysics = intvalue;
}

G4int PadCentralData::GetVerbosePhysics() {
  return VerbosePhysics;
}

void PadCentralData::SetEventWriteOut(G4int intvalue) {
  EventWriteOut = intvalue;
}

G4int PadCentralData::GetEventWriteOut() {
  return EventWriteOut;
}

void PadCentralData::SetStartingPosition(vector<double> position) {
  StartingPosition = position;
}

vector<double> PadCentralData::GetStartingPosition() {
  return StartingPosition;
}

void PadCentralData::SetStartingMomentum(vector<double> momentum) {
  StartingMomentum = momentum;
}

vector<double> PadCentralData::GetStartingMomentum() {
  return StartingMomentum;
}

void PadCentralData::SetMaxFilmStep(G4double value) {
  MaxFilmStep = value;
}

G4double PadCentralData::GetMaxFilmStep() {
  return MaxFilmStep;
}
