#include "PadPrimaryGeneratorAction.hh"
#include "PadAnalysisManager.hh"
#include "PadDetectorConstruction.hh"
#include "PadCentralData.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"

#include "globals.hh"
#include "Randomize.hh"
#include <math.h>


PadPrimaryGeneratorAction::PadPrimaryGeneratorAction(PadCentralData *dataObject,
						       PadAnalysisManager* analysisObject):
  dataPointer(dataObject), analysisPointer(analysisObject)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* electron = particleTable->FindParticle("e-");

  particleGun->SetParticleDefinition(electron);
}


PadPrimaryGeneratorAction::~PadPrimaryGeneratorAction()
{
  delete particleGun;
}

void PadPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double MaxRadius = dataPointer->GetMaxRadius();
  G4double MeanDepth = dataPointer->GetMeanDepth();
  G4double SigmaDepth = dataPointer->GetSigmaDepth();
  G4double ThetaRotation = dataPointer->GetThetaRotation();
  G4double PhiRotation = dataPointer->GetPhiRotation();

  G4double eTemp1X = 0, eTemp1Y = 0, eTemp1Z = 0, Phi = 0;
  G4double eTemp2X = 0, eTemp2Y = 0, eTemp2Z = 0;
  G4double SinTheta = 0, CosTheta = 0, Radius = 0, Radius2 = 0;
  G4double eMomentumX = 0, eMomentumY = 0, eMomentumZ = 0;
  G4double ePositionX = 0, ePositionY = 0, ePositionZ = 0;

  G4int const NumberConvElecLines = dataPointer->GetNumberOfConvElecLines();
  G4int i;
  
  vector<double> ConvElecWeights = dataPointer->GetConvElecWeights();
  vector<double> ConvElecEnergies = dataPointer->GetConvElecEnergies();
  G4double ConvElecAbsoluteWeight = dataPointer->GetConvElecAbsoluteWeight();

  vector<double> StartingPosition = dataPointer->GetStartingPosition();
  vector<double> StartingMomentum = dataPointer->GetStartingMomentum();
 
  G4double eEnergy = 0.;
  
  // If a fixed starting position is given then StartingPosition[0]=1
  // otherwise it has to be taken random within the implantation profile
  
  
  eTemp1Z = -(CLHEP::RandGauss::shoot() * SigmaDepth + MeanDepth);
  if (eTemp1Z > 0.) eTemp1Z = 0.;
  Phi = 2. * pi * CLHEP::RandFlat::shoot();
  Radius2 = MaxRadius * MaxRadius * CLHEP::RandFlat::shoot();
  Radius = sqrt(Radius2);
  eTemp1X = Radius * cos(Phi);
  eTemp1Y = Radius * sin(Phi);
  
  // TODO check if values need to be converted to radian
  eTemp2X = eTemp1X * cos(PhiRotation) - eTemp1Y * sin(PhiRotation);
  eTemp2Y = eTemp1X * sin(PhiRotation) + eTemp1Y * cos(PhiRotation);
  eTemp2Z = eTemp1Z;
  
  ePositionX = eTemp2X * cos(ThetaRotation) - eTemp2Z * sin(ThetaRotation); 
  ePositionY = eTemp2Y;
  ePositionZ = eTemp2X * sin(ThetaRotation) + eTemp2Z * cos(ThetaRotation);
  
  if (StartingPosition[0] == 1 ) {
    ePositionX = StartingPosition[1];
    ePositionY = StartingPosition[2];
    ePositionZ = StartingPosition[3];
  };
  
  particleGun->SetParticlePosition(G4ThreeVector(ePositionX,ePositionY,ePositionZ));
  
  G4int elecType = dataPointer->GetElecType();
  G4double Random;

  switch (elecType) {
  
  case 0: // beta and conversion electrons
    // Decide first to generate a beta particle or CE electron

    Random = CLHEP::RandFlat::shoot(1.+ConvElecAbsoluteWeight);
    if (Random <= 1) { // beta particle
      eEnergy = (analysisPointer->GenerateRandomBeta())*keV;
    } else { // conversion electron
      Random = CLHEP::RandFlat::shoot();
      for (i = NumberConvElecLines; i >= 0; i--) {
	if (Random < ConvElecWeights[i]) 
	  {eEnergy = ConvElecEnergies[i];} 
      };
    }
    break;
  
  case 1: // beta electrons

    eEnergy = (analysisPointer->GenerateRandomBeta())*keV;
    break;

  case 2: // conversion electrons

    G4double RandomNumber = CLHEP::RandFlat::shoot();
    for (i = NumberConvElecLines; i >= 0; i--) {
      if (RandomNumber < ConvElecWeights[i]) 
	{eEnergy = ConvElecEnergies[i];} 
    };
    break;

  };

  analysisPointer->BeginEnergy(eEnergy/keV);

  particleGun->SetParticleEnergy(eEnergy);
  //std::cout<<"total E" <<eEnergy<<std::endl;
  // If StartingMomentum[0]=1 then a fixed momentum was given in the input
  // file, otherwise take random momentum direction.

  CosTheta = 2. *  CLHEP::RandFlat::shoot() - 1.;
  SinTheta = sqrt(1. - CosTheta * CosTheta);
  Phi = 2. * pi * CLHEP::RandFlat::shoot();

  eMomentumX = SinTheta * cos(Phi);
  eMomentumY = SinTheta * sin(Phi);
  eMomentumZ = CosTheta;

  if (StartingMomentum[0] == 1) {
    eMomentumX = StartingMomentum[1];
    eMomentumY = StartingMomentum[2];
    eMomentumZ = StartingMomentum[3];
  };

  particleGun->SetParticleMomentumDirection(G4ThreeVector(eMomentumX,eMomentumY,eMomentumZ));
    
  particleGun->GeneratePrimaryVertex(anEvent);
}

