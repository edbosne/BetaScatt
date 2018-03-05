#include "PadPhysicsList.hh"
#include "globals.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include <iomanip>
#include "PadCentralData.hh"

PadPhysicsList::PadPhysicsList(): G4VUserPhysicsList()
{
  PadCentralData* dataPointer = PadCentralData::GetDataPointer();
  SetVerboseLevel(dataPointer->GetVerbosePhysics());
}

PadPhysicsList::~PadPhysicsList()
{}

void PadPhysicsList::ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void PadPhysicsList::ConstructProcess()
{

  //Register processes transportation and EM interactions
  AddTransportation();
  ConstructEM();
}

#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyRayleigh.hh"

#include "G4MultipleScattering.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyIonisation.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"
//add in geant4.7
//#include "G4StepLimiter.hh"

void PadPhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      
      theLEPhotoElectric   = new G4LowEnergyPhotoElectric();
      theLECompton         = new G4LowEnergyCompton();
      theLEGammaConversion = new G4LowEnergyGammaConversion();
      theLERayleigh        = new G4LowEnergyRayleigh();
      
      pmanager->AddDiscreteProcess(theLEPhotoElectric);
      pmanager->AddDiscreteProcess(theLECompton);
      pmanager->AddDiscreteProcess(theLERayleigh);
      pmanager->AddDiscreteProcess(theLEGammaConversion);

    }
    else if (particleName == "e-") {

      theeminusLEIonisation = new G4LowEnergyIonisation();
      theeminusLEBremsstrahlung = new G4LowEnergyBremsstrahlung();
      theeminusMultipleScattering = new G4MultipleScattering();
      //add in geant4.7
      //theeminusStepLimiter = new G4StepLimiter();

      pmanager->AddProcess(theeminusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeminusLEIonisation,-1,2,2);
      //never use this form of LEIonisation: e- don't stop
      //pmanager->AddProcess(theLEIonisation,-1,-1,2);
      pmanager->AddProcess(theeminusLEBremsstrahlung,-1,-1,3);
      //add in geant4.7
      //pmanager->AddProcess(theeminusStepLimiter,-1,-1,4);
    }
    else if (particleName == "e+") {
      
      theeplusIonisation = new G4eIonisation();
      theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusMultipleScattering = new G4MultipleScattering();
      theeplusAnnihilation = new G4eplusAnnihilation();

      pmanager->AddProcess(theeplusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeplusIonisation,-1,2,2);
      pmanager->AddProcess(theeplusBremsstrahlung,-1,-1,3);
      pmanager->AddProcess(theeplusAnnihilation,0,-1,4);
    }

  }
}

#include "G4UserSpecialCuts.hh"

void PadPhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 
  //

  PadCentralData* dataPointer = PadCentralData::GetDataPointer();
  G4double ParticleCut = dataPointer->GetDefaultCutValue();
  if (ParticleCut > 0) this->SetDefaultCutValue(ParticleCut);

  SetCutsWithDefault();
     
  if (verboseLevel>0) DumpCutValuesTable();
}
