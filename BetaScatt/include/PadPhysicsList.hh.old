#ifndef PadPhysicsList_h
#define PadPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4LowEnergyPhotoElectric;
class G4LowEnergyCompton;
class G4LowEnergyGammaConversion;
class G4LowEnergyRayleigh;

class G4MultipleScattering;
class G4LowEnergyBremsstrahlung;
class G4LowEnergyIonisation;

class G4eBremsstrahlung;
class G4eIonisation;
class G4eplusAnnihilation;
//add in geant4.7
//class G4StepLimiter;

class PadPhysicsList: public G4VUserPhysicsList
{
public:
  PadPhysicsList();
  ~PadPhysicsList();
  
protected:
  void ConstructParticle();
  void ConstructProcess();
  
  void SetCuts();
  
protected:
  void ConstructEM();
  
private:
  G4LowEnergyPhotoElectric* theLEPhotoElectric;
  G4LowEnergyCompton* theLECompton;
  G4LowEnergyGammaConversion* theLEGammaConversion;
  G4LowEnergyRayleigh* theLERayleigh;
  
  G4MultipleScattering* theeminusMultipleScattering;
  G4LowEnergyBremsstrahlung* theeminusLEBremsstrahlung;
  G4LowEnergyIonisation* theeminusLEIonisation;
  //add in geant4.7
  //G4StepLimiter* theeminusStepLimiter;
  
  G4MultipleScattering* theeplusMultipleScattering;
  G4eBremsstrahlung* theeplusBremsstrahlung;
  G4eIonisation* theeplusIonisation;
  G4eplusAnnihilation* theeplusAnnihilation;  
};

#endif

