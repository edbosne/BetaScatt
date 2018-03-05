#ifndef PadPrimaryGeneratorAction_h
#define PadPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class PadAnalysisManager;
class G4ParticleGun;
class G4Event;
class PadCentralData;

class PadPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction 
{
public:
  PadPrimaryGeneratorAction(PadCentralData*, PadAnalysisManager*);
  ~PadPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun *particleGun;
  PadCentralData *dataPointer;
  PadAnalysisManager *analysisPointer;
};

#endif
