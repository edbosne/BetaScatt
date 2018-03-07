#include "PadPSEnergyDeposit.hh"


PadPSEnergyDeposit::PadPSEnergyDeposit(G4String name, G4int depth)
 : G4PSEnergyDeposit(name,depth)
{}

PadPSEnergyDeposit::PadPSEnergyDeposit(G4String name, const G4String& unit, G4int depth)
 : G4PSEnergyDeposit(name, unit, depth)
{}

PadPSEnergyDeposit::~PadPSEnergyDeposit()
{}

G4int PadPSEnergyDeposit::GetIndex(G4Step* step)
{
    G4StepPoint* point = step->GetPreStepPoint();
    G4TouchableHandle handle = point->GetTouchableHandle();
    G4int X = handle->GetCopyNumber(0); // GetReplicaNumber
    G4int Y = handle->GetCopyNumber(1);
    //G4cout << "copy number X Y " << X << "," << Y << G4endl;
    return X + Y*1000;
}
