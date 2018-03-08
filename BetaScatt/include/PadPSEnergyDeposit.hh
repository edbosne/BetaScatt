#ifndef PADPSENERGYDEPOSIT_HH_
#define PADPSENERGYDEPOSIT_HH_ 1

#include "G4PSEnergyDeposit.hh"

class G4Step;


class PadPSEnergyDeposit : public G4PSEnergyDeposit
{

public:
    PadPSEnergyDeposit(G4String name, G4int depth=0); // default unit
    PadPSEnergyDeposit(G4String name, const G4String& unit, G4int depth=0);
    virtual ~PadPSEnergyDeposit();

protected:
	virtual G4int GetIndex(G4Step* step);

};

#endif /* PADPSENERGYDEPOSIT_HH_ */
