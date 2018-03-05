#ifndef PadDetectorConstruction_h
#define PadDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"


class G4Box;
class G4Tubs;
class G4Cons;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PadCentralData;

class PadDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  PadDetectorConstruction();
  ~PadDetectorConstruction();
  
public:
  
  G4VPhysicalVolume* Construct();
  
private:
  
  void DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();
  PadCentralData* dataPointer;
  
private:
  G4Box*             solidWorld;
  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;
  
  G4Box*             solidVirtualBox;
  G4LogicalVolume*   logicVirtualBox;
  G4VPhysicalVolume* physiVirtualBox;
  
  G4Box*             solidTaPart1;
  G4LogicalVolume*   logicTaPart1;
  G4VPhysicalVolume* physiTaPart1;

  G4Box*             solidTaPart2;
  G4LogicalVolume*   logicTaPart2;
  G4VPhysicalVolume* physiTaPart2;

  G4Box*             solidFilm;
  G4LogicalVolume*   logicFilm;
  G4VPhysicalVolume* physiFilm;

  G4Box*             solidSubstrate;
  G4LogicalVolume*   logicSubstrate;
  G4VPhysicalVolume* physiSubstrate;
  
  G4Box*             solidSampleHolder;
  G4LogicalVolume*   logicSampleHolder;
  G4VPhysicalVolume* physiSampleHolder;
  
  G4Tubs*            solidShield;
  G4LogicalVolume*   logicShield;
  G4VPhysicalVolume* physiShield;

  G4Box*            solidShield2;
  G4LogicalVolume*   logicShield2;
  G4VPhysicalVolume* physiShield2;

  G4Tubs*            solidTube1;
  G4LogicalVolume*   logicTube1;
  G4VPhysicalVolume* physiTube1;
  
  G4Tubs*            solidTube2;
  G4LogicalVolume*   logicTube2;
  G4VPhysicalVolume* physiTube2;
  
  G4Tubs*            solidTube3;
  G4LogicalVolume*   logicTube3;
  G4VPhysicalVolume* physiTube3;

  G4Tubs*            solidValve1;
  G4LogicalVolume*   logicValve1;
  G4VPhysicalVolume* physiValve1;
  
  G4Tubs*            solidValve2;
  G4LogicalVolume*   logicValve2;
  G4VPhysicalVolume* physiValve2;

  G4Tubs*            solidTube4;
  G4LogicalVolume*   logicTube4;
  G4VPhysicalVolume* physiTube4;

  G4Tubs*            solidTube5;
  G4LogicalVolume*   logicTube5;
  G4VPhysicalVolume* physiTube5;

  G4Tubs*            solidTube6;
  G4LogicalVolume*   logicTube6;
  G4VPhysicalVolume* physiTube6;

  G4Tubs*            solidTube7;
  G4LogicalVolume*   logicTube7;
  G4VPhysicalVolume* physiTube7;

  G4Tubs*            solidTube8;
  G4LogicalVolume*   logicTube8;
  G4VPhysicalVolume* physiTube8;

  G4Tubs*            solidTube9;
  G4LogicalVolume*   logicTube9;
  G4VPhysicalVolume* physiTube9;

  G4Tubs*            solidTube10;
  G4LogicalVolume*   logicTube10;
  G4VPhysicalVolume* physiTube10;

  G4Cons*            solidCone;
  G4LogicalVolume*   logicCone;
  G4VPhysicalVolume* physiCone;

  G4Cons*            solidCone2;
  G4LogicalVolume*   logicCone2;
  G4VPhysicalVolume* physiCone2;

  G4Tubs*            solidCollimator;
  G4LogicalVolume*   logicCollimator;
  G4VPhysicalVolume* physiCollimator;

  G4Box*             solidDetector1;
  G4LogicalVolume*   logicDetector1;
  G4VPhysicalVolume* physiDetector1;
  
  G4Box*             solidDetector2;
  G4LogicalVolume*   logicDetector2;
  G4VPhysicalVolume* physiDetector2;
  
  G4Box*             solidDetector3;
  G4LogicalVolume*   logicDetector3;
  G4VPhysicalVolume* physiDetector3;
  
  G4Box*             solidVacuumbox;
  G4LogicalVolume*   logicVacuumbox;
  G4VPhysicalVolume* physiVacuumbox;

  G4Box*             solidShieldplate;
  G4LogicalVolume*   logicShieldplate;
  G4VPhysicalVolume* physiShieldplate;

  G4Box*             solidShieldplateA;
  G4LogicalVolume*   logicShieldplateA;
  G4VPhysicalVolume* physiShieldplateA;

  G4Box*             solidShieldplateB;
  G4LogicalVolume*   logicShieldplateB;
  G4VPhysicalVolume* physiShieldplateB;

  G4Box*             solidShieldplateC;
  G4LogicalVolume*   logicShieldplateC;
  G4VPhysicalVolume* physiShieldplateC;

  G4Box*             solidShieldplateD;
  G4LogicalVolume*   logicShieldplateD;
  G4VPhysicalVolume* physiShieldplateD;

};

#endif

