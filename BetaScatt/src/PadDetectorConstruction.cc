#include "PadDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "PadCentralData.hh"
#include "G4ios.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include <fstream>

using namespace std;

PadDetectorConstruction::PadDetectorConstruction()
  :  solidWorld(0), logicWorld(0), physiWorld(0), 
     solidVirtualBox(0), logicVirtualBox(0), physiVirtualBox(0),
     solidFilm(0), logicFilm(0), physiFilm(0),
     solidSubstrate(0), logicSubstrate(0), physiSubstrate(0),
     solidSampleHolder(0), logicSampleHolder(0), physiSampleHolder(0),
     solidShield(0), logicShield(0), physiShield(0),
     solidShield2(0), logicShield2(0), physiShield2(0),
     solidTube1(0), logicTube1(0), physiTube1(0),
     solidTube2(0), logicTube2(0), physiTube2(0),
     solidTube3(0), logicTube3(0), physiTube3(0),
     solidValve1(0), logicValve1(0), physiValve1(0),
     solidValve2(0), logicValve2(0), physiValve2(0),
     solidTube4(0), logicTube4(0), physiTube4(0),
     solidTube5(0), logicTube5(0), physiTube5(0),
     solidTube6(0), logicTube6(0), physiTube6(0),
     solidTube7(0), logicTube7(0), physiTube7(0),
     solidTube8(0), logicTube8(0), physiTube8(0),
     solidTube9(0), logicTube9(0), physiTube9(0),
     solidTube10(0), logicTube10(0), physiTube10(0),
     solidCone(0), logicCone(0), physiCone(0),
     solidCone2(0), logicCone2(0), physiCone2(0),
     solidCollimator(0), logicCollimator(0), physiCollimator(0),
     solidDetector1(0), logicDetector1(0), physiDetector1(0),
     solidDetector2(0), logicDetector2(0), physiDetector2(0),
     solidDetector3(0), logicDetector3(0), physiDetector3(0),
	 solidVacuumbox(0), logicVacuumbox(0), physiVacuumbox(0)
{
dataPointer = PadCentralData::GetDataPointer();
}    

PadDetectorConstruction::~PadDetectorConstruction()
{;}

G4VPhysicalVolume* PadDetectorConstruction::Construct()
{
  fstream& logFile = dataPointer->GetOutputStream();

  //------------------------------------------------------------
  //                         Construct Materials
  //------------------------------------------------------------

  G4double a; // atomic mass
  G4double z; // atomic number
  G4double density;
  G4String name, symbol;
  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  a = 1.01*g/mole;
  z = 1.;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z, a);

  a = 9.012182*g/mole;
  z = 4.;
  G4Element* elBe = new G4Element(name="Beryllium", symbol="Be", z, a);

  a = 12.01*g/mole;
  z = 6.;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z, a);

  a = 14.00674*g/mole;
  z = 7.;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z, a);

  a = 15.999*g/mole;
  z = 8.;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z, a);

  a = 26.982*g/mole;
  z = 13.;
  G4Element* elAl = new G4Element(name="Aluminum", symbol="Al", z, a);

  a = 28.086*g/mole;
  z = 14.;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z, a);

  a = 30.973762*g/mole;
  z = 15.;
  G4Element* elP = new G4Element(name="Phosphorus", symbol="P", z, a);

  a = 47.867*g/mole;
  z = 22.;
  G4Element* elTi = new G4Element(name="Titanium", symbol="Ti", z, a);

  a = 51.996*g/mole;
  z = 24.;
  G4Element* elCr = new G4Element(name="Chromium", symbol="Cr", z, a);

  a = 54.938*g/mole;
  z = 25.;
  G4Element* elMn = new G4Element(name="Manganese", symbol="Mn", z, a);

  a = 55.847*g/mole;
  z = 26.;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z, a);

  a = 58.693*g/mole;
  z = 28.;
  G4Element* elNi = new G4Element(name="Nickel", symbol="Ni", z, a);

  a = 63.546*g/mole;
  z = 29.;
  G4Element* elCu = new G4Element(name="Copper", symbol="Cu", z, a);

  a = 65.390*g/mole;
  z = 30.;
  G4Element* elZn = new G4Element(name="Zinc", symbol="Zn", z, a);

  a = 69.723*g/mole;
  z = 31.;
  G4Element* elGa = new G4Element(name="Gallium", symbol="Ga", z, a);

  a = 72.61*g/mole;
  z = 32.;
  G4Element* elGe = new G4Element(name="Germanium", symbol="Ge", z, a);

  a = 74.9216*g/mole;
  z = 33.;
  G4Element* elAs = new G4Element(name="Arsenic", symbol="As", z, a);

  a = 87.62*g/mole;
  z = 38.;
  G4Element* elSr = new G4Element(name="Strontium", symbol="Sr", z, a);

  a = 95.941*g/mole;
  z = 42.;
  G4Element* elMo = new G4Element(name="Molybdenum", symbol="Mo", z, a);

  a = 107.87*g/mole;
  z = 47.;
  G4Element* elAg = new G4Element(name="Silver", symbol="Ag", z, a);

  a = 114.818*g/mole;
  z = 49;
  G4Element* elIn = new G4Element(name="Indium", symbol="In", z, a);

  a = 180.9479*g/mole;
  z = 73.;
  G4Element* elTa = new G4Element(name="Tantalum", symbol="Ta", z, a);

  a = 196.967*g/mole;
  z = 79.;
  G4Element* elAu = new G4Element(name="Gold", symbol="Au", z, a);


  density = 2.700*g/cm3;
  G4Material* Al = new G4Material(name="Aluminum", density, ncomponents=1);
  Al->AddElement(elAl, natoms=1);

  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", density, ncomponents=1);
  Si->AddElement(elSi, natoms=1);

  density = 7.874*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", density, ncomponents=1);
  Fe->AddElement(elFe, natoms=1);

  density = 8.920*g/cm3;
  G4Material* Cu = new G4Material(name="Copper", density, ncomponents=1);
  Cu->AddElement(elCu, natoms=1);

  density = 10.280*g/cm3;
  G4Material* Mo = new G4Material(name="Molybdenum", density, ncomponents=1);
  Mo->AddElement(elMo, natoms=1);
 
  density = 10.490*g/cm3;
  G4Material* Ag = new G4Material(name="Silver", density, ncomponents=1);
  Ag->AddElement(elAg, natoms=1);

  density = 16.650*g/cm3;
  G4Material* Ta = new G4Material(name="Tantalum", density, ncomponents=1);
  Ta->AddElement(elTa, natoms=1);

  density = 19.300*g/cm3;
  G4Material* Au = new G4Material(name="Gold", density, ncomponents=1);
  Au->AddElement(elAu, natoms=1);

  density = 8.06*g/cm3;
  G4Material* Steel = new G4Material(name="Steel", density, ncomponents=6);
  Steel->AddElement(elC, fractionmass=0.001);
  Steel->AddElement(elSi, fractionmass=0.007);
  Steel->AddElement(elCr, fractionmass=0.18);
  Steel->AddElement(elMn, fractionmass=0.01);
  Steel->AddElement(elFe, fractionmass=0.712);
  Steel->AddElement(elNi, fractionmass=0.09);

  density = 1.370*g/cm3;
  G4Material* Mylar = new G4Material(name="Mylar", density, ncomponents=6);
  Mylar->AddElement(elH, fractionmass=0.0419);
  Mylar->AddElement(elC, fractionmass=0.6250);
  Mylar->AddElement(elO, fractionmass=0.3330);

  G4double const mbar = 1.e-03 * bar;
  temperature = 293.*kelvin;
  pressure = 1.e-06*mbar;
  density = 1.29e-03*(g/cm3) * pressure/(1.*atmosphere);
  G4Material* Vacuum = new G4Material(name="Vacuum", density, ncomponents=2,
                                   kStateGas, temperature, pressure);
  Vacuum->AddElement(elN, fractionmass=0.7);
  Vacuum->AddElement(elO, fractionmass=0.3);

  density = 1.*g/cm3;
  G4Material* Water = new G4Material(name="Water", density, ncomponents=2);
  Water->AddElement(elH, natoms=2);
  Water->AddElement(elO, natoms=1);

  density = 6.81*g/cm3;
  G4Material* InN = new G4Material(name="Indium Nitride", density, ncomponents=2);
  InN->AddElement(elIn, natoms=1);
  InN->AddElement(elN, natoms=1);

  density = 6.15*g/cm3;
  G4Material* GaN = new G4Material(name="Gallium Nitride", density, ncomponents=2);
  GaN->AddElement(elGa, natoms=1);
  GaN->AddElement(elN, natoms=1);

  density = 4.81*g/cm3;
  G4Material * InP = new G4Material (name="Indium Phosphide", density, ncomponents=2);
  InP->AddElement(elIn, natoms=1);
  InP->AddElement(elP, natoms=1);

  density = 3.23*g/cm3;
  G4Material* AlN = new G4Material(name="Aluminum Nitride", density, ncomponents=2);
  AlN->AddElement(elAl, natoms=1);
  AlN->AddElement(elN, natoms=1);

  density = 5.657*g/cm3;
  G4Material* ZnO = new G4Material(name="Zinc Oxide", density, ncomponents=2);
  ZnO->AddElement(elZn, natoms=1);
  ZnO->AddElement(elO, natoms=1);

  density = 3.98*g/cm3;
  G4Material* Sapphire = new G4Material(name="Sapphire", density, ncomponents=2);
  Sapphire->AddElement(elO, natoms=3);
  Sapphire->AddElement(elAl, natoms=2);

  density = 5.323*g/cm3;
  G4Material* Ge = new G4Material(name="Germanium", density, ncomponents=1);
  Ge->AddElement(elGe, natoms=1);

  density = 5.117*g/cm3;
  G4Material* SrTiO3 = new G4Material(name="Strontium titanate", density, ncomponents=3);
  SrTiO3->AddElement(elSr, natoms=1);
  SrTiO3->AddElement(elTi, natoms=1);
  SrTiO3->AddElement(elO, natoms=3);
  
  density = 8.96*g/cm3;
  G4Material* CuBe = new G4Material(name="Copper Beryllium", density, ncomponents=2);
  CuBe->AddElement(elCu, fractionmass=0.97);
  CuBe->AddElement(elBe, fractionmass=0.03);

  density = 5.316*g/cm3;
  G4Material* GaAs = new G4Material(name="Gallium Arsenide", density, ncomponents=2);
  GaAs->AddElement(elGa, natoms=1);
  GaAs->AddElement(elAs, natoms=1);

  density = 3.210*g/cm3;
  G4Material* SiC = new G4Material(name="Silicon Carbide", density, ncomponents=2);
  SiC->AddElement(elSi, natoms=1);
  SiC->AddElement(elC, natoms=1);


  //------------------------------------------------------------
  //           Choose Material for Film and Substrate
  //------------------------------------------------------------

  G4String FilmMatString = dataPointer->GetFilmMaterial();
  G4Material* FilmMaterial;

  if (FilmMatString == "GaN") {
    FilmMaterial = GaN;
  } else if (FilmMatString == "Si") {
    FilmMaterial = Si;
  } else if (FilmMatString == "ZnO") {
    FilmMaterial = ZnO;
  } else if (FilmMatString == "AlN") {
    FilmMaterial = AlN;
  } else if (FilmMatString == "InN") {
    FilmMaterial = InN;
  } else if (FilmMatString == "Ge") {
    FilmMaterial = Ge;
  } else if (FilmMatString == "Sapphire") {
    FilmMaterial = Sapphire;
  } else if (FilmMatString == "GaAs") {
    FilmMaterial = GaAs;
  } else if (FilmMatString == "InP"){
    FilmMaterial = InP;
  } else if (FilmMatString == "SrTiO3"){
    FilmMaterial = SrTiO3;
  } else if (FilmMatString == "SiC"){
    FilmMaterial = SiC;
  } else {
    logFile << "Material " << FilmMatString << " not implemented yet." 
	    << G4endl << "Ask responsable person (Bart De Vries) " 
	    << "to add this material" << G4endl<< "For the moment "
	    << "falling back to standard Si" << G4endl;
    FilmMaterial = Si;
  };


  G4String SubstrateMatString = dataPointer->GetSubstrateMaterial();
  G4Material* SubstrateMaterial;

  if (SubstrateMatString == "Sapphire") {
    SubstrateMaterial = Sapphire;
  } else if (SubstrateMatString == "ZnO") {
    SubstrateMaterial = ZnO;
  } else if (SubstrateMatString == "Si") {
    SubstrateMaterial = Si;
  } else if (SubstrateMatString == "Ge") {
    SubstrateMaterial = Ge;
  } else if (SubstrateMatString == "GaAs") {
    SubstrateMaterial = GaAs;
  } else if (SubstrateMatString == "InP"){
    SubstrateMaterial = InP;
  } else if (SubstrateMatString == "GaN"){
    SubstrateMaterial = GaN;
  } else if (SubstrateMatString == "AlN"){
    SubstrateMaterial = AlN;
  } else if (SubstrateMatString == "InN"){
    SubstrateMaterial = InN;
  }else {
    logFile << "Material " << SubstrateMatString << " not implemented yet." 
	    << G4endl << "Ask responsable person (Stefan Decoster) " 
	    << "to add this material" << G4endl << "For the moment "
	    << "falling back to standard Si" << G4endl;
    SubstrateMaterial = Si;
  };


  //------------------------------------------------------------
  //         Attributes needed for color in visualization
  //------------------------------------------------------------

  G4VisAttributes* whiteVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
  //G4VisAttributes* greyVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  //G4VisAttributes* blackVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.));
  G4VisAttributes* redVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
  G4VisAttributes* greenVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
  G4VisAttributes* blueVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
  G4VisAttributes* cyanVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.));
  G4VisAttributes* magentaVisAtt = new G4VisAttributes(G4Colour(1.,0.,1.));
  G4VisAttributes* yellowVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));



  //------------------------------------------------------------
  //                         Construct Volumes
  //------------------------------------------------------------


  G4double HalfX, HalfY, HalfZ;
  G4double InnerRadius, OuterRadius, StartingAngle, SegmentAngle;
  G4double InnerRadius1, OuterRadius1, InnerRadius2, OuterRadius2;
  G4ThreeVector Position;
  
  G4double theta = dataPointer->GetThetaRotation();
  G4double phi = dataPointer->GetPhiRotation();

  G4double FilmThickness = dataPointer->GetFilmThickness();
  G4double SampleThickness = dataPointer->GetSampleThickness();

  

  //--------------------
  // Virtual Box Rotation Matrix
  //--------------------

  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateY(theta);
  rot->rotateZ(-phi);
  logFile << "Rotation matrix for sample holder (defined by angles theta and phi):"
	  << G4endl
	  << rot->xx() << " " << rot->xy() << " " << rot->xz() << G4endl
	  << rot->yx() << " " << rot->yy() << " " << rot->yz() << G4endl 
	  << rot->zx() << " " << rot->zy() << " " << rot->zz() << G4endl 
	  << G4endl;


  //--------------------
  // World
  //--------------------

  HalfX = 17.*0.5*cm;
  HalfY = 17.*0.5*cm;
  HalfZ = 140.*0.5*cm;

  Position = G4ThreeVector(0.,0.,0.);
  
  solidWorld = new G4Box("world",HalfX,HalfY,HalfZ);
  logicWorld = new G4LogicalVolume(solidWorld, Vacuum,
				   "World", 0, 0, 0);
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  physiWorld = new G4PVPlacement(0, Position, "World",
				 logicWorld, 0, false, 0);

  //--------------------
  // Virtual Sample Box
  //--------------------
  
  HalfX = 5.*0.5*cm;
  HalfY = 5.*0.5*cm;
  HalfZ = 1.*0.5*cm;
  
  Position = G4ThreeVector(0.,0.,0.);
  
  solidVirtualBox = new G4Box("virtual box",HalfX,HalfY,HalfZ);
  logicVirtualBox = new G4LogicalVolume(solidVirtualBox, Vacuum,
					"VirtualBox", 0, 0, 0);
  logicVirtualBox->SetVisAttributes(G4VisAttributes::Invisible);
  physiVirtualBox = new G4PVPlacement(rot, Position, "VirtualBox",
				      logicVirtualBox, physiWorld, false, 0);


				        
  //--------------------
  // Ta Holder Parts 1
  //--------------------
  
  HalfX = 1.*0.5*cm;
  HalfY = 0.5*0.5*cm;
  HalfZ = 0.02*0.5*cm;

  G4RotationMatrix* rotTa1 = new G4RotationMatrix();
  rotTa1->rotateX(-30.*deg);
  
  G4double heightTa = 0.25*cm * sin(30.*deg)+ 0.01*cm * cos(30.*deg);
  Position = G4ThreeVector(0.,0.4*cm,heightTa);
  
  solidTaPart1 = new G4Box("tapart1",HalfX,HalfY,HalfZ);
  logicTaPart1 = new G4LogicalVolume(solidTaPart1, Ta,
					"TaHolderPart1", 0, 0, 0);
  logicTaPart1->SetVisAttributes(whiteVisAtt);
  physiTaPart1 = new G4PVPlacement(rotTa1, Position, "TaHolderPart1",
				      logicTaPart1, physiVirtualBox, false, 0);
  
  //--------------------
  // Ta Holder Parts 2
  //--------------------
  
  HalfX = 1.*0.5*cm;
  HalfY = 0.5*0.5*cm;
  HalfZ = 0.02*0.5*cm;

  G4RotationMatrix* rotTa2 = new G4RotationMatrix();
  rotTa2->rotateX(30.*deg);
  
  Position = G4ThreeVector(0.,-0.4*cm,heightTa);
  
  solidTaPart2 = new G4Box("tapart2",HalfX,HalfY,HalfZ);
  logicTaPart2 = new G4LogicalVolume(solidTaPart2, Ta,
					"TaHolderPart2", 0, 0, 0);
  logicTaPart2->SetVisAttributes(whiteVisAtt);
  physiTaPart2 = new G4PVPlacement(rotTa2, Position, "TaHolderPart2",
				      logicTaPart2, physiVirtualBox, false, 0);
  

  //--------------------
  // Thin film
  //--------------------
  
  HalfX = 0.5*0.5*cm;
  HalfY = 0.5*0.5*cm;
  HalfZ = FilmThickness*0.5;
  
  Position = G4ThreeVector(0.,0.,(-FilmThickness)*0.5);
  
  solidFilm = new G4Box("film",HalfX,HalfY,HalfZ);
  logicFilm = new G4LogicalVolume(solidFilm, FilmMaterial,
				  "Film", 0, 0, 0);
  logicFilm->SetVisAttributes(redVisAtt);

  // If MaxFilmStep < 0 then no special step cut is applied otherwise 
  // the maximum step is limited to MaxFilmStep
  G4double maxStep = dataPointer->GetMaxFilmStep();
  if (maxStep > 0) logicFilm->SetUserLimits(new G4UserLimits(maxStep));

  physiFilm = new G4PVPlacement(0, Position, "Film",
				logicFilm, physiVirtualBox, false, 0);
  
  
  //--------------------
  // Substrate
  //--------------------
  
  HalfX = 0.5*0.5*cm;
  HalfY = 0.5*0.5*cm;
  HalfZ = (SampleThickness-FilmThickness)*0.5;
  
  Position = G4ThreeVector(0.,0.,-(SampleThickness+FilmThickness)*0.5);
  
  solidSubstrate = new G4Box("substrate",HalfX,HalfY,HalfZ);
  logicSubstrate = new G4LogicalVolume(solidSubstrate, SubstrateMaterial,
				       "Substrate", 0, 0, 0);
  logicSubstrate->SetVisAttributes(blueVisAtt);
  physiSubstrate = new G4PVPlacement(0, Position, "Substrate",
				     logicSubstrate, physiVirtualBox, false, 0);
  
  
  //--------------------
  // Mo sample holder
  //--------------------
  
  HalfX = 1.58*0.5*cm;
  HalfY = 4.*0.5*cm;
  HalfZ = 0.1*0.5*cm;
  
  Position = G4ThreeVector(0.,0.,(-0.05*cm-SampleThickness));
  
  solidSampleHolder = new G4Box("sampleholder",HalfX,HalfY,HalfZ);
  logicSampleHolder = new G4LogicalVolume(solidSampleHolder, Mo,
					  "SampleHolder", 0, 0, 0);
  logicSampleHolder->SetVisAttributes(blueVisAtt);
  physiSampleHolder = new G4PVPlacement(0, Position, "SampleHolder",
					logicSampleHolder, physiVirtualBox, false, 0);

  
  //----------------------------------------------------------
  //                  PAD1 Setup only
  //----------------------------------------------------------

  
  if (dataPointer->GetPadSetup() == 1) {
    
        
    //---------------------
    // Vacuum tube 1 (pad1)
    //---------------------
    
    InnerRadius = 14.63*0.5*cm;
    OuterRadius = (14.63*0.5 + 0.304)*cm; 
    HalfZ = 9.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,12.6*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(cyanVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);
    
    
    
    //---------------------
    // Vacuum tube 2 (pad1)
    //---------------------
    
    InnerRadius = 5.89*0.5*cm;
    OuterRadius = (14.63*0.5 + 0.304)*cm; 
    HalfZ = 2.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,18.3*cm);
    
    solidTube2 = new G4Tubs("tube2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube2 = new G4LogicalVolume(solidTube2, Steel,
				     "Tube2", 0, 0, 0);
    logicTube2->SetVisAttributes(magentaVisAtt);
    physiTube2 = new G4PVPlacement(0, Position, "Tube2",
				   logicTube2, physiWorld, false, 0);
    
    
    
    //---------------------
    // Vacuum tube 3 (pad1)
    //---------------------
    
    InnerRadius = 5.89*0.5*cm;
    OuterRadius = (5.89*0.5 + 0.163)*cm; 
    HalfZ = 2.5*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,20.65*cm);
    
    solidTube3 = new G4Tubs("tube3",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube3 = new G4LogicalVolume(solidTube3, Steel,
				     "Tube3", 0, 0, 0);
    logicTube3->SetVisAttributes(yellowVisAtt);
    physiTube3 = new G4PVPlacement(0, Position, "Tube3",
				   logicTube3, physiWorld, false, 0);
    
    
    
    //---------------------
    // Al collimator (pad1)
    //---------------------
    
    InnerRadius = 5.89*0.5*cm;
    OuterRadius = (5.89*0.5 + 4.)*cm; 
    HalfZ = 2.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,25.85*cm);
    
    solidCollimator = new G4Tubs("collimator",InnerRadius,OuterRadius,HalfZ,
				 StartingAngle,SegmentAngle);
    logicCollimator = new G4LogicalVolume(solidCollimator, Al,
					  "Collimator", 0, 0, 0);
    logicCollimator->SetVisAttributes(whiteVisAtt);
    physiCollimator = new G4PVPlacement(0, Position, "Collimator",
					logicCollimator, physiWorld, false, 0);
    
    
    
    //-----------------------------
    // Si detector (22 x 22) (pad1)
    //-----------------------------
    
    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,28.625*cm);
    
    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);
    
    
    
    //---------------------------
    // Si detector (3 x 3) (pad1)
    //---------------------------
    
    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);
    
    
    
    //---------------------------
    // Si detector (1 x 1) (pad1)
    //---------------------------
    
    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);
    

  //----------------------------------------------------------
  //                  PAD2 Setup only
  //----------------------------------------------------------


  } else if (dataPointer->GetPadSetup() == 2) {

 
    //-------------------------------
    // Al contamination shield (pad2)
    //-------------------------------
    
    InnerRadius = 3.5*0.5*cm;
    OuterRadius = 10*0.5*cm; 
    HalfZ = 0.02*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,5*cm);
    
    solidShield = new G4Tubs("shield",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicShield = new G4LogicalVolume(solidShield, Al,
				      "Shield", 0, 0, 0);
    logicShield->SetVisAttributes(whiteVisAtt);
    physiShield = new G4PVPlacement(0, Position, "Shield",
				    logicShield, physiWorld, false, 0);
    
    //------------------------
    // Vacuum tube cone (pad2)
    //------------------------
    
    G4RotationMatrix* rotcone = new G4RotationMatrix();
    rotcone->rotateY(18*deg);

    InnerRadius1 = 10*0.5*cm;
    OuterRadius1 = (10*0.5 + 0.163)*cm;
    InnerRadius2 = 15*0.5*cm;
    OuterRadius2 = (15*0.5 + 0.163)*cm;
    HalfZ = 5.*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(-2.63*cm,0.,8.08*cm);
    
    solidCone = new G4Cons("cone",InnerRadius1,OuterRadius1,
			    InnerRadius2,OuterRadius2,HalfZ,
			    StartingAngle,SegmentAngle);
    logicCone = new G4LogicalVolume(solidCone, Steel,
				     "Cone", 0, 0, 0);
    logicCone->SetVisAttributes(cyanVisAtt);
    physiCone = new G4PVPlacement(rotcone, Position, "Cone",
				   logicCone, physiWorld, false, 0);
    
  
    //---------------------
    // Vacuum tube 1 (pad2)
    //---------------------
    
    InnerRadius = 6.3*0.5*cm;
    OuterRadius = (6.3*0.5 + 0.163)*cm; 
    HalfZ = 10*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,16.4*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(magentaVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);


    //---------------------
    // Vacuum tube 2 (pad2)
    //---------------------
    
    InnerRadius = 6.*0.5*cm;
    OuterRadius = (6.*0.5 + 0.313)*cm; 
    HalfZ = 6.*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,24.4*cm);
    
    solidTube2 = new G4Tubs("tube2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube2 = new G4LogicalVolume(solidTube2, Steel,
				     "Tube2", 0, 0, 0);
    logicTube2->SetVisAttributes(yellowVisAtt);
    physiTube2 = new G4PVPlacement(0, Position, "Tube2",
				   logicTube2, physiWorld, false, 0);
    
 
    
    //---------------------
    // Al collimator (pad2)
    //---------------------
    
    InnerRadius = 6*0.5*cm;
    OuterRadius = (6*0.5 + 2.)*cm; 
    HalfZ = 1.5*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,28.15*cm);
    
    solidCollimator = new G4Tubs("collimator",InnerRadius,OuterRadius,HalfZ,
				 StartingAngle,SegmentAngle);
    logicCollimator = new G4LogicalVolume(solidCollimator, Al,
					  "Collimator", 0, 0, 0);
    logicCollimator->SetVisAttributes(whiteVisAtt);
    physiCollimator = new G4PVPlacement(0, Position, "Collimator",
					logicCollimator, physiWorld, false, 0);
    
   


    //-----------------------------
    // Si detector (22 x 22) (pad2)
    //-----------------------------
    
    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,30.2*cm);
    
    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);
    
    
    
    //---------------------------
    // Si detector (3 x 3) (pad2)
    //---------------------------
    
    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);
    
    
    
    //---------------------------
    // Si detector (1 x 1) (pad2)
    //---------------------------
    
    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);
    
 

  //----------------------------------------------------------
  //                  PAD3 Setup only (PAD3 is old PAD4)
  //----------------------------------------------------------


  } else if (dataPointer->GetPadSetup() == 3) {
  
    //-------------------------------
    // Al contamination shield (pad3)
    //-------------------------------
    
    InnerRadius = 3.5*0.5*cm;
    OuterRadius = 10*0.5*cm; 
    HalfZ = 0.02*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,5*cm);
    
    solidShield = new G4Tubs("shield",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicShield = new G4LogicalVolume(solidShield, Al,
				      "Shield", 0, 0, 0);
    logicShield->SetVisAttributes(whiteVisAtt);
    physiShield = new G4PVPlacement(0, Position, "Shield",
				    logicShield, physiWorld, false, 0);
    

    //---------------------
    // Vacuum tube 1 (pad3)
    //---------------------
    
    InnerRadius = 10*0.5*cm;
    OuterRadius = (10*0.5 + 0.163)*cm; 
    HalfZ = 8.5*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,9.25*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(cyanVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);
    
    
    
    //---------------------
    // Vacuum tube 2 (pad3)
    //---------------------
    
    InnerRadius = 7*0.5*cm;
    OuterRadius = (10*0.5 + 0.163)*cm; 
    HalfZ = 2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,14.5*cm);
    
    solidTube2 = new G4Tubs("tube2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube2 = new G4LogicalVolume(solidTube2, Steel,
				     "Tube2", 0, 0, 0);
    logicTube2->SetVisAttributes(magentaVisAtt);
    physiTube2 = new G4PVPlacement(0, Position, "Tube2",
				   logicTube2, physiWorld, false, 0);
    
    
    
    //---------------------
    // Vacuum tube 3 (pad3)
    //---------------------
    
    InnerRadius = 7*0.5*cm;
    OuterRadius = (7*0.5 + 0.163)*cm; 
    HalfZ = 3.7*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,17.35*cm);
    
    solidTube3 = new G4Tubs("tube3",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube3 = new G4LogicalVolume(solidTube3, Steel,
				     "Tube3", 0, 0, 0);
    logicTube3->SetVisAttributes(yellowVisAtt);
    physiTube3 = new G4PVPlacement(0, Position, "Tube3",
				   logicTube3, physiWorld, false, 0);
    
    
    
    //------------------
    // Al valve 1 (pad3)
    //------------------
    
    InnerRadius = 6.5*0.5*cm;
    OuterRadius = 7*0.5*cm; 
    HalfZ = 0.7*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,19.55*cm);
    
    solidValve1 = new G4Tubs("valve1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve1 = new G4LogicalVolume(solidValve1, Al,
				      "Valve1", 0, 0, 0);
    logicValve1->SetVisAttributes(whiteVisAtt);
    physiValve1 = new G4PVPlacement(0, Position, "Valve1",
				    logicValve1, physiWorld, false, 0);
    
 
    //-----------------
    // Al valve 2(pad3)
    //-----------------
    
    InnerRadius = 6.5*0.5*cm;
    OuterRadius = 7*0.5*cm; 
    HalfZ = 0.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,22.1*cm);
    
    solidValve2 = new G4Tubs("valve2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve2 = new G4LogicalVolume(solidValve2, Al,
				      "Valve2", 0, 0, 0);
    logicValve2->SetVisAttributes(whiteVisAtt);
    physiValve2 = new G4PVPlacement(0, Position, "Valve2",
				    logicValve2, physiWorld, false, 0);
    

   //---------------------
    // Vacuum tube 4 (pad3)
    //---------------------
    
    InnerRadius = 6.8*0.5*cm;
    OuterRadius = (6.8*0.5 + 0.163)*cm; 
    HalfZ = 5.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,25.2*cm);
    
    solidTube4 = new G4Tubs("tube4",InnerRadius,OuterRadius,HalfZ,
	        	    StartingAngle,SegmentAngle);
    logicTube4 = new G4LogicalVolume(solidTube4, Steel,
				     "Tube4", 0, 0, 0);
    logicTube4->SetVisAttributes(cyanVisAtt);
    physiTube4 = new G4PVPlacement(0, Position, "Tube4",
				   logicTube4, physiWorld, false, 0);
   
  
    //---------------------
    // Al collimator (pad3)
    //---------------------
    
    InnerRadius = 6*0.5*cm;
    OuterRadius = (6.8*0.5 + 4.)*cm; 
    HalfZ = 1*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,28.5*cm);
    
    solidCollimator = new G4Tubs("collimator",InnerRadius,OuterRadius,HalfZ,
				 StartingAngle,SegmentAngle);
    logicCollimator = new G4LogicalVolume(solidCollimator, Al,
					  "Collimator", 0, 0, 0);
    logicCollimator->SetVisAttributes(magentaVisAtt);
    physiCollimator = new G4PVPlacement(0, Position, "Collimator",
					logicCollimator, physiWorld, false, 0);
    
    
    
    //-----------------------------
    // Si detector (22 x 22) (pad3)
    //-----------------------------
    
    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,30.895*cm);
    
    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);
    
    
    
    //---------------------------
    // Si detector (3 x 3) (pad3)
    //---------------------------
    
    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);
    
    
    
    //---------------------------
    // Si detector (1 x 1) (pad3)
    //---------------------------
    
    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);
    
    
/*
//COMMENTED BY LA
  //----------------------------------------------------------
  //                  PAD4 Setup only
  //----------------------------------------------------------
    
  }   else if (dataPointer->GetPadSetup() == 4) {  
  
    //-------------------------------
    // Al contamination shield (pad4)
    //-------------------------------
    
    InnerRadius = 3.5*0.5*cm;
    OuterRadius = 10.0*0.5*cm; 
    HalfZ = 0.02*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,5*cm);
    
    solidShield = new G4Tubs("shield",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicShield = new G4LogicalVolume(solidShield, Al,
				      "Shield", 0, 0, 0);
    logicShield->SetVisAttributes(whiteVisAtt);
    physiShield = new G4PVPlacement(0, Position, "Shield",
				    logicShield, physiWorld, false, 0);
    
    
    //---------------------
    // Vacuum tube 1 (pad4)
    //---------------------
    
    InnerRadius = 10.0*0.5*cm;
    OuterRadius = (10.0*0.5 + 0.163)*cm; 
    HalfZ = 19.9*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,14.95*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(cyanVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);
    
    
    //------------------
    // Al valve 1 (pad4)
    //------------------
    
    InnerRadius = 9.5*0.5*cm;
    OuterRadius = 11*0.5*cm; 
    HalfZ = 0.7*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,25.25*cm);
    
    solidValve1 = new G4Tubs("valve1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve1 = new G4LogicalVolume(solidValve1, Al,
				      "Valve1", 0, 0, 0);
    logicValve1->SetVisAttributes(whiteVisAtt);
    physiValve1 = new G4PVPlacement(0, Position, "Valve1",
				    logicValve1, physiWorld, false, 0);
    
 
    //-----------------
    // Al valve 2(pad4)
    //-----------------
    
    InnerRadius = 9.5*0.5*cm;
    OuterRadius = 11*0.5*cm; 
    HalfZ = 0.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,27.8*cm);
    
    solidValve2 = new G4Tubs("valve2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve2 = new G4LogicalVolume(solidValve2, Al,
				      "Valve2", 0, 0, 0);
    logicValve2->SetVisAttributes(whiteVisAtt);
    physiValve2 = new G4PVPlacement(0, Position, "Valve2",
				    logicValve2, physiWorld, false, 0);
    

    //---------------------
    // Al collimator (pad4)
    //---------------------
    
    InnerRadius = 10*0.5*cm;
    OuterRadius = (10*0.5 + 2.)*cm; 
    HalfZ = 1.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,28.7*cm);
    
    solidCollimator = new G4Tubs("collimator",InnerRadius,OuterRadius,HalfZ,
				 StartingAngle,SegmentAngle);
    logicCollimator = new G4LogicalVolume(solidCollimator, Al,
					  "Collimator", 0, 0, 0);
    logicCollimator->SetVisAttributes(magentaVisAtt);
    physiCollimator = new G4PVPlacement(0, Position, "Collimator",
					logicCollimator, physiWorld, false, 0);
    
    
    

    //-----------------------------
    // Si detector (22 x 22) (pad4)
    //-----------------------------
    
    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,30.895*cm);
    
    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);
    
    
    
    //---------------------------
    // Si detector (3 x 3) (pad4)
    //---------------------------
    
    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);
    
    
    
    //---------------------------
    // Si detector (1 x 1) (pad4)
    //---------------------------
    
    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);

*/
    
  //----------------------------------------------------------
  //                  PAD5 Setup only
  //----------------------------------------------------------
    
  }   else if (dataPointer->GetPadSetup() == 5) {  
  
    //-------------------------------
    // Cu contamination shield (pad5)
    //-------------------------------
    
    InnerRadius = 2.4*0.5*cm;
    OuterRadius = 10.0*0.5*cm; 
    HalfZ = 0.05*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,5.*cm);
    
    solidShield = new G4Tubs("shield",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicShield = new G4LogicalVolume(solidShield, Cu,
				      "Shield", 0, 0, 0);
    logicShield->SetVisAttributes(whiteVisAtt);
    physiShield = new G4PVPlacement(0, Position, "Shield",
				    logicShield, physiWorld, false, 0);
    
    
    //---------------------
    // Vacuum tube 1 (pad5)
    //---------------------
    
    InnerRadius = 6.67*0.5*cm;
    OuterRadius = (6.67*0.5 + 2.)*cm; 
    HalfZ = 3.5*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,8.75*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(cyanVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);
    
    //---------------------
    // Vacuum tube 2 (pad5)
    //---------------------
    
    InnerRadius = 6.0*0.5*cm;
    OuterRadius = (6.0*0.5 + 2.)*cm; 
    HalfZ = 10.06*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,15.58*cm);
    
    solidTube2 = new G4Tubs("tube2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube2 = new G4LogicalVolume(solidTube2, Steel,
				     "Tube2", 0, 0, 0);
    logicTube2->SetVisAttributes(cyanVisAtt);
    physiTube2 = new G4PVPlacement(0, Position, "Tube2",
				   logicTube2, physiWorld, false, 0);
    
    //---------------------
    // Vacuum tube 3 (pad5)
    //---------------------
    
    InnerRadius = 6.5*0.5*cm;
    OuterRadius = (6.5*0.5 + 2.)*cm; 
    HalfZ = 3.35*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,22.285*cm);
    
    solidTube3 = new G4Tubs("tube3",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube3 = new G4LogicalVolume(solidTube3, Steel,
				     "Tube3", 0, 0, 0);
    logicTube3->SetVisAttributes(cyanVisAtt);
    physiTube3 = new G4PVPlacement(0, Position, "Tube3",
				   logicTube3, physiWorld, false, 0);
        
    //------------------
    // Al valve 1 (pad5)
    //------------------
    
    InnerRadius = 6.5*0.5*cm;
    OuterRadius = (6.5*0.5 + 2.)*cm; 
    HalfZ = 0.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,24.26*cm);
    
    solidValve1 = new G4Tubs("valve1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve1 = new G4LogicalVolume(solidValve1, Al,
				      "Valve1", 0, 0, 0);
    logicValve1->SetVisAttributes(whiteVisAtt);
    physiValve1 = new G4PVPlacement(0, Position, "Valve1",
				    logicValve1, physiWorld, false, 0);
    
 
    //-----------------
    // Al valve 2(pad5)
    //-----------------
    
    InnerRadius = 6.5*0.5*cm;
    OuterRadius = (6.5*0.5 + 2.)*cm; 
    HalfZ = 0.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,26.86*cm);
    
    solidValve2 = new G4Tubs("valve2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve2 = new G4LogicalVolume(solidValve2, Al,
				      "Valve2", 0, 0, 0);
    logicValve2->SetVisAttributes(whiteVisAtt);
    physiValve2 = new G4PVPlacement(0, Position, "Valve2",
				    logicValve2, physiWorld, false, 0);
    
    //---------------------
    // Vacuum tube 4 (pad5)
    //---------------------
    
    InnerRadius = 6.5*0.5*cm;
    OuterRadius = (6.5*0.5 + 2.)*cm; 
    HalfZ = 4.76*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,29.54*cm);
    
    solidTube4 = new G4Tubs("tube4",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube4 = new G4LogicalVolume(solidTube4, Steel,
				     "Tube4", 0, 0, 0);
    logicTube4->SetVisAttributes(cyanVisAtt);
    physiTube4 = new G4PVPlacement(0, Position, "Tube4",
				   logicTube4, physiWorld, false, 0);
        
    //---------------------
    // Vacuum tube 5 (pad5)
    //---------------------

    InnerRadius = 6.2*0.5*cm;
    OuterRadius = (6.2*0.5 + 2.)*cm; 
    HalfZ = 0.9*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    Position = G4ThreeVector(0.,0.,32.37*cm);

    solidTube5 = new G4Tubs("tube5",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube5 = new G4LogicalVolume(solidTube5, Steel,
			     "Tube5", 0, 0, 0);
    logicTube5->SetVisAttributes(cyanVisAtt);
    physiTube5 = new G4PVPlacement(0, Position, "Tube5",
     		   logicTube5, physiWorld, false, 0);
        
    //---------------------
    // Al collimator (pad5)
    //---------------------
    
    InnerRadius = 6.2*0.5*cm;
    OuterRadius = (6.2*0.5 + 2.)*cm; 
    HalfZ = 0.5*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,33.07*cm);
    
    solidCollimator = new G4Tubs("collimator",InnerRadius,OuterRadius,HalfZ,
				 StartingAngle,SegmentAngle);
    logicCollimator = new G4LogicalVolume(solidCollimator, Al,
					  "Collimator", 0, 0, 0);
    logicCollimator->SetVisAttributes(magentaVisAtt);
    physiCollimator = new G4PVPlacement(0, Position, "Collimator",
					logicCollimator, physiWorld, false, 0);
    
    
    //-----------------------------
    // Si detector (22 x 22) (pad5)
    //-----------------------------
    
    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,35.81*cm);
    
    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);
    
    
    
    //---------------------------
    // Si detector (3 x 3) (pad5)
    //---------------------------
    
    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);
    
    
    
    //---------------------------
    // Si detector (1 x 1) (pad5)
    //---------------------------
    
    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.);
    
    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);

  //----------------------------------------------------------
  //              PAD6 Timepix-quad and Timepix Setup
  //----------------------------------------------------------
    
  }   else if (dataPointer->GetPadSetup() == 6   ||
		  	   dataPointer->GetPadSetup() == 512 ||
			   dataPointer->GetPadSetup() == 256) {
	// Setups in the EC-SLI chamber are equivalent for all  3 detectors

    //---------------------
    // Vacuum tube 1 (pad6)
    //---------------------
    
    InnerRadius = 5.2*0.5*cm;
    OuterRadius = 7.*cm; 
    HalfZ = 0.4*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,15.6*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(whiteVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 2 (pad6)
    //---------------------
    
    InnerRadius = 5.2*0.5*cm;
    OuterRadius = 5.8*0.5*cm; 
    HalfZ = 0.4*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,16.0*cm);
    
    solidTube2 = new G4Tubs("tube2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube2 = new G4LogicalVolume(solidTube2, Steel,
				     "Tube2", 0, 0, 0);
    logicTube2->SetVisAttributes(cyanVisAtt);
    physiTube2 = new G4PVPlacement(0, Position, "Tube2",
				   logicTube2, physiWorld, false, 0);

    //------------------------
    // Vacuum tube cone (pad6)
    //------------------------

    InnerRadius1 = 5.2*0.5*cm;
    OuterRadius1 = 5.8*0.5*cm;
    InnerRadius2 = 8.*0.5*cm;
    OuterRadius2 = 8.6*0.5*cm;
    HalfZ = 5.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0,0.,19.*cm);
    
    solidCone = new G4Cons("cone",InnerRadius1,OuterRadius1,
			    InnerRadius2,OuterRadius2,HalfZ,
			    StartingAngle,SegmentAngle);
    logicCone = new G4LogicalVolume(solidCone, Steel,
				     "Cone", 0, 0, 0);
    logicCone->SetVisAttributes(cyanVisAtt);
    physiCone = new G4PVPlacement(0, Position, "Cone",
				   logicCone, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 3 (pad6)
    //---------------------
    
    InnerRadius = 8.*0.5*cm;
    OuterRadius = 8.6*0.5*cm; 
    HalfZ = 1.4*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,22.5*cm);
    
    solidTube3 = new G4Tubs("tube3",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube3 = new G4LogicalVolume(solidTube3, Steel,
				     "Tube3", 0, 0, 0);
    logicTube3->SetVisAttributes(cyanVisAtt);
    physiTube3 = new G4PVPlacement(0, Position, "Tube3",
				   logicTube3, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 4 (pad6)
    //---------------------
    
    InnerRadius = 8.*0.5*cm;
    OuterRadius = 12.*0.5*cm; 
    HalfZ = 2.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,24.5*cm);
    
    solidTube4 = new G4Tubs("tube4",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube4 = new G4LogicalVolume(solidTube4, Steel,
				     "Tube4", 0, 0, 0);
    logicTube4->SetVisAttributes(cyanVisAtt);
    physiTube4 = new G4PVPlacement(0, Position, "Tube4",
				   logicTube4, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 5 (pad6)
    //---------------------
    
    InnerRadius = 8.*0.5*cm;
    OuterRadius = 12.*0.5*cm; 
    HalfZ = 1.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,29.6*cm);
    
    solidTube5 = new G4Tubs("tube5",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube5 = new G4LogicalVolume(solidTube5, Steel,
				     "Tube5", 0, 0, 0);
    logicTube5->SetVisAttributes(greenVisAtt);
    physiTube5 = new G4PVPlacement(0, Position, "Tube5",
				   logicTube5, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 6: nozzle tip (pad6)
    //---------------------
    
    InnerRadius = 0.8*0.5*cm;
    OuterRadius = 1.08*0.5*cm; 
    HalfZ = 1.28*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    G4RotationMatrix* rottube6 = new G4RotationMatrix();
    rottube6->rotateY(343*deg);
    
    Position = G4ThreeVector(1.474*cm,0.,4.820*cm); //fully screwed
//    Position = G4ThreeVector(1.403*cm,0.,4.8*cm); //5mm unscrewed
    
    solidTube6 = new G4Tubs("tube6",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube6 = new G4LogicalVolume(solidTube6, CuBe,
				     "Tube6", 0, 0, 0);
    logicTube6->SetVisAttributes(yellowVisAtt);
    physiTube6 = new G4PVPlacement(rottube6, Position, "Tube6",
				   logicTube6, physiWorld, false, 0);
    
    //---------------------
    // Vacuum cone 2: nozzle cone (pad6)
    //---------------------
        
    InnerRadius1 = 0.8*0.5*cm;
    OuterRadius1 = 1.08*0.5*cm;
    InnerRadius2 = 1.4*0.5*cm;
    OuterRadius2 = 2.*0.5*cm;
    HalfZ = 1.59*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    G4RotationMatrix* rotcone2 = new G4RotationMatrix();
    rotcone2->rotateY(343*deg);
    
    Position = G4ThreeVector(1.895*cm,0.,6.197*cm);
    
    solidCone2 = new G4Cons("cone2",InnerRadius1,OuterRadius1,
			    InnerRadius2,OuterRadius2,HalfZ,
			    StartingAngle,SegmentAngle);
    logicCone2 = new G4LogicalVolume(solidCone2, CuBe,
				     "Cone2", 0, 0, 0);
    logicCone2->SetVisAttributes(yellowVisAtt);
    physiCone2 = new G4PVPlacement(rotcone2, Position, "Cone2",
				   logicCone2, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 7: nozzle tube2 (pad6)
    //---------------------
    
    InnerRadius = 1.4*0.5*cm;
    OuterRadius = 2.*0.5*cm;
    HalfZ = 1.29*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    G4RotationMatrix* rottube7 = new G4RotationMatrix();
    rottube7->rotateY(343*deg);
    
    Position = G4ThreeVector(2.316*cm,0.,7.574*cm);
    
    solidTube7 = new G4Tubs("tube7",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube7 = new G4LogicalVolume(solidTube7, CuBe,
				     "Tube7", 0, 0, 0);
    logicTube7->SetVisAttributes(yellowVisAtt);
    physiTube7 = new G4PVPlacement(rottube7, Position, "Tube7",
				   logicTube7, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 8: nozzle tube3 (pad6)
    //---------------------
    
    InnerRadius = 1.4*0.5*cm;
    OuterRadius = 2.*0.5*cm;
    HalfZ = 8.7*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    G4RotationMatrix* rottube8 = new G4RotationMatrix();
    rottube8->rotateY(343*deg);
    
    Position = G4ThreeVector(3.775*cm,0.,12.346*cm);
    
    solidTube8 = new G4Tubs("tube8",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube8 = new G4LogicalVolume(solidTube8, Steel,
				     "Tube8", 0, 0, 0);
    logicTube8->SetVisAttributes(cyanVisAtt);
    physiTube8 = new G4PVPlacement(rottube8, Position, "Tube8",
				   logicTube8, physiWorld, false, 0);

    //-----------------------------
    // Si detector (22 x 22) (pad6)
    //-----------------------------
    
    if (dataPointer->GetPadSetup() == 6) {

		HalfX = 2.86*0.5*cm;
		HalfY = 2.86*0.5*cm;
		HalfZ = 0.05*0.5*cm;

		Position = G4ThreeVector(0.,0.,31.6*cm);

		solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
		logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
						 "Detector1", 0, 0, 0);
		logicDetector1->SetVisAttributes(greenVisAtt);
		physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
						   logicDetector1, physiWorld, false, 0);



		//---------------------------
		// Si detector (3 x 3) (pad6)
		//---------------------------

		HalfX = 0.39*0.5*cm;
		HalfY = 0.39*0.5*cm;
		HalfZ = 0.05*0.5*cm;

		Position = G4ThreeVector(0.,0.,0.*cm);

		solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
		logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
						 "Detector2", 0, 0, 0);
		logicDetector2->SetVisAttributes(greenVisAtt);
		physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
						   logicDetector2, physiDetector1, false, 0);



		//---------------------------
		// Si detector (1 x 1) (pad6)
		//---------------------------

		HalfX = 0.13*0.5*cm;
		HalfY = 0.13*0.5*cm;
		HalfZ = 0.05*0.5*cm;

		Position = G4ThreeVector(0.,0.,0.*cm);

		solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
		logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
						 "Detector3", 0, 0, 0);
		logicDetector3->SetVisAttributes(greenVisAtt);
		physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
						   logicDetector3, physiDetector2, false, 0);

	//----------------------------------------------------------
	//                  Timepix quad and 256 Setup
	//----------------------------------------------------------
    } else if (dataPointer->GetPadSetup() == 512 ||
    		   dataPointer->GetPadSetup() == 256) {

    	G4int fNofPixels = 0;
    	if (dataPointer->GetPadSetup() == 512){
    		fNofPixels = 512 + 4; //4 to take into account the larger central pixels
    	} else if (dataPointer->GetPadSetup() == 256)
    		fNofPixels = 256;

    	// Full detector
    	G4double detector_sizeZ = 300*um;
		G4double pixel_sizeXY = 55*um;
		G4double detector_sizeXY = fNofPixels * pixel_sizeXY;

		Position = G4ThreeVector(0.,0.,31.6*cm);

		solidDetector1 = new G4Box("Detector1",
				0.5*detector_sizeXY, 0.5*detector_sizeXY, 0.5*detector_sizeZ);
		logicDetector1 = new G4LogicalVolume(solidDetector1, Si, "Detector1");
		logicDetector1->SetVisAttributes(greenVisAtt);

		// Detector line

		solidDetector2 = new G4Box("Detector2",
				0.5*detector_sizeXY, 0.5*pixel_sizeXY, 0.5*detector_sizeZ);
		logicDetector2 = new G4LogicalVolume(solidDetector2, Si, "Detector2");

		// Detector pixel

		solidDetector3 = new G4Box("Detector3",
				0.5*pixel_sizeXY, 0.5*pixel_sizeXY, 0.5*detector_sizeZ);
		logicDetector3 = new G4LogicalVolume(solidDetector3, Si, "Detector3");

		//Pixel Replicas to form line

		new G4PVReplica("DetectorPixelReplica",
				logicDetector3, logicDetector2,
				kXAxis, fNofPixels, pixel_sizeXY);

		//Line replicas to form detector

		new G4PVReplica("DetectorLineReplica",
				logicDetector2, logicDetector1,
				kYAxis, fNofPixels, pixel_sizeXY);

		physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
										   logicDetector1, physiWorld, false, 0);

		// User Limits
		// Sets a max step length in the tracker region, with G4StepLimiter

		G4double maxStep = 1*um;
		logicDetector3->SetUserLimits(new G4UserLimits(maxStep));

    } // end of tpx but not (tpx + pad)

//BY LA
  //----------------------------------------------------------
  //                  PAD4 Setup + c saw structure
  //----------------------------------------------------------

  }   else if (dataPointer->GetPadSetup() == 4) {

    //-------------------------------
    // Al contamination shield (pad4 + saw)
    //-------------------------------

    InnerRadius = 3.5*0.5*cm;
    OuterRadius = 10.0*0.5*cm;
    HalfZ = 0.02*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    Position = G4ThreeVector(0.,0.,5*cm);

    solidShield = new G4Tubs("shield",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicShield = new G4LogicalVolume(solidShield, Al,
				      "Shield", 0, 0, 0);
    logicShield->SetVisAttributes(whiteVisAtt);
    physiShield = new G4PVPlacement(0, Position, "Shield",
				    logicShield, physiWorld, false, 0);


    //---------------------
    // Vacuum tube 1 (pad4 + saw)
    //---------------------

    InnerRadius = 10.0*0.5*cm;
    OuterRadius = (10.0*0.5 + 0.163)*cm;
    HalfZ = 9.4*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    Position = G4ThreeVector(0.,0.,9.675*cm);

    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(cyanVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);

    //---------------------
    // Vacuum tube 2 with saw atructure (pad4 + saw)
    //---------------------
	//TO TEST OUTSIDE THE INSIDE SAW STRUCTURE
/*	G4double zPlanes[] = 		{0.00*cm,     1.90*cm,     2.00*cm,     2.50*cm,     3.19*cm,     3.19*cm,     3.88*cm,     3.88*cm,     4.57*cm,     4.57*cm,     5.26*cm,     5.26*cm,     5.95*cm,     5.95*cm,     6.64*cm,     6.64*cm,     7.33*cm,     7.33*cm,     8.02*cm,     8.02*cm,     8.71*cm,     8.71*cm,     9.40*cm,     9.40*cm,     10.5*cm};
	G4double OuterRadiousA[] = 	{10.3*0.5*cm, 10.3*0.5*cm, 10.0*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.0*0.5*cm};
	G4double InnerRadiousA[] = 	{5.2*0.5*cm, 5.2*0.5*cm, 0.6*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 1.1*0.5*cm, 0.6*0.5*cm};*/


	G4double zPlanes[] = 		{0.00*cm,     1.90*cm,     2.00*cm,     2.50*cm,     3.19*cm,     3.19*cm,     3.88*cm,     3.88*cm,     4.57*cm,     4.57*cm,     5.26*cm,     5.26*cm,     5.95*cm,     5.95*cm,     6.64*cm,     6.64*cm,     7.33*cm,     7.33*cm,     8.02*cm,     8.02*cm,     8.71*cm,     8.71*cm,     9.40*cm,     9.40*cm,     10.5*cm};
	G4double InnerRadiousA[] = 	{10.3*0.5*cm, 10.3*0.5*cm, 10.0*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.5*0.5*cm, 10.0*0.5*cm, 10.0*0.5*cm};
	G4double OuterRadiousA[] = 	{15.2*0.5*cm, 15.2*0.5*cm, 10.6*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 11.1*0.5*cm, 10.6*0.5*cm};

	G4int numSegments = 25;

	StartingAngle = 0.*deg;
	SegmentAngle = 360.*deg;
	Position = G4ThreeVector(0.,0.,14.375*cm);
	G4Polycone * solidTubes2 = new G4Polycone ("tube2", StartingAngle, SegmentAngle, numSegments, zPlanes, InnerRadiousA, OuterRadiousA);
	logicTube2 = new G4LogicalVolume(solidTubes2, Steel, "Tube2", 0, 0, 0);
	logicTube2->SetVisAttributes(yellowVisAtt);
	physiTube2 = new G4PVPlacement(0, Position, "Tube2", logicTube2, physiWorld, false, 0);

    //------------------
    // Al valve 1 (pad4 + saw)
    //------------------

    InnerRadius = 9.5*0.5*cm;
    OuterRadius = 11*0.5*cm;
    HalfZ = 0.7*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    Position = G4ThreeVector(0.,0.,25.25*cm);

    solidValve1 = new G4Tubs("valve1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve1 = new G4LogicalVolume(solidValve1, Al,
				      "Valve1", 0, 0, 0);
    logicValve1->SetVisAttributes(whiteVisAtt);
    physiValve1 = new G4PVPlacement(0, Position, "Valve1",
				    logicValve1, physiWorld, false, 0);


    //-----------------
    // Al valve 2(pad4 + saw)
    //-----------------

    InnerRadius = 9.5*0.5*cm;
    OuterRadius = 11*0.5*cm;
    HalfZ = 0.6*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    Position = G4ThreeVector(0.,0.,27.8*cm);

    solidValve2 = new G4Tubs("valve2",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicValve2 = new G4LogicalVolume(solidValve2, Al,
				      "Valve2", 0, 0, 0);
    logicValve2->SetVisAttributes(whiteVisAtt);
    physiValve2 = new G4PVPlacement(0, Position, "Valve2",
				    logicValve2, physiWorld, false, 0);


    //---------------------
    // Al collimator (pad4 + saw)
    //---------------------

    InnerRadius = 10*0.5*cm;
    OuterRadius = (10*0.5 + 2.)*cm;
    HalfZ = 1.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    Position = G4ThreeVector(0.,0.,28.7*cm);

    solidCollimator = new G4Tubs("collimator",InnerRadius,OuterRadius,HalfZ,
				 StartingAngle,SegmentAngle);
    logicCollimator = new G4LogicalVolume(solidCollimator, Al,
					  "Collimator", 0, 0, 0);
    logicCollimator->SetVisAttributes(magentaVisAtt);
    physiCollimator = new G4PVPlacement(0, Position, "Collimator",
					logicCollimator, physiWorld, false, 0);




    //-----------------------------
    // Si detector (22 x 22) (pad4 + saw)
    //-----------------------------

    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;

    Position = G4ThreeVector(0.,0.,30.895*cm);

    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);



    //---------------------------
    // Si detector (3 x 3) (pad4 + saw)
    //---------------------------

    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;

    Position = G4ThreeVector(0.,0.,0.);

    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);



    //---------------------------
    // Si detector (1 x 1) (pad4 + saw)
    //---------------------------

    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;

    Position = G4ThreeVector(0.,0.,0.);

    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);


//end BY LA

//BY LA
  //----------------------------------------------------------
  //                  PAD6 Setup with longer distance Pad7
  //----------------------------------------------------------
    
  }   else if (dataPointer->GetPadSetup() == 7) {  
	//TA COLLIMATOR RADIUS
	G4double colR = 4.2*0.5*cm;
//	G4double colR = 6.0*0.5*cm;
//	G4double colR = 4.6*0.5*cm;
//	G4double colR = 3.8*0.5*cm;
//	G4double colR = 3.4*0.5*cm;
//	G4double colR = 3.0*0.5*cm;
//	G4double colR = 2.6*0.5*cm;

    //---------------------
    // Vacuum tube 1 : chamber wall (pad6 + long)
    //---------------------

        
    InnerRadius = 5.2*0.5*cm;
    OuterRadius = 7.*cm; 
    HalfZ = 0.4*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,15.6*cm);
    
    solidTube1 = new G4Tubs("tube1",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube1 = new G4LogicalVolume(solidTube1, Steel,
				     "Tube1", 0, 0, 0);
    logicTube1->SetVisAttributes(whiteVisAtt);
    physiTube1 = new G4PVPlacement(0, Position, "Tube1",
				   logicTube1, physiWorld, false, 0);

   //----------------------
   // Vacuum Tube 2 : steel Polygon until main tube (pad6 + long)
    //---------------------

	G4double zPlanes2[] = 		{0.00*cm,     0.40*cm,     0.40*cm,    6.37*cm,     6.37*cm,     7.40*cm,     7.40*cm,      14.4*cm,     14.4*cm,     14.7*cm};
	//WITH COLLIMATOR
	G4double InnerRadiousA2[] = 	{05.2*0.5*cm, 05.2*0.5*cm, 05.2*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, colR+0.6*cm,  colR+0.6*cm};
	//WITHOUT COLLIMATOR
//	G4double InnerRadiousA2[] = 	{05.2*0.5*cm, 05.2*0.5*cm, 05.2*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 07.9*0.5*cm, 7.9*0.5*cm,  7.9*0.5*cm};
	G4double OuterRadiousA2[] = 	{05.8*0.5*cm, 05.8*0.5*cm, 05.8*0.5*cm, 08.6*0.5*cm, 08.6*0.5*cm, 08.6*0.5*cm, 12.0*0.5*cm, 12.0*0.5*cm, 12.0*0.5*cm, 12.0*0.5*cm};
	G4int numSegments = 10;

	StartingAngle = 0.*deg;
	SegmentAngle = 360.*deg;
	
	Position = G4ThreeVector(0.,0.,15.46*cm);
	
	G4Polycone * solidTubes2 = new G4Polycone ("tube2", StartingAngle, 
			SegmentAngle, numSegments, zPlanes2, InnerRadiousA2, OuterRadiousA2);
	logicTube2 = new G4LogicalVolume(solidTubes2, Steel, "Tube2", 0, 0, 0);
	logicTube2->SetVisAttributes(cyanVisAtt);
	physiTube2 = new G4PVPlacement(0, Position, "Tube2", logicTube2, physiWorld, false, 0);

   //----------------------
    // Vacuum tube 3 : Ta colimator(pad6 + long)
    //---------------------

    InnerRadius = colR;
    OuterRadius = 7.64*0.5*cm; 
    HalfZ = 0.1*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,29.61*cm);
    
    solidTube3 = new G4Tubs("tube3",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube3 = new G4LogicalVolume(solidTube3, Ta,
				     "Tube3", 0, 0, 0);
    logicTube3->SetVisAttributes(whiteVisAtt);
    physiTube3 = new G4PVPlacement(0, Position, "Tube3",
				   logicTube3, physiWorld, false, 0);

   //----------------------
    // Vacuum tube 4 : Cu colimator(pad6 + long)
    //---------------------

    InnerRadius = colR+0.7*0.5*cm;
    OuterRadius = 7.64*0.5*cm; 
    HalfZ = 0.2*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,29.76*cm);
    
    solidTube4 = new G4Tubs("tube4",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube4 = new G4LogicalVolume(solidTube4, Cu,
				     "Tube3", 0, 0, 0);
    logicTube4->SetVisAttributes(yellowVisAtt);
    physiTube4 = new G4PVPlacement(0, Position, "Tube4",
				   logicTube4, physiWorld, false, 0);
	
   //----------------------
    // Vacuum tube 5 with vaccum entrance (pad6 + long)
    //---------------------
    //new long tube
    InnerRadius = 8.49*0.5*cm;
    OuterRadius = 8.69*0.5*cm; 
    HalfZ = 26.9*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,(29.6+26.9*0.5)*cm);

    //Main tube body
    G4Tubs * tubeTmp1 = new G4Tubs ("tubeTmp1", InnerRadius, OuterRadius, HalfZ, StartingAngle, SegmentAngle);
    //Tube to be subtracted to the main
    G4Tubs * tubeTmp2 = new G4Tubs ("tubeTmp2", 0., 2.74*0.5*cm, 5.*0.5*cm, StartingAngle, SegmentAngle);
    //vacuum tube
    G4Tubs * tubeTmp3 = new G4Tubs ("tubeTmp3", 2.54*0.5*cm, 2.74*0.5*cm, 3.22*0.5*cm, StartingAngle, SegmentAngle);
    //tube to be subtracted from the vacuum tube
//    G4Tubs * tubeTmp4 = new G4Tubs ("tubeTmp4", 0., 8.49*0.5*cm, 4.0*0.5*cm, StartingAngle, SegmentAngle);
    
    //orientation of the subtraction tube in relation to the main tube
    G4RotationMatrix* rottubeTmp2 = new G4RotationMatrix();
    rottubeTmp2->rotateY(90*deg);
    G4ThreeVector positionTmp2 (-(8.49+3.0-0.22)*0.5*cm, 0, 0);
//     G4RotationMatrix* rottubeTmp4 = new G4RotationMatrix();
//     rottubeTmp4->rotateY(0*deg);    
//     G4ThreeVector positionTmp4 (0,0,0);
    
    G4SubtractionSolid * tubeSub1 = new G4SubtractionSolid ("tubeSub1", tubeTmp1, tubeTmp2, rottubeTmp2, positionTmp2);
    G4UnionSolid * tube5 = new G4UnionSolid ("tube5", tubeSub1, tubeTmp3, rottubeTmp2, positionTmp2);
//    G4SubtractionSolid * tube5 = new G4SubtractionSolid ("tube5", tubeSum1, tubeTmp4, 0, 0);

    logicTube5 = new G4LogicalVolume(tube5, Steel,
				     "Tube5", 0, 0, 0);
    logicTube5->SetVisAttributes(greenVisAtt);
    physiTube5 = new G4PVPlacement(0, Position, "Tube5",
				   logicTube5, physiWorld, false, 0);

   //----------------------
    // Vacuum tube 6 flange tube->detector (pad6 + long)
    //---------------------
    //
    InnerRadius = 8.*0.5*cm;
    OuterRadius = 12.*0.5*cm; 
    HalfZ = 1.4*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,56.6*cm);
    
    solidTube6 = new G4Tubs("tube6",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube6 = new G4LogicalVolume(solidTube6, Steel,
				     "Tube6", 0, 0, 0);
    logicTube6->SetVisAttributes(greenVisAtt);
    physiTube6 = new G4PVPlacement(0, Position, "Tube6",
				   logicTube6, physiWorld, false, 0);
  
   //----------------------
    // Vacuum tube 7 Al Heat shield suport (pad6 + long)
    //---------------------
    //
    InnerRadius = 7.94*0.5*cm;
    OuterRadius = 9.0*0.5*cm; 
    HalfZ = 2.53*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,57.165*cm);
    
    solidTube7 = new G4Tubs("tube7",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube7 = new G4LogicalVolume(solidTube7, Al,
				     "Tube7", 0, 0, 0);
    logicTube7->SetVisAttributes(yellowVisAtt);
    physiTube7 = new G4PVPlacement(0, Position, "Tube7",
				   logicTube7, physiWorld, false, 0);

   //----------------------
    // Vacuum tube 10 Mylar shield (pad6 + long)
    //---------------------
    InnerRadius = 0.;
    OuterRadius = 7.94*0.5*cm; 
    HalfZ = 0.0002*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;
    
    Position = G4ThreeVector(0.,0.,58.43*cm);
    
    solidTube10 = new G4Tubs("tube10",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube10 = new G4LogicalVolume(solidTube10, Al,
				     "Tube10", 0, 0, 0);
    logicTube10->SetVisAttributes(yellowVisAtt);
    physiTube10 = new G4PVPlacement(0, Position, "Tube10",
				   logicTube10, physiWorld, false, 0);

   //----------------------
   // Vacuum Tube 8: nozzle (pad6 + long)
    //---------------------
 
	G4double zPlanes8[] = 		{0.00*cm,     1.28*cm,     1.28*cm,     2.87*cm,     2.87*cm,     4.16*cm};
	G4double InnerRadiousA8[] = 	{00.8*0.5*cm, 00.8*0.5*cm, 00.8*0.5*cm, 01.4*0.5*cm, 01.4*0.5*cm, 01.4*0.5*cm};
	G4double OuterRadiousA8[] = 	{1.08*0.5*cm, 1.08*0.5*cm, 1.08*0.5*cm, 02.0*0.5*cm, 02.0*0.5*cm, 02.0*0.5*cm};
	numSegments = 6;

	StartingAngle = 0.*deg;
	SegmentAngle = 360.*deg;
	
	G4RotationMatrix* rottube8 = new G4RotationMatrix();
	rottube8->rotateY(343*deg);
	
	Position = G4ThreeVector((1.474-0.187)*cm,0.,(4.820-0.612)*cm); //fully screwed
//	Position = G4ThreeVector(1.403*cm,0.,4.8*cm); //5mm unscrewed
    	
	G4Polycone * solidTubes8 = new G4Polycone ("tube8", StartingAngle, SegmentAngle, numSegments, zPlanes8, InnerRadiousA8, OuterRadiousA8);
	
	logicTube8 = new G4LogicalVolume(solidTubes8, CuBe, "Tube8", 0, 0, 0);
	
	logicTube8->SetVisAttributes(yellowVisAtt);
	physiTube8 = new G4PVPlacement(rottube8, Position, "Tube8", logicTube8, physiWorld, false, 0);
    
    //---------------------
    // Vacuum tube 9: nozzle tube (pad6 + long)
    //---------------------
    
    InnerRadius = 1.4*0.5*cm;
    OuterRadius = 2.*0.5*cm;
    HalfZ = 8.7*0.5*cm;
    StartingAngle = 0.*deg;
    SegmentAngle = 360.*deg;

    G4RotationMatrix* rottube9 = new G4RotationMatrix();
    rottube9->rotateY(343*deg);
    
    Position = G4ThreeVector(3.775*cm,0.,12.346*cm);
    
    solidTube9 = new G4Tubs("tube9",InnerRadius,OuterRadius,HalfZ,
			    StartingAngle,SegmentAngle);
    logicTube9 = new G4LogicalVolume(solidTube9, Steel,
				     "Tube9", 0, 0, 0);
    logicTube9->SetVisAttributes(cyanVisAtt);
    physiTube9 = new G4PVPlacement(rottube9, Position, "Tube9",
				   logicTube9, physiWorld, false, 0);

    //-----------------------------
    // Si detector (22 x 22) (pad6 + long)
    //-----------------------------
    
    HalfX = 2.86*0.5*cm;
    HalfY = 2.86*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,60.2*cm);
    
    solidDetector1 = new G4Box("detector1",HalfX,HalfY,HalfZ);
    logicDetector1 = new G4LogicalVolume(solidDetector1, Si,
					 "Detector1", 0, 0, 0);
    logicDetector1->SetVisAttributes(greenVisAtt);
    physiDetector1 = new G4PVPlacement(0, Position, "Detector1",
				       logicDetector1, physiWorld, false, 0);
    
    
    
    //---------------------------
    // Si detector (3 x 3) (pad6 + long)
    //---------------------------
    
    HalfX = 0.39*0.5*cm;
    HalfY = 0.39*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.*cm);
    
    solidDetector2 = new G4Box("detector2",HalfX,HalfY,HalfZ);
    logicDetector2 = new G4LogicalVolume(solidDetector2, Si,
					 "Detector2", 0, 0, 0);
    logicDetector2->SetVisAttributes(greenVisAtt);
    physiDetector2 = new G4PVPlacement(0, Position, "Detector2",
				       logicDetector2, physiDetector1, false, 0);
    
    
    
    //---------------------------
    // Si detector (1 x 1) (pad6 + long)
    //---------------------------
    
    HalfX = 0.13*0.5*cm;
    HalfY = 0.13*0.5*cm;
    HalfZ = 0.05*0.5*cm;
    
    Position = G4ThreeVector(0.,0.,0.*cm);
    
    solidDetector3 = new G4Box("detector3",HalfX,HalfY,HalfZ);
    logicDetector3 = new G4LogicalVolume(solidDetector3, Si,
					 "Detector3", 0, 0, 0);
    logicDetector3->SetVisAttributes(greenVisAtt);
    physiDetector3 = new G4PVPlacement(0, Position, "Detector3",
				       logicDetector3, physiDetector2, false, 0);

//end BY LA
  };

  return physiWorld;
}

void PadDetectorConstruction::ConstructSDandField()
{
	if (dataPointer->GetPadSetup() == 512 ||
		dataPointer->GetPadSetup() == 256) {
		//
		// Scorers
		//

		// declare Absorber as a MultiFunctionalDetector scorer
		//
		G4MultiFunctionalDetector* PixDetector
		= new G4MultiFunctionalDetector("Pixel");
		G4SDManager::GetSDMpointer()->AddNewDetector(PixDetector);

		G4VPrimitiveScorer* primitive;
		primitive = new PadPSEnergyDeposit("Edep");

		PixDetector->RegisterPrimitive(primitive);

		G4SDManager* SDman = G4SDManager::GetSDMpointer();
		SDman->AddNewDetector( PixDetector );

		logicDetector3->SetSensitiveDetector(PixDetector);

		//
		// Create global magnetic field messenger.
		// Uniform magnetic field is then created automatically if
		// the field value is not zero.
		//
		G4ThreeVector fieldValue = G4ThreeVector();
		fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
		fMagFieldMessenger->SetVerboseLevel(1);
		// register the field messenger for deleting
		G4AutoDelete::Register(fMagFieldMessenger);
	}
}
