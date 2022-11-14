//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN03DetectorConstruction.cc,v 1.19 2003/11/25 14:23:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03DetectorConstruction.hh"
#include "ExN03DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Orb.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ParticleTypes.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "G4NistManager.hh"

//#include "G4ParticleTable.hh"

//#include "Randomize.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorConstruction::ExN03DetectorConstruction():
	AbsorberMaterial(0),
	//GapMaterial(0),
	defaultMaterial(0),
	solidWorld(0),logicWorld(0),physiWorld(0),
	solidAbsorber(0),logicAbsorber(0),physiAbsorber(0),
	detectorMessenger(0),
	stepLimit(0),
	maxStep(0.0),
	sensitiveDetector(0),
	is_valid(false),
	allParticles(false),
	outParticle(0)
{
  
  // default parameter values of the calorimeter


 //diameters the of cell, cell core and nanoparticle
   
	//G4int			Nnp;	

//  ComputeCalorParameters();
  
  // materials
	G4NistManager* man = G4NistManager::Instance();
	man->SetVerbose(1);

  DefineMaterials();
 // SetAbsorberMaterial("EnrichedMagnetite");
 // SetGapMaterial("Water");
  
  // create commands for interactive definition of the calorimeter
  detectorMessenger = new ExN03DetectorMessenger(this);
	
  G4cout<<"Detector object created."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorConstruction::~ExN03DetectorConstruction()
{
	if (detectorMessenger)
		delete detectorMessenger;

	if (stepLimit)
		delete stepLimit;

	//sensitiveDetector is destroyed automatically
		
	G4cout<<"Detector object dectroyed."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN03DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::DefineMaterials()
{ 
//This function illustrates the possible ways to define materials
 
	G4String name, symbol;
	G4double a, z, density;      //a=mass of a mole; 
                             //z=mean number of protons;  
														 
 G4int ncomponents, natoms;
 G4double abundance;

//
// define Elements
//

G4Element* H  = new G4Element("Hydrogen", symbol="H",  z= 1., a=  1.01*g/mole);
//G4Element* C  = new G4Element("Carbon"  , symbol="C",  z= 6., a= 12.01*g/mole);
//G4Element* N  = new G4Element("Nitrogen", symbol="N",  z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  , symbol="O",  z= 8., a= 15.999*g/mole);
G4Element* Fe = new G4Element("Iron",     symbol="Fe", z=26., a= 55.845*g/mole);
//G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);
//G4Element* Cu = new G4Element("Copper", symbol="Cu"  , z= 29., a= 63.55*g/mole);
//G4Element* Mo = new G4Element("Molybdenum", symbol="Mo"  , z= 42., a= 95.94*g/mole);

//
// define an Element from isotopes, by relative abundance 
//

//
// define simple materials
//

//new G4Material("Beryllium", z= 4., a=  9.012*g/mole, density=  1.848*g/cm3);
//new G4Material("Aluminium", z=13., a=  26.98*g/mole, density=  2.699*g/cm3);
//new G4Material("Chromium" , z=24., a=  52.00*g/mole, density=  7.19*g/cm3);
//new G4Material("Iron"     , z=26., a=  55.845*g/mole, density=  7.874*g/cm3); //X-Ray data booklet
//new G4Material("Copper"   , z=29., a=  63.55*g/mole, density=  8.96*g/cm3);
//new G4Material("Zinc"     , z=30., a=  65.39*g/mole, density=  7.13*g/cm3);
//new G4Material("Zirconium", z=40., a=  91.22*g/mole, density=  6.45*g/cm3);
//new G4Material("Niobium"  , z=41., a=  92.91*g/mole, density=  8.57*g/cm3);
//new G4Material("Molybdenum",z=42., a=  95.94*g/mole, density= 10.22*g/cm3);
//new G4Material("Rhodium"  , z=45., a= 102.91*g/mole, density= 12.40*g/cm3);
//new G4Material("Silver"   , z=47., a= 107.87*g/mole, density= 10.50*g/cm3);
//new G4Material("Tungsten" , z=74., a= 183.85*g/mole, density= 19.35*g/cm3);
//new G4Material("Gold"     , z=79., a= 196.97*g/mole, density= 19.3*g/cm3);
//new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

//
// define a material from elements.   case 1: chemical molecule
//

	G4Material* H2O = new G4Material(
		name="Water",
		density= 1.000*g/cm3,
		ncomponents=2);

	H2O->AddElement(H, natoms=2);
	H2O->AddElement(O, natoms=1);
	// overwrite computed meanExcitationEnergy with ICRU recommended value 
	H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);


	//Natural magnetite:
	G4Material* Magnetite = new G4Material(
		name="Magnetite",
		density = 5.17*g/cm3,
		ncomponents=2);

	Magnetite->AddElement(Fe, natoms=3);
	Magnetite->AddElement(O,  natoms=4);

	//natural magnetite:
	G4Isotope* isoFe54 = new G4Isotope("Fe56", 26, 54, a=53,9396105*g/mole);
	G4Isotope* isoFe56 = new G4Isotope("Fe56", 26, 56, a=55.9349375*g/mole);
	G4Isotope* isoFe57 = new G4Isotope("Fe57", 26, 57, a=56.9353940*g/mole);
	G4Isotope* isoFe58 = new G4Isotope("Fe58", 26, 58, a=57,9332756*g/mole);
	G4Element* NaturalFe = new G4Element("natural Fe", "Fe", ncomponents=4);

	//Natural Fe: 5.845% of 54Fe, 91.754% of 56Fe, 2.119% of 57Fe and 0.282% of 58Fe
	NaturalFe->AddIsotope(isoFe54, abundance=5.845*perCent);
	NaturalFe->AddIsotope(isoFe56, abundance=91.754*perCent);
	NaturalFe->AddIsotope(isoFe57, abundance=2.119*perCent);
	NaturalFe->AddIsotope(isoFe58, abundance=0.282*perCent);

	G4Material* NaturalMagnetite = new G4Material(
		name="NaturalMagnetite",
		density = 5.21*g/cm3, //5.17 * ((57*0.6 + 56*0.4)*3 + 16*4) / (56*3 + 16*4) = 5.21
		ncomponents=2);

	NaturalMagnetite->AddElement(NaturalFe, natoms=3);
	NaturalMagnetite->AddElement(O,            natoms=4);

	//Completely enriched magnetite:
	G4Element* Fe57 = new G4Element("Fe57", "Fe", ncomponents=1);
	Fe57->AddIsotope(isoFe57, abundance=100.*perCent);

	G4Material* Magnetite57 = new G4Material(
		name="Magnetite57",
		density = 5.24*g/cm3, //5.17 * (57*3 + 16*4) / (56*3 + 16*4) = 5.24
		ncomponents=2);

	Magnetite57->AddElement(Fe57, natoms=3);
	Magnetite57->AddElement(O,    natoms=4);

/*
	G4Material* Air = new G4Material(
		"Air",
		density= 1.290*mg/cm3,
		ncomponents=2);

	Air->AddElement(N, fractionmass=79.0*perCent);
	Air->AddElement(O, fractionmass=21.0*perCent);
*/

//
// examples of vacuum
//

G4Material* Galactic =  new G4Material( "Galactic", z=1., a=1.01*g/mole, density=universe_mean_density, kStateGas, 3.e-18*pascal, 2.73*kelvin);

 //default materials of the World
 defaultMaterial = H2O;

//	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*G4double* ExN03DetectorConstruction::ArrangeNanoparticles (const G4int N )
{
	G4int k=10; //sphere divisions
	G4double rtf [3][10000];

}*/

G4VPhysicalVolume* ExN03DetectorConstruction::ConstructCalorimeter()
{
//.G4cout<<"###################### ExN03DetectorConstruction::ConstructCalorimeter... #######################"<<G4endl;


  G4double WorldSize = 20000 * nm;
  G4double CellSize = 15000 * nm ;			
  G4double CoreCellSize = 1000* nm;
  G4double MeanNPSize = 15* nm;
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
 // ComputeCalorParameters();


  G4VisAttributes* wiredVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.1));
  wiredVisAtt->SetForceAuxEdgeVisible(true);
  wiredVisAtt->SetDaughtersInvisible(false);
  wiredVisAtt->SetForceSolid(false);
  wiredVisAtt->SetForceWireframe(true);
  wiredVisAtt->SetLineWidth(1);

  G4VisAttributes* solidVisAtt= new G4VisAttributes(true, G4Colour(0.5,0.5,0.5));
  solidVisAtt->SetForceAuxEdgeVisible(true);
  solidVisAtt->SetDaughtersInvisible(false);
  solidVisAtt->SetForceWireframe(false);
  solidVisAtt->SetForceSolid(true);

  //     
  // World
  //
  //G4GeometryManager::GetInstance()->SetWorldMaximumExtent(CalorSize/2.0);
  //******************

	solidWorld = new G4Box( "World", WorldSize/2.0, WorldSize/2.0, WorldSize/2.0);	

	if (sensitiveDetector) {
		G4cout << "NOTICE: reusing sensitive detector." << G4endl;
	} else {
		sensitiveDetector = new VoxelSD("voxel", WorldSize);
		G4cout << "NOTICE: have created sensitive detector." << G4endl;
	}

	sensitiveDetector -> Activate(false);
	G4cout << "NOTICE: sensitive detector is not activated, use \"/hits/activate voxel\"."<< G4endl;

	logicWorld = new G4LogicalVolume(
		solidWorld,	      //its solid
		GapMaterial,	    //its material (not defaultMaterial)
		"World",          //its name
		0,                //field manager
		sensitiveDetector //sensitive detector
	);

	physiWorld = new G4PVPlacement(
		0,			//no rotation
		G4ThreeVector(),	//at (0,0,0)
		logicWorld,		//its logical volume				 
		"World",		//its name
		 0,			//its mother  volume
		 false,			//no boolean operation
		 0);			//copy number
  
  logicWorld->SetVisAttributes(wiredVisAtt);

	//Limit step for collecting histograms and energy deposit output.
	if (maxStep > 0.0) {
		if (stepLimit) {
			stepLimit->SetMaxAllowedStep(maxStep);
			G4cout << "NOTICE: have reconfigured existing step limits object." << G4endl;
		} else {
			stepLimit = new G4UserLimits(maxStep);
			G4cout << "NOTICE: have created new step limits object." << G4endl; 
		}
		logicWorld->SetUserLimits(stepLimit);
		G4cout << "NOTICE: step limit now is "<<G4BestUnit(maxStep,"Length")<< G4endl;
	}

  //                               
  // Cell
  //
	CellSize= 15000* nm;


  //experimentType == "bulk"; G4cout<< "experiment type  "<< experimentType<<G4endl;

		G4Orb* solidCell = new G4Orb("Cell", d_cell / 2.0);
		G4LogicalVolume* logicCell = new G4LogicalVolume(solidCell,defaultMaterial,"Cell_log");
		G4VPhysicalVolume* physCell = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),  logicCell,"Cell",logicWorld,false,0);   
																										  
// nanoparticles																										    
		
		G4int N=1000;																								 
		G4Orb* solidNP = new G4Orb("NP", d_core / 2.0);
		G4LogicalVolume* logicNP = new G4LogicalVolume(solidNP,AbsorberMaterial,"NP");
		G4VPhysicalVolume* physNP = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),logicNP,"NP_phys",logicCell, false,0); 
																										  


		//G4double orb_steplimit=0.25 * nm;
		
		//just to avoid compillation warnings:
		if (!physNP) {
			G4cerr<<"ERROR while creating object."<<G4endl;
		}

		logicNP->SetVisAttributes(solidVisAtt);
		
	
    
  //PrintCalorParameters();     
  
  //                                        
  // Visualization attributes
  //

  
  //
  //always return the physical World
  //
	is_valid = true;

	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	G4cout<<"Have constructed detector for experiment type "<<experimentType<<G4endl;
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void ExN03DetectorConstruction::PrintCalorParameters()
{
	if (experimentType == "orb") {
		G4cout << "Orb diameter: "<<AbsorberThickness/mm << "mm, material:"<< AbsorberMaterial->GetName()<<G4endl;
	} else if (experimentType == "bulk") {
		//G4cout << "Bulk sample, no secondary target (world only)."<<G4endl;
	} else {
		G4cout << "Experiment type: "<<experimentType<<G4endl;
	}
	G4cout << "World material: " << GapMaterial->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name

	// may use NIST materals, such as G4_WATER:
	G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

	// user-defined materials only:
	//G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);

  if (pttoMaterial) {
		AbsorberMaterial = pttoMaterial;
	} else {
		G4cerr<<"ERROR: Invalid absorber material: "<<materialChoice<<G4endl;
		exit(3);
	}

	is_valid = false; //detector required reconstruction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name

	// may use NIST materals, such as G4_WATER:
	G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

	// user-defined materials only:
	// G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);

  if (pttoMaterial) {
		GapMaterial = pttoMaterial;
	} else {
		G4cerr<<"ERROR: Invalid gap material: "<<materialChoice<<G4endl;
		exit(3);
	}
	is_valid = false; //detector required reconstruction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
	is_valid = false; //detector required reconstruction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::SetCalorSize(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  CalorSize = val;
	is_valid = false; //detector required reconstruction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::SetExperimentType(G4String newtype)
{
	experimentType = newtype;
	is_valid = false; //detector required reconstruction
}

*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN03DetectorConstruction::SetOutParticle(const G4String& partName)
{
	allParticles = false;

	if (partName == "none") {
		outParticle = 0;
	} else if (partName == "all") {
		outParticle = 0;
		allParticles = true;
	} else if (partName == "gamma") {
		outParticle = G4Gamma::Gamma();
	} else if (partName == "e-") {
		outParticle = G4Electron::Electron();
	} else if (partName == "proton") {
		outParticle = G4Proton::Proton();
	} else {
		G4cerr<<"ERROR: invalid particle name "<<partName<<G4endl;
		exit(2);
	}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void ExN03DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
