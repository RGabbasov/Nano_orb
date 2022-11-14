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
// $Id: ExN03DetectorConstruction.hh,v 1.6 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN03DetectorConstruction_h
#define ExN03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UserLimits.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4Orb.hh"

#include "VoxelSD.hh"

class G4Box;
class G4Orb;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ExN03DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
	ExN03DetectorConstruction();
	~ExN03DetectorConstruction();

public:
     
   //  void SetAbsorberMaterial (G4String);     
   //  void SetAbsorberThickness(G4double);     

  //   void SetGapMaterial (G4String);     
  //void SetGapThickness(G4double);
     
     void SetWorldSize(G4double);          
//.     void SetNbOfLayers (G4int);   
      
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  //   void SetExperimentType(G4String newtype);
     
  public:
  
   //  void PrintCalorParameters(); 
                    
	//.     G4double GetWorldSizeX()           {return WorldSizeX;}; 
	//.     G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
//.     G4double GetCalorThickness()       {return CalorThickness;}; 
    // G4double GetCalorSize()          {return WorldSize;};
     
     //G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     //G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     
     //G4Material* GetGapMaterial()       {return GapMaterial;};
//.     G4double    GetGapThickness()      {return GapThickness;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
//.     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
//.     const G4VPhysicalVolume* GetGap()        {return physiGap;};
                 
  private:
     
	G4Material*        AbsorberMaterial;
	G4double           AbsorberThickness;
     
	//G4Material*        GapMaterial;
     
	G4double           World_Size;
	//new possible parameters
	G4double		CellSize;			//diameters the of cell, cell core and nanoparticle
	G4double		CoreCellSize;
	G4double		MeanNPSize;
	G4int			Nnp;							//number of nanoparticles in cell

	//G4double d_cell,d_core,d_NP; 
	
	G4Material*        defaultMaterial;
  
	G4Box*             solidWorld;    //pointer to the solid World 
	G4LogicalVolume*   logicWorld;    //pointer to the logical World
	G4VPhysicalVolume* physiWorld;    //pointer to the physical World

	G4Orb*             solidCell; //pointer to the solid Absorber Cell Volume
	G4LogicalVolume*   logicCell; //pointer to the logical Absorber  Cell Volume
	G4VPhysicalVolume* physiCell; //pointer to the physical Absorber  Cell Volume
     
	ExN03DetectorMessenger* detectorMessenger;  //pointer to the Messenger
	G4UserLimits* stepLimit;                   // pointer to user step limits
	G4double      maxStep;

	VoxelSD* sensitiveDetector;
     
	//G4String experimentType;
	G4bool is_valid;
	G4bool allParticles;
	G4ParticleDefinition* outParticle;

	

public:
	G4bool IsValid() const {return is_valid;}
	void SetOutParticle(const G4String& partName);
	const G4ParticleDefinition* GetOutParticleDefinition() {return outParticle;}
	G4bool AllOutParticles() const {return allParticles;}

	G4double GetMaxStep() const {return maxStep;}
	void SetMaxStep(G4double val) {maxStep = val; is_valid = false;}


private:
    
	void DefineMaterials();
	//void ComputeCalorParameters();
	//G4VPhysicalVolume* ConstructCalorimeter();
	//G4double* ArrangeNanoparticles (const G4int );
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*inline void ExN03DetectorConstruction::ComputeCalorParameters()
{
}*/

//-----------------------------------------------------------------------------



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

