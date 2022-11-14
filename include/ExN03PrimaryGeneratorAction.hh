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
// $Id: ExN03PrimaryGeneratorAction.hh,v 1.6 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN03PrimaryGeneratorAction_h
#define ExN03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class ExN03DetectorConstruction;
class ExN03PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    ExN03PrimaryGeneratorAction(ExN03DetectorConstruction*);    
   ~ExN03PrimaryGeneratorAction();

public:
	void GeneratePrimaries(G4Event*);
	void SetRndmFlag(G4String val) {rndmFlag = (val == "on");}
	G4double GetParticleEnergy() const;
	G4String GetParticleName() const;

	G4double GetBeamDiameter()  {return beamDiameter;}
	void SetBeamDiameter(G4double val) {beamDiameter = val;}

private:
	G4ParticleGun *particleGun;         //pointer a to G4  class
	
	//G4GeneralParticleSource* ParticleGun;
	ExN03DetectorConstruction*    ExN03Detector;  //pointer to the geometry
	
	ExN03PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
	G4bool                          rndmFlag;	  //flag for a rndm impact point
	G4double                        beamDiameter;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


