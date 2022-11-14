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
// $Id: ExN03PrimaryGeneratorMessenger.cc,v 1.8 2002/12/16 16:37:27 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03PrimaryGeneratorMessenger.hh"

#include "ExN03PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PrimaryGeneratorMessenger::ExN03PrimaryGeneratorMessenger(
                                          ExN03PrimaryGeneratorAction* ExN03Gun)
:ExN03Action(ExN03Gun)
{
  gunDir = new G4UIdirectory("/N03/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  RndmCmd = new G4UIcmdWithAString("/N03/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot particles in random directions.");
  RndmCmd->SetGuidance("  Choice : on, off (default)");
  RndmCmd->SetParameterName("gunRndm",false);
	//  RndmCmd->SetDefaultValue("off");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BeamDiamCmd = new G4UIcmdWithADoubleAndUnit("/N03/gun/beamDiam",this);
  BeamDiamCmd->SetGuidance("Set beam diameter (default: zero).");
  BeamDiamCmd->SetParameterName("beamDiam",true);
	BeamDiamCmd->SetDefaultValue(0.0);
  BeamDiamCmd->SetRange("beamDiam>=0.");
  BeamDiamCmd->SetUnitCategory("Length");
  BeamDiamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PrimaryGeneratorMessenger::~ExN03PrimaryGeneratorMessenger()
{
  delete gunDir;
  delete RndmCmd;
  delete BeamDiamCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
	if (command == RndmCmd){
  	ExN03Action->SetRndmFlag(newValue);
  	if (newValue=="on") {
			G4cout << "Will shoot particles in random directions."<<G4endl;
		} else {
			G4cout << "Will shoot particles as parallel beam."<<G4endl;
		}

	} else if (command == BeamDiamCmd) {
		ExN03Action->SetBeamDiameter(BeamDiamCmd->GetNewDoubleValue(newValue));
		G4cout << "Set beam diameter to "<<newValue<<G4endl;

	} else {
		G4cerr << "ERROR: Invalid command !" << G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

