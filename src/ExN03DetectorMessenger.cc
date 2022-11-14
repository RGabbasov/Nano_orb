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
// $Id: ExN03DetectorMessenger.cc,v 1.9 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03DetectorMessenger.hh"

#include "ExN03DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorMessenger::ExN03DetectorMessenger(
                                           ExN03DetectorConstruction* ExN03Det)
:ExN03Detector(ExN03Det)
{ 
  N03Dir = new G4UIdirectory("/N03/");
  N03Dir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/N03/det/");
  detDir->SetGuidance("detector control");
       
  AbsMaterCmd = new G4UIcmdWithAString("/N03/det/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("absMat",false);
  AbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set thickness of the absorber");
  AbsThickCmd->SetParameterName("absThick",false);
  AbsThickCmd->SetRange("absThick>0.");
  AbsThickCmd->SetUnitCategory("Length");
  AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  GapMaterCmd = new G4UIcmdWithAString("/N03/det/setGapMat",this);
  GapMaterCmd->SetGuidance("Select Material of the Gap.");
  GapMaterCmd->SetParameterName("gapMat",false);
  GapMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  OutParticle = new G4UIcmdWithAString("/N03/det/setOutParticle",this);
  OutParticle->SetGuidance("Select output particles to count (none to disable).");
  OutParticle->SetGuidance("  Choice : none (default), gamma, e-, all");
  OutParticle->SetParameterName("outParticle",true);
  OutParticle->SetDefaultValue("none");
  OutParticle->SetCandidates("none gamma e- all");
  OutParticle->AvailableForStates(G4State_PreInit,G4State_Idle);

  ExperimentType = new G4UIcmdWithAString("/N03/det/setExperimentType",this);
  ExperimentType->SetGuidance("Select experiment type.");
  ExperimentType->SetGuidance("  Choice : orb (default), bulk, test");
  ExperimentType->SetParameterName("experimentType",true);
  ExperimentType->SetDefaultValue("foil");
  ExperimentType->SetCandidates("orb bulk test");
  ExperimentType->AvailableForStates(G4State_PreInit,G4State_Idle);
    
/*  GapThickCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setGapThick",this);
  GapThickCmd->SetGuidance("Set Thickness of the Gap");
  GapThickCmd->SetParameterName("Size",false);
  GapThickCmd->SetRange("Size>=0.");
  GapThickCmd->SetUnitCategory("Length");  
  GapThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
  SizeCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setSize",this);
  SizeCmd->SetGuidance("Set tranverse size of the calorimeter");
  SizeCmd->SetParameterName("calorSize",false);
  SizeCmd->SetRange("calorSize>0.");
  SizeCmd->SetUnitCategory("Length");    
  SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
//.  NbLayersCmd = new G4UIcmdWithAnInteger("/N03/det/setNbOfLayers",this);
//.  NbLayersCmd->SetGuidance("Set number of layers.");
//.  NbLayersCmd->SetParameterName("NbLayers",false);
//.  NbLayersCmd->SetRange("NbLayers>0 && NbLayers<500");
//.  NbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/N03/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  StepLimitCmd = new G4UIcmdWithADoubleAndUnit("/N03/det/setStepLimit",this);
  StepLimitCmd->SetGuidance("Set step limit, zero to disable step limiting (the default value).");
  StepLimitCmd->SetParameterName("stepLimit",false);
	//.	StepLimitCmd->SetDefaultValue(0.0*mm);
  StepLimitCmd->SetRange("stepLimit>=0.");
  StepLimitCmd->SetUnitCategory("Length");    
  StepLimitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorMessenger::~ExN03DetectorMessenger()
{
//.  delete NbLayersCmd;
  delete AbsMaterCmd; 
  delete GapMaterCmd; 
  delete OutParticle;
	delete ExperimentType;
  delete AbsThickCmd; 
//.  delete GapThickCmd;
  delete SizeCmd;   
  delete UpdateCmd;
  delete StepLimitCmd;
  delete detDir;
  delete N03Dir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == AbsMaterCmd ) { 
		ExN03Detector->SetAbsorberMaterial(newValue);
		//?ExN03Detector->UpdateGeometry();
		
	} else if( command == AbsThickCmd ) { 
		ExN03Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));
		//?ExN03Detector->UpdateGeometry();

	} else if( command == GapMaterCmd ) { 
		ExN03Detector->SetGapMaterial(newValue);
//?		ExN03Detector->UpdateGeometry(); 

	} else if( command == OutParticle ) {
		ExN03Detector->SetOutParticle(newValue);

	} else if( command == ExperimentType ) { 
		ExN03Detector->SetExperimentType(newValue);
		//		ExN03Detector->UpdateGeometry();
	
	} else if( command == SizeCmd ) { 
		ExN03Detector->SetCalorSize(SizeCmd->GetNewDoubleValue(newValue));
//?		ExN03Detector->UpdateGeometry();

	} else if(command == StepLimitCmd) { 
		ExN03Detector->SetMaxStep(StepLimitCmd->GetNewDoubleValue(newValue));

	} else if( command == UpdateCmd) { 
		ExN03Detector->UpdateGeometry();

	} else {
		G4cerr << "ERROR: Invalid command !" << G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
