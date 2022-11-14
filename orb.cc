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
// $Id: exampleN03.cc,v 1.19 2003/09/15 15:38:08 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"


#ifdef G4UI_USE_XM
# include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
# include "G4VisExecutive.hh"
#endif


#include "Randomize.hh"

#include "ExN03DetectorConstruction.hh"
#include "ExN03PhysicsList.hh"
#include "ExN03PrimaryGeneratorAction.hh"
#include "ExN03RunAction.hh"
#include "ExN03EventAction.hh"
#include "ExN03SteppingAction.hh"
#include "ExN03SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
	// random engine
	CLHEP::Ranlux64Engine randomEngine;
  CLHEP::HepRandom::setTheEngine(&randomEngine);

  G4long seed = time(NULL); //G4long time_t
	if (seed == G4long(((time_t) -1))) {
		G4cerr << "WARNING: Failed to get time for random seed." << G4endl;
	}

  CLHEP::HepRandom::setTheSeed(seed);
	G4cout<<"Set seed " <<seed<<G4endl;

  // Verbose output class
  G4VSteppingVerbose::SetInstance(new ExN03SteppingVerbose);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  ExN03DetectorConstruction* detector = new ExN03DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new ExN03PhysicsList);
  
	G4UIsession* session = 0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else           
# ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
# else
      session = new G4UIterminal();
# endif
#endif
    }
  
#ifdef G4VIS_USE
  // visualization manager
 // G4VisManager* visManager = new ExN03VisManager;
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new ExN03PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new ExN03RunAction);
  ExN03EventAction* eventaction = new ExN03EventAction;
  runManager->SetUserAction(eventaction);
  runManager->SetUserAction(new ExN03SteppingAction(detector, eventaction));
  
  //Initialize G4 kernel
	//?runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
//      UI->ApplyCommand("/control/execute vis.mac");    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute gui.mac");
#endif
#ifdef G4VIS_USE
      UI->ApplyCommand("/control/execute vis.mac");     
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
