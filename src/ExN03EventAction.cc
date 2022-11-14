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
// $Id: ExN03EventAction.cc,v 1.23 2003/11/12 16:15:48 johna Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03EventAction.hh"
#include "ExN03EventActionMessenger.hh"
#include "ExN03DetectorConstruction.hh"


//.#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03EventAction::ExN03EventAction():
	drawFlag("none"),
	printModulo(1),
	eventMessenger(0)
{
  eventMessenger = new ExN03EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03EventAction::~ExN03EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03EventAction::BeginOfEventAction(const G4Event* evt)
{  
	if (evt == NULL) {
		G4cerr << "NULL event" << G4endl;
	}


	G4int evtNb = evt->GetEventID();

	if (evtNb % printModulo == 0) { 
		fprintf(stderr, ".");
		// G4cout << G4endl << "---> Begin of event: " << evtNb << G4endl;
		//		HepRandom::showEngineStatus();
 	}
  dE_el=0;
  dE_gamma=0;
  dE_prim_gamma=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03EventAction::EndOfEventAction(const G4Event* evt)
{
  // extract the trajectories and draw them

  // You can get a default drawing without this code by using, e.g.,
  // /vis/scene/add/trajectories 1000
  // The code here adds sophistication under control of drawFlag.
	//drawFlag="all";
	if (drawFlag != "none") {
		G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();

		if (pVisManager) {
			G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
			G4int n_trajectories = 0;

			if (trajectoryContainer) {
				n_trajectories = trajectoryContainer->entries();
			}

			for (G4int i=0; i<n_trajectories; i++) { 
				G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[i]);
				if (drawFlag == "all")
					trj->DrawTrajectory();
				else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.0))
					trj->DrawTrajectory();
				else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.0))
					trj->DrawTrajectory();
			}
		}
	}

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillH1(7, dE_el/eV);
	analysisManager->FillH1(8, dE_gamma/eV);
	analysisManager->FillH1(9, dE_prim_gamma/eV);
	
	/*if (dE_el+dE_gamma+dE_prim_gamma)
	{
	G4cout <<"elec_dep "<<  dE_el/eV<<G4endl;
	G4cout <<"gamma_dep "<<  dE_gamma/keV<<G4endl;
	G4cout <<"primary_gamma_dep "<<  dE_prim_gamma/keV<<G4endl;
	}*/
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
