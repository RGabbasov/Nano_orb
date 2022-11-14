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
// $Id: ExN03SteppingAction.cc,v 1.8 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03SteppingAction.hh"

#include "ExN03DetectorConstruction.hh"
#include "ExN03EventAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "ExN03PrimaryGeneratorAction.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

////#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03SteppingAction::ExN03SteppingAction(ExN03DetectorConstruction* det,
                                         ExN03EventAction* evt)
:detector(det), eventaction(evt), verboseLevel(0)
{
	eventaction=evt;
//	verboseLevel = 3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03SteppingAction::~ExN03SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4Track* track = aStep->GetTrack();
	G4String name = track->GetParticleDefinition()->GetParticleName();
	//Record exiting particles.
	if (track->GetNextVolume() == 0) {
		if (name == "e-") {
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
			const G4double en = track->GetKineticEnergy() / eV;
			analysisManager->FillH1(1, en);
			analysisManager->FillH1(2, en);

			if (verboseLevel >= 2)
				G4cout << "added e- to histogram, energy " << en << " eV." << G4endl;

		} else if (name == "gamma") {
			if (! track->GetDynamicParticle()->GetPrimaryParticle()) {
				G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
				const G4double en = track->GetKineticEnergy() / eV;
				analysisManager->FillH1(3, en);
				analysisManager->FillH1(4, en);

				if (verboseLevel >= 2)
					G4cout << "adding secondary gamma to histogram, energy " << en << " eV." << G4endl;
			} else {
				if (verboseLevel >= 3) {
					const G4double en = track->GetKineticEnergy() / eV;
					G4cout << "primary gamma exit, energy " << en << " eV." << G4endl;
				}
			}
		} else {
			G4cerr << "WARNING: skipping unknown particle: "<<name<<G4endl;
			return;
		}

		if (detector->AllOutParticles() || 
				(detector->GetOutParticleDefinition() && 
				 detector->GetOutParticleDefinition() == track->GetDefinition())) {
			
			const G4ThreeVector& pos = track->GetPosition();
			
			//project point to orb boundary (orb center is at 0,0,0)
			//todo: fixme: project along momentum direction
			//AbsorberThickness is diameter
			const double radius = detector->GetAbsorberThickness() / 2.0;
			
			G4ThreeVector p = pos * (radius / (pos.r() * nm)); // pos will be in nanometers
			const G4ThreeVector& md = track->GetMomentumDirection(); //unit vector
			const double en = track->GetKineticEnergy();
			//G4cout <<"energy "<<en<<G4endl;
			//Will exclude the primary beam.
			//get primary beam energy:
			double pen;
			const ExN03PrimaryGeneratorAction* pga = dynamic_cast<const ExN03PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

			if (pga) {
				pen = pga->GetParticleEnergy();
			} else {
				G4cerr << "ERROR: Failed to get primary particle energy." << G4endl;
				exit(3);
			}

			if (md.x() != 1.0 || md.y() != 0.0 || md.z() != 0.0 || en != pen) {
				G4cout<<"exit at " <<p.x()<<'\t'<<p.y()<<'\t'<<p.z()<<" nm, "<<md.x()<<'\t'<<md.y()<<'\t'<<md.z()<<" direction\t" << G4BestUnit(en, "Energy");
				if (detector->AllOutParticles()) {
					G4cout<<'\t'<<(track->GetDefinition()->GetParticleName());
				}
				G4cout<<G4endl;
			}
		}
	} 
	
	//accumulating energy deposits inside world
	//************************************
	G4String PhVolume = track->GetVolume()->GetName(); 
	if (PhVolume != "Orb")
		{
		G4double dep=aStep->GetTotalEnergyDeposit();
		if (name == "e-")
			eventaction->AddElectronEnergyDeposit(dep);
		else if (name == "gamma")	//secondary gamma
			if (! track->GetDynamicParticle()->GetPrimaryParticle()) 
				eventaction->AddGammaSecondaryEnergyDeposit(dep);
			else					//primary  gamma
				eventaction->AddGammaPrimaryEnergyDeposit(dep);
		}
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
