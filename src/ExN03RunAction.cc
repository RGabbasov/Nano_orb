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
// $Id: ExN03RunAction.cc,v 1.15 2003/11/25 16:50:13 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03RunAction.hh"
#include "Analysis.hh"

//#include "ExN03PrimaryGeneratorAction.hh"
#include "ExN03DetectorConstruction.hh"

#include "G4Run.hh"
//#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"

#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03RunAction::ExN03RunAction()
{
 G4cout<<"###################### ExN03RunAction #######################"<<G4endl;

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
	analysisManager->SetFirstHistoId(1);

	//1D histograms.
  analysisManager->CreateH1("h1","Secondary electron energy, eV (range 0-15 keV)", ENERGY_BINS, 0.0, 15e+3);
  analysisManager->CreateH1("h2","Secondary electron energy, eV (range 0-1 keV)", ENERGY_BINS, 0.0, 1e+3);
  analysisManager->CreateH1("h3","Secondary gamma energy, eV (range 0-15 keV)", ENERGY_BINS, 0.0, 15e+3);
  analysisManager->CreateH1("h4","Secondary gamma energy, eV (range 0-1 keV)", ENERGY_BINS, 0.0, 1e+3);

  analysisManager->CreateH1("h5","Radial dose by electrons, eV", HALF_BOX_BINS, 0.0, 1.0);
  analysisManager->CreateH1("h6","Radial dose by secondary gammas, eV", HALF_BOX_BINS, 0.0, 1.0);
  //**************
   analysisManager->CreateH1("h7","Secondary electron energy deposit, eV (range 0-15 keV)", ENERGY_BINS, 0.0, 15e+3);
   analysisManager->CreateH1("h8","Secondary electron energy deposit, eV (range 0-15 keV)", ENERGY_BINS, 0.0, 15e+3);
   analysisManager->CreateH1("h9","Primary gamma energy deposit, eV (range 0-15 keV)", ENERGY_BINS, 0.0, 15e+3);
   //*****************

	for (G4int i=1; i<=6; ++i) {
		analysisManager->SetH1Ascii(i, true); //write this histogram to ascii file
		analysisManager->SetH1XAxisTitle(i, "Energy, eV");
		analysisManager->SetH1YAxisTitle(i, "Count");
	}

	//2D histograms.
  analysisManager->CreateH2("hh1","Dose by electrons, eV",
		HALF_BOX_BINS*2, -1.0, 1.0,  //normalized coordinate along beam
		HALF_BOX_BINS,    0.0, 1.0); //along radius

  analysisManager->CreateH2("hh2","Dose by secondary gammas, eV",
		HALF_BOX_BINS*2, -1.0, 1.0,  //normalized coordinate along beam
		HALF_BOX_BINS,    0.0, 1.0); //along radius

  analysisManager->CreateH2("hh3","Dose by primary gammas, eV",
		HALF_BOX_BINS*2, -1.0, 1.0,  //normalized coordinate along beam
		HALF_BOX_BINS,    0.0, 1.0); //along radius

	for (G4int i=1; i<=3; ++i) {
		analysisManager->SetH2Ascii(i, false); //do not write this histogram to ascii file
		analysisManager->SetH2XAxisTitle(i, "X");
		analysisManager->SetH2YAxisTitle(i, "R");
		analysisManager->SetH2ZAxisTitle(i, "Dose, eV");
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03RunAction::~ExN03RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
//.G4cout<<"###################### ExN03RunAction::BeginOfRunAction... #######################"<<G4endl;

	G4cout<<"Run "<<aRun->GetRunID()<<" started."<<G4endl;

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "hist";
  analysisManager->OpenFile(fileName);

#ifdef ACTIVITY_LED
		system("xset led 3");
#endif

  //inform the runManager to save random number seed
  //!G4RunManager::GetRunManager()->SetRandomNumberStore(true);
//.G4cout<<"###################### ...ExN03RunAction::BeginOfRunAction #######################"<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03RunAction::EndOfRunAction(const G4Run* aRun)
{
	G4cout<<"Run "<<aRun->GetRunID()<<" completed."<<G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

#ifdef ACTIVITY_LED
		system("xset -led 3");
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
