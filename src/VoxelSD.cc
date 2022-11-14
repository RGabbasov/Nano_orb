//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: VoxelSD.cc,v 1.1 2007/11/16 14:29:33 kmura Exp $
// $Name: geant4-09-02 $
//
// ====================================================================
//   VoxelSD.cc
//
//                                         2007 Q
// ====================================================================

#include "VoxelSD.hh"
#include "Analysis.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


const double pi  = 3.14159265358979323846;

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////
VoxelSD::VoxelSD(const G4String& name, G4double sz)
  : G4VSensitiveDetector(name), calorSize(sz)
//////////////////////////////////////
{
  G4SDManager::GetSDMpointer()-> AddNewDetector(this);
}


///////////////////
VoxelSD::~VoxelSD()
///////////////////
{
}

///////////////////////////////////////////////////////////////
G4bool VoxelSD::ProcessHits(G4Step* astep, G4TouchableHistory*)
///////////////////////////////////////////////////////////////
{
	G4double edep = astep -> GetTotalEnergyDeposit() / eV;

	if (edep == 0.0)
		return false;

	const G4Track* track = astep->GetTrack();

	//Select histogram.
	G4String name = track->GetDefinition()->GetParticleName();
	G4int h1id = 0;
	G4int h2id;

	if (name == "e-") {
		h1id = 5;
		h2id = 1;
	} else if (name == "gamma") {
		if (! track->GetDynamicParticle()->GetPrimaryParticle()) {
			h1id = 6;
			h2id = 2; //secondary
		} else {
			h2id = 3; //primary
		}
	} else {
		G4cerr << "WARNING: skipping energy deposit by unknown particle "<<name<<G4endl;
		return true;
	}


	G4ThreeVector x0 = astep->GetPreStepPoint()->GetPosition();
	G4ThreeVector x1 = astep->GetPostStepPoint()->GetPosition();
	G4ThreeVector p = G4UniformRand()*(x1-x0) + x0; // position sampling

	//Make normalized coordinate -1..+1
	G4double x = p.x() / (0.5 * calorSize);

	//Normalized distance from beam axis, useful range 0..1
	G4double r = sqrt(p.y()*p.y() + p.z()*p.z()) / (0.5 * calorSize);

	//Normalized distance from particle center, useful range 0..1
	G4double R = p.mag() / (0.5 * calorSize);

	if (x < -1.000001 || x > 1.000001 || r > sqrt(2.0) + 1e-6 || R > sqrt(3.0) + 1e-6)
		G4cerr << "WARNING: energy deposit beyond detector."<<G4endl;

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	//Scale to eV per box volume
	//ring numbers: 1..HALF_BOX_BINS
	//ring width = ring height = 1/HALF_BOX_BINS
	//world volume = 8 (!)

	const G4double world_volume = 8;
	const G4double d = 1.0 / HALF_BOX_BINS; //size of histogram bin

	//Sperical average.
	if (h1id) {
		const G4int N = 1 + R*HALF_BOX_BINS; //number of spherical shell
		const G4double shellvolume = (4.0/3.0)*pi*(N*N*N-(N-1)*(N-1)*(N-1))*d*d*d;
		analysisManager->FillH1(h1id, R, edep/shellvolume);
	}

	//Cylindric average.
	const G4int nring = 1 + r*HALF_BOX_BINS;
	//N^2 - (N-1)^2 = 2N - 1
	const G4double ringarea = (2*nring - 1) * pi * d * d;
	const G4double ringvolume = ringarea * d;
	analysisManager->FillH2(h2id, x, r, edep/ringvolume);

	/*
	G4cout << "energy deposit at " << p.x()/nm << '\t' << p.y()/nm << '\t' << p.z()/nm << " nm\t"
				 << edep << " eV\tby " << name << G4endl;
	*/

  return true;
}

////////////////////////////////////////////
void VoxelSD::Initialize(G4HCofThisEvent*)
////////////////////////////////////////////
{
}


////////////////////////////////////////////
void VoxelSD::EndOfEvent(G4HCofThisEvent*)
////////////////////////////////////////////
{
}

