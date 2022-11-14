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
// $Id: ExN03PrimaryGeneratorAction.cc,v 1.7 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03PrimaryGeneratorAction.hh"

#include "ExN03DetectorConstruction.hh"
#include "ExN03PrimaryGeneratorMessenger.hh"

//#include "G4GeneralParticleSource.hh"
//#include "G4SPSAngDistribution.hh"
//#include "G4SPSEneDistribution.hh"
//#include "G4SPSPosDistribution.hh"



#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PrimaryGeneratorAction::ExN03PrimaryGeneratorAction(ExN03DetectorConstruction* ExN03DC):
	ExN03Detector(ExN03DC),
	rndmFlag(false),
	beamDiameter(0.0)
{
//.G4cout<<"############ ExN03PrimaryGeneratorAction... #################"<<G4endl;
//?  const G4double pi = 3.14159265358979323846;
//?  G4int n_particle = 1;

  particleGun =  new G4ParticleGun();
  //ParticleGun= new G4GeneralParticleSource();
  //create a messenger for this class
  gunMessenger = new ExN03PrimaryGeneratorMessenger(this);

  // default particle kinematic

  //G4SPSEneDistribution* EnergyDist = ParticleGun->GetCurrentSource()->GetEneDist();
  //G4SPSAngDistribution* AngularDist = ParticleGun->GetCurrentSource()->GetAngDist();
  //G4SPSPosDistribution* PositionDist = ParticleGun->GetCurrentSource()->GetPosDist();

  //particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma"); //see also below
  particleGun->SetParticleDefinition(particle);

  //position 
  //PositionDist->SetPosDisType("Point");
  G4double position = -0.4999999999 * ExN03Detector->GetCalorSize(); //0.45?
  particleGun->SetParticlePosition(G4ThreeVector(position, 0.0*cm, 0.0*cm));

 //angular dist 
 /* AngularDist->SetAngDistType("Beam");
  AngularDist->SetMaxTheta(pi/2);
  AngularDist->SetMinTheta(pi/2);
  AngularDist->SetMaxPhi(pi/2);
  AngularDist->SetMinPhi(pi/2);*/
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));

  //energy dist
  //EnergyDist->SetEnergyDisType("Mono");
  //EnergyDist->SetMonoEnergy((14.414+0.00000000002)*keV);
  


//.G4cout<<"############ ...ExN03PrimaryGeneratorAction #################"<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PrimaryGeneratorAction::~ExN03PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//-----------------------------------------------------------------------------
//Returns randomly oriented unit vector.

void RandDirection(double& x, double& y, double& z) {
	double len;

	do {
		x = 2 * (G4UniformRand() - 0.5); // -1 .. 1
		y = 2 * (G4UniformRand() - 0.5);
		z = 2 * (G4UniformRand() - 0.5);
		len = sqrt(x*x + y*y + z*z);
	} while (len <= 0.1 || len >= 1.0);

	x /= len;
	y /= len;
	z /= len;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//.G4cout<<"############ ExN03PrimaryGeneratorAction::GeneratePrimaries... #################"<<G4endl;
  //this function is called at the begining of event
  // 

	if (! ExN03Detector->IsValid()) {
		G4cerr << "ERROR: detector is not valid." << G4endl;
		exit(3);
	}

	G4double x(0.0), y(0.0), z(0.0); //particle position
	G4double vx(1.0), vy(0.0), vz(0.0); //particle direction

	if (rndmFlag) {
		RandDirection(vx, vy, vz); //random direction, otherwise 1,0,0
	} else {
		x = -0.4999999999 * ExN03Detector->GetCalorSize();
	}

	if (beamDiameter > 0.0) {
		do {
			y = beamDiameter * (G4UniformRand() - 0.5);
			z = beamDiameter * (G4UniformRand() - 0.5);
		} while (y*y + z*z >= (beamDiameter/2.0) * (beamDiameter/2.0));
	}

  particleGun->SetParticlePosition(G4ThreeVector(x, y,z));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(vx, vy, vz));

  G4double ResEnergy=14.414*keV;
  G4double width = 1*0.000001*eV;  //10^-10 keV
  G4double energy= G4RandGauss::shoot(ResEnergy, width);
  particleGun->SetParticleEnergy(energy);

  particleGun->GeneratePrimaryVertex(anEvent);
//.G4cout<<"############ ...ExN03PrimaryGeneratorAction::GeneratePrimaries #################"<<G4endl;
}

//-----------------------------------------------------------------------------

G4double ExN03PrimaryGeneratorAction::GetParticleEnergy() const
{
	return particleGun->GetParticleEnergy();
}

//-----------------------------------------------------------------------------

G4String ExN03PrimaryGeneratorAction::GetParticleName() const
{
	return particleGun->GetParticleDefinition()->GetParticleName();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

