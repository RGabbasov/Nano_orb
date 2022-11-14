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
// $Id: ExN03PhysicsList.cc,v 1.17 2003/10/24 12:34:15 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN03PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4ProcessManager.hh"

#include "G4ios.hh"
#include "G4StepLimiter.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmProcessOptions.hh"



//#include "G4UserSpecialCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PhysicsList::ExN03PhysicsList():
	G4VUserPhysicsList(),
	fEmPhysicsList(NULL),
	fMessenger(NULL),
	activateMossbauer(true)
{
	defaultCutValue = 0.1*mm;
	//	SetVerboseLevel(1);

	fEmPhysicsList = new G4EmLivermorePhysics();
	fMessenger = new PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PhysicsList::~ExN03PhysicsList()
{
	if (fMessenger)
		delete fMessenger;

	if (fEmPhysicsList)
		delete fEmPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PhysicsList::ConstructParticle()
{
//.G4cout<<"###################### ExN03PhysicsList::ConstructParticle... #######################"<<G4endl;

  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

	if (fEmPhysicsList) {
		fEmPhysicsList->ConstructParticle();
	}
//.G4cout<<"###################### ...ExN03PhysicsList::ConstructParticle #######################"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PhysicsList::ConstructProcess()
{
//.G4cout<<"###################### ExN03PhysicsList::ConstructProcess... #######################"<<G4endl;

  AddTransportation();

	if (fEmPhysicsList) {
		fEmPhysicsList->ConstructProcess();
	}

	G4EmProcessOptions emOptions;
	emOptions.SetFluo(true); // To activate deexcitation processes and fluorescence
	emOptions.SetAuger(true); // To activate Auger effect if deexcitation is activated
	emOptions.SetPIXE(true); // To activate Particle Induced X-Ray Emission (PIXE)
	emOptions.SetDeexcitationActiveRegion("World", true, true, true); // mandatory to activate Auger emission

	if (activateMossbauer)
		AddMossbauer();

	AddStepMax();

//.G4cout<<"###################### ...ExN03PhysicsList::ConstructProcess #######################"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PhysicsList::AddStepMax()
{
  theParticleIterator->reset();

  while ((*theParticleIterator)()) {
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		pmanager->AddDiscreteProcess(new G4StepLimiter());
  }
}

//*******************************************************
#include "MossbauerScattering.hh"

void ExN03PhysicsList::AddMossbauer()
{ 
 MossbauerScattering* moss= new MossbauerScattering("MossbauerScattering");
 moss->ConstructSpectrum();
 moss->SetFe57abundance(1.00); 
 theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
   
    if (moss->IsApplicable(*particle)) {
	 pmanager->AddDiscreteProcess (moss);
     pmanager->SetProcessOrderingToLast(moss,idxPostStep);
	G4cout <<" Mossbauer Scattering proceess registered"<< G4endl;
  }
 }
}
//**************************************************************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PhysicsList::SetCuts()
{
	
   if (verboseLevel >0) {
    G4cout << "ExN03PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
	

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //

//    SetCutValue(defaultCutValue, "gamma");
//    SetCutValue(defaultCutValue, "e-");
//    SetCutValue(defaultCutValue, "e+");
    
  SetCutsWithDefault();

	//recomended 250*eV
//	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(20*eV, 1*GeV);

    //if (verboseLevel>0) {
     //   DumpCutValuesTable();}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
