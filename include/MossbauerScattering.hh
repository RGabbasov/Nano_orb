// ************************************************************************
//
// ZepIII-V03 (ZepIIIScintillation)
//
// Henrique Araujo, Imperial College London
// 
// ************************************************************************

#ifndef MossbauerScattering_h
#define MossbauerScattering_h 1

//#include "G4VRestDiscreteProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "globals.hh"

class MossbauerScattering : public G4VDiscreteProcess {

  public:

    MossbauerScattering(const G4String& processName = "MossbauerScattering");

    ~MossbauerScattering();	

public:

//	G4VParticleChange* AtRestDoIt ( G4Track& aTrack,  G4Step& aStep);

	G4double GetMeanFreePath(const G4Track& aTrack,G4double,G4ForceCondition*);

	void ConstructSpectrum ();
	void SetFe57abundance(G4double abund);

	//G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*);

	//G4double GetCrossSection(const G4DynamicParticle& aPart );

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

    //?inline
		G4bool IsApplicable(const G4ParticleDefinition& partDef);


  // inlines
 // public : 

  private:
  G4int CalculateChannel (G4double energy);
  G4double GetAmlplitude (G4int n);
  G4String ResonancePhysVolume;
  G4double ResonanceEnergy, EnergyResonanceRange, SpectrumVelocityRange;
  G4double ElectronConversionProb;
  G4double KShellProb, LShellProb, MShellProb;
  G4double clight,pi;
  G4double exit_t;
  G4double v_st,n0;
  G4double abund57;
  G4int N;
  G4double spectr[256];
};


#endif 
