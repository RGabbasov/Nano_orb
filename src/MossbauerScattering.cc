

#include "MossbauerScattering.hh"

#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh" 
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include <cmath>
#include "G4Gamma.hh"
#include "G4Electron.hh"


MossbauerScattering::MossbauerScattering
   (const G4String& processName) : G4VDiscreteProcess(processName) {

  G4cout << GetProcessName() << " is created " << G4endl;
 pi=3.141592;
 ResonanceEnergy         = 14.414*keV;  // Mossbauer level for Fe57
 ResonancePhysVolume	 = "NP_phys";		// Physical volume for nanoparticle;
 //SpectrumVelocityRange   =10;			// max distance from resonance level in velocity, mm/s 
 			
// EnergyResonanceRange	 =ResonanceEnergy*(SpectrumVelocityRange/clight);  //max distance from resonance level in energy;
 ElectronConversionProb  =0.896;        //electron conversion probality integral and separetely to K,L,M shells
 KShellProb				 =0.802;  
 LShellProb				 =0.082;
 MShellProb				 =0.012;
 exit_t					 =100*ns;		//lifetime of exitated nucleon Fe57

 //Fe57 isotope abundance in  asnoparticles;
 abund57=0.022;

 //spectrum calibration parameters
  G4double v_st256=0.094810*mm/s;  
  
  n0=127.290;							//zero channel
  N=256;								//number of channels
  v_st=v_st256*256/(N*1.0);				//channel price in velocity;

	//to output created particles
	SetVerboseLevel(0);

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MossbauerScattering::~MossbauerScattering() {
}


// ..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..
// 00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..ooOO

G4int MossbauerScattering::CalculateChannel (G4double energy){ 

	clight=300000000000* mm/s;		// light speed in mm/s
	G4double v=(energy-ResonanceEnergy)/ResonanceEnergy*clight;
	G4double nR=(v/v_st)+n0; 
	G4int n=std:: floor (nR);
	//G4cout <<" velocity "<<v/ (mm/s)<<"   "<<v_st/ (mm/s)<<G4endl;
	if ((nR-n)<(n+1-nR))  return n;
	else return n+1;
}

G4double MossbauerScattering::GetAmlplitude (G4int n) {
	if ((n>=0) && (n<N))
	return spectr [n];
	else return 0;
}

void MossbauerScattering::ConstructSpectrum () {
 
#include <iostream>
#include <fstream>
	
	std :: ifstream inFile;
	inFile.open("SpectrumY.txt");
	G4String xin;
	
	for (G4int i=0;i<N;i++)
	{
	 inFile >> xin;
	 spectr [i]=atof(xin);
	 //G4cout <<i<<"  d "<<xin<<"    "<< spectr [i]<<G4endl;
	}
	
	//for (G4int i=0;i<N;i++)
	 //spectr [i]= 0.99+(G4UniformRand()/100.0);  //random amlitude. temporarily. needs import from spectrum data file




// G4double spe [3]={ 1.0,1.0,1.0};
	 

 G4cout <<"spectrum for "<< N <<" channels constructed"<<G4endl; 
}
void MossbauerScattering::SetFe57abundance(G4double abund)
{	abund57=abund;	}


G4bool MossbauerScattering::IsApplicable(const G4ParticleDefinition& partDef)
      {
        return &partDef == G4Gamma::GammaDefinition();
      }
// ..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..
// 00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..ooOO

G4double MossbauerScattering::GetMeanFreePath
    (const G4Track& aTrack, G4double, G4ForceCondition*)
 {

  //constants definition

  //G4double cros0 = 2500000000*millibarn; 
  //?G4double cros0 = 25 ; //10^-18 cm^2
  //?G4double fm=0.276;	//mossbauer effect factor
  G4double A;  //amplitude effect
  G4double leng =DBL_MAX;  //interaction length

  //physical sample parameters (which has Mossbauer spectra)
  G4double enr2=0.022; // Fe57partial abundance (natural)
  G4double d2=0.01*mm;  //sample thickness
  G4double n2=5.4/(enr2*1000.0); // Fe (nonisotop) concentration in phys. sample (10^23 cm^-3)
  
  //model sample parameters
  G4double enr1=abund57; // Fe57partial abundance 
  //?G4double d1=10*nm;
  G4double n1=0.43;// iron nucleus concentaration (10^23 cm^-3)
   
//dynamics parameters definition
  const G4DynamicParticle* theParticle = aTrack.GetDynamicParticle();
  G4String theParticleName = theParticle->GetDefinition()->GetParticleName();
  G4String thePhysVolume = aTrack.GetVolume()->GetName();
  G4double energy=theParticle->GetKineticEnergy();
  G4int n=CalculateChannel (energy);
	  //G4bool reson=std::fabs(energy-ResonanceEnergy);
 
  A=GetAmlplitude (n);
 // G4double cros= cros0 *(1-A)*(n1*d1*enr1)/(n2*d2*enr2);
  if (A>0)
   leng= d2* (n2*enr2)/(n1*enr1)/A;
 // G4cout <<"leng "<< leng/mm<<G4endl;
  //_______________________________________________

  if(( theParticleName== "gamma") && ( thePhysVolume==ResonancePhysVolume) && ((n>=0) && (n<N))) { 
	//G4cout << "channel "<<n<<" ampl  "<<A<<G4endl;
		return leng;	}		//gamma in Fe57 volume and kinetic energy close to nuclear resonance level
  else {return DBL_MAX; }		//out of resonance
 
}

// ..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..
// 00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..oo00oo..ooOO


G4VParticleChange* MossbauerScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {

  aParticleChange.Initialize(aTrack);
  const G4DynamicParticle* theParticle = aTrack.GetDynamicParticle();
 // G4String theParticleType = theParticle->GetDefinition()->GetParticleType();
 // G4String theParticleName = theParticle->GetDefinition()->GetParticleName();
 // G4String thePhysVolume = aTrack.GetVolume()->GetName();
    G4double energy= theParticle->GetKineticEnergy();
  // step information
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double      t0 = pPreStepPoint->GetGlobalTime();
  G4double rand = G4UniformRand();
  G4double delta =aStep.GetStepLength() * rand;
  G4double deltaTime = delta /
      ((pPreStepPoint->GetVelocity()+pPostStepPoint->GetVelocity())/2.);
  deltaTime = deltaTime + exit_t * std::exp (- G4UniformRand() );
  //_______________________________________________
  // secondaries
  
  G4double aSecondaryTime = t0 + deltaTime;
  G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();
  aParticleChange.SetNumberOfSecondaries(1); 

  G4double cost,sint, phi, sinp, cosp;  //angles;
  G4double probchannel = G4UniformRand(); //probability calculation;
  G4double en=energy;
 
    if (probchannel<ElectronConversionProb) { 
	
	//electron conversion channel 
	 if (probchannel<KShellProb	)						//K shell works
		en=energy-(7.3*keV);
	 else if (probchannel-KShellProb<LShellProb)		// L shell works
	 { en=energy-((13.45+0.25*G4UniformRand())*keV);  }
	 else if (probchannel-KShellProb-LShellProb<MShellProb)
	 { en=energy-((14.31+0.099*G4UniformRand())*keV);	} // M shell works
	 //create electron
	 cost = 1. - 2.*G4UniformRand();
     sint = std::sqrt((1.-cost)*(1.+cost));
     phi = 2*pi*G4UniformRand();
     sinp = std::sin(phi);
     cosp = std::cos(phi);
	// random momentum direction
	G4ParticleMomentum electronMomentum(sint*cosp, sint*sinp, cost);
	// generate a new electron
	G4DynamicParticle* convElectron = new G4DynamicParticle (G4Electron::Electron(), electronMomentum);
	convElectron->SetKineticEnergy(en );
	// Generate new G4Track
	G4Track* aSecondaryTrack = new G4Track(convElectron,aSecondaryTime,aSecondaryPosition);
    aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    aParticleChange.AddSecondary(aSecondaryTrack);

	if (GetVerboseLevel() >= 1)
		G4cout <<"electron created "<<aSecondaryPosition << "  pos, energy "<<G4BestUnit(en, "Energy")<<G4endl;

	} else {								   
	// gamma irradiation channel
    // random momentum direction
    cost = 1. - 2.*G4UniformRand();
     sint = std::sqrt((1.-cost)*(1.+cost));
     phi = 2*pi*G4UniformRand();
     sinp = std::sin(phi);
     cosp = std::cos(phi);
     G4ParticleMomentum gammaMomentum(sint*cosp, sint*sinp, cost);
    // polarization
    G4ThreeVector photonPol(cost*cosp, cost*sinp, -sint);
    G4ThreeVector perp = gammaMomentum.cross(photonPol);
    phi = 2*pi*G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);
    photonPol = cosp * photonPol + sinp * perp;
    photonPol = photonPol.unit();
    // generate a new gamma
    G4DynamicParticle* aIrrGamma = new G4DynamicParticle (G4Gamma::Gamma(), gammaMomentum);
    aIrrGamma->SetKineticEnergy(energy );
    aIrrGamma->SetPolarization (photonPol.x(), photonPol.y(), photonPol.z());
    // Generate new G4Track
    G4Track* aSecondaryTrack = new G4Track(aIrrGamma,aSecondaryTime,aSecondaryPosition);
    aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    aParticleChange.AddSecondary(aSecondaryTrack);

		if (GetVerboseLevel() >= 1)
			G4cout <<"gamma created "<< aSecondaryPosition << "  pos, energy "<<G4BestUnit(en, "Energy")<<G4endl;

   }
	//kill initial gamma
	//G4Step* st= aParticleChange.UpdateStepForPostStep( & aStep);
	//aParticleChange.

	//G4VParticleChange::SetStatus ( fStopAndKill );
	aParticleChange. ProposeTrackStatus ( fStopAndKill );



  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


