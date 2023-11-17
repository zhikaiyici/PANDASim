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
// $Id: PANDASimSteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file PANDASimSteppingAction.cc
/// \brief Implementation of the PANDASimSteppingAction class

#include "PANDASimRun.hh"
#include "PANDASimSteppingAction.hh"
#include "PANDASimDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"

#include "UserDataInput.hh"
#include "Spline.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

PANDASimSteppingAction::PANDASimSteppingAction(/*PANDASimEventAction* eventAction, PANDASimRunAction* runAction*/)
	: G4UserSteppingAction(),
	/*fRunAction(runAction),*/ fEventAction(nullptr),
	fScoringVolume(nullptr), fPhotoelectricScoringVolume(nullptr), fGdFilmScoringVolume(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PANDASimSteppingAction::~PANDASimSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PANDASimSteppingAction::UserSteppingAction(const G4Step* step)
{
	if (!fScoringVolume) {
		const PANDASimDetectorConstruction* detectorConstruction
			= static_cast<const PANDASimDetectorConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fScoringVolume = detectorConstruction->GetScoringVolume();
	}

	if (!fGdFilmScoringVolume) {
		const PANDASimDetectorConstruction* detectorConstruction
			= static_cast<const PANDASimDetectorConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fGdFilmScoringVolume = detectorConstruction->GetGdFilmScoringVolume();
	}

	if (!fPhotoelectricScoringVolume) {
		const PANDASimDetectorConstruction* detectorConstruction
			= static_cast<const PANDASimDetectorConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fPhotoelectricScoringVolume = detectorConstruction->GetPhotoelectricScoringVolume();
	}

	if (!fEventAction)
		fEventAction = static_cast<PANDASimEventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());

	const G4StepPoint* postStepPoint = step->GetPostStepPoint();
	//const G4TouchableHandle& postTouch = postStepPoint->GetTouchableHandle();
	//const G4VPhysicalVolume* physPostVolume = postTouch->GetVolume();
	//if (physPostVolume == nullptr) return;
	//const G4LogicalVolume* logicPostVolume = physPostVolume->GetLogicalVolume();

	// get volume of the current step
	const G4StepPoint* preStepPoint = step->GetPreStepPoint();
	const G4TouchableHandle& preTouch = preStepPoint->GetTouchableHandle();
	const G4LogicalVolume* logicPreVolume = preTouch->GetVolume()->GetLogicalVolume();

	//G4String logicVolName = logicVolume->GetName();
	//G4String physVolumeName = physVolume->GetName();
	//G4cout << "logicVolName:" << logicVolName << G4endl;
	//G4cout << "physVolumeName: " << physVolumeName << G4endl;
	////getchar();

	// check if we are in scoring volume
	if (logicPreVolume != fScoringVolume && 
		logicPreVolume != fGdFilmScoringVolume &&
		logicPreVolume != fPhotoelectricScoringVolume ) return;

	G4Track* theTrack = step->GetTrack();
	const G4ParticleDefinition* particleDefinition = theTrack->GetDynamicParticle()->GetParticleDefinition();
	const G4String& strPrtclName = particleDefinition->GetParticleName();

	if (logicPreVolume == fGdFilmScoringVolume)
	{
		if (strPrtclName == "opticalphoton") return;

		const G4int parentID = theTrack->GetParentID();
		const G4TrackStatus trackStatus = theTrack->GetTrackStatus();
		// 判断是否为中子
		if (strPrtclName == "neutron" && parentID == 0)
		{
			const G4String processName = postStepPoint->GetProcessDefinedStep()->GetProcessName();
			//const G4String processName_ = preStepPoint->GetProcessDefinedStep()->GetProcessName();
			//G4cout << "processName: " << processName << G4endl;
			if (processName == "nCapture")
			{
				const G4int moduleRowReplicaNumber = preTouch->GetReplicaNumber(2);
				const G4int moduleRepliaNumber = preTouch->GetReplicaNumber(1);
				const G4double capTimeGd = postStepPoint->GetGlobalTime() / us;
				const G4double capTimeGd_ = theTrack->GetGlobalTime() / us;

				fEventAction->SetDelayFlagGd(true);
				fEventAction->CapTimeGd(moduleRepliaNumber, moduleRowReplicaNumber, capTimeGd);
				//G4cout << "processName: " << processName << G4endl;
				//G4cout << "processName pre: " << processName_ << G4endl;
				//G4cout << "capTime: " << capTimeGd << G4endl;
				//G4cout << "capTime_track: " << capTimeGd_ << G4endl;
				//getchar();

				//PANDASimRun* fPANDASimRun
				//	= static_cast<PANDASimRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
				//fPANDASimRun->PushBackCaptureTimeGd(capTimeGd);
			}
		}
	}

	if (logicPreVolume == fScoringVolume)
	{
		if (strPrtclName == "opticalphoton") return;

		const G4int parentID = theTrack->GetParentID();
		//const G4TrackStatus trackStatus = theTrack->GetTrackStatus();
		const G4String processName = postStepPoint->GetProcessDefinedStep()->GetProcessName();
		// 判断是否为中子 
		if (strPrtclName == "neutron" && parentID == 0)
		{
			// G4cout << "processName: " << processName << G4endl;
			if (processName == "nCapture")
			{
				const G4int moduleRowReplicaNumber = preTouch->GetReplicaNumber(2);
				const G4int moduleRepliaNumber = preTouch->GetReplicaNumber(1);
				const G4double capTimeH = postStepPoint->GetGlobalTime() / us;

				fEventAction->SetDelayFlagH(true);
				fEventAction->CapTimeH(moduleRepliaNumber, moduleRowReplicaNumber, capTimeH);
				//G4cout << "capTime: " << capTime << G4endl;
				//getchar();

				//PANDASimRun* fPANDASimRun
				//	= static_cast<PANDASimRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
				//fPANDASimRun->PushBackCaptureTimeH(capTimeH);
			}
		}
		else if (strPrtclName == "mu+" || strPrtclName == "mu-")
		{
			if (processName == "Decay")
			{
				const G4int moduleRowReplicaNumber = preTouch->GetReplicaNumber(2);
				const G4int moduleRepliaNumber = preTouch->GetReplicaNumber(1);
				//G4cout << "mouleRowReplicaNumber: " << moduleRowReplicaNumber 
				//	   << ". moduleRepliaNumber:" << moduleRepliaNumber << G4endl;
				//getchar();

				const G4double muDecayTime = postStepPoint->GetGlobalTime() / us;
				PANDASimRun* fPANDASimRun
					= static_cast<PANDASimRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
				fPANDASimRun->PushBackDecayTimeMu(muDecayTime);
				G4cout << "muDecayTime: " << muDecayTime << G4endl;
				G4cout << "strPrtclName: " << strPrtclName << G4endl;
				//getchar();
			}
		}

		// collect energy deposited in this step
		const G4double edepStep = step->GetTotalEnergyDeposit() / MeV;

		if (edepStep != 0.)
		{
			const G4int moduleRowReplicaNumber = preTouch->GetReplicaNumber(2);
			const G4int moduleRepliaNumber = preTouch->GetReplicaNumber(1);
			//G4cout << "mouleRowReplicaNumber: " << moduleRowReplicaNumber 
			//	   << ". moduleRepliaNumber:" << moduleRepliaNumber << G4endl;
			//getchar();

			G4bool delayFlagH = fEventAction->GetDelayFlagH();
			G4bool delayFlagGd = fEventAction->GetDelayFlagGd();
			G4bool decayFlagMu = fEventAction->GetDecayFlagMu();
			if (delayFlagH) // captured by H
				// delay signal
				fEventAction->AddEnergyDepositDelayH(moduleRepliaNumber, moduleRowReplicaNumber, edepStep);
			else if (delayFlagGd) // captured by Gd
				// delay signal
				fEventAction->AddEnergyDepositDelayGd(moduleRepliaNumber, moduleRowReplicaNumber, edepStep);
			else if (decayFlagMu)
				fEventAction->AddEnergyDepositDecayMu(moduleRepliaNumber, moduleRowReplicaNumber, edepStep);
			else
			{
				// prompt signal
				fEventAction->AddEnergyDeposit(moduleRepliaNumber, moduleRowReplicaNumber, edepStep);
				fEventAction->AddEdep(edepStep);
			}
		}
	}

	if (logicPreVolume == fPhotoelectricScoringVolume)
	{
		if (strPrtclName != "opticalphoton") return;
		//{
		//	G4cout << "not optical photon: " << strPrtclName << G4endl;
		//}

		const G4int moduleRowReplicaNumber = preTouch->GetReplicaNumber(4);
		const G4int moduleRepliaNumber = preTouch->GetReplicaNumber(3);
		const G4int PMTCopyNumber = preTouch->GetCopyNumber(2);

		G4double photonEnergy[] =
		{ 1.8338 * eV, 1.8419 * eV, 1.8527 * eV, 1.8721 * eV, 1.8890 * eV, 1.9172 * eV, 1.9530 * eV, 1.9800 * eV,
		  2.0022 * eV, 2.0413 * eV, 2.0845 * eV, 2.1479 * eV, 2.2163 * eV, 2.2922 * eV, 2.4194 * eV, 2.5563 * eV,
		  2.7037 * eV, 2.8891 * eV, 3.0268 * eV, 3.1703 * eV, 3.3728 * eV, 3.6556 * eV, 3.9353 * eV, 4.0806 * eV,
		  4.2007 * eV, 4.2506 * eV };
		const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

		G4double quatumnEfficiency[] =
		{ 0.0005, 0.0006, 0.0009, 0.0013, 0.0021, 0.0034, 0.0068, 0.0093, 0.0129, 0.0184, 0.0289, 0.0436,
		  0.0624, 0.0903, 0.1354, 0.1785, 0.2165, 0.2461, 0.2530, 0.2460, 0.2268, 0.1802, 0.1222, 0.0847,
		  0.0510, 0.0387 };
		assert(sizeof(quatumnEfficiency) == sizeof(photonEnergy));

		try
		{
			static G4ThreadLocal SplineSpace::SplineInterface* sp = new SplineSpace::Spline(photonEnergy, quatumnEfficiency, nEntries);
			const G4double photonKE = preStepPoint->GetKineticEnergy();
			G4double pheNum;
			sp->SinglePointInterp(photonKE, pheNum);
			//G4cout << "x =" << photonKE / eV << ": " << peNum << G4endl;
			if (pheNum != 0.)
			{
				G4bool delayFlagH = fEventAction->GetDelayFlagH();
				G4bool delayFlagGd = fEventAction->GetDelayFlagGd();
				G4bool decayFlagMu = fEventAction->GetDecayFlagMu();

				if (delayFlagH)
					fEventAction->AddModuleCalPhDelayH(moduleRepliaNumber, moduleRowReplicaNumber, PMTCopyNumber, pheNum);
				else if (delayFlagGd)
					fEventAction->AddModuleCalPhDelayGd(moduleRepliaNumber, moduleRowReplicaNumber, PMTCopyNumber, pheNum);
				else if(decayFlagMu)
					fEventAction->AddModuleCalPhDecayMu(moduleRepliaNumber, moduleRowReplicaNumber, PMTCopyNumber, pheNum);
				else
					fEventAction->AddModuleCalPh(moduleRepliaNumber, moduleRowReplicaNumber, PMTCopyNumber, pheNum);
			}
			theTrack->SetTrackStatus(fStopAndKill);
		}
		catch (SplineSpace::SplineFailure sf)
		{
			G4cout << sf.GetMessage() << G4endl;
		}

		/*G4OpBoundaryProcessStatus boundaryStatus = Undefined;
		static G4ThreadLocal G4OpBoundaryProcess* opticalBoundary = nullptr;
		//find the boundary process only once
		if (!opticalBoundary)
		{
			G4ProcessManager* processManager = particleDefinition->GetProcessManager();
			G4ProcessVector* processVector = processManager->GetProcessList();
			G4int nProcesses = processManager->GetProcessListLength();
			for (G4int i = 0; i < nProcesses; i++) {
				if ((*processVector)[i]->GetProcessName() == "OpBoundary")
				{
					opticalBoundary = (G4OpBoundaryProcess*)(*processVector)[i];
					break;
				}
			}
		}
		if (opticalBoundary)
		{
			boundaryStatus = opticalBoundary->GetStatus();
			if (boundaryStatus == Absorption)
			{
				//fEventAction->AddAllAbPh(PMTCopyNumber);
				fEventAction->AddModuleAbPh(moduleRowReplicaNumber, moduleRepliaNumber, PMTCopyNumber);
				//G4cout << "Absorption" << G4endl;
			}
			else if (boundaryStatus == Detection)
			{
				//fEventAction->AddAllDtPh(PMTCopyNumber);
				fEventAction->AddModuleDtPh(moduleRowReplicaNumber, moduleRepliaNumber, PMTCopyNumber);
				//G4cout << "Detection" << G4endl;
			}
		}*/
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
