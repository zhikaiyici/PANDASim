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
// $Id: PANDASimPrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PANDASimPrimaryGeneratorAction.cc
/// \brief Implementation of the PANDASimPrimaryGeneratorAction class

#include "PANDASimPrimaryGeneratorAction.hh"
#include "Spline.hh"

#include "CalculateReferencePoints.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

PANDASimPrimaryGeneratorAction::PANDASimPrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fParticle(nullptr), fPositron(nullptr), fNeutron(nullptr),
	fMuonN(nullptr), fMuonP(nullptr), fProton(nullptr), fGamma(nullptr),
	fParticleGun(nullptr), fParticleGunP(nullptr),
	arraySize(0), neutrinoPosition(0),
	scinitillatorXHalfLength(0.), scinitillatorYHalfLength(0.), scinitillatorZHalfLength(0.),
	containerXHalfLength(0.), containerYHalfLength(0.), containerZHalfLength(0.),
	distanceBetweenModules(0.) , sourceType("NEUTRINO"), sourcePosition("CENTER"), referencePoints(0)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);
	fParticleGunP = new G4ParticleGun(n_particle);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	fPositron = particleTable->FindParticle(particleName = "e+");
	fNeutron = particleTable->FindParticle(particleName = "neutron");
	fMuonN = particleTable->FindParticle(particleName = "mu-");
	fMuonP = particleTable->FindParticle(particleName = "mu+");
	fGamma = particleTable->FindParticle(particleName = "gamma");
	fProton = particleTable->FindParticle(particleName = "proton");

	G4LogicalVolumeStore* logicVolStroe = G4LogicalVolumeStore::GetInstance();
	G4LogicalVolume* logicPlasticScintillator = logicVolStroe->GetVolume("PlasticScintillatorLV");
	G4Box* scinitillatorBox = dynamic_cast<G4Box*>(logicPlasticScintillator->GetSolid());
	scinitillatorXHalfLength = scinitillatorBox->GetXHalfLength();
	scinitillatorYHalfLength = scinitillatorBox->GetYHalfLength();
	scinitillatorZHalfLength = scinitillatorBox->GetZHalfLength();

	auto containerLV = logicVolStroe->GetVolume("ContainerLV");
	G4Box* containerBox = dynamic_cast<G4Box*>(containerLV->GetSolid());
	containerXHalfLength = containerBox->GetXHalfLength();
	containerYHalfLength = containerBox->GetYHalfLength();
	containerZHalfLength = containerBox->GetZHalfLength();

	sourceType = UserDataInput::GetSoureType();
	sourcePosition = UserDataInput::GetSourePosition();
	neutrinoPosition = UserDataInput::GetPositionOfNeutrino();
	arraySize = UserDataInput::GetSizeOfArray();
	distanceBetweenModules = UserDataInput::GetDistanceBetweenModules();
	neutrinoPercentage = UserDataInput::GetNeutrinoPercentage();
	neutronEnergy = UserDataInput::GetNeutronEnergy();
	neutronCDFSpectrum = UserDataInput::GetNeutronCDFSpectrum();
	positronEnergy = UserDataInput::GetPositronEnergy();
	positronCDFSpectrum = UserDataInput::GetPositronCDFSpectrum();

	CalculateReferencePoints refP;
	referencePoints = refP.GetRefrencePoints();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PANDASimPrimaryGeneratorAction::~PANDASimPrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fParticleGunP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PANDASimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//this function is called at the begining of ecah event
	//

	G4ThreeVector positionVector, directionVector = G4ThreeVector();
	G4double primaryParticleEnergy = 0. * MeV;

	if (sourceType != "NEUTRINO")
	{
		if (sourceType == "Am-Be")
		{
			fParticle = fNeutron;
			primaryParticleEnergy = EnergySampling(neutronCDFSpectrum, neutronCDFSpectrum);

			G4double phi = twopi * G4UniformRand();
			G4double costheta = -1. * G4UniformRand();
			G4double sintheta = sqrt(1 - pow(costheta, 2));
			directionVector = G4ThreeVector(sintheta * cos(phi), sintheta * sin(phi), costheta);

			fParticleGun->SetParticleDefinition(fParticle);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			fParticleGun->SetParticleMomentumDirection(directionVector);
		}
		else if (sourceType == "Cs137g")
		{
			fParticle = fGamma;
			primaryParticleEnergy = 0.66166 * MeV;
			G4double costheta = 2. * G4UniformRand() - 1.;
			G4double sintheta = sqrt(1 - pow(costheta, 2));
			G4double phi = twopi * G4UniformRand();
			directionVector = G4ThreeVector(sintheta * cos(phi), sintheta * sin(phi), costheta);

			fParticleGun->SetParticleDefinition(fParticle);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			fParticleGun->SetParticleMomentumDirection(directionVector);
		}
		else if (sourceType == "Co60g")
		{
			fParticle = fGamma;
			if (G4UniformRand() > 0.5)
				primaryParticleEnergy = 1.1732 * MeV;
			else
				primaryParticleEnergy = 1.3325 * MeV;
			G4double costheta = 2. * G4UniformRand() - 1.;
			G4double sintheta = sqrt(1 - pow(costheta, 2));
			G4double phi = twopi * G4UniformRand();
			directionVector = G4ThreeVector(sintheta * cos(phi), sintheta * sin(phi), costheta);

			fParticleGun->SetParticleDefinition(fParticle);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			fParticleGun->SetParticleMomentumDirection(directionVector);
		}

		if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino())
		{
			G4IonTable* ionTable = G4IonTable::GetIonTable();
			if (sourceType == "Co60")
			{
				fParticle = ionTable->GetIon(27, 60);
				primaryParticleEnergy = 0. * MeV;
				directionVector = G4ThreeVector();
			}
			else if (sourceType == "Cs137")
			{
				fParticle = ionTable->GetIon(55, 137);
				primaryParticleEnergy = 0. * MeV;
				directionVector = G4ThreeVector();
			}
			else if (sourceType == "Na22")
			{
				fParticle = ionTable->GetIon(11, 22);
				primaryParticleEnergy = 0. * MeV;
				directionVector = G4ThreeVector();
			}
			fParticleGun->SetParticleDefinition(fParticle);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			fParticleGun->SetParticleMomentumDirection(directionVector);
			fParticleGun->SetParticleCharge(0. * eplus);
		}

		if (sourcePosition == "LEFT")
		{
			positionVector = G4ThreeVector(-(scinitillatorXHalfLength - distanceBetweenModules), 0., containerZHalfLength);
		}
		else if (sourcePosition == "RIGHT")
		{
			positionVector = G4ThreeVector(scinitillatorXHalfLength - distanceBetweenModules, 0., containerZHalfLength);
		}
		else
		{
			positionVector = G4ThreeVector(0., 0., containerZHalfLength);
		}
		fParticleGun->SetParticlePosition(positionVector);
		fParticleGun->GeneratePrimaryVertex(anEvent);
	}
	else
	{
		G4double random = G4UniformRand();
		if (random < neutrinoPercentage)
		{
			fParticleGun->SetParticleDefinition(fNeutron);
			primaryParticleEnergy = EnergySampling(neutronEnergy, neutronCDFSpectrum);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			Sampling(positionVector, directionVector);
			fParticleGun->SetParticlePosition(positionVector);
			fParticleGun->SetParticleMomentumDirection(directionVector);
			fParticleGun->GeneratePrimaryVertex(anEvent);

			fParticleGun->SetParticleDefinition(fPositron);
			primaryParticleEnergy = EnergySampling(positronEnergy, positronCDFSpectrum);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			fParticleGun->SetParticlePosition(positionVector);
			Sampling(positionVector, directionVector);
			fParticleGun->SetParticleMomentumDirection(directionVector);
			fParticleGun->GeneratePrimaryVertex(anEvent);
		}
		else
		{
			G4double randmu = G4UniformRand();
			if (randmu < 5. / 11) // mu+/mu- = 1.2
			{
				fParticle = fMuonN;
			}
			else
			{
				fParticle = fMuonP;
			}
			primaryParticleEnergy = (2. * G4UniformRand() + 3.) * GeV;// 4. * GeV;
			fParticleGun->SetParticleDefinition(fParticle);
			fParticleGun->SetParticleEnergy(primaryParticleEnergy);
			Sampling(positionVector, directionVector);
			fParticleGun->SetParticlePosition(positionVector);
			fParticleGun->SetParticleMomentumDirection(directionVector);
			fParticleGun->GeneratePrimaryVertex(anEvent);
		}
	}
}


G4double PANDASimPrimaryGeneratorAction::EnergySampling(vector<G4double> energy, vector<G4double> cdfSpectrum)
{
	G4double random = G4UniformRand();
	G4int i = 0;
	while (random > cdfSpectrum[i])
	{
		++i;
	}
	return energy[i] * MeV;
	/*while (random < cdfSpectrum[0] || random > *(cdfSpectrum.end() - 1))
	{
		random = G4UniformRand();
		i++;
		if (i > 1e6)
		{
			G4ExceptionDescription msg;
			msg << "Cannot sample primaryParticleEnergy. 1e6 times samplings have been tried." << G4endl;
			msg << "Break the while loop now.";
			G4Exception("PANDASimPrimaryGeneratorAction::EnergySampling()", 
				"while (random < cdfSpectrum[0] || random > *(cdfSpectrum.end() - 1))", JustWarning, msg);
			break;
		}
	}
	try
	{
		SplineSpace::SplineInterface* sp = new SplineSpace::Spline(&cdfSpectrum[0], &energy[0], energy.size());
		sp->SinglePointInterp(random, primaryParticleEnergy);
		//G4cout << "random =: " << primaryParticleEnergy << G4endl;
	}
	catch (SplineSpace::SplineFailure sf)
	{
		G4cout << sf.GetMessage() << G4endl;
		//primaryParticleEnergy = 15.75 * keV;
	}
	return primaryParticleEnergy;*/
}

void PANDASimPrimaryGeneratorAction::Sampling(G4ThreeVector& positionVector, G4ThreeVector& directionVector)
{
	G4double sourcePositionX = 0.;
	G4double sourcePositionY = 0.;
	G4double sourcePositionZ = 0.;// 2. * scinitillatorZHalfLength * G4UniformRand() - scinitillatorZHalfLength;

	G4double theta = pi * G4UniformRand();
	G4double phi = twopi * G4UniformRand();

	auto thisParticle = fParticleGun->GetParticleDefinition();
	//auto thisParticleP = fParticleGunP->GetParticleDefinition();
	if (thisParticle == fNeutron || thisParticle == fPositron)
	{
		G4int randi;
		if (neutrinoPosition < arraySize * arraySize)
		{
			randi = neutrinoPosition;
		}
		else
		{
			randi = (G4int)(arraySize * arraySize * G4UniformRand());
		}
		G4double minZ = referencePoints[randi][0] - scinitillatorZHalfLength;
		G4double maxZ = referencePoints[randi][0] + scinitillatorZHalfLength;
		G4double minY = referencePoints[randi][1] - scinitillatorYHalfLength;
		G4double maxY = referencePoints[randi][1] + scinitillatorYHalfLength;

		//G4cout << minX << ", " << maxZ << "; " << maxX << ", " << maxZ << G4endl;
		//G4cout << minX << ", " << minZ << "; " << maxX << ", " << minZ << G4endl;
		//getchar();

		sourcePositionX = 2. * scinitillatorXHalfLength * G4UniformRand() - scinitillatorXHalfLength;
		sourcePositionY = minY + (maxY - minY) * G4UniformRand();
		sourcePositionZ = minZ + (maxZ - minZ) * G4UniformRand();

		G4double costheta = 2. * G4UniformRand() - 1.;
		G4double sintheta = sqrt(1 - pow(costheta, 2));
		directionVector = G4ThreeVector(sintheta * cos(phi), sintheta * sin(phi), costheta);
	}
	else if (thisParticle == fMuonP || thisParticle == fMuonN)
	{
		sourcePositionX = 2. * containerXHalfLength * G4UniformRand() - containerXHalfLength;
		sourcePositionY = 2. * containerYHalfLength * G4UniformRand() - containerYHalfLength;
		sourcePositionZ = containerZHalfLength;
		 
		// for mu- / mu+, f(theta) = 4/pi * cos^2(theta) (pi/2, pi),替换抽样，抽样效率1/2
		G4double random = G4UniformRand();
		G4double candidateTheta = 0.5 * pi * (G4UniformRand() + 1);
		G4double cos2theta = pow(cos(candidateTheta), 2);
		G4int i = 0;
		while (random > cos2theta)
		{
			candidateTheta = 0.5 * pi * (G4UniformRand() + 1);
			cos2theta = pow(cos(candidateTheta), 2);
			i++;
			if (i > 1e6)
			{
				G4ExceptionDescription msg;
				msg << "Cannot sample theta. 1e6 times samplings have been tried." << G4endl;
				msg << "Break the while loop now.";
				G4Exception("PANDASimPrimaryGeneratorAction::Sampling()", "while (random > cos2theta)", JustWarning, msg);
				break;
			}
		}
		theta = candidateTheta;
		directionVector = G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
	}
	positionVector = G4ThreeVector(sourcePositionX, sourcePositionY, sourcePositionZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
