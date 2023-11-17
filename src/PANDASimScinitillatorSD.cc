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
//
/// \file PANDASimScinitillatorSD.cc
/// \brief Implementation of the PANDASimScinitillatorSD class

#include "PANDASimScinitillatorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PANDASimScinitillatorSD::PANDASimScinitillatorSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofArray)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
    fArraySize(nofArray)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PANDASimScinitillatorSD::~PANDASimScinitillatorSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PANDASimScinitillatorSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new PANDASimScinitillatorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits and fill  hits with zero energy deposition
  for (G4int i = 0; i < fArraySize; i++ ) 
  {
      for (G4int j = 0; j < fArraySize; j++)
      {
          fHitsCollection->insert(new PANDASimScinitillatorHit());
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PANDASimScinitillatorSD::ProcessHits(G4Step* step,
    G4TouchableHistory*)
{
    // energy deposit
    G4double edep = step->GetTotalEnergyDeposit();

    if (edep == 0.) return false;

    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();

    size_t depth = touchable->GetHistory()->GetDepth();

    // Get calorimeter cell id 
    G4int layerNumber = touchable->GetReplicaNumber(depth);

    // Get hit accounting data for this cell
    auto hit = (*fHitsCollection)[layerNumber];
    if (!hit) {
        G4ExceptionDescription msg;
        msg << "Cannot access hit " << layerNumber;
        G4Exception("PANDASimScinitillatorSD::ProcessHits()",
            "MyCode0004", FatalException, msg);
    }

    // Get hit for total accounting
    auto hitTotal
        = (*fHitsCollection)[fHitsCollection->entries() - 1];

    // Add values
    hit->Add(edep);
    hitTotal->Add(edep);

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PANDASimScinitillatorSD::EndOfEvent(G4HCofThisEvent*)
{
    if (verboseLevel > 1) {
        auto nofHits = fHitsCollection->entries();
        G4cout
            << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
        for (std::size_t i = 0; i < nofHits; ++i) (*fHitsCollection)[i]->Print();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
