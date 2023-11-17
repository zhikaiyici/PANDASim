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
/// \file PANDASimScinitillatorHit.hh
/// \brief Definition of the PANDASimScinitillatorHit class

#ifndef PANDASimScinitillatorHit_h
#define PANDASimScinitillatorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class PANDASimScinitillatorHit : public G4VHit
{
  public:
    PANDASimScinitillatorHit();
    PANDASimScinitillatorHit(const PANDASimScinitillatorHit&);
    virtual ~PANDASimScinitillatorHit();

    // operators
    const PANDASimScinitillatorHit& operator=(const PANDASimScinitillatorHit&);
    G4bool operator==(const PANDASimScinitillatorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de);

    // get methods
    G4double GetEdep() const;
      
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using PANDASimScinitillatorHitsCollection = G4THitsCollection<PANDASimScinitillatorHit>;

extern G4ThreadLocal G4Allocator<PANDASimScinitillatorHit>* PANDASimScinitillatorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PANDASimScinitillatorHit::operator new(size_t)
{
  if (!PANDASimScinitillatorHitAllocator) {
    PANDASimScinitillatorHitAllocator = new G4Allocator<PANDASimScinitillatorHit>;
  }
  void *hit;
  hit = (void *) PANDASimScinitillatorHitAllocator->MallocSingle();
  return hit;
}

inline void PANDASimScinitillatorHit::operator delete(void *hit)
{
  if (!PANDASimScinitillatorHitAllocator) {
    PANDASimScinitillatorHitAllocator = new G4Allocator<PANDASimScinitillatorHit>;
  }
  PANDASimScinitillatorHitAllocator->FreeSingle((PANDASimScinitillatorHit*) hit);
}

inline void PANDASimScinitillatorHit::Add(G4double de) {
  fEdep += de; 
}

inline G4double PANDASimScinitillatorHit::GetEdep() const { 
  return fEdep; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
