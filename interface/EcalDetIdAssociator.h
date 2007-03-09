#ifndef HTrackAssociator_HEcalDetIdAssociator_h
#define HTrackAssociator_HEcalDetIdAssociator_h 1
// -*- C++ -*-
//
// Package:    HTrackAssociator
// Class:      HEcalDetIdAssociator
// 
/*

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dmytro Kovalskyi
// Modified for ECAL+HCAL by Michal Szleper
//

#include "Calibration/Tools/interface/CaloDetIdAssociator.h"

class HEcalDetIdAssociator: public HCaloDetIdAssociator{
 public:
   HEcalDetIdAssociator():HCaloDetIdAssociator(180,150,0.04){};
 protected:
   virtual std::set<DetId> getASetOfValidDetIds(){
      std::set<DetId> setOfValidIds;
      std::vector<DetId> vectOfValidIds = geometry_->getValidDetIds(DetId::Ecal, 1);//EB
      for(std::vector<DetId>::const_iterator it = vectOfValidIds.begin(); it != vectOfValidIds.end(); ++it)
         setOfValidIds.insert(*it);

      vectOfValidIds.clear();
      vectOfValidIds = geometry_->getValidDetIds(DetId::Ecal, 2);//EE
      for(std::vector<DetId>::const_iterator it = vectOfValidIds.begin(); it != vectOfValidIds.end(); ++it)
         setOfValidIds.insert(*it);

      return setOfValidIds;

   }
};
#endif
