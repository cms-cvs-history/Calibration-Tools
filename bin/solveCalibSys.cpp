#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "Calibration/Tools/interface/calibXMLwriter.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibTools.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibMapEcal.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalBarrel.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalEndcap.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Calibration/Tools/interface/matrixSaver.h"
#include "Calibration/EcalCalibAlgos/interface/BlockSolver.h"
#include "boost/shared_ptr.hpp"

//#include "Calibration/EcalAlCaRecoProducers/interface/trivialParser.h"
//#include "Calibration/EcalAlCaRecoProducers/bin/trivialParser.h"

#include "TH2.h"
#include "TProfile.h"
#include "TH1.h"
#include "TFile.h"

#include "CLHEP/Matrix/GenMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#define PI_GRECO 3.14159265

std::auto_ptr<CLHEP::HepMatrix>
sumMatrices (std::vector<std::string> matrices, 
             std::vector<double> weights)
{
/* 
1. get matrix size
2. create mtr
3. read and fill
4. return
*/

//  std::auto_ptr<CLHEP::HepMatrix> = new CLHEP::HepMatrix () ;

  matrixSaver leggo ;

  std::vector<std::string>::const_iterator mtrIt = matrices.begin () ;    
  CLHEP::HepMatrix * chi2Mtr = dynamic_cast<CLHEP::HepMatrix *> (
    leggo.getMatrix (*mtrIt)) ;
  chi2Mtr->num_row () ;
  std::auto_ptr<CLHEP::HepMatrix> sumMatrix (new CLHEP::HepMatrix (
    chi2Mtr->num_row (), chi2Mtr->num_row (), 0)) ;
  (*sumMatrix) += (*chi2Mtr) ;
  delete chi2Mtr ;

  //PG loop on matrices
  for (++mtrIt ;
       mtrIt != matrices.end () ;
       ++mtrIt)
    {
      CLHEP::HepMatrix * chi2Mtr = dynamic_cast<CLHEP::HepMatrix *> (
        leggo.getMatrix (*mtrIt)) ;
      (*sumMatrix) += (*chi2Mtr) ;
      delete chi2Mtr ;
    } //PG loop on matrices
    
  return sumMatrix ;
}              


// ------------------------------------------------------------------------


inline int etaShifter (const int etaOld) 
   {
     if (etaOld < 0) return etaOld + 85 ;
     else if (etaOld > 0) return etaOld + 84 ;
   }


// ------------------------------------------------------------------------


int main (int argc, char* argv[]) 
{

  if(argc != 2) {
    std::cout << "Usage: edmParameterSetDump <cfgfile>" << std::endl;
  }
  std::string fileName (argv[1]) ;
  boost::shared_ptr<edm::ProcessDesc> processDesc = edm::readConfigFile (fileName) ;

  boost::shared_ptr<edm::ParameterSet> pSet = processDesc->getProcessPSet () ;
  std::cout << pSet->dump () << std::endl;

  //PG get the names of the matrices and vectors
  typedef std::vector<std::string> namelist ;
  namelist chi2Matrices = pSet->getParameter<namelist> ("chi2Matrices") ;
  namelist chi2Vectors = pSet->getParameter<namelist> ("chi2Vectors") ;

  if (chi2Matrices.size () != chi2Vectors.size ())
    {
      std::cerr << "size of \"chi2Matrices\" and \"chi2Vectors\" differ" << std::endl ;
      return 1 ; 
    }

  std::vector<double> chi2Weights = pSet->getParameter<std::vector<double> > ("chi2Weights") ;
  if (chi2Matrices.size () != chi2Weights.size ())
    {
      std::cerr << "size of \"chi2Matrices\" and \"chi2Weights\" differ" << std::endl ;
      return 1 ; 
    }



  if (argc < 3) return 1 ;
  std::string chi2MtrFile = argv[0] ; 
  std::string chi2VtrFile = argv[1] ; 
  std::string cfgFile = argv[2] ;

  matrixSaver leggo ;
  
  CLHEP::HepMatrix * chi2Mtr = dynamic_cast<CLHEP::HepMatrix *> (
    leggo.getMatrix (chi2MtrFile)) ;
  CLHEP::HepVector * chi2Vtr = dynamic_cast<CLHEP::HepVector *> (
    leggo.getMatrix (chi2VtrFile)) ;

  double min = 0.5 ; //FIXME
  double max = 1.5 ; //FIXME
  bool usingBlockSolver = false ; //FIXME
  int region = 0 ; //FIXME
    
  CLHEP::HepVector result = CLHEP::solve (*chi2Mtr,*chi2Vtr) ;
  if (result.normsq () < min * chi2Mtr->num_row () ||
      result.normsq () > max * chi2Mtr->num_row ()) 
    {
    if (usingBlockSolver)  
      {
        edm::LogWarning ("IML") << "using  blocSlover " << std::endl ;
        BlockSolver() (*chi2Mtr,*chi2Vtr,result) ;
      }
    else 
      {
        edm::LogWarning ("IML") <<"coeff out of range " <<std::endl;
        for (int i = 0 ; i < chi2Vtr->num_row () ; ++i)
              result[i] = 1. ;
      }
    }
        
  unsigned int numberOfElements = chi2Mtr->num_row () ;
  std::map<unsigned int, float> coefficients ;
  for (unsigned int i=0; i < numberOfElements; ++i)
    coefficients[i] = result[i] ;

  if (region == 0) //PG EB
    {
//      int index = 
    
    }
  else if (region = -1) //PG EB-
    {}
  else //PG EB+ 
    {}    
    
  //FIXME salva la mappa
  return 0 ;


}
