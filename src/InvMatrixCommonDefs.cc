/**   
    $Date: 2007/11/17 17:57:07 $
    $Revision: 1.1 $
    $Id: InvMatrixCommonDefs.cc,v 1.1 2007/11/17 17:57:07 govoni Exp $ 
    \author $Author: govoni $
*/

#include "Calibration/Tools/interface/InvMatrixCommonDefs.h"

int
ecalIM::uniqueIndex (int eta, int phi)
{
  return eta * SCMaxPhi + phi ;
}
