/**
 * @class ITkrMcPatCandTool
 *
 * @brief Interface to the tool for returning information on MC/Recon Pattern Recognition tables 
 *        Basically, the algorithm TkrBuildMcRelationsAlg will create a series of relational tables
 *        which relate Monte Carlo track information to the recon. This tool will help return 
 *        specific information with respect to those relations for the pattern recognition.
 *
 * @author Tracy Usher
 */

#ifndef ITkrMcPatCandTOOL_H
#define ITkrMcPatCandTOOL_H

#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_ITkrMcPatCandTool("ITkrMcPatCandTool", 2 , 0);

class ITkrMcPatCandTool : virtual public IAlgTool 
{
 public:
  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ITkrMcPatCandTool; }

  /// @brief Return the number of Monte Carlo tracks
  virtual int getNumMcTracks()=0;

};
#endif
