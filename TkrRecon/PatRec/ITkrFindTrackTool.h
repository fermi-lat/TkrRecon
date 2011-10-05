/**
 * @class ITkrFindTrackTool
 *
 * @brief Interface to the track fitting tools. 
 * Basically, provides an interface to various single track fit 
 * classes. Currently there is but one, but this allows for future
 * ideas/expansion/etc.
 *
 * @author Tracy Usher
 */

#ifndef ITKRFINDTRACKTOOL_H
#define ITKRFINDTRACKTOOL_H

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IProperty.h"

static const InterfaceID IID_ITkrFindTrackTool("ITkrFindTrackTool", 7111 , 0);

class ITkrFindTrackTool : virtual public IAlgTool
{
 public:
  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ITkrFindTrackTool; }

  /// @brief Perform first pass track finding tasks (implementation dependent)
  virtual StatusCode firstPass()=0;

  /// @brief Perform second pass track finding tasks (implementation dependent)
  virtual StatusCode secondPass()=0;
};
#endif
