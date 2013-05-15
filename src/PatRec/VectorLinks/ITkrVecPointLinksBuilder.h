/** @file ITkrVecPointLinksBuilder.h
 * @class ITkrVecPointLinksBuilder
 *
 * @brief Interface class to the TkrVecPointsLink builder
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/VectorLinks/ITkrVecPointLinksBuilder.h,v 1.2 2012/04/25 04:54:36 heather Exp $
 *
*/

#ifndef __ITkrVecPointsLinkBuilder_H
#define __ITkrVecPointsLinkBuilder_H 1


#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/IProperty.h"

static const InterfaceID IID_ITkrVecPointsLinkBuilder("ITkrVecPointsLinkBuilder", 7111 , 0);

#include "Event/Recon/TkrRecon/TkrVecPointsLinkInfo.h"

class ITkrVecPointsLinkBuilder : virtual public IAlgTool
{
 public:
  /// Retrieve interface ID
  static const InterfaceID& interfaceID() { return IID_ITkrVecPointsLinkBuilder; }

  /// @brief Retrieve pointer to links collection which contains only single layer links
  virtual Event::TkrVecPointsLinkInfo* getSingleLayerLinks(const Point&  refPoint, 
                                                           const Vector& refAxis,
                                                           double        refError = 1000.,
                                                           double        energy   = 30.) = 0;

  /// @brief Retrieve pointer to links collection for all possible links
  virtual Event::TkrVecPointsLinkInfo* getAllLayerLinks(const Point&  refPoint, 
                                                        const Vector& refAxis,
                                                        double        refError = 1000.,
                                                        double        angError = M_PI,
                                                        double        energy   = 30.)    = 0;
};

#endif
