/**
 * @class TkrLinkAndTreeFitTool
 *
 * @brief Implements a Gaudi Tool for performing a track fit. This version uses 
 *        candidate tracks from the Link and Tree Pattern Recognition which are then fit
 *        with KalFitTrack 
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrLinkAndTreeFitTool.h,v 1.2 2002/08/29 19:18:47 usher Exp $
 */

#ifndef TKRLINKANDTREEFITTOOL_H
#define TKRLINKANDTREEFITTOOL_H

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "TkrRecon/Track/ITkrFitTool.h"
#include "TkrUtil/ITkrGeometrySvc.h"

class TkrLinkAndTreeFitTool : public AlgTool, virtual public ITkrFitTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrLinkAndTreeFitTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrLinkAndTreeFitTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode doTrackFit(Event::TkrPatCand* patCand);

private:
    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc* pTkrGeoSvc;

    /// Pointer to the Gaudi data provider service
    DataSvc*        pDataSvc;
};

#endif
