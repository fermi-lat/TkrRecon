/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/Track/TkrTrackEnergyTool.h,v 1.2 2004/09/08 15:32:45 usher Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"

//#include "TkrRecon/GaudiAlg/TkrTrackFitAlg.h"
#include "src/Track/TkrControl.h"
#include "src/TrackFit/KalFitTrack/KalFitter.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "TkrUtil/ITkrQueryClustersTool.h"

static const InterfaceID IID_TkrTrackEnergyTool("TkrTrackEnergyTool", 1 , 0);

class TkrTrackEnergyTool : public AlgTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrTrackEnergyTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrTrackEnergyTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode SetTrackEnergies(double totalEnergy);

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_TkrTrackEnergyTool; }

private:
    /// Internal methods
    inline Point  getPosAtZ(const Event::TkrPatCand* track, double deltaZ)const
                {return track->getPosition() + track->getDirection() * deltaZ;} 

    double getTotalEnergy(Event::TkrPatCand* track, double CalEnergy);

    /// Pointer to the local Tracker geometry service
    ITkrGeometrySvc*       m_tkrGeo;

    /// Pointer to the cluster tool
    ITkrQueryClustersTool* m_clusTool;

    TkrControl*            m_control;

    /// Pointer to the Gaudi data provider service
    DataSvc*         m_dataSvc;
};
