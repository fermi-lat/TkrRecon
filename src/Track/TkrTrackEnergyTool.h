/**
 * @class TkrTrackEnergyTool
 *
 * @brief Implements a Gaudi Tool for setting the candidate track energies before 
 *        the track fit
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrTrackEnergyTool.h,v 1.3 2004/10/12 19:03:38 lsrea Exp $
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"

#include "src/Track/TkrControl.h"
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
    //inline Point  getPosAtZ(const Event::TkrPatCand* track, double deltaZ)const
    //            {return track->getPosition() + track->getDirection() * deltaZ;} 

    //double getTotalEnergy(Event::TkrPatCand* track, double CalEnergy);

    /// Pointer to the local Tracker geometry service
    //ITkrGeometrySvc*       m_tkrGeom;

    /// Pointer to the cluster tool
    //ITkrQueryClustersTool* m_clusTool;

    //TkrControl*            m_control;

    /// Pointer to the Gaudi data provider service
    //DataSvc*         m_dataSvc;
};
