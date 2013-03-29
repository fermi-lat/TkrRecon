/**
 * @class ITkrFitTool
 *
 * @brief Interface to the track fitting tools. 
 * Basically, provides an interface to various single track fit 
 * classes. Currently there is but one, but this allows for future
 * ideas/expansion/etc.
 *
 * @author Tracy Usher
 */

#ifndef ITKRFITTOOL_H
#define ITKRFITTOOL_H

#include "GaudiKernel/IAlgTool.h"

namespace Event
{
    class TkrTrack;
    class TkrTrackHit;
}

static const InterfaceID IID_ITkrFitTool("ITkrFitTool", 7111 , 0);

class ITkrFitTool : virtual public IAlgTool 
{
 public:
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrFitTool; }

    /// @brief Given a pattern track, perform the track fit
    virtual StatusCode doTrackFit(Event::TkrTrack* patCand)=0;

    /// @brief Given a pattern track, perform the track re-fit
    virtual StatusCode doTrackReFit(Event::TkrTrack* patCand)=0;

    /// @brief Method to set (get) type of hit energy loss for a track
    virtual void         setHitEnergyLoss(const std::string& energyLossType) = 0; 
    virtual std::string  getHitEnergyLoss() = 0; 

    /// @brief Method to set (get) the particle type assumed for a track
    virtual void         setParticleType(const std::string ParticleName) = 0;
    virtual std::string  getParticleType() = 0;

    /// @brief Method to set method for determing cluster errors in fit
    virtual void       setClusErrCompType(const std::string& clusErrorType) = 0; 

    /// @brief Method to set multiple scattering matrix computation
    virtual void       setMultipleScatter(const bool doMultScatComp, const std::string& procNoiseType) = 0; 

    /// @brief Method to set Kalman Filter projection matrix type
    virtual void       setProjectionMatrix(const bool measOnly) = 0; 

    /// @brief This method does a "Filter" fit only - sets Filtered Chi-Square in track
    virtual void       doFilterFit(Event::TkrTrack& track) = 0;

    /// @brief Same as above but allows for "kinks" in the track
    virtual void       doFilterFitWithKinks(Event::TkrTrack& track) = 0;

    /// @brief This method does a "Smoother" fit only - sets Smoothed Chi-Square in track
    virtual void       doSmootherFit(Event::TkrTrack& track) = 0;

    /// @brief This method runs the filter for the next hit
    virtual double     doFilterStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& filterHit)
                        = 0; //{return 0.;}

    /// @brief This method runs the smoother for the next hit
    virtual double     doSmoothStep(Event::TkrTrackHit& referenceHit, Event::TkrTrackHit& smoothHit)
                        = 0; //{return 0.;}

};
#endif
