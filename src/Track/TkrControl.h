/**
 * @class TkrControl
 *
 * @brief A singleton class which contains control parameters 
 * for the tracking reconstruction.

 * Class is intended to be instantiated first by TkrInitSvc which 
 * has the ability to change parameters via job options parameters.
 * Once instantiated, used by pattern recognition and track fitting. 
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrControl.h,v 1.6 2002/09/05 13:48:49 burnett Exp $
 */


#ifndef TKRCONTROL_H
#define TKRCONTROL_H

//############################################
class TkrControl
//############################################
{
public:

    /// Define as a singleton object
    static TkrControl* getPtr();

    /// Retrieve values
    const int    getMaxCandidates()      {return m_maxCandidates;     }
    const int    getMinTermHitCount()    {return m_minTermHitCount;   }
    const double getFEneParticle()       {return m_fEneParticle;      }
    const double getSigmaCut()           {return m_sigmaCut;          }
    const double getMinEnergy()          {return m_minEnergy;         }
    const int    getMaxConsecutiveGaps() {return m_maxConsecutiveGaps;}
    const int    getMinSegmentHits()     {return m_minSegmentHits;    }
    const double getMaxChisqCut()        {return m_maxChiSqCut;       }
    const double getIniErrSlope()        {return m_iniErrorSlope;     }
    const double getIniErrPosition()     {return m_iniErrorPosition;  }
    const bool   getPlaneEnergies()      {return m_planeEnergies;     }


    /// Allow for control variables to be set at initialization
    void setMaxCandidates(  int    maxCand)   {m_maxCandidates     = maxCand; }
    void setMinTermHitCount(int    minCount)  {m_minTermHitCount   = minCount;}
    void setFEneParticle(   double enePart)   {m_fEneParticle      = enePart; }
    void setSigmaCut(       double sigmaCut)  {m_sigmaCut          = sigmaCut;}
    void setMinEnergy(     double minEnergy) {m_minEnergy         = minEnergy;}
    void setMaxConsGaps(    int    maxGaps)   {m_maxConsecutiveGaps = maxGaps;}
    void setMinSegmentHits( int    minHits)   {m_minSegmentHits    = minHits; }
    void setMaxChisqCut(    double maxChi)    {m_maxChiSqCut       = maxChi;  }
    void setIniErrSlope(    double errSlp)    {m_iniErrorSlope     = errSlp;  }
    void setIniErrPos(      double errPos)    {m_iniErrorPosition  = errPos;  }
    void setPlaneEnergies(  bool   enePlane)  {m_planeEnergies     = enePlane;}

private:
    /// private Constructor
    TkrControl();
 

    /// Pointer to the singleton object
    static TkrControl* m_this;

    /// Data members
    int    m_maxCandidates;      // Max number of Candidates 
    int    m_minTermHitCount;    // Number of hits to terminate Combo PR

    double m_fEneParticle;       // Fraction of Cal energy to use in PR.

    double m_sigmaCut;           // PR search window (in  sigmas)
    double m_minEnergy;              // Min tracking energy (MeV)

    int    m_maxConsecutiveGaps; // Max consecutive Gaps - Stop
    int    m_minSegmentHits;     // Min number of hits for segment
    double m_maxChiSqCut;        // Max allow PR Chisq. 
    double m_iniErrorSlope;      // First Hit error in Kalman: 10 deg 
    double m_iniErrorPosition;   // First Hit error in Kalman: .1 mm

    bool   m_planeEnergies;      // Decrease particle energies by exp(-rad_len)
};

#endif
