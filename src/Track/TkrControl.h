/**
 * @class TkrControl
 *
 * @brief A singleton class which contains control parameters for the tracking reconstruction.
 *        Class is intended to be instantiated first by TkrInitSvc which has the ability to 
 *        change parameters via job options parameters. Once instantiated, used by pattern 
 *        recognition and track fitting. 
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Track/TkrControl.h,v 1.5 2002/09/05 04:15:50 lsrea Exp $
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
    const int    getMaxCandidates()      {return m_MaxCandidates;     }
    const int    getMinTermHitCount()    {return m_MinTermHitCount;   }
    const double getFEneParticle()       {return m_FEneParticle;      }
    const double getSigmaCut()           {return m_SigmaCut;          }
    const double getMinEnergy()          {return m_MinEnergy;         }
    const int    getMaxConsecutiveGaps() {return m_MaxConsecutiveGaps;}
    const int    getMinSegmentHits()     {return m_MinSegmentHits;    }
    const double getMaxChisqCut()        {return m_MaxChiSqCut;       }
    const double getIniErrSlope()        {return m_IniErrorSlope;     }
    const double getIniErrPosition()     {return m_IniErrorPosition;  }
    const bool   getPlaneEnergies()      {return m_PlaneEnergies;     }


    /// Allow for control variables to be set at initialization
    void setMaxCandidates(  int    maxCand)   {m_MaxCandidates      = maxCand;  }
    void setMinTermHitCount(int    minCount)  {m_MinTermHitCount    = minCount; }
    void setFEneParticle(   double enePart)   {m_FEneParticle       = enePart;  }
    void setSigmaCut(       double sigmaCut)  {m_SigmaCut           = sigmaCut; }
    void setMinEnergy(      double minEnergy) {m_MinEnergy          = minEnergy;}
    void setMaxConsGaps(    int    maxGaps)   {m_MaxConsecutiveGaps = maxGaps;  }
    void setMinSegmentHits( int    minHits)   {m_MinSegmentHits     = minHits;  }
    void setMaxChisqCut(    double maxChi)    {m_MaxChiSqCut        = maxChi;   }
    void setIniErrSlope(    double errSlp)    {m_IniErrorSlope      = errSlp;   }
    void setIniErrPos(      double errPos)    {m_IniErrorPosition   = errPos;   }
    void setPlaneEnergies(  bool   enePlane)  {m_PlaneEnergies      = enePlane; }

private:
    /// private Constructor
    TkrControl();
 

    /// Pointer to the singleton object
    static TkrControl* m_this;

    /// Data members
    int    m_MaxCandidates;          // Max number of Candidates 
    int    m_MinTermHitCount;        // Number of hits to terminate Combo PR

    double m_FEneParticle;           // Fraction of Cal energy to use in PR.

    double m_SigmaCut;               // PR search window (in  sigmas)
    double m_MinEnergy;              // Min tracking energy (MeV)

    int    m_MaxConsecutiveGaps;     // Max consecutive Gaps - Stop
    int    m_MinSegmentHits;         // Min number of hits for segment
    double m_MaxChiSqCut;            // Max allow PR Chisq. 
    double m_IniErrorSlope;          // First Hit error in Kalman: 10 deg 
    double m_IniErrorPosition;       // First Hit error in Kalman: .1 mm

    bool   m_PlaneEnergies;          // Decrease particle energies by exp(-rad_len)

};

#endif
