
#ifndef TKRCONTROL_H
#define TKRCONTROL_H

// Temporary (I think...)
#include "src/TrackFit/KalFitTrack/GFcontrol.h"

//############################################
class TkrControl
//############################################
{
public:

    // Define as a singleton object
    static TkrControl* GetPtr();

    // Retrieve values
    inline int    GetMaxCandidates()      const {return m_MaxCandidates;     }
    inline int    GetMinTermHitCount()    const {return m_MinTermHitCount;   }
    inline double GetFEneParticle()       const {return m_FEneParticle;      }
    inline double GetSigmaCut()           const {return m_SigmaCut;          }
    inline double GetMinEnergy()          const {return m_MinEnergy;         }
    inline int    GetMaxConsecutiveGaps() const {return m_MaxConsecutiveGaps;}
    inline int    GetMinSegmentHits()     const {return m_MinSegmentHits;    }
    inline double GetMaxChisqCut()        const {return m_MaxChiSqCut;       }
    inline double GetIniErrSlope()        const {return m_IniErrorSlope;     }
    inline double GetIniErrPosition()     const {return m_IniErrorPosition;  }
    inline bool   GetPlaneEnergies()      const {return m_PlaneEnergies;     }

    // Allow for control variables to be set at initialization
    void SetMaxCandidates(  int    maxCand)   {m_MaxCandidates      = maxCand;   GFcontrol::maxCandidates     = maxCand;   }
    void SetMinTermHitCount(int    minCount)  {m_MinTermHitCount    = minCount;  GFcontrol::minTermHitCount    = minCount; }
    void SetFEneParticle(   double enePart)   {m_FEneParticle       = enePart;   GFcontrol::FEneParticle       = enePart;  }
    void SetSigmaCut(       double sigmaCut)  {m_SigmaCut           = sigmaCut;  GFcontrol::sigmaCut           = sigmaCut; }
    void SetMinEnergy(      double minEnergy) {m_MinEnergy          = minEnergy; GFcontrol::minEnergy          = minEnergy;}
    void SetMaxConsGaps(    int    maxGaps)   {m_MaxConsecutiveGaps = maxGaps;   GFcontrol::maxConsecutiveGaps = maxGaps;  }
    void SetMinSegmentHits( int    minHits)   {m_MinSegmentHits     = minHits;   GFcontrol::minSegmentHits     = minHits;  }
    void SetMaxChisqCut(    double maxChi)    {m_MaxChiSqCut        = maxChi;    GFcontrol::maxChisqCut        = maxChi;   }
    void SetIniErrSlope(    double errSlp)    {m_IniErrorSlope      = errSlp;    GFcontrol::iniErrorSlope      = errSlp;   }
    void SetIniErrPos(      double errPos)    {m_IniErrorPosition   = errPos;    GFcontrol::iniErrorPosition   = errPos;   }
    void SetPlaneEnergies(  bool   enePlane)  {m_PlaneEnergies      = enePlane;  GFcontrol::planeEnergies      = enePlane; }

private:
    // Constructor and destructor are private
    TkrControl();
   ~TkrControl() {}

    // Pointer to the singleton object
    static TkrControl* m_this;

    // Data members
    int    m_MaxCandidates;
    int    m_MinTermHitCount;

    double m_FEneParticle;

    double m_SigmaCut;
    double m_MinEnergy;

    int    m_MaxConsecutiveGaps;
    int    m_MinSegmentHits;
    double m_MaxChiSqCut;
    double m_IniErrorSlope;
    double m_IniErrorPosition;

    bool   m_PlaneEnergies;
};

#endif
