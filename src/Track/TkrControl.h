
#ifndef TKRCONTROL_H
#define TKRCONTROL_H

//############################################
class TkrControl
//############################################
{
public:

    // Define as a singleton object
    static TkrControl* getPtr();

    // Retrieve values

    const inline int    getMaxCandidates()      {return m_MaxCandidates;     }
    const inline int    getMinTermHitCount()    {return m_MinTermHitCount;   }
    const inline double getFEneParticle()       {return m_FEneParticle;      }
    const inline double getSigmaCut()           {return m_SigmaCut;          }
    const inline double getMinEnergy()          {return m_MinEnergy;         }
    const inline int    getMaxConsecutiveGaps() {return m_MaxConsecutiveGaps;}
    const inline int    getMinSegmentHits()     {return m_MinSegmentHits;    }
    const inline double getMaxChisqCut()        {return m_MaxChiSqCut;       }
    const inline double getIniErrSlope()        {return m_IniErrorSlope;     }
    const inline double getIniErrPosition()     {return m_IniErrorPosition;  }
    const inline bool   getPlaneEnergies()      {return m_PlaneEnergies;     }


    // Allow for control variables to be set at initialization
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
