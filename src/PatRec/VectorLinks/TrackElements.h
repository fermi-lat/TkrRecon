/**
 * @class TrackElements
 *
 * @brief Implementation of a particular "matrix" class for the generic Kalman Filter fitter. 
 *        This version based on CLHEP HepMatrix. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/src/PatRec/VectorLinks/TrackElements.h,v 1.1 2005/05/26 20:33:07 usher Exp $
 */

#ifndef TrackElements_h
#define TrackElements_h

#include "VecPointsLink.h"

class TrackElements 
{
public:
    enum StatusBits {UNUSED   = 0x0001,
                     USED     = 0x0002,
                     BESTTRK  = 0x0004,
                     NOTBEST  = 0x0008};

    // Constructors
    TrackElements(int numLinks, double rmsAngle, VecPointsLink* firstLink) : 
                  m_statusBits(UNUSED), 
                  m_numLinks(numLinks), 
                  m_rmsAngle(rmsAngle),
                  m_firstLink(firstLink)
                  {}

   ~TrackElements() {}

    unsigned int   getStatusBits() const {return m_statusBits;}
    int            getNumLinks()   const {return m_numLinks;}
    double         getRmsAngle()   const {return m_rmsAngle;}
    VecPointsLink* getFirstLink()  const {return m_firstLink;}

    bool         isUnUsed()      const {return (m_statusBits & UNUSED) == UNUSED;}

    const bool operator<(const TrackElements& right)  const;
    const bool operator==(const TrackElements& right) const;

    void setStatusBits(unsigned int bits) {m_statusBits = bits;}
    void setUsed()                        {m_statusBits |= USED; m_statusBits &= ~UNUSED;}

private:
    unsigned int   m_statusBits;
    int            m_numLinks;
    double         m_rmsAngle;
    VecPointsLink* m_firstLink;
};

const bool TrackElements::operator<(const TrackElements& right) const
{
    if (m_numLinks > 16 && right.getNumLinks() > 16)     // 18 is max
    {
        // calculate weight factors
        double leftWeight  = m_rmsAngle / (10. * m_numLinks);
        double rightWeight = right.getRmsAngle() / (10. * right.getNumLinks());

        if (leftWeight < rightWeight) return true;
    }
    else
    {
        int lenDiff = abs(m_numLinks - right.getNumLinks());

        // Divide between close in length and not close
        if (lenDiff > 0)
        {
            // Longer track wins
            if (m_numLinks > right.getNumLinks()) return true;
        }
        else
        {
            // straighter wins
            //if (m_numLinks == right.getNumLinks()) 
            //{
                if (m_rmsAngle < right.getRmsAngle()) return true;
            //}
        }
    }

    return false;
}

const bool TrackElements::operator==(const TrackElements& right) const
{
    // same length
    if (m_numLinks != right.getNumLinks()) return false;
    
    // and same angle
    if (m_rmsAngle != right.getRmsAngle()) return false;

    return true;
}

typedef std::vector<TrackElements>  TrackElementsVec;
typedef std::vector<TrackElements*> TrackElementsPtrVec;

#include "StdRelTable/RelTable.h"
#include "VecPoint.h"

typedef TkrRecon::RelTable<TrackElements, VecPoint>               TrackElemToPointsTab;
typedef TkrRecon::Relation<TrackElements, VecPoint>               TrackElemToPointsRel;
typedef ObjectList< TkrRecon::Relation<TrackElements, VecPoint> > TrackElemToPointsTabList;

#endif

