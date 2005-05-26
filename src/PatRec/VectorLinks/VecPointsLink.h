/**
 * @class VecPointLinks
 *
 * @brief Implementation of a particular "matrix" class for the generic Kalman Filter fitter. 
 *        This version based on CLHEP HepMatrix. 
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/TrackFit/KalmanFilterFit/VecPointLinks.h,v 1.4 2004/09/23 21:30:29 usher Exp $
 */

#ifndef VecPointsLink_h
#define VecPointsLink_h

#include "VecPoint.h"
//#include "GaudiKernel/DataObject.h"

class VecPointsLink //: virtual public DataObject
{
public:
    
    enum StatusBits {ASSOCIATED = 0x0001,
                     FIRSTLINK  = 0x0002,
                     SAMETOWER  = 0x0100};

    // Constructors
    VecPointsLink(const VecPoint* firstPoint, const VecPoint* secondPoint, double ang);

    ~VecPointsLink() {}

    void setUnAssociated()              {m_statusBits      &= ~ASSOCIATED;}
    void setAssociated()                {m_statusBits      |=  ASSOCIATED;}
    void setNotFirstLink()              {m_statusBits      &= ~FIRSTLINK;}
    void setFirstLink()                 {m_statusBits      |=  FIRSTLINK;}
    void setSameTower()                 {m_statusBits      |=  SAMETOWER;}
    void setMaxScatAngle(double ang)    {m_maxScatAngle     =  ang;}
    void setAngleToNextLink(double ang) {m_angleToNextLink  =  ang;}

    const bool      associated()         const {return (m_statusBits & ASSOCIATED) == ASSOCIATED;}
    const bool      firstLink()          const {return (m_statusBits & FIRSTLINK)  == FIRSTLINK;}
    const bool      sameTower()          const {return (m_statusBits & SAMETOWER)  == SAMETOWER;}
    const Point&    getPosition()        const {return m_position;}
    const Vector&   getVector()          const {return m_vector;}
    const VecPoint* getFirstVecPoint()   const {return m_firstVecPoint;}
    const VecPoint* getSecondVecPoint()  const {return m_secondVecPoint;}
    const double    getAngleToNextLink() const {return m_angleToNextLink;}
    const double    getMaxScatAngle()    const {return m_maxScatAngle;}

    // This used to match links
    const bool      matchFirst( const VecPointsLink& linkToMatch) {return linkToMatch.getSecondVecPoint() == m_firstVecPoint;}
    const bool      matchSecond(const VecPointsLink& linkToMatch) {return linkToMatch.getFirstVecPoint()  == m_secondVecPoint;}

    double          angleToNextLink(const VecPointsLink& linkToMatch);

    
    const bool operator<(const VecPointsLink* right) const {return m_position.z() < right->getPosition().z();}

private:
    const VecPoint* m_firstVecPoint;
    const VecPoint* m_secondVecPoint;
    
    unsigned int    m_statusBits;
    Point           m_position;
    Vector          m_vector;

    // Calculated expected maximum scattering angle
    double          m_maxScatAngle;

    // warning! a very volatile variable for specific use in the link attachment stage!
    double          m_angleToNextLink;
};

VecPointsLink::VecPointsLink(const VecPoint* firstPoint, const VecPoint* secondPoint, double ang) : 
                             m_firstVecPoint(firstPoint), 
                             m_secondVecPoint(secondPoint), 
                             m_statusBits(0), 
                             m_maxScatAngle(ang),
                             m_angleToNextLink(0.)
{
    Ray rayBetweenPoints = firstPoint->getRayTo(secondPoint);
    
    if (firstPoint->getXCluster()->tower() == secondPoint->getXCluster()->tower())
        m_statusBits = SAMETOWER;

    m_position = rayBetweenPoints.position();
    m_vector   = rayBetweenPoints.direction();
}

double VecPointsLink::angleToNextLink(const VecPointsLink& linkToMatch)
{
    double cosAng = m_vector.dot(linkToMatch.getVector());

    m_angleToNextLink = acos(cosAng);

    return m_angleToNextLink;
}

typedef std::vector<VecPointsLink>  VecPointsLinkVec;
typedef std::vector<VecPointsLink*> VecPointsLinkPtrVec;

#endif

