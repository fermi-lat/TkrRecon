/**
 * @class McEventStructure
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/MonteCarlo/McEventStructure.h,v 1.1 2003/08/04 19:57:40 usher Exp $
 */
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/MonteCarlo/McParticle.h"

#ifndef McEventStructure_h
#define McEventStructure_h

namespace Event {

typedef SmartRefVector<Event::McParticle> McParticleRefVec; 
typedef SmartRef<Event::McParticle>       McParticleRef;

class McEventStructure : public DataObject
{
public:
    //! Define some bits to help classify the event
    enum ClassificationBits{  
        NOPRIMARY = 1 ,    //! No primary particle found (can't happen?)
        CHARGED   = 1<<1,  //! Primary particle is charged
        NEUTRAL   = 1<<2,  //! Primary particle is neutral
        GAMMA     = 1<<3,  //! Primary is a gamma
        CONVERT   = 1<<4,  //! Secondaries from gamma conversion
        BREM      = 1<<5,  //! Secondaries from Bremstrahlung
        COMPT     = 1<<6,  //! Secondaries from Compton scatter
        PHOT      = 1<<7,  //! Secondaries from Photoelectric effect (?)
        OTHER     = 1<<8,  //! Secondaries from some other process 
        RUNBIT    = 1<<12  //! Indicates code was called (does nothing)
    };

    /// Standard Gaudi Tool interface constructor
    McEventStructure(IDataProviderSvc* dataSvc, IParticlePropertySvc* partPropSvc);
   ~McEventStructure();

    /// Retrieve results
    const unsigned long                     getClassificationBits() const {return m_classification;}
    
    const Event::McParticleRef              getPrimaryParticle()    const {return m_primary;}

    const int                               getNumSecondaries()     const {return m_secondaries.size();}
    Event::McParticleRefVec::const_iterator beginSecondaries()      const {return m_secondaries.begin();}
    Event::McParticleRefVec::const_iterator endSecondaries()        const {return m_secondaries.end();}

    const int                               getNumAssociated()      const {return m_associated.size();}
    Event::McParticleRefVec::const_iterator beginAssociated()       const {return m_associated.begin();}
    Event::McParticleRefVec::const_iterator endAssociated()         const {return m_associated.end();}

    Event::McParticleRefVec                 getTrackVector();

private:
    bool isPrimaryDaughter(const Event::McParticle* mcPart);

    /// Bit-field for classification
    unsigned long             m_classification;

    Event::McParticleRef      m_primary;
    Event::McParticleRefVec   m_secondaries;
    Event::McParticleRefVec   m_associated;
};

};

#endif