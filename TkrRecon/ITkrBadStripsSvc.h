#ifndef __ITKRBADSTRIPSSVC_H
#define __ITKRBADSTRIPSSVC_H 1

#include "GaudiKernel/IInterface.h"
#include "idents/GlastAxis.h"

#include <string>
#include <vector>
#include <fstream>

typedef std::vector<int> v_strips;
typedef v_strips::const_iterator v_strips_it;

//----------------------------------------------
//
//   TkrBadStripsSvc
//
//	 Tracker BadStrips Service. Supplies the bad strips 
//   for use by SiClustersAlg, for example
//----------------------------------------------
//             Leon Rochester, 3-June-2001
//----------------------------------------------

static const InterfaceID IID_ITkrBadStripsSvc(907, 1 , 0); 

class ITkrBadStripsSvc : public virtual IInterface
{
public:

    //! Constructor of this form must be provided

    static const InterfaceID& interfaceID() { return IID_ITkrBadStripsSvc; }
   
    //virtual StatusCode initialize() = 0;
    //virtual StatusCode finalize() = 0;
    

    virtual int getIndex(const int tower, const int layer, const idents::GlastAxis::axis axis) = 0;
    virtual v_strips* getBadStrips(const int tower, const int layer, 
        const idents::GlastAxis::axis axis) = 0;
    virtual v_strips* getBadStrips(const int index)= 0;
    virtual bool isBadStrip(const int tower, const int layer, 
        const idents::GlastAxis::axis axis, const int strip) = 0;
    virtual bool isBadStrip(const v_strips* v, const int strip) = 0;
    virtual bool isTaggedBad(const int taggedStrip) = 0;
    virtual int tagBad(const int strip) = 0;
    virtual int tagGood(const int strip) = 0;
    virtual int untag(const int strip) = 0;

 };

#endif