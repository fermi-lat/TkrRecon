
#ifndef TKRBADSTRIPSSVC_H
#define TKRBADSTRIPSSVC_H 

/** 
 * @class BadStripsSvc
 *
 * @brief Maintains lists of bad strips, and provides access methods 
 *
 * First version 3-Jun-2001
 *
 * The bad strips are kept in ascii files, which are read in under
 * the control of the jobOptions file. In the ascii files, strips are
 * marked as hot or dead, but in memory, strips are only bad.
 *
 * The service creates an array of vectors.  The singly indexed array
 * corresponds to a doubly indexed array by tower and layer.
 * 
 * The original design was a vector of vectors, but this was abandoned
 * because some of the code failed on unix.
 * 
 * The use of the bad strips in the clustering algorithm, for example, depends
 * on mixing good and bad strips and still being able to sort in asending strip order.
 *
 * To this end strips are tagged good or bad by left-shifting the strip number by 1
 * and ORing with 0 (good) or 1 (bad). This mechanism is internal to the service.
 *
 * This can be extended to a more comprehensive
 * tagging arrangement if necessary. Another possiblity would be to create
 * a vector of pairs.  This might be more transparent in the long run.
 *
 * @author Leon Rochester
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/Services/TkrBadStripsSvc.h,v 1.3 2002/08/31 17:51:39 lsrea Exp $
 */

#include "GaudiKernel/Service.h"

#include "TkrRecon/ITkrBadStripsSvc.h"
#include "TkrRecon/ITkrGeometrySvc.h"

#include <string>
#include <vector>

class TkrBadStripsSvc : public Service,
        virtual public ITkrBadStripsSvc
{
public:

    /// Constructor of this form must be provided
    TkrBadStripsSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~TkrBadStripsSvc() {}
    
    //. reads in bad strips and creates in-memory vectors
	StatusCode initialize();
    StatusCode finalize();

	/// converts from (tower, layer, view) to index into array
    int getIndex(const int tower, const int layer, const idents::GlastAxis::axis);
	/// returns a pointer to a vector of bad strips for a given (tower, layer and view)
    v_strips* getBadStrips(const int tower, const int layer, const idents::GlastAxis::axis);
	/// returns a pointer to a vector of bad strips for a given array index
    v_strips* getBadStrips(const int index);
	/// returns true if the strip (tower, layer, view, strip) is bad
    bool isBadStrip(const int tower, const int layer, const idents::GlastAxis::axis, const int strip);
    /// returns true if the given strip is found in the vector pointed to by v_strips
	bool isBadStrip(const v_strips* v, const int strip);
    /// returns true if the tagged strip taggedStrop is tagged bad
	bool isTaggedBad(const int taggedStrip);
	/// tags a strip bad
    int tagBad(const int strip);
	/// tags a strip good
    int tagGood(const int strip);
	/// untags a strip
    int untag(const int strip);
       
    /// queryInterface - required for a service
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);
    /// required for a service
    static const InterfaceID& interfaceID() { return ITkrBadStripsSvc::interfaceID(); }
    /// returns the service type
    const IID& type() const;
   
private:

    //void makeCol(const int size); // left over from original attempt

	/// reads bad strips from file file
    void readFromFile(std::ifstream* file);
	/// adds a strip to a badstrip vector
    void addStrip(v_strips* v, const int strip);

    /// pointer to the geometry service
	ITkrGeometrySvc* pTkrGeom;

    /// File name for constants
	std::string m_badStripsFile;  

	/// number of towers
    int m_ntowers;
	/// number of bilayers
	int m_nlayers;
	/// number of views
    int m_nviews;


    /// array to hold bad strips vectors  [ <= 16*18*2 ]   
	v_strips m_stripsCol[576];
};


#endif // TKRBADSTRIPSSVC_H
