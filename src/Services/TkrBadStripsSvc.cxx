
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/Services/TkrBadStripsSvc.h"
#include <fstream>
#include <algorithm>
#include <iostream>

#include "xml/IFile.h"

static const SvcFactory<TkrBadStripsSvc> s_factory;
const ISvcFactory& TkrBadStripsSvcFactory = s_factory;


// Service parameters which can be set at run time must be declared.
// This should be done in the constructor.

TkrBadStripsSvc::TkrBadStripsSvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the file to get data from   
    declareProperty("badStripsFile", m_badStripsFile);
    
    return;	
}

StatusCode TkrBadStripsSvc::initialize()
{
    // Purpose: reads in bad strips and constructs in-memory bad strip vectors
    // Inputs:  None
    // Outputs: Status code (Success/Failure)
    
    MsgStream log(msgSvc(), name());
    log.setLevel(MSG::DEBUG);
    StatusCode sc = StatusCode::SUCCESS;
    
    Service::initialize();
    
    m_badStripsFile = "";
    
    setProperties();
    
    // commented code in this routine was the original attempt to implement
    //  the bad strips as a vector of vectors.
    
    // int size = 0;
    // makeCol(size); //make sure that Collection is sensibly initialized
    
    // If there is no bad strips file, service will do nothing
    if (m_badStripsFile=="") {        
        log << MSG::INFO << "No bad strips file was requested." << endreq;
        log << MSG::INFO << "  No strip filtering will be done." << endreq;
        return sc;
    }
    
    // this method resolves environmental variables in the file name
    xml::IFile::extractEnvVar(&m_badStripsFile);    
    log << MSG::INFO << "Input file for bad strips: " << m_badStripsFile << endreq;
    
    // open bad strips file
    std::ifstream file;
    file.open( m_badStripsFile.c_str());
    
    if (!file) {
        log << MSG::ERROR << "  File not found: check jobOptions." << endreq;
        return StatusCode::FAILURE;
    }
    
    //int size = m_ntowers*m_nlayers*m_nviews;    
    //makeCol(size);
    
    readFromFile(&file);
    
    file.close();
    
    // log << MSG::DEBUG<< "m_stripsCol has " << 
    //      m_stripsCol.size() << " elements" << endreq;
    
    return sc;
}

StatusCode TkrBadStripsSvc::finalize()
{
    return StatusCode::SUCCESS;
}

/*void TkrBadStripsSvc::makeCol(const int size)
{
//    m_stripsCol.assign(size);
//  return;
}
*/

void TkrBadStripsSvc::readFromFile(std::ifstream* file)
{    
    // Purpose: read bad strips from file and make in-memory vectors
    // Inputs:  File name
    // Outputs: None
    // Dependencies: None
    // Caveats: None
    
    bool read = true;           // for testing
    bool makestrips = true;     // for testing
    
    int nStrips = 0;
    std::string junk;
    
    // format of file:
    //
    // -1 at beginning of line is a comment
    //
    // -2 at beginning of line sets the tag (default is 1)
    // 
    // for each layer with bad strips:
    // tower# sequence#  [strip#] [strip#] ...  -1
    // all whitespace is ignored, so data for a layer may span lines
    //
    // 
    
    while(read && !file->eof()) {
        int tower;
        int sequence;
        
        int tag = 1;
        
        *file >> tower ;
        if(tower==-1) { // comment line, just skip 
            std::getline(*file, junk);
            continue;
        }
        if (tower==-2) { // set the tag
            *file >> tag;
            std::getline(*file, junk);
            continue;
        }
        
        if (file->eof()) break;
        *file >> sequence;
        
        // kludge until the geometry supplies this info
        // converts layer (0...35) to bilayer and view
        int layer = sequence/2;
        int element = (sequence+3)%4;
        int view = element/2;
        
        v_strips* v;
        if (makestrips) v = getBadStrips(tower, layer, static_cast<idents::GlastAxis::axis>(view) );
        int strip = -1;
        *file >>  strip;
        while (strip>=0) {
            if (makestrips) {
                addStrip(v, tagBad(strip, tag));
                *file >> strip;
                nStrips++;
            }
            
            // sort strips in ascending order after each line is read in
            if (makestrips) std::sort(v->begin(), v->end());
            
        }  
        return;
    }
}

int TkrBadStripsSvc::getIndex(const int tower, const int layer, const idents::GlastAxis::axis axis) 
{
    // Purpose:  calculate index into array of vectors
    // Inputs:   tower, bilayer, axis
    // Outputs:  index
    
    int view; 
    // this is to decouple the store from the current definition of axes
    // not that it will ever change
    if (axis==idents::GlastAxis::X)  {view = 0;}
    else if (axis==idents::GlastAxis::Y) {view = 1;}
    else {return 0;}
    // for now, hardwired to be as large as will ever by needed
    return view + NVIEWS*(layer + NLAYERS*tower);
}

void TkrBadStripsSvc::addStrip(v_strips* v, const int taggedStrip) 
{
    // Purpose: add a bad strip to the list, already tagged bad
    // Inputs:  strip number
    // Outputs: None
    
    v->push_back(taggedStrip);
    return;
}

v_strips* TkrBadStripsSvc::getBadStrips(const int tower, const int layer, const idents::GlastAxis::axis axis)
{
    // Purpose:  return pointer to a bad strip vector
    // Inputs:   tower, layer, axis
    // Outputs:  pointer to that vector
    
    int index = getIndex(tower, layer, axis);
    
    return getBadStrips(index);
}


v_strips* TkrBadStripsSvc::getBadStrips(const int index)
{
    // Purpose:  return pointer to a bad strip vector
    // Inputs:   index
    // Outputs:  pointer to that vector
    
    // int ind = (m_stripsCol.size()==0) ? m_stripsCol.size() : index;
    int ind = (index<0 || index>575) ? 0 : index;
    return &m_stripsCol[ind];
}


int TkrBadStripsSvc::tagBad(const int strip, const int tag) 
{
    // Purpose: tag a strip bad
    
    return (strip | tag << tagShift);
}


int TkrBadStripsSvc::tagGood(const int strip) 
{
    // Purpose: tag a strip good
    
    return strip;
}


int TkrBadStripsSvc::untag(const int strip) 
{
    // Purpose:  untag a tagged strip
    return (strip & stripMask);
}

int TkrBadStripsSvc::tagField(const int strip) 
{
    // Purpose: recover the tag field from the tagged strip
    // Input:   tagged strip
    // Output:  tag field
    return (strip>>tagShift);
}


bool TkrBadStripsSvc::isTaggedBad(const int taggedStrip) 
{
    // Purpose: return status of tagged strip
    // Input:   tagged strip
    // Output:  true if strip is tagged bad
    
    return (tagField(taggedStrip)>0);
}


bool TkrBadStripsSvc::isBadStrip(const int tower, const int layer, 
                                 const idents::GlastAxis::axis axis, const int strip) 
{
    // Purpose: determine if a given strip is bad
    // Inputs:  tower, bilayer, axis, strip#
    // Output:  true if strip is in the list ( that is, is bad)
    
    v_strips* v = getBadStrips(tower, layer, axis);
    return isBadStrip(v, strip);
}

bool TkrBadStripsSvc::isBadStrip(const v_strips* v, const int strip)
{
    // Purpose: determine if a given strip is bad
    // Inputs:  index, strip#
    // Output:  true if strip is in the list
    
    bool isBad = false;
    v_strips_it it;
    for (it=v->begin(); it!=v->end(); it++) {
        if ( untag(*it)==strip ) {
            isBad = true;
            break;
        }
    }
    return isBad;
    
    // this was the orginal code, before the tag was variable
    // might be useful again some day
    //v_strips_it it = std::find(v->begin(), v->end(), tagBad(strip));
    //return (it!=v->end());
}

int TkrBadStripsSvc::swapForSort( const int strip)
{
    // Purpose: swaps high and low order stuff, so tagged strips can be sorted by strip number
    // Inputs:  strip, tagged or untagged (normal or swapped)
    // Output:  same strip (swapped or normal)
    
    return ((strip&stripMask)<< tagShift) | (strip>>tagShift); 
}

// queryInterface

StatusCode  TkrBadStripsSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ITkrBadStripsSvc == riid) {
        *ppvIF = dynamic_cast<ITkrBadStripsSvc*> (this);
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
    return StatusCode::SUCCESS;
}

// access the type of this service

const IID&  TkrBadStripsSvc::type () const {
    return IID_ITkrBadStripsSvc;
}
