
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "TkrRecon/TkrBadStripsSvc.h"
#include <fstream>
#include <algorithm>
#include <iostream>

#include "xml/IFile.h"

static const SvcFactory<TkrBadStripsSvc> s_factory;
const ISvcFactory& TkrBadStripsSvcFactory = s_factory;

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

TkrBadStripsSvc::TkrBadStripsSvc(const std::string& name, ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{
    //Name of the file to get data from   
    declareProperty("badStripsFile", m_badStripsFile);

    return;	
}


StatusCode TkrBadStripsSvc::initialize()
{
    bool debug = false;

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
       
    Service::initialize();

	m_badStripsFile = "";

    setProperties();

    int size = 0;
    makeCol(size); //make sure that Collection is sensibly initialized

    if (m_badStripsFile=="") {        
        log << MSG::INFO << "No bad strips file was requested." << endreq;
        log << MSG::INFO << "  No strip filtering will be done." << endreq;
        return sc;
    }

	//test 1
	log << MSG::DEBUG << "Test 1"<< endreq;

    xml::IFile::extractEnvVar(&m_badStripsFile);    
    log << MSG::INFO << "Input file for bad strips: " << m_badStripsFile << endreq;

	//test 2
	log << MSG::DEBUG << "Test 2"<< endreq;

    std::ifstream file;
	file.open( m_badStripsFile.c_str());
	
    if (!file) {
        log << MSG::ERROR << "  File not found: check jobOptions." << endreq;
        return StatusCode::FAILURE;
    }
	

    sc = service("TkrGeometrySvc", pTkrGeom, true);
    
    m_ntowers = pTkrGeom->numYTowers()*pTkrGeom->numYTowers();

    m_nlayers = pTkrGeom->numLayers();
    m_nviews  = pTkrGeom->numViews();

    //This is to test that TkrGeometrySvc is initialized... is there a better way?
	
    if ((m_ntowers<1) || (m_ntowers>25) || (m_nlayers<1) 
        || (m_nlayers>20) || (m_nviews != 2)) { 
        log << MSG::ERROR << "TkrGeometrySvc must precede BadStripsSvc"
            << " in the jobOptions file." << endreq;
        return StatusCode::FAILURE;
    }//
	
    size = m_ntowers*m_nlayers*m_nviews;    
    makeCol(size);
    
    readFromFile(&file);
	
	file.close();
	
	log << MSG::DEBUG<< "m_stripsCol has " << 
		m_stripsCol.size() << " elements" << endreq;
	for (int i = 0; i < m_stripsCol.size(); i++) {
		log << "Element " << i << ": " ;
		v_strips v = m_stripsCol[i];
		log << MSG::DEBUG << v.size() << " strips" << endreq;
		for (int j = 0; j < v.size(); j++) {
			log << v[j] << " " ;
		}
		log << MSG::DEBUG << endreq;
	}
	
           
    return sc;
}


StatusCode TkrBadStripsSvc::finalize()
{
    return StatusCode::SUCCESS;
}

void TkrBadStripsSvc::makeCol(const int size)
{
    m_stripsCol.assign(size);
    return;
}


void TkrBadStripsSvc::readFromFile(std::ifstream* file)
{    
    bool debug = false;
	bool read = true;
	bool makestrips = true;

    MsgStream log(msgSvc(), name());
    
    int nStrips = 0;
	std::string junk;

    while(read && !file->eof()) {
        int tower;
		int sequence;

        *file >> tower ;
		if(tower==-1) { // comment line, just skip 
			std::getline(*file, junk);
			continue;
		}
        if (file->eof()) break;
        *file >> sequence;
        // kludge until the geometry supplies this info
		int layer = sequence/2;
		int element = (sequence+3)%4;
		int view = element/2;

        v_strips* v;
		if (makestrips) v = getBadStrips(tower, layer, view);
        int strip = -1;
        *file >>  strip;
        while (strip>=0) {
            if (makestrips) addStrip(v, strip);
            *file >> strip;
            nStrips++;
        }

        if (makestrips) std::sort(v->begin(), v->end());
        
    }
    log << MSG::INFO << nStrips << " bad strips read from file" << endreq;
   
    return;
}


int TkrBadStripsSvc::getIndex(const int tower, const int layer, const int view) 
{
    return view + m_nviews*(layer + m_nlayers*tower);
}


void TkrBadStripsSvc::addStrip(v_strips* v, const int strip) {
    int tagged_strip = tagBad(strip);
    v->push_back(tagged_strip);
    int size = v->size();
    int capa = v->capacity();
    return;
}


v_strips* TkrBadStripsSvc::getBadStrips(const int tower, const int layer, const int view)
{
    int index = getIndex(tower, layer, view);

    return getBadStrips(index);
}


v_strips* TkrBadStripsSvc::getBadStrips(const int index)
{
    int ind = (m_stripsCol.size()==0) ? m_stripsCol.size() : index;
    return &m_stripsCol[ind];
}


int TkrBadStripsSvc::tagBad(const int strip) 
{
    return (strip << 1) | 1;
}


int TkrBadStripsSvc::tagGood(const int strip) 
{
    return (strip << 1);
}


int TkrBadStripsSvc::untag(const int strip) 
{
    return (strip >> 1);
}


bool TkrBadStripsSvc::isTaggedBad(const int taggedStrip) {
    return ((taggedStrip & 1) != 0);
}


bool TkrBadStripsSvc::isBadStrip(const int tower, const int layer, 
                                 const int view, const int strip) 
{
    v_strips* v = getBadStrips(tower, layer, view);
    return isBadStrip(v, strip);
}

bool TkrBadStripsSvc::isBadStrip(const v_strips* v, const int strip)
{
    return std::find(v->begin(), v->end(), tagBad(strip));
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
    //addRef(); 
    return StatusCode::SUCCESS;
}

// access the type of this service

const IID&  TkrBadStripsSvc::type () const {
    return IID_ITkrBadStripsSvc;
}
