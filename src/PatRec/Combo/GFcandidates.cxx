


#include "src/PatRec/Combo/GFcandidates.h"
#include "src/PatRec/Combo/GFpoints.h"
#include "TkrRecon/Track/GFtutor.h"


GFcandidates::GFcandidates(enum GFcandidates::type t, double ene,
                           double cut,
                           Point Pend, Point Pini):m_type(t),
                           m_eneCandidate(ene),m_Pend(Pend),m_Pini(Pini)
{
    ini();   
    bool ok= findSeedCandidates(m_candidates,m_seedtype, m_eneCandidate);  
}

//###########################################################
void GFcandidates::clear()
//###########################################################
{
    m_candidates.clear();
}

//###########################################################
GFdata GFcandidates::GFconstructor(enum GFcandidates::type type, double ene,
                                   int ilayer,const Ray testRay)
//###########################################################
{
    GFdata data;
    
    if (type == GFcandidates::PARTICLE) {
        GFparticle* _par = new GFparticle(GFcontrol::sigmaCut,	
            ene, ilayer, testRay);
        if (!_par->empty() && _par->accept()) data = _par->getGFdata();
        delete _par;
    } else if (type == GFcandidates::GAMMA) {
        GFgamma* _gamma = new GFgamma(GFcontrol::FEne, GFcontrol::sigmaCut,
            ene, ilayer, testRay);
        if (!_gamma->empty() && _gamma->accept()) {
            data = _gamma->getGFdata();
        }
        delete _gamma;
    } else if (type == GFcandidates::TRACK) {
        GFtrack* _track = new GFtrack(GFcontrol::sigmaCut,
            ene, ilayer, testRay);
        if (!_track->empty() && _track->accept()) data = _track->getGFdata();
        delete _track;
    } else if (type == GFcandidates::PAIR) {
        GFpair* _pair = new GFpair(GFcontrol::FEne, TkrCluster::X, GFcontrol::sigmaCut,
            ene, ilayer, testRay);
        if (!_pair->empty()) {
            if (_pair->accept()) {
                data = _pair->getBest()->getGFdata();
            }
        }
        delete _pair;
    }
    
    return data;
}
//-----------  Private drivers  ----------------------------- 

//###########################################################
void GFcandidates::ini()
//###########################################################
{
    clear();
    if (m_type == GAMMA) m_seedtype = PAIR;
    if (m_type == PARTICLE) m_seedtype = TRACK;
}

//###########################################################
bool  GFcandidates::findCandidates(std::vector<GFdata>& candidates,
                                   const GFdata& candidate, 
                                   double ene,
                                   enum GFcandidates::type typ)
//###########################################################
{
    bool ok = false;
    int naccepted = 0;
    
    int iniLayer  = candidate.firstLayer();
    int lastLayer = candidate.firstLayer();
    
    if (lastLayer - iniLayer >= GFcontrol::maxConsecutiveGaps) return ok;
    
    Ray testRay = candidate.ray();
    
    for (int ilayer = iniLayer ; ilayer <= lastLayer; ilayer++) {
        
        GFdata candidateGFdata = GFconstructor(typ, ene, ilayer, testRay);
        if (candidateGFdata.Q() > GFcontrol::minQ) {
            ok = true;
            naccepted++;
            incorporate(candidates, candidateGFdata);
        }
    }
    
    return ok;
}
//------------- Utilities -----------------------------------

//###########################################################
bool GFcandidates::findSeedCandidates(std::vector<GFdata>& candidates,GFcandidates::type typ, 
                                      double ene)
//###########################################################
{
    bool OK = false;
//    for (int iplane = 0 ; iplane < GFtutor::numPlanes() - 2; iplane++){
    // Temporary for dev. work....
    for (int iplane = 0 ; iplane < GFtutor::numPlanes() - 2; iplane++){
        bool ok = findSeedCandidates(candidates, typ, ene, iplane);
        OK = OK || ok;
    }
    return OK;
}

//###########################################################
bool GFcandidates::findSeedCandidates(std::vector<GFdata>& candidates, 
                                      GFcandidates::type typ,
                                      double ene, int ilayer, int itower)
//###########################################################
{
    //unused:	int nconstructed = 0;
    int naccepted = 0;
   
    // Create space point loops and check for hits
    GFpoints first_Hit(ilayer);
    if(first_Hit.finished()) return false;
    GFpoints secnd_Hit(ilayer+1);
    if(secnd_Hit.finished()) return false;

    while(!first_Hit.finished()) {
         Point x1(first_Hit.getSpacePoint());
         if(first_Hit.finished()) break;

         while(!secnd_Hit.finished()) {
             Point x2(secnd_Hit.getSpacePoint());
             if(secnd_Hit.finished()) continue;
       
        
            Vector VDir(x2.x()-x1.x(),x2.y()-x1.y(),x2.z()-x1.z());
            Ray testRay = Ray(x1, VDir.unit());
            if(fabs(testRay.direction().z()) < .19) continue; 

            //Flag which layer testRay starts in... 
            if(first_Hit.x_Layer()) testRay.setFlag(0);  
            else                    testRay.setFlag(1);

            GFdata candidateGFdata = GFconstructor(typ, ene, ilayer, testRay);
        
            if (candidateGFdata.Q() > GFcontrol::minQ) {
                naccepted++;
                incorporate(candidates, candidateGFdata);
            }
        }
    }
    return (naccepted > 0);
    
}
//###########################################################
void GFcandidates::incorporate(std::vector<GFdata>& pDatalist, const GFdata pData)
//###########################################################
{
    bool ienter = false;
    for (unsigned int i =0; i < pDatalist.size(); i++) {
        if (pData.Q() >= pDatalist[i].Q()) {
            pDatalist.insert(&pDatalist[i],pData);
            ienter = true;
        }
        if (ienter) break;
    }
    if (!ienter) pDatalist.push_back(pData);
    if (pDatalist.size()>GFcontrol::maxCandidates) pDatalist.pop_back();
    
}

//########################################################
Point GFcandidates::createPend(int ilayer, const Point& Pini)
//########################################################
{
    /*
    double weight = GFcontrol::minEnergy/(3.*m_eneCandidate); // weight;
    
    Point PCal = m_Pend;
    
    if (m_eneCandidate < GFcontrol::minEnergy) weight = 1.;
    if (PCal.mag() == 0) weight = 1.;
    double side = GFtutor::trayWidth();
    Point PTrk = GFtutor::_DATA -> meanHitInside(axis, ilayer+1,0.5*side, Pini);
    
    if (PTrk.mag() == 0.) {
        PTrk = GFtutor::_DATA -> meanHitInside(axis, ilayer+2,0.5*side, Pini);
    }
    
    if (PTrk.mag() == 0.) weight = 0.;
    double x = (1.-weight)*PCal.x()+weight*PTrk.x();
    double y = (1.-weight)*PCal.y()+weight*PTrk.y();
    double z = (1.-weight)*PCal.z()+weight*PTrk.z();
    
    Point PRef(x,y,z);
    
    if (PRef.mag() != 0.) {
        if (axis == TkrCluster::X ) y = Pini.y();
        else x = Pini.x();
    } else {
        x = Pini.x();
        y = Pini.y();
        z = Pini.z()-GFtutor::trayGap();
    }
    PRef = Point(x,y,z);
*/    
    double x = Pini.x();
    double y = Pini.y();
    double z = 10.; 
    Point PRef(x,y,z); 
    return PRef;
}
/*
//###########################################################
bool GFcandidates::findCandidates()
//###########################################################
{
    bool OK = false;
    
    m_candidates.clear();
    int ix = 0;
    for (; ix < m_Xcandidates.size(); ix++) {
        for (int iy = 0; iy < m_Ycandidates.size(); iy++) {
            GFdata Xcandidate = m_Xcandidates[ix];
            GFdata Ycandidate = m_Ycandidates[iy];
            bool ok = false;
            ok = findCandidates(m_candidates, 
                Xcandidate, Ycandidate, m_eneCandidate, m_type);
            OK = OK || ok;
        } // Y candidates
    } // X candidates
    
    if (OK)  return OK;
    
    // Force a PairFit with no veto
    bool save_veto = GFtutor::CUT_veto;
    std::vector<GFdata> candidates;
    candidates.clear();
    for (ix = 0 ; ix < m_Xcandidates.size(); ix++) {
        GFdata Xcandidate = m_Xcandidates[ix];
        candidates.clear();
        GFtutor::CUT_veto = false;
        
        findSeedCandidates(candidates, m_seedtype,m_eneCandidate,
            Xcandidate.firstLayer(),Xcandidate.tower());
        for (int iy =0 ; iy < candidates.size(); iy++) {
            GFtutor::CUT_veto = save_veto;
            GFdata Ycandidate = candidates[iy];
            bool ok = findCandidates(m_candidates, 
                Xcandidate,Ycandidate, m_eneCandidate, m_type);
            OK = OK || ok;
        } // Y candidates
    } // X candidates
    
    candidates.clear();
    for (int iy = 0 ; iy < m_Ycandidates.size(); iy++) {
        GFdata Ycandidate = m_Ycandidates[iy];
        candidates.clear();
        GFtutor::CUT_veto = false;
        
        findSeedCandidates(candidates, m_seedtype,m_eneCandidate,
            Ycandidate.firstLayer(),Ycandidate.tower());
        for (ix =0 ; ix < candidates.size(); ix++) {
            GFtutor::CUT_veto = save_veto;
            GFdata Xcandidate = candidates[ix];
            bool ok = findCandidates(m_candidates, 
                Xcandidate,Ycandidate, m_eneCandidate, m_type);
            OK = OK || ok;
        } // Y candidates
    } // X candidates
    GFtutor::CUT_veto = save_veto;
  
    return OK;
}
*/