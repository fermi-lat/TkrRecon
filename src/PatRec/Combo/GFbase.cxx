
#include "src/PatRec/Combo/GFbase.h"
#include "TkrRecon/Track/GFtutor.h"

//##########################################
GFbase::GFbase(double sigmaCut, double ene, int ist, const Ray& testRay): 
m_sigmaCut(sigmaCut), 
m_iniEnergy(ene),
m_iniLayer(ist) 
//##########################################
{
    // control defauls
    m_alive = true;
    
    // input data
    if (m_iniEnergy < GFcontrol::minEnergy) m_iniEnergy = GFcontrol::minEnergy;
    m_inVertex = testRay.position();
    m_inDirection = testRay.direction();
    
}       
//########################################################
void GFbase::doit()                              
//########################################################
{
    int kplane = m_iniLayer;

    for ( ; -1 < kplane && kplane < GFtutor::numPlanes(); kplane++) {
        step(kplane);
                anastep(kplane);
        if (!alive()) {
            break;
        }
    }
}