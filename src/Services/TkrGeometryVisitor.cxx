#include "TkrRecon/Services/TkrGeometryVisitor.h"
#include <iostream>
#include <algorithm>
using namespace std;

TkrGeometryVisitor::TkrGeometryVisitor() : m_mode("propagate") {}


IGeometry::VisitorRet
TkrGeometryVisitor::pushShape(ShapeType /*s*/, 
                              const UintVector& idvec, 
                              string name, 
                              string /*material*/, 
                              const DoubleVector& params, 
                              VolumeType /*type*/)
{
    // Purpose:  returns at each node of the geometry
    // Parameters:  s -- type of shape
    //              id -- vector of inspecified ints
    //              name 
    //              material
    //              params   -- vector of the 6 transformation parameters
    //                          followed by 3 or more dimensions, 
    //                          depending on th shape
    
    // This is a toy example, which illustrates the features of a visitor
    //
    // This visitor looks for silicon layers in tower number 9, 
    // and prints out the parameters associated with them.   
    //
    // For each volume, idvec[0] is that volume's piece of 
    // the volume identifier. For the tray, the view is stored in idvec[1].
    //  
    // For a box, there are 9 parameters. the first 3 are the rotation 
    // from the local frame to the mother frame, the next three the 
    // translation, and the last three, the (half) dimensions of the box.
    
    int i;
    
    // each call is further down the chain.  We will encounter the 
    // row of y-towers first
    
    if(name=="oneCAL") {
        cout << " encountered a CAL " << endl;
        return AbortSubtree;
    }
    
    
    if (name=="towerRow") {
        m_towerY = idvec[0];
    }
    
    //then the x tower
    
    if (name=="oneTower")   {
        m_towerX = idvec[0];
        //now we can make up a TowerID
        m_tower = idents::TowerId(m_towerX, m_towerY);
        if(m_tower.id()==9) {
            cout << "Tower " << m_tower.id() << endl;
        } else {
            return AbortSubtree;
        }   
    }
    
    //these are all the tray names.
    if (name=="traySuper"||name=="trayTop"||name=="trayBot"
        ||name=="trayReg"||name=="trayNoConv") {
        m_tray = idvec[0]; m_view = idvec[1];
        for (i=0;i<9;i++){
            m_param[i] = params[i];
        }
        
        if(m_tower.id()==9) {
            cout << "     Tray " << m_tray << " " << name << " view " 
                << m_view << endl;
            cout << "        params " ;
            for (i=0;i<9;i++) {
                cout << m_param[i] << " " ;
            }
            cout << endl;
        }
    }
    
    // and the silicon
    if (name.find("SiLayerBox")!=std::string::npos) {
        m_botTop = idvec[0];
        m_layer = m_tray - 1 + m_botTop;
        
        if (m_tower.id()==9) {
            cout << "           layer " << m_layer << " topBot " 
                << m_botTop << endl;
        }
    }
    
    return More;
}

void TkrGeometryVisitor::popShape()
{
    return;
}

