
#include "TkrRecon/GFdata.h"

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

//##########################################
void GFdata::ini()
//##########################################
{
    m_vertex = Point(0.,0.,0.);
    m_direction = Vector(0.,0.,0.);
    m_RCenergy=0.;
    m_quality=-100.;
    m_firstLayer=-1;
    m_nhits = -1;
    m_itower =-1;
}
//##########################################
bool GFdata::empty() const
//##########################################
{
    bool empty = false;
    if (m_firstLayer < 0) empty = true;
    return empty;
}

//########################################################
void GFdata::writeOut(std::ostream& out) const
//########################################################
{
    out << " --- GFdata::writeOut --- " << "\n";
    out << " Quality       = " << Q() << "\n";
    out << " Vertex        = " << vertex().x() << " " <<	vertex().y() << " " << vertex().z() << "\n";
    out << " Direction     = " << direction().x() << " " << direction().y() << " " << direction().z() << "\n";
    out << " RCenergy      = " << RCenergy() << "\n";
    out << " first Layer   = " << firstLayer() << "\n";
    out << " Tower         = " << tower() << "\n";
    out << " num Hits      = " << nhits() << "\n";
}

//#########################################################################
Point GFdata::doVertex(const Ray& r1, const Ray& r2)
//#########################################################################
{
    Point Vertex(0.,0.,0.);
    
    double z1 = r1.position().z();
    double z2 = r2.position().z();
    double cosz1 = r1.direction().z();
    double cosz2 = r2.direction().z();
    if (cosz2 == 0. || cosz1 == 0.) return Vertex;
    
    // comon position z2ref = z1
    Vector Vdir = r2.direction();
    Point Porigin = r2.position((z1-z2)/cosz2);
    
    Ray r2ref(Porigin,Vdir);
    
    double deltaX = r1.position().x()-r2ref.position().x();
    double deltaY = r1.position().y()-r2ref.position().y();
    double deltaSlopeX = (r2ref.direction().x()/cosz2) - (r1.direction().x()/cosz1);
    double deltaSlopeY = (r2ref.direction().y()/cosz2) - (r1.direction().y()/cosz1);
    
    bool xview = true;
    bool yview = true;
    double deltaZX = 0.;
    double deltaZY = 0.;
    if (deltaSlopeX == 0.) xview = false;
    if (deltaSlopeY == 0.) yview = false;
    if (xview) deltaZX = deltaX/deltaSlopeX;
    if (yview) deltaZY = deltaY/deltaSlopeY;
    
    if (xview && !yview) Vertex = r1.position(deltaZX/cosz1);
    else if (!xview && yview) Vertex = r1.position(deltaZY/cosz1);
    else if (xview && yview) {
	Vertex	= (z1 >= z2? r1.position(0.) : r2.position(0.));
	// We assume r1 = X ray and r2 = Y ray;
	if (Vertex.x() == 0.) { 
	    double x0 = r1.position((z2-z1)/cosz1).x();
	    Vertex = Point( x0, r2.position().y(),r2.position().z());
	}
	if (Vertex.y() == 0.) {
	    double y0 = r2.position((z1-z2)/cosz2).y();
	    Vertex = Point(r1.position().x(), y0 , r1.position().z());
	}
    }
    return Vertex;
}

//#########################################################################
Vector GFdata::doDirection(const Vector& xdir, const Vector& ydir)
//#########################################################################
{
    // we assume xdir goes in the X direction and ydir goes in the Y direction
    double slopeX = xdir.x()/xdir.z();
    double slopeY = ydir.y()/ydir.z();
    
    double factor = -1;
    Vector dir = Vector(factor*slopeX,factor*slopeY,factor).unit();
    
    return dir;
}
