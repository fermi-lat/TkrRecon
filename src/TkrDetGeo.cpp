#include "src/tkrDetGeo.h"

//#########################################################
void TkrGeo::setGeo(Point p, Point s) 
//#########################################################
{
	setPosition(p);
	setSize(s);
}
//#########################################################
void TkrMatGeo::setMaterial(std::string name, double d) 
//#########################################################
{
	m_material = name;
	m_X0 = d;
}
//#########################################################
double TkrMatGeo::radLen() 
//#########################################################
{
	return 2.*size().z()/X0();
}
/*
//##########################################
void tkrDetGeo::writeOut() const
//##########################################
{
	if (messageManager::instance()->acceptLevel("DEBUG_MEDIUM")) return;
	std::ostream& out = messageManager::instance()->out();
	out << name() <<" " << material() <<" "<<  m_layer << " " << m_axis << " " << m_id << " " ;
	out << position().x() << " " << position().y() <<" "<< position().z() << " ";
	out << 2.*size().x() << " " << 2.*size().y() <<" "<< 2.*size().z();
	out << "\n";

}
*/
//##########################################
void tkrDetGeo::draw(gui::DisplayRep& v) const
//##########################################
{
	v.setColor("red");

	double sx = size().x();
	double sy = size().y();
	double sz = size().z();
	double x0 = position().x();
	double y0 = position().y();
	double z0 = position().z();
	int i =0;
	for (i = 0; i<2; i++) {
		double fz = ((i==0)?-1:1);
		v.moveTo(Point(x0-sx,y0-sy,z0+fz*sz));
		v.lineTo(Point(x0+sx,y0-sy,z0+fz*sz));
		v.lineTo(Point(x0+sx,y0+sy,z0+fz*sz));
		v.lineTo(Point(x0-sx,y0+sy,z0+fz*sz));
		v.lineTo(Point(x0-sx,y0-sy,z0+fz*sz));
	}
	for (i = 0; i<4; i++) {
		double fx = ( i < 2?-1:1);
		double fy = ( (i == 0 || i == 2)?-1:1);
		v.moveTo(Point(x0+fx*sx,y0+fy*sy,z0-sz));
		v.lineTo(Point(x0+fx*sx,y0+fy*sy,z0+sz));
	}
}
