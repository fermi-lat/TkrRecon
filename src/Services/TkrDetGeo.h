#ifndef __tkrDetGeo_H
#define __tkrDetGeo_H 1

#include <vector>
#include <string>
#include "TkrRecon/Services/TkrAxis.h"
#include "geometry/Point.h"
#include "gui/DisplayRep.h"

//##############################################
class TkrGeo
//##############################################
{
public:
	// construct
	TkrGeo(Point p = Point(0.,0.,0.), Point s=Point(0.,0.,0))
		: m_position(p),m_size(s) {}
	void setGeo(Point p, Point s);
	void setPosition(Point p)   {m_position = p;};
	void setSize(Point s)       {m_size = s;};

	// access
	Point position() const  {return m_position;}
	Point size() const      {return m_size;}

private:

	Point m_position;
	Point m_size;
};

//##############################################
class TkrMatGeo : public TkrGeo
//##############################################
{
public:

	// construct
	TkrMatGeo(){}
	TkrMatGeo(Point p, Point s):m_material("vacuum"),m_X0(0), TkrGeo(p,s){}
	TkrMatGeo(std::string name, double d, Point p, Point s):m_material("vacuum"),
		m_X0(0), TkrGeo(p,s){}
	void setMaterial(std::string name, double d);

	// access
	std::string material() const  {return m_material;}
	double X0() const        {return m_X0;}
	double radLen();

private:

	std::string   m_material;
	double        m_X0;
};
//##############################################
class tkrDetGeo : public TkrMatGeo, public TkrAxis
//##############################################
{
public:

	friend class trackerGeo;

public:

	// construct
	tkrDetGeo(int ilayer, axis a, int id, Point p, Point s):m_layer(ilayer),
		m_axis(a), m_id(id), TkrMatGeo(p,s) {}
	tkrDetGeo(int ilayer, axis a, int id): m_layer(ilayer),
		m_axis(a), m_id(id), TkrMatGeo() {}
	~tkrDetGeo() {};
	void setName(std::string n)  {m_name = n;}

	// access
//	type getType() const        {return m_type;}
	int layer()            const {return m_layer;}
	axis getAxis()         const {return m_axis;}
	int id()               const {return m_id;}
	std::string name()     const {return m_name;}

	// operations
	static tkrDetGeo::axis makeAxis(int i) {return (i == 0? tkrDetGeo::X : tkrDetGeo::Y);}
	static int makeAxis(tkrDetGeo::axis a) {return (int) a;}
	// void writeOut() const;
	void draw(gui::DisplayRep& v) const;

private:

	std::string m_name;
	int m_layer;
	axis m_axis;
	int m_id;

};

#endif
