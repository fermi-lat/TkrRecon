
#ifndef TKRGEOMETRYVISITOR_H
#define TKRGEOMETRYVISITOR_H 

/**  
* @class TkrGeometryVisitor
*
* @brief A sample Geometry visitor for the TKR
*
* Just a toy example
*
* @author Leon Rochester 
*
* $Header$
*/

#include "GlastSvc/GlastDetSvc/IGeometry.h"
#include "idents/VolumeIdentifier.h"
#include "idents/TowerId.h"
  
class TkrGeometryVisitor: public IGeometry
{
public:

	TkrGeometryVisitor();

	~TkrGeometryVisitor() {}
   
	/// Standard interface to the detModel
    virtual IGeometry::VisitorRet pushShape(ShapeType s, const UintVector& id, std::string name, 
			 std::string material, const DoubleVector& params, VolumeType type);
  
    /// called to signal end of nesting 
    virtual void popShape();

    /// Need a setMode in case someone wants other than default 
    virtual void setMode(std::string pmode) {m_mode = pmode;}
    /// Implements getMode for IGeometry interface
    virtual std::string getMode() {return m_mode;}

private:

    /// Store the relevant information from each call

	idents::TowerId m_tower;
	int m_towerX;
	int m_towerY;
	int m_tray;
	int m_botTop;
	int m_view;
	int m_layer;

	double m_param[9];

	/// mode for traversing geometry
    std::string m_mode;


};

#endif  // TKRGEOMETRYVISITOR_H