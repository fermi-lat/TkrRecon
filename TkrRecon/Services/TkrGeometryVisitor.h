#include "GlastSvc/GlastDetSvc/IGeometry.h"
#include "idents/VolumeIdentifier.h"
#include "idents/TowerId.h"


class TkrGeometryVisitor: public IGeometry
{
public:

	TkrGeometryVisitor();

	~TkrGeometryVisitor() {}

/**  
* @class TkrGeometryVisitor
*
* @brief For now not used, but the code works, and utilizes the all the features of a visitor.
*
* @author Leon Rochester 
*   
* @param s type of the shape
* @param id vector of unsigned ints (maybe null)
* @param name
* @param material
* @param params vector with the six transformation parameters, 
* followed by 3 or so dimensions depending on the shape.
*/
   
	virtual IGeometry::VisitorRet pushShape(ShapeType s, const UintVector& id, std::string name, 
			 std::string material, const DoubleVector& params, VolumeType type);
  
    /// called to signal end of nesting 
    virtual void popShape();

    /// Need a setMode in case someone wants other than default 
    /// choice mode to implement getMode for IGeometry interface
    virtual void setMode(std::string pmode) {m_mode = pmode;}
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

	/// choice mode for traversing geometry
    std::string m_mode;


};