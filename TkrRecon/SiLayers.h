
#ifndef __SILAYERS_H
#define __SILAYERS_H 1

// #include <stdlib.h>
#include <iostream>
#include <vector>
#include "Gaudi/Kernel/DataObject.h"

extern const CLID& CLID_SiLayers;

//----------------------------------------------
//
//   SiLayers
//
//   Transient Storage Clases
//----------------------------------------------
//   It contains the Raw tracker data
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
/*!
SiLayer class, element container of the RAW silicon data (organized by layer)
*/
//##############################################
class SiLayer
//##############################################
{
public:

	//! constructor with layer and view and ToT information
	SiLayer(int layer, int iview, int ToT = -1);
	//! destructor, clear the list of strips addresses
	virtual ~SiLayer() {clear();}
	
	//! set ToT value
	inline void setToT(int ToT) {m_ToT = ToT;}
	//! add strips address to the list of strips in the SiLayer
	void addStrip(int istrip);

	//! returns ToT value
	inline int ToT() const   {return m_ToT;}
	//! returns layer id
	inline int layer() const {return m_layer;} 
	//! returns view id (default X=0, but crosscheck!)
	inline int view() const  {return m_view;}

	//! returns the number of strips in the SiLayer
	inline int nstrips()      {return m_stripList.size();}
	//! returns the strip id (address) of the strip in position i on the list
	inline int idstrip(int i) {return m_stripList[i];}

	//! clear the list of strips
	void clear()              {m_stripList.clear();}
	//! writes out the information of the SiLayer
	void writeOut() const;

private:

	//! returns true if the strip id is valid (TB particular class due to different configurations)
	// bool validStrip(int istrip);

private:

	//! view 
	int m_view;
	//! layer
	int m_layer;
	//! ToT value
	int m_ToT;

	//! list of strips addresses
	std::vector<int> m_stripList;

};

//const static CLID CLID_SiLayers = 250;
/*!
SiLayers class, container class of the SiLayer elements and service of SiLayer
*/
//##############################################
class SiLayers : public DataObject
//##############################################
{
public:

	//! constructor -ini the list
	SiLayers() {ini();}
	//! destructor - clear the list (delete pointers to elements)
	virtual ~SiLayers() {clear();}

	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_SiLayers;}
	virtual const CLID& clID() const {return classID();}

	//! add a SiLayer pointer into the list
	void add(SiLayer* sd)   {m_SiLayersList.push_back(sd);}

	//! clear the data
	virtual void clear();
	//! empty method
	virtual void make() {} 
	
	//! returns the total number of strips 
	int numStrips();
	//! returns a pointer to the SiLayer with view and plane 
	SiLayer* getSiLayer(int ilayer, int iview);

	//! returns the total number of strips
	int num()             const {return m_SiLayersList.size();}
	//! returns a SiLayer located in position i of the list
	SiLayer* Layer(int i) const {return m_SiLayersList[i];}

	//! writes out the information of the SiLayers
	void writeOut() const;

protected:

	//! initialize the list
	virtual void ini();

private:

	//! vector of pointer to SiLayer objects
	std::vector<SiLayer*>  m_SiLayersList;

};

#endif
