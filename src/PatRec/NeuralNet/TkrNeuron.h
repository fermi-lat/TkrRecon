/**
* @class TkrNeuron
*
* @brief Container class for Neurons to be used in the Neural Net PR routine. 
* 
* This class contains information about possible subtracks (neurons)
* which are used to determine track candidates by the Nerual Net PR
* routine for the Kalman filter.
*
* last modified 3/02
*   
* @authors b. allgood and w. atwood
*
* @todo impliment a real equality operator.
* @todo decide whether to keep 'position' global.
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/NeuralNet/TkrNeuron.h,v 1.1 2002/04/01 19:22:37 allgood Exp $
*/

#ifndef __TKRNEURON_H
#define __TKRNEURON_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "src/PatRec/Utilities/TkrPoint.h"
#include <vector>
#include <assert.h>

// forward declarations
class TkrNeuralNet;

// global definitions
/// top is defined as being the end of the neuron with the lower layer number.
typedef enum {top,bottom} position;

class TkrNeuron
{
public:

	// constructors
	TkrNeuron(TkrPoint* pnt0, TkrPoint* pnt1, float act);

	// destructors
	~TkrNeuron() {}

	/**
	* @class Synapse
	* 
	* @brief Synapse is an inner class of TkrNeuron which defines the edges in a neural net graph.
	* 
	* Each neuron will have two lists of synapses, each associated with the 
    * connections to one end of the line segments.  To create an instance one
    * must pass it at least a pointer to a neuron.
	*
	*/
	class Synapse
	{
	public:

		// constructor
		Synapse(TkrNeuron* neuron, float wieght = 0.0):
          m_neuron(neuron), m_wieght(wieght) {}

		// destructor
		virtual ~Synapse() {}

		/** @name access methods
		*/
		//@{
		TkrNeuron* getNeuron() const { return m_neuron; }
		float getWieght()      const { return m_wieght; }
		//}@

		/** @name manipulation methods
		*/
		//@{
		void setWieght(float wieght) { m_wieght = wieght;}
		//}@

	private:

		// data members

		/// pointer to neuron which shares the connection
		TkrNeuron* m_neuron;
		/// wieght or cost of having both neurons active
		float      m_wieght;

	};

	// definitions
	typedef std::vector<Synapse*> SynapseList;

	/** @name access methods
	*/
	//@{	
	// in the following 'position' is an enum which has values of 'top' or 'bottom'.
	Point  getPnt(position pos)   const	{
        return pos==top ? m_pnt0->getPoint() : m_pnt1->getPoint();}
	int    getLayer(position pos) const	{
        return pos==top ? m_pnt0->getLayer() : m_pnt1->getLayer();}
	int    getTower(position pos) const	{
        return pos==top ? m_pnt0->getTower() : m_pnt1->getTower();}
	int	   getIdX(position pos)   const	{
        return pos==top ? m_pnt0->getIdX() : m_pnt1->getIdX();}
	int	   getIdY(position pos)	  const {
        return pos==top ? m_pnt0->getIdY() : m_pnt1->getIdY();}
	int    getLayerDiff()		  const {return m_l;}
	Vector getDirection()		  const {return m_direction;}
	float  getActivity()          const {return m_activity;}
	float  getBias()              const {return m_bias;}

	// access methods for synapses
	unsigned int numSynapse(position pos) const 
	{ return pos==top ? m_synapseList0.size():m_synapseList1.size(); }

	// getNextNeuron pre: the neuron pointed to by the synapse is != NULL.
	TkrNeuron*   getNextNeuron(position pos, unsigned int syn) const;		
	float		 getWieght(position pos, unsigned int syn) const;
	//@}

	/** @name manipulation methods
	*/
	//@{
	void addConnection(float (TkrNeuralNet::*func)(TkrNeuron*,TkrNeuron*), 
		               TkrNeuralNet* const net, position pos, TkrNeuron* neuron);
	void setActivity(float act) { 
        if((act >= 0.0) && (act <= 1.0)) m_activity = act;}
	void setBias(float bias)    { m_bias = bias;}
	//@}

	/** @name other methods
	*/
	//@{
	TkrNeuron& operator=(const TkrNeuron& neuron);  // assignment operator
	bool operator==(const TkrNeuron& neuron) const; // equality operator
	//@}

private:

	// data members

	/// hits and info about hits
	TkrPoint *m_pnt0, *m_pnt1;
	/// level of activity (a float betweeen 0 and 1)
	float       m_activity;
	/// l = layer0 - layer1  (layer1 is always greater than layer0)
	int         m_l;
	/// bias parameter to be used in some min. finding alg.
	float       m_bias;
	/// unit vector pointing from pnt0 to pnt1
	Vector      m_direction;
	/// List of synapses at pnt0
	SynapseList m_synapseList0; 
	/// List of synapses at pnt1
	SynapseList m_synapseList1; 
};

#endif // __TKRNEURON_H
