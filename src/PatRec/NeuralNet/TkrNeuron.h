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
* $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/PatRec/NeuralNet/TkrNeuron.h,v 1.2 2004/09/08 15:32:44 usher Exp $
*/

#ifndef __TKRNEURON_H
#define __TKRNEURON_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "src/Utilities/TkrPoint.h"
#include <vector>
#include <assert.h>
#include <iostream>

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
        Synapse(TkrNeuron* neuron, float weight = 0.0):
          m_neuron(neuron), m_weight(weight) {}

        // destructor
        virtual ~Synapse() {}

        /** @name access methods
        */
        //@{
        TkrNeuron* getNeuron() const { return m_neuron; }
        float getWeight()      const { return m_weight; }
        //}@

        /** @name manipulation methods
        */
        //@{
        void setWeight(float weight) { m_weight = weight;}
        //}@

    private:

        // data members

        /// pointer to neuron which shares the connection
        TkrNeuron* m_neuron;
        /// weight or cost of having both neurons active
        float      m_weight;

    };

    // definitions
    typedef std::vector<Synapse*> SynapseList;

    /** @name access methods
    */
    //@{    
    // in the following 'position' is an enum which has values of 'top' or 'bottom'.
    Point  getPnt(position pos)   const {
        return pos==top ? m_pnt0->getPoint() : m_pnt1->getPoint();}

    int    getLayer(position pos) const {
        return pos==top ? m_pnt0->getLayer() : m_pnt1->getLayer();}

    int    getTower(position pos) const {
        return pos==top ? m_pnt0->getTower() : m_pnt1->getTower();}

    int    getIdX(position pos)   const {
        return pos==top ? m_pnt0->getIdX() : m_pnt1->getIdX();}

    int    getIdY(position pos)   const {
        return pos==top ? m_pnt0->getIdY() : m_pnt1->getIdY();}

    int    getLayerDiff()         const {return m_l;}
    Vector getDirection()         const {return m_direction;}
    float  getActivity()          const {return m_activity;}
    float  getBias()              const {return m_bias;}

    // access methods for synapses
    unsigned int numSynapse(position pos) const 
      { return pos==top ? m_synapseList0.size():m_synapseList1.size(); }

    // getNextNeuron pre: the neuron pointed to by the synapse is != NULL.
    TkrNeuron*   getNextNeuron(position pos, unsigned int syn) const;       
    float        getWeight(position pos, unsigned int syn) const;
    //@}

    /** @name manipulation methods
    */
    //@{
    void addConnection(position, TkrNeuron*, double, double);

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

    /**
     * update()
     *
     * @param temp is the tempurature of the system.
     * 
     * This function updates a given neurons activity level based on the update
     * rule.  Neurons which are incoming or outgoing reinforce the neuron.
     * Neurons which are competing detract from the neurons activity.
     *
     * There currently many global parameters used in this function.  I would
     * like to change it so that these parameters are passed from the joboptions
     * file.
     */
    void update(float, double, double, double);


    friend std::ostream & operator<<(std::ostream& os, const TkrNeuron& n);

private:

    /**
     * calcWeight()
     *
     * This function returns the weight value for the two neurons' connection.
     * The formula is:
     * \f[
     * T_{ijk} = \frac{\cos^{\lambda}\psi_{ijk}}{d^{\mu}_{ij} + d^{\mu}_{jk}}
     * \f]
     */
     float calcWeight(TkrNeuron*,double,double);

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

typedef std::vector<TkrNeuron> TkrNeuronList;
typedef std::vector<TkrNeuron>::const_iterator neuron_const_iterator;

class TopToBottom
{
public:
  bool operator()(const TkrNeuron& N0, const TkrNeuron& N1)
    {
      double z0 = N0.getPnt(top).z();
      double z1 = N1.getPnt(top).z();
      return (z0>z1);
    }
};

class BottomToTop
{
public:
  bool operator()(const TkrNeuron& N0, const TkrNeuron& N1)
    {
      double z0 = N0.getPnt(bottom).z();
      double z1 = N1.getPnt(bottom).z();
      return (z0<z1);
    }
};


#endif // __TKRNEURON_H
