/**
* @class TkrNeuralNet
*
* @brief The TkrNeuralNet class searches for track candidates using a neural net.
*
* Based on the work of: 
* Stimpfl-Abele & Garrido Comp. Phys. Commun. 64 (1991) 46-50
*
* last modified 3/02
* 
* @authors b. allgood and w. atwood
*
* @todo Impliment graph searching function to better pick out candidate tracks.
* @todo Move global parameter definitions to joboptions file
*
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/NeuralNet/TkrNeuralNet.h,v 1.4 2002/05/10 21:53:57 usher Exp $
*/

#ifndef __TKR_NEURALNET_H
#define __TKR_NEURALNET_H

#include "Event/Recon/TkrRecon/TkrPatCandCol.h"
#include "src/TrackFit/KalFitTrack/KalFitTrack.h"
#include "src/PatRec/NeuralNet/TkrNeuron.h"
#include "TkrRecon/Track/TkrPoints.h"
#include "TkrRecon/Track/TkrPoint.h"
#include "TkrRecon/Track/TkrBase.h"
#include <vector>

using namespace Event;

class TkrNeuralNet : public TkrPatCandCol 
{
public:

	// definitions
	typedef TkrBase Candidate;
    typedef std::vector<Candidate> CandidateList; 
    typedef std::vector<Candidate>::const_iterator const_iterator;
	typedef std::vector<TkrNeuron> TkrNeuronList;
    typedef std::vector<TkrNeuron>::const_iterator neuron_const_iterator;
	typedef std::vector<TkrPoint> TkrPointList;


	// constructor
	TkrNeuralNet(ITkrGeometrySvc* pTkrGeo, TkrClusterCol* pClusters,
		double calEne = 0., Point calHit = Point(0.,0.,0.));

	// destructor
	~TkrNeuralNet() {}

	/** @name access methods
	*/
	//@{
	unsigned int numNeurons()         const {return m_numNeurons;}
    const CandidateList& candidates() const {return m_candidates;}
    const_iterator begin()            const {return m_candidates.begin();}
    const_iterator end()              const {return m_candidates.end();}
    int numCandidates()               const {return (int) m_candidates.size();}
    const TkrNeuronList& neurons()    const {return m_neuronList;}
    //@}

private:

	/** @name inner methods
    */
    //@{

    /**
     * generateNeurons()
     * 
     * This function first gathers the cluster pairs and stores them in 
     * m_pointList.  It then constructs neurons and adds them to m_neuronList
     * if they meet the criteria of layer seperation and max pitch angle.
     * Each neuron is assigned a random bias using the standard rand() C++ call.
     * We may want to have this use randSvc() from Gaudi?  This returns
     * the total number of neurons created.
     */
	unsigned int generateNeurons();

    /**
     * buildNet()
     * 
     * This function constructs the synapse lists for all of the neurons.  The
     * synapse lists contains pointers to other neurons sharing the same end
     * points and the wieghts associated with that connection.
     */
	void buildNet();

    /**
     * relax()
     *
     * This function evolves the neural net to the ground state.  This is done
     * by shuffling the neurons and then updating the neurons using the update 
     * function (explained below).  This is repeated until the normalized 
     * difference in activity of the last two update rounds is less then some
     * convergence value.
     */
	void relax();

    /**
     * buildCand()
     *
     * This is function is based on one taken from the TkrCombo PatRec routine.
     * This funciton takes all neurons with an activity above 0.9 and places
     * them in m_candidates.  It then proceed to do some preliminary fitting
     * of the track.  This will be changed soon.
     */
    void buildCand();

    /**
     * calcWieght()
     *
     * This function returns the weight value for the two neurons' connection.
     * The formula is:
     * \f[ 
     * T_{ijk} = \frac{\cos^{\lambda}\psi_{ijk}}{d^{\mu}_{ij} + d^{\mu}_{jk}}
     * \f]
     */
	float calcWieght(TkrNeuron* neuron1, TkrNeuron* neuron2);

    /**
     * updateNeurons()
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
    void updateNeuron(TkrNeuron* neuron,float temp);
    //@}

	//data members

    /// list of candidate tracks to be passed to the Kalman fit.
	CandidateList  m_candidates;
    /// list of tracks to be used with Kalman fit.
	TkrFitCol      m_tracks;
    /// list of all neurons
	TkrNeuronList  m_neuronList;
    /// list of all (x,y,z) points (not all are used in neurons)
	TkrPointList   m_pointList;
    /// number of neurons in m_neuronList
	unsigned int   m_numNeurons;
    /// number of points in m_pointList
	unsigned int   m_numPoints;
    /// position of cal hit
    Point          m_Pcal;
    /// energy for the event
    double         m_energy;

};

#endif // __TKR_NEURLNET_H