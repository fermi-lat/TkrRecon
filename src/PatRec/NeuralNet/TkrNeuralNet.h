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
* @todo replace std::rand in generateNeurons() with call to Gaudi Service
*
* $Header: /nfs/slac/g/glast/ground/cvs/users/TkrGroup/TkrRecon/src/PatRec/NeuralNet/TkrNeuralNet.h,v 1.2 2004/09/08 15:32:44 usher Exp $
*/

#ifndef __TKR_NEURALNET_H
#define __TKR_NEURALNET_H

#include "TkrUtil/ITkrQueryClustersTool.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "src/PatRec/NeuralNet/TkrNeuron.h"
#include "src/Utilities/TkrBase.h"
#include "GaudiKernel/DataObject.h"
#include <vector>
#include <map>

class TkrNeuralNet : public DataObject 
{
 public:

  struct ltstr
  {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) < 0;
    }
  };
  
  
  // constructor
  TkrNeuralNet(Event::TkrClusterCol* pClusters, ITkrQueryClustersTool* clusTool,
	       std::map<const char*,double, ltstr>& params,
	       double calEne = 0., Point calHit = Point(0.,0.,0.));
  
  // destructor
  ~TkrNeuralNet() {}
  
  /** @name access methods
   */
  //@{
  unsigned int numNeurons()         const {return m_numNeurons;}
  const std::vector<TkrBase>& candidates() const {return m_candidates;}
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

  
 private:
  
  //data members
  
  /// list of candidate tracks to be passed to the Kalman fit.
  std::vector<TkrBase>  m_candidates;
    
  /// list of all neurons
  TkrNeuronList  m_neuronList;
      
  /// list of all (x,y,z) points (not all are used in neurons)
  TkrPointList   m_pointList;

  /// number of neurons in m_neuronList
  unsigned int   m_numNeurons;
	  
  /// position of cal hit
  Point          m_Pcal;

  /// energy for the event
  double         m_energy;
	      
  Event::TkrClusterCol*   m_clusters;

  ITkrQueryClustersTool*  m_clusTool;

  std::map<const char*,double, ltstr>& m_params;
};

#endif // __TKR_NEURLNET_H
