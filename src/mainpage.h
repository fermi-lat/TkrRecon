// Mainpage for doxygen

/** @mainpage TkrRecon
 *
@section introduction Introduction
 *
TkrRecon contains the code that reconstructs tracks from
hit strips in the TKR, and then reconstructs vertices from the tracks. 

It's organized as a series of Gaudi algorithms that act
successively, taking input from, and sending output to,
the Transient Data Store. Starting from TkrDigi objects, it generates the
TDS objects TkrCluster, then TkrPatCand, TkrFitTrack and finally, TkrVertex.

The major driver algorithms are written to support multiple 
implementations of a particular reconstruction. The interface between 
the algorithm and the implementation appears as a Gaudi tool. 
For example, TkrFindAlg currently has three implementations
in various stages of development:
TkrComboPatRec, TkrLinkAndTree and TkrNeuralNet, with Gaudi tool 
interfaces ComboFindTrackTool, LinkAndTreeFindTrackTool and 
NeuralNetFindTrackTool, respectively. 
 *
In order to simplify the interface to the casual user, the Tracker
reconstruction chain is controlled by a top level algorithm, TkrReconAlg, 
which treats the major driver algorithms as Gaudi subalgorithms. 
 *
The diagram below attempts to illustrate this procedure:
 *
 \image html Tkr_diagram.gif
 *
 *
Here is an example of a reconstructed gamma. The green x's represent the
cluster positions reconstructed from the hit strips. After this, the blue
lines represent 
the fit tracks found from the pattern recognition. The yellow line is the 
resultant direction from the vertexing algorithm.
 *
 *
 \image html Event1.gif
 *
@section algorithms Major algorithms

The algorithms themselves don't contain the "algorithms." They merely handle
the Gaudi manipulations and then pass the work on to other classes, which
in turn, operate on the TDS classes. The main algorithms in the chain are:
 *
TkrClusterAlg: The work is done by TkrMakeClusters, which groups hits into 
clusters. Takes as input the digitized hits in TkrDigi, and fills 
TkrClusterCol with TkrCluster objects. Another class, TkrQueryClusters, has
some methods for doing some detailed calculations using the clusters.
 *
TkrFindAlg: The work is currently done by one of ComboFindTrackTool, 
LinkAndTreeFindTrackTool or NeuralNetFindTrackTool which call, 
respectively,  TkrComboPatRec, TkrLinkAndTree, or TkrNeuralNet. 
These find "all" possible track candidates, taking as input the list of 
clusters and outputting TkrPatCand objects (collected in the TDS object 
TkrCandidateCol).
 *
TkrTrackFitAlg: The work of this algorithm is currently done by one of 
TkrComboFitTool, TkrLinkAndTreeFitTool or TkrNeuralNetFitTool. These tools 
take the candidate tracks from their respective pattern recognition 
analogues and perform a Kalman Filter track fit on them through 
the use of the KalFitTrack class. This allows for the slight differences 
which may exist in the output of a particular pattern recognition 
implementation. Successfully fit tracks, appearing to the TDS as 
TkrFitTrack objects, are then stored in the TDS TkrFitTrackCol collection. 
 *
TkrVertexAlg: The work of this algorithm is currently done by one of 
ComboVtxTool, VtxKalFitTool or VtxSingleTrkTool. The first of these 
follows the TkrComboVtxRecon method of finding and reconstructing 
vertices from the TkrFitTrack objects. VtxKalFitTool finds and reconstructs
vertices in the case of more than one TkrFitTrack objects using a Kalman 
Filter vertexing algorithm. For more information, follow this 
<A HREF="../../doc/KalVtx/KalVtx.html"> link</A>.

VtxSingleTrkTool returns a vertex for the 
case of single tracks (fairly frequent). All methods return TkrVertex 
objects which are stored in the TDS TkrVertexCol collection. 


TkrReconAlg: The main algorithm controlling the above reconstruction sequence. 
 *
@section access Access to the "interesting" quantities
 *
The abstract class TkrRecInfo provides the user of TkrRecon TDS objects 
with a common interface to the basic information about TkrPatCand, 
TkrFitTrack and TkrVertex objects. The class lives in Event, 
so that other subsystems may access this information without having 
to depend on TkrRecon.
 *
The information includes the track (or vertex) parameters and 
error matrices at both ends of the track (or one "end" for the vertex), 
the energy and quality.

@section services Services

TkrGeometrySvc: an interface to GlastDetSvc, it supplies 
the geometry constants. 
 *
TkrBadStripsSvc: manages the list of bad strips. 
This service currently interacts with TkrClusterAlg, but ultimately will 
also be used by other algorithms.
 *
TkrInitSvc: Used to initialize the TkrControl object which provides 
parameters for controlling the TkrRecon reconstruction process. 
 *
@section information More Information
 *
For more information, please see the Tracker Recon Group 
<A HREF="http://www-glast.slac.stanford.edu/software/TKR/NewTracker/TrkRecon.htm"> 
website. </A>
 *
@authors Tracy Usher, Leon Rochester (SLAC) <br>
Bill Atwood, Brandon Allgood(UC Santa Cruz) <br>
Michael Kuss, Johann Cohen-Tanugi (INFN Pisa) <br>
emeritus: Jose Angel Hernando Morata (Santiago de Compostela)
<hr>
@section requirements CMT Requirements
@verbinclude requirements
<hr>
@section notes Release notes
release.notes
<hr>
@todo
Develop mechanism to handle constants required by code. Current plan 
is to put the constants in an xml file, and to make them properties 
as well, so that they can be over-written from the jobOptions.
 *
Consolidate references to geometry - remove from GFtutor and other places
 *
Root out last hardwired numbers
 *
Doxygenate remainder of code
 */
