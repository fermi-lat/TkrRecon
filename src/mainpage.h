// Mainpage for doxygen

/** @mainpage TkrRecon
 *
 * @section introduction Introduction
 *
 * TkrRecon contains the code that reconstructs tracks from
 * hit strips in the TKR.
 * 
 * It's organized a a series of algorithms that act
 * successively, taking input from and sending output to
 * the Transient Data Store.
 * 
 * The major driver algorithms are written to support multiple 
 * implementations. 
 * For example, TkrFindAlg currently has two implementations:
 * ComboPR and LinkAndTreePR.
 *
 * A diagram of the program flow can be seen
 * <A href= "../images/Tkr_diagram.gif>here</A>.
 * 
 * @section algorithms Major algorithms
 * 
 * The algorithms are:
 *
 * TkrClusterAlg: Groups hits into clusters. Takes as input the digitized hits in TkrDigi,
 * and fills TkrClusters
 *
 * TkrFindAlg: Finds possible tracks. Takes the list of clusters as input and fills TkrCandidates. 
 *
 * TkrReconAlg: (bad name! It's called TkrFitAlg in the diagram.) Fits the candidate tracks and fills TkrTracks.
 *
 * TkrVertexAlg: Tries to find vertices among the TkrTracks, to produce TkrVertices. (not yet implemented)
 *
 * @section services Services
 * 
 * TkrGeometrySvc: supplies the geometry constants 
 *
 * TkrBadStripsSvc: manages the list of bad strips.
 * this service currently interacts with TkrClusterAlg, but 
 * ultimately will also be used by other algorithms.
 *
 * TkrInitSvc: allows for loading of code required by the
 * various specific implementations of the algorithms.
 * This service will probably no longer be necessary after
 * the algorithms are taught about sub-algorithms.
 *
 * @authors Tracy Usher, Bill Atwood, Leon Rochester
 * <hr>
 * @section requirements CMT Requirements
 * @verbinclude requirements
 * <hr>
 * @section notes Release notes
 * release.notes
 * <hr>
 * @todo
 * 
 */
