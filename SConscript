# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/SConscript,v 1.8 2008/12/08 21:40:11 ecephas Exp $ 
# Authors: Leon Rochester <lsrea@slac.stanford.edu>, Tracy Usher <usher@slac.stanford.edu>
# Version: TkrRecon-10-16-05
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('TkrReconLib', depsOnly = 1)
TkrRecon = libEnv.SharedLibrary('TkrRecon', listFiles(['src/Dll/*.cxx','src/GaudiAlg/*.cxx',
                                                       'src/Display/*.cxx','src/Services/*.cxx',
                                                       'src/Cluster/*.cxx','src/Filter/*.cxx',
                                                       'src/PatRec/*.cxx','src/PatRec/LinkAndTree/*.cxx',
                                                       'src/PatRec/Combo/*.cxx','src/PatRec/NeuralNet/*.cxx',
                                                       'src/PatRec/VectorLinks/*.cxx',
                                                       'src/PatRec/VectorLinks/StdRelTable/*.cxx',
                                                       'src/PatRec/MonteCarlo/*.cxx','src/PatRec/KalFitTrack/*.cxx'
                                                       'src/PatRec/Utilities/*.cxx','src/Track/*.cxx','src/TrackFit/*.cxx',
                                                       'src/TrackFit/KalmanFilter/*.cxx','src/TrackFit/KalmanFilterUtils/*.cxx',
                                                       'src/TrackFit/KalmanFilterFit/*.cxx','src/TrackFit/KalmanFilterFit/FitMatrices/*.cxx',
                                                       'src/TrackFit/KalmanFilterFit/HitErrors/*.cxx','src/TrackFit/KalmanFilterFit/TrackEnergy/*.cxx',
                                                       'src/TrackFit/KalFitTrack/*.cxx','src/TrackFit/LineFit2D/*.cxx','src/TrackFit/LineFit3D/*.cxx',
                                                       'src/Utilities/*.cxx','src/Vertex/*.cxx',
                                                       'src/Vertex/Combo/*.cxx',
                                                       'src/Vertex/DocaVtx/*.cxx']))
progEnv.Tool('TkrReconLib')
test_TkrRecon = progEnv.GaudiProgram('test_TkrRecon', listFiles(['src/test/*.cxx']), test = 1)

progEnv.Tool('registerObjects', package = 'TkrRecon', libraries = [TkrRecon], testApps = [test_TkrRecon],
	includes = listFiles(['TkrRecon/*'], recursive = 1))



