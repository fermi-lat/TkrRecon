# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/TkrRecon/TkrReconLib.py,v 1.1 2008/08/15 21:42:37 ecephas Exp $
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrRecon'])
    env.Tool('CalUtilLib')
    env.Tool('CalibDataLib')
    env.Tool('GlastSvcLib')
    env.Tool('guiLib')
def exists(env):
    return 1;
