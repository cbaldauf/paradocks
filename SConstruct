#manual configuration
#compile error check code
ERROR=1
#compile debug code
DEBUG=0
#compile debug code for more cases
DEBUG2=0

if ERROR:
	cflags="-DERROR "

if DEBUG:
	cflags=cflags+"-DDEBUG "

if DEBUG2:
	cflags=cflags+"-DDEBUG2 "
	#no DEBUG2 without DEBUG
	if not DEBUG:
		cflags=cflags+"-DDEBUG "

import os
env = Environment(ENV = {'PATH' : os.environ['PATH']}) #, tools = ['default', 'doxygen'], toolpath = ['doc'])

#if you need to adjust your include and lib path to non standard places put them in this line sepperated by spaces
#for example:
env.Append(CPPPATH=Split("""#"""))
env.Append(LIBPATH=Split("""#build"""))

#CFLAGS for complilaton
env.Append(CCFLAGS = cflags+'-O2 -fomit-frame-pointer -pipe')

if not env.GetOption('clean'):
	conf = Configure(env)
	if not conf.CheckLib('boost_system'):
		print 'Can not find libboost_system from Boost C++ Libraries!'
		Exit(1)
	env = conf.Finish()


# build documentation
#env.Doxygen('doc/Doxyfile')

# build the molecule graph lib
SConscript(['mgraph/SConscript'],variant_dir='mgraph/build', duplicate=0, exports=("""env"""))

SConscript(['math/SConscript'], variant_dir='math/build', dublicate=0, exports=("""env"""))

# buld fitness function common base class
SConscript(['fitness/SConscript'],variant_dir='fitness/build', duplicate=0, exports=("""env"""))
# build xscore lib
SConscript(['fitness/pscore/SConscript'],variant_dir='fitness/pscore/build', duplicate=0, exports=("""env"""))
# build pmf04 potential lib
SConscript(['fitness/pmf04/SConscript'],variant_dir='fitness/pmf04/build', duplicate=0, exports=("""env"""))

# build optimizer common base class
SConscript(['optimizer/SConscript'],variant_dir='optimizer/build', duplicate=0, exports=("""env"""))
SConscript(['optimizer/pso/SConscript'],variant_dir='optimizer/pso/build', duplicate=0, exports="""env""")
# build random lib                                                                                                                                                                      
SConscript(['optimizer/random/SConscript'],variant_dir='optimizer/random/build', duplicate=0, exports="""env""")                                                                        
# build random hill climbing lib                                                                                                                                                        
SConscript(['optimizer/random_hill_climbing/SConscript'],variant_dir='optimizer/random_hill_climbing/build', duplicate=0, exports="""env""")
SConscript(['optimizer/nelder_mead/SConscript'],variant_dir='optimizer/nelder_mead/build', duplicate=0, exports="""env""")
SConscript(['optimizer/hbee/SConscript'],variant_dir='optimizer/hbee/build', duplicate=0, exports="""env""")

# build the framework and the program
SConscript(['framework/SConscript'], variant_dir='build', duplicate=0, exports=("""env"""))

SConscript(['rmsd/SConscript'], variant_dir='rmsd/build', duplicate=0, exports=("""env"""))
