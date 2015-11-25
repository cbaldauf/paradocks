all:
	@scons -Q -j2

clean:
	@scons -Q --clean
	@rm -f `find . -name "*~"`
	@rm -fr .sconf_temp
	@rm -f `find . -name ".sconsign*"`
	@rm -f config.log		
	@rm -fr `find . -name build`
