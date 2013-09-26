AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = planet.dox planet.page txt_common fig_common

DISTCLEANFILES  = planet.dox
DISTCLEANFILES += Makefile

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

