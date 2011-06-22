add_library(gss GenSetBase.C GenSetBox2d.C GenSetMin.C GenSetStd.C OptGSS.C ../globals.C)

INSTALL(TARGETS gss
	ARCHIVE DESTINATION lib/optpp
)
