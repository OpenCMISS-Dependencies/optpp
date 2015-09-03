add_library(gss GenSetBase.C GenSetBox2d.C GenSetMin.C GenSetStd.C OptGSS.C ../globals.C)

INSTALL(TARGETS gss
    EXPORT optpp-config
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include/optpp
)
