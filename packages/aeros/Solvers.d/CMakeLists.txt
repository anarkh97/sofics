ADD_LIBRARY(solver
        Solver.h
        DSCsolver.C
        MumpsCore.C
        Rbm.C
        MappedAssembledSolver.C
        SolverCore.C
        SpoolesCore.C
        GoldfarbIdnani.C
        SolverFactoryCore.C
		SolverFactory.C
        MultiDomainSolver.h MultiDomainSolver.C
		CholmodSolver.cpp CholmodSolver.h
		PardisoSolver.h PardisoSolver.cpp
		)
if(USE_CHOLMOD)
	message("Looking for CHOLMOD")
	find_package(Cholmod)
	message("Cholmod include: ${CHOLMOD_INCLUDE_DIR}")
	message("Cholmod libraries: ${CHOLMOD_LIBRARIES}")
	if (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
		message("Will use CholMod and ${CHOLMOD_FOUND}")
		find_library(RT_LIBRARY rt)
		if(RT_LIBRARY)
			list(APPEND CHOLMOD_LIBRARIES ${RT_LIBRARY})
		endif(RT_LIBRARY)
		# pass whether we use choldmod or not to the base CMakeLists.
		set (CHOLMOD_FOUND ${CHOLMOD_FOUND} PARENT_SCOPE)
		target_sources(solver
				PRIVATE
				CholmodImp.h CholmodImp.cpp
				)
		target_compile_definitions(solver PUBLIC WITH_CHOLMOD=1)
		target_include_directories(solver PRIVATE ${CHOLMOD_INCLUDE_DIR})
		target_link_libraries(solver PUBLIC ${CHOLMOD_LIBRARIES})
	endif (CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARIES)
endif(USE_CHOLMOD)

FIND_PACKAGE(BLAS REQUIRED)
if(USE_PARDISO)
	message("Looking for MKL Pardiso")
	if(BLAS_mkl_core_LIBRARY)
		message("We probably have Pardiso")
		set( PARDISO_FOUND true PARENT_SCOPE )
		target_compile_definitions(solver PRIVATE WITH_PARDISO)
		target_include_directories(solver PRIVATE /opt/intel/mkl/include)
	else(BLAS_mkl_core_LIBRARY)
		set( PARDISO_FOUND false PARENT_SCOPE )
		message("We probably do not have Pardiso")
	endif(BLAS_mkl_core_LIBRARY)
endif(USE_PARDISO)


TARGET_LINK_LIBRARIES(solver PUBLIC Eigen3::Eigen)

#SET(INC_DIR ..)

#INCLUDE_DIRECTORIES( ${INC_DIR} )

