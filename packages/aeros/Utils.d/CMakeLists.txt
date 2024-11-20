ADD_LIBRARY(Util
		Connectivity.C
		dofset.C
		BlockAlloc.C
		DistHelper.C
		MFTT.C
		BinaryOutputFile.C
		CompositeInfo.C
		pstress.C
		b2Rotation.C
		list.C
		BinaryResultFile.C
		NodeSpaceArray.C
		MathUtils.C
		dbg_alloca.C
		SolverInfo.C
		GlobalToLocalMap.C
		Conwep.d/BlastLoading.C
		Scotch.C
		AutoTimer.h
		BinFileHandler.h
		)
# TODO Find out how to make the .h public without trouble for the users of Util in other directories.
TARGET_SOURCES(Util PRIVATE ConnectivityT.h ConnectivityTImpl.h)
TARGET_LINK_LIBRARIES(Util Eigen3::Eigen)
