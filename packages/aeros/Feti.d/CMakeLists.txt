ADD_LIBRARY(Feti
        CoarseSetCore.C
        CoarseSet.C
        CoarseSet.h
        CGOrthoSet.C CGOrthoSet.h
        DistrVectorCore.C
        DistrVector.h
        DistrVectorSet.h
        FetiCore.C
        FetiOpCore.C
        FetiOpControler.h
        FetiDP.C
        Feti.C
        Feti.h
        FetiInfo.h
        FetiOp.h
        FetiOp.C
        NLFetiCore.C
        CCtSolver.d/BlockCCt.C
        CCtSolver.d/CCtSolver.h
        CCtSolver.d/GlobalCCt.h
        CCtSolver.d/GlobalCCt.C
        CCtSolver.d/SuperBlockCCt.h
        CCtSolver.d/SuperBlockCCt.C
        FetiSub.h
		FETISolver.cpp FETISolver.h
		FetiSub.cpp
		Operators.cpp Operators.h
		FetiBaseSub.C
		FetiBaseClass.cpp FetiBaseClass.h
		CornerSelector.cpp CornerSelector.h
		MultiLevelFetiDP.cpp
		DecomposedSystem.cpp DecomposedSystem.h
		SetOfSubs.cpp SetOfSubs.h
		)

TARGET_LINK_LIBRARIES(Feti Eigen3::Eigen)

SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
