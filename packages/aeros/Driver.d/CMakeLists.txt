ADD_LIBRARY(Driver DecDomainCore.C
		GeoSourceCore.C HData.C Mpc.C BaseSub.C
		MultiFront.C Dynam.C TempDynam.C
		NLDynam.C NLStatic.C Static.C
		SubDomainCore.C Domain.C SComm.C
		Eigen.C MPIComm.C SubDomainFactory.C
		PolygonSet.C Sower.C CornerMaker.C
		jacobi.f pade.f zpade.f
		BinaryOutput.C
		DomainPita.C
		GeoSourceInstance.C OpMakeInstance.C
		ControlLawInfo.C BinaryOutputInstance.C
		EFrameDataInstance.C CommunicatorCore.cpp
		DecDomain.h
		DecDomainImpl.h
		LOpsImpl.h
		SubDomain.C)
TARGET_SOURCES(Driver PRIVATE GeoSource.h)
TARGET_LINK_LIBRARIES(Driver Eigen3::Eigen)
