ADD_LIBRARY( paral MultiDomainBase.C MDDynam.C MDNLStatic.C MDNLDynam.C MDTemp.C MDModalBase.C MDModal.C
		DomainGroupTask.C GenMS.h)
TARGET_LINK_LIBRARIES(paral Eigen3::Eigen)
