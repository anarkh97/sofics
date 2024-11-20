ADD_LIBRARY(Timers Timing.C StaticTimers.C DistTimer.C
		GetTime.C NLTimers.C FetiDPtimers.C )
TARGET_LINK_LIBRARIES(Timers Eigen3::Eigen)
