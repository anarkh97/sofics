ADD_LIBRARY(Sfem Sfem.C cijk.C chaos.C)
TARGET_LINK_LIBRARIES(Sfem Eigen3::Eigen)
