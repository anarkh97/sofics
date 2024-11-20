ADD_LIBRARY(Material Material.cpp LinearElasticMaterial.cpp 
            IsotropicLinearElasticJ2PlasticMaterial.cpp
            IsotropicLinearElasticJ2PlasticPlaneStressMaterial.cpp
            KorkolisKyriakidesPlaneStressMaterial.cpp
            KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding.cpp 
            KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2.cpp)
TARGET_LINK_LIBRARIES(Material Eigen3::Eigen)
