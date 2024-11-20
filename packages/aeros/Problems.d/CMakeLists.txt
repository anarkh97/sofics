add_library (Problems EigenDescr.C  
       SingleDomainBase.C
       DynamDescr.C  
       NonLinStatic.C
       CondDescr.C   
       NonLinDynam.C 
       TempDescr.C   
       ModalBase.C   
       DEMProblem.C)
TARGET_LINK_LIBRARIES(Problems Eigen3::Eigen)
