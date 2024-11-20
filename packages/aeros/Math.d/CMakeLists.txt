ADD_LIBRARY(Math spsmvp.f
        sptmv.f cspsmvp.f cdspsmvp.f
        matrix.C matrix.h
        FullSquareMatrix.h FullSquareMatrix.C
        FullRectMatrix.C
        IntFullM.C
        SparseData.C
        matrixFull.C BigMatrixCore.C BLKSparseMatrixCore.C
        CuCSparseCore.C DBSparseMatrixCore.C
        FullRectMatrix.C NBSparseMatrixCore.C TTensor.C
        EiSparseMatrixCore.C
        SparseMatrix.h
        CuCSparse.h
        BLKSparseMatrix.h
        DistVector.h
        DiagMatrix.h
        Vector.h
Skyline.d/forbackr7.f
Skyline.d/BlockSkyCore.C  Skyline.d/forbackr7ns.f
Skyline.d/forbackr8.f
Skyline.d/forbackr8ns.f
Skyline.d/ffulpvm.f       
Skyline.d/forbackc1.f    
Skyline.d/forbackc2.f   
Skyline.d/forbackc3.f     
#sgi only Skyline.d/pfactorC.f
Skyline.d/forbackc4.f 
#sgi only Skyline.d/pfactor.f
Skyline.d/forbackr1.f
Skyline.d/forbackr1ns.f   Skyline.d/SkyData.C
Skyline.d/forbackr2.f     
Skyline.d/forbackr2ns.f   Skyline.d/SkyMatrixCore.C
Skyline.d/forbackr3.f     Skyline.d/skymul.f
Skyline.d/forbackr3ns.f   Skyline.d/slacol.f
Skyline.d/forbackr4.f     Skyline.d/svbu4cb.f
Skyline.d/forbackr4ns.f   Skyline.d/svbu4gmb.f
Skyline.d/forbackr5.f     Skyline.d/svbu4rb.f
Skyline.d/forbackr5ns.f   Skyline.d/svbu4rb.save.f
Skyline.d/forbackr6.f     Skyline.d/utility.C
Skyline.d/forbackr6ns.f
ComplexFactorSparse.d/zassmb.f    ComplexFactorSparse.d/zblkslv.f
ComplexFactorSparse.d/zblkl2.f    ComplexFactorSparse.d/zblkslvp.f
ComplexFactorSparse.d/zblkldl.f   ComplexFactorSparse.d/zgecp.f
ComplexFactorSparse.d/zblkns.f    ComplexFactorSparse.d/zgers.f
ComplexFactorSparse.d/zblkslv2.f  ComplexFactorSparse.d/zmmpyi.f
ComplexFactorSparse.d/zblkslv3.f  ComplexFactorSparse.d/zpotf2_ldl.f
ComplexFactorSparse.d/zblkslv4.f  ComplexFactorSparse.d/zzeromat.f
FactorSparse.d/addMat.f      FactorSparse.d/epost2.f  FactorSparse.d/invinv.f
FactorSparse.d/addOne.f      FactorSparse.d/etordr.f  FactorSparse.d/ldindx.f
FactorSparse.d/assmb.f       FactorSparse.d/etpost.f  FactorSparse.d/lstats.f
FactorSparse.d/betree.f      FactorSparse.d/etree.f   FactorSparse.d/ltfrep.f
FactorSparse.d/bfinit.f      FactorSparse.d/fcnthn.f  FactorSparse.d/main.f
FactorSparse.d/blkl2.f       FactorSparse.d/fntsiz.f  FactorSparse.d/mmdelm.f
FactorSparse.d/blkldl.f      FactorSparse.d/fsup1.f   FactorSparse.d/mmdint.f
FactorSparse.d/blkns.f       FactorSparse.d/fsup2.f   FactorSparse.d/mmdnum.f
FactorSparse.d/genmmd.f  FactorSparse.d/mmdupd.f
FactorSparse.d/blkslv.f      FactorSparse.d/getadj.f  FactorSparse.d/mmpyi.f
FactorSparse.d/blkslvp.f     FactorSparse.d/getnrm.f  FactorSparse.d/ordmmd.f
FactorSparse.d/btree2.f      FactorSparse.d/getres.f  FactorSparse.d/ordnat.f
FactorSparse.d/chkns.f       FactorSparse.d/getrhs.f  FactorSparse.d/sfinit.f
FactorSparse.d/chordr.f      FactorSparse.d/gtimer.f  FactorSparse.d/symfc2.f
FactorSparse.d/dgecp.f       FactorSparse.d/igathr.f  FactorSparse.d/symfct.f
FactorSparse.d/dgers.f       FactorSparse.d/inpdof.f  FactorSparse.d/zero.f
FactorSparse.d/dpotf2_ldl.f  FactorSparse.d/inpmat.f
FactorSparse.d/dpotrf_ldl.f  FactorSparse.d/inpnv.f
SCMatrix.d/fpdlapiv.f        SCMatrix.d/fpslapiv.f
SCMatrix.d/fpdlapv2.f        SCMatrix.d/fpslapv2.f    SCMatrix.d/ilasrt2.f
SCMatrix.d/SCBaseMatrix.C
SCMatrix.d/SCIntMatrix.C
SCMatrix.d/SCDoubleMatrix.C BLAS.h
        Blockwise/BlockDiagonalMatrix.cpp Blockwise/BlockDiagonalMatrix.h
        Blockwise/BlockProduct.cpp Blockwise/BlockProduct.h )
TARGET_LINK_LIBRARIES(Math Eigen3::Eigen)
