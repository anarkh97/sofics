LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE     := FemLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -DCREATE_DSO -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC 
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../main.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := AcmeLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Acme.d/enforcement/ContactEnfModel.C \
	../Acme.d/enforcement/ContactEnforcement.C \
	../Acme.d/enforcement/ContactEnforcementData.C \
	../Acme.d/enforcement/ContactEnforcementUserQuery.C \
	../Acme.d/enforcement/ContactEnfUserQuery.C \
	../Acme.d/enforcement/ContactEnfZoltan.C \
	../Acme.d/enforcement/ContactGapRemoval.C \
	../Acme.d/enforcement/ContactMPCs.C \
	../Acme.d/enforcement/ContactQSEnforcement.C \
	../Acme.d/enforcement/ContactTDAdhesion.C \
	../Acme.d/enforcement/ContactTDAreaWeld.C \
	../Acme.d/enforcement/ContactTDCohesiveZone.C \
	../Acme.d/enforcement/ContactTDConstantFriction.C \
	../Acme.d/enforcement/ContactTDEnfModel.C \
	../Acme.d/enforcement/ContactTDEnforcement.C \
	../Acme.d/enforcement/ContactTDEnfPenalty.C \
	../Acme.d/enforcement/ContactTDFaceFaceEnf.C \
	../Acme.d/enforcement/ContactTDFrictionless.C \
	../Acme.d/enforcement/ContactTDJunction.C \
	../Acme.d/enforcement/ContactTDPressureDependent.C \
	../Acme.d/enforcement/ContactTDPVDependent.C \
	../Acme.d/enforcement/ContactTDShared.C \
	../Acme.d/enforcement/ContactTDSpotWeld.C \
	../Acme.d/enforcement/ContactTDSpringWeld.C \
	../Acme.d/enforcement/ContactTDThreaded.C \
	../Acme.d/enforcement/ContactTDTied.C \
	../Acme.d/enforcement/ContactTDUser.C \
	../Acme.d/enforcement/ContactTDUserStubs.f \
	../Acme.d/enforcement/ContactTDVelocityDependent.C \
	../Acme.d/enforcement/ContactTiedKinematics.C \
	../Acme.d/enforcement/ContactVolumeTransfer.C \
	../Acme.d/enforcement/Enforcement_Interface.C \
	../Acme.d/enforcement/NodeEntityInteractionList.C \
	../Acme.d/search/ContactAnalyticCylinderInside.C \
	../Acme.d/search/ContactAnalyticCylinderOutside.C \
	../Acme.d/search/ContactAnalyticPlane.C \
	../Acme.d/search/ContactAnalyticSphere.C \
	../Acme.d/search/ContactAnalyticSurface.C \
	../Acme.d/search/ContactAsymComm.C \
	../Acme.d/search/ContactBlockEntityList.C \
	../Acme.d/search/ContactBoundingBox.C \
	../Acme.d/search/ContactBoundingBoxHierarchy.C \
	../Acme.d/search/ContactBoundingBoxHierarchy_Int.C \
	../Acme.d/search/ContactCommBuffer.C \
	../Acme.d/search/ContactCommList.C \
	../Acme.d/search/Contact_Communication.C \
	../Acme.d/search/ContactDoublyLinkedList.C \
	../Acme.d/search/ContactEdgeBlock.C \
	../Acme.d/search/ContactEdgeInstance.C \
	../Acme.d/search/ContactElementBlock.C \
	../Acme.d/search/ContactElementInstance.C \
	../Acme.d/search/ContactElementElementInteraction.C \
	../Acme.d/search/ContactEntityDataHash.C \
	../Acme.d/search/ContactErrors.C \
	../Acme.d/search/ContactExodusOutput.C \
	../Acme.d/search/ContactExodusResults.C \
	../Acme.d/search/ContactFaceBlock.C \
	../Acme.d/search/ContactFaceInstance.C \
	../Acme.d/search/ContactFaceCoverageInteraction.C \
	../Acme.d/search/ContactFaceFaceInteraction.C \
	../Acme.d/search/ContactFaceFaceSearchInstance.C \
	../Acme.d/search/ContactFixedSizeAllocator.C \
	../Acme.d/search/ContactGlobalSearch.C \
	../Acme.d/search/ContactHexElementL8Instance.C \
	../Acme.d/search/ContactHexOverlap.C \
	../Acme.d/search/ContactHostGlobalID.C \
	../Acme.d/search/ContactInteractionEntity.C \
	../Acme.d/search/ContactLineEdgeL2Instance.C \
	../Acme.d/search/ContactLineEdgeQ3.C \
	../Acme.d/search/ContactLineFaceL2.C \
	../Acme.d/search/ContactLineFaceQ3.C \
	../Acme.d/search/ContactNodeBlock.C \
	../Acme.d/search/ContactNodeInstance.C \
	../Acme.d/search/ContactNodeEntityInteraction.C \
	../Acme.d/search/ContactNodeFaceInteraction.C \
	../Acme.d/search/ContactNodeNodeInteraction.C \
	../Acme.d/search/ContactNodeSurfaceInteraction.C \
	../Acme.d/search/ContactParOStream.C \
	../Acme.d/search/ContactPolygon.C \
	../Acme.d/search/ContactProcessFaceCoverage.C \
	../Acme.d/search/ContactQuadFaceL4Instance.C \
	../Acme.d/search/ContactQuadFaceQ8.C \
	../Acme.d/search/ContactQuadFaceQ9.C \
	../Acme.d/search/ContactRangeSearch.C \
	../Acme.d/search/ContactRestart.C \
	../Acme.d/search/ContactScratchManager.C \
	../Acme.d/search/ContactSearch.C \
	../Acme.d/search/ContactSearchData.C \
	../Acme.d/search/ContactSecondaryDecomposition.C \
	../Acme.d/search/ContactSequentialAllocator.C \
	../Acme.d/search/ContactShellHandler.C \
	../Acme.d/search/ContactShellNode.C \
	../Acme.d/search/ContactShellQuadFaceL4Instance.C \
	../Acme.d/search/ContactShellTriFaceL3Instance.C \
	../Acme.d/search/contact_sorting.f \
	../Acme.d/search/ContactSurfaceGeometry.C \
	../Acme.d/search/ContactSymComm.C \
	../Acme.d/search/ContactTable.C \
	../Acme.d/search/ContactTimer.C \
	../Acme.d/search/ContactTopLevelSearch.C \
	../Acme.d/search/ContactTopology.C \
	../Acme.d/search/ContactTopologyEntityInstance.C \
	../Acme.d/search/ContactTopologyEntityHash.C \
	../Acme.d/search/ContactTopologyEntityList.C \
	../Acme.d/search/ContactTrackedSearch.C \
	../Acme.d/search/ContactTracking.C \
	../Acme.d/search/ContactTriFaceL3Instance.C \
	../Acme.d/search/ContactTriFaceQ6.C \
	../Acme.d/search/ContactUpdateTopology.C \
	../Acme.d/search/ContactVariableRestart.C \
	../Acme.d/search/ContactWedgeElementL6Instance.C \
	../Acme.d/search/contact_zapotec.c \
	../Acme.d/search/ContactZoltan.C \
	../Acme.d/search/ContactZoltanComm.C \
	../Acme.d/search/ContactZoltanCommUtils.C \
	../Acme.d/search/ContactZoltanID.C \
	../Acme.d/search/CString.C \
	../Acme.d/search/Search_Interface.C \
	../Acme.d/search/search_methods.f \
	../Acme.d/search/Zoltan_Interface.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := CommLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Comm.d/CommunicatorCore.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := CorotationalLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Corotational.d/Corotator.C \
	../Corotational.d/BeamCorotator.C \
	../Corotational.d/inc_rottensor.C \
	../Corotational.d/Shell3Corotator.C \
	../Corotational.d/BarCorotator.C \
	../Corotational.d/BarThermalCorotator.C \
	../Corotational.d/rightmult_rotvar.C \
	../Corotational.d/mat_to_vec.C \
	../Corotational.d/pseudorot_2var.C \
	../Corotational.d/GeomState.C \
	../Corotational.d/leftmult_rotvar.C \
	../Corotational.d/pseudorot_var.C \
	../Corotational.d/vec_to_quat.C \
	../Corotational.d/form_rottensor.C \
	../Corotational.d/mat_to_quat.C \
	../Corotational.d/fitalg3_3nodnew.C \
	../Corotational.d/vec_to_mat.C \
	../Corotational.d/quat_to_mat.C \
	../Corotational.d/crossprod.C \
	../Corotational.d/tran_fsl.C \
	../Corotational.d/mat_mult_mat.C \
	../Corotational.d/GeomNLSolver.C \
	../Corotational.d/orthonorm3.C \
	../Corotational.d/DistrGeomState.C \
	../Corotational.d/normalize.C \
	../Corotational.d/tran_force.C \
	../Corotational.d/tran_stiff.C \
	../Corotational.d/tran_rvec.C \
	../Corotational.d/tran_veloc.C \
	../Corotational.d/tran_accel.C \
	../Corotational.d/SpringCorotator.C \
	../Corotational.d/TetCorotator.C \
	../Corotational.d/BrickCorotator.C \
	../Corotational.d/SuperCorotator.C \
	../Corotational.d/PentaCorotator.C \
	../Corotational.d/BarFCorotator.C \
	../Corotational.d/QuadThermalCorotator.C \
	../Corotational.d/TriangleThermalCorotator.C \
	../Corotational.d/PhantomCorotator.C \
	../Corotational.d/MatNLCorotator.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := DecLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Dec.d/Decomp.d/Greedy.C \
	../Dec.d/Decomp.d/DecOpti.C \
	../Dec.d/Geom.d/IPoints.C \
	../Dec.d/Geom.d/Mesh.C \
	../Dec.d/Geom.d/PolyObj.C \
	../Dec.d/dec.C \
	../Dec.d/ElemMFCheck.C \
	../Dec.d/AddE.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := DriverLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Driver.d/DecDomainCore.C \
	../Driver.d/GeoSourceCore.C \
	../Driver.d/HData.C \
	../Driver.d/Mpc.C \
	../Driver.d/BaseSub.C \
	../Driver.d/MultiFront.C \
	../Driver.d/Dynam.C \
	../Driver.d/TempDynam.C \
	../Driver.d/NLDynam.C \
	../Driver.d/NLStatic.C \
	../Driver.d/Static.C \
	../Driver.d/SubDomainCore.C \
	../Driver.d/Domain.C \
	../Driver.d/SComm.C \
	../Driver.d/Eigen.C \
	../Driver.d/MPIComm.C \
	../Driver.d/SubDomainFactory.C \
	../Driver.d/PolygonSet.C \
	../Driver.d/Sower.C \
	../Driver.d/CornerMaker.C \
	../Driver.d/BinaryOutput.C \
	../Driver.d/DomainPita.C \
	../Driver.d/jacobi.f \
	../Driver.d/pade.f \
	../Driver.d/zpade.f

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := ElementLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Element.d/AddElem.C \
	../Element.d/Coordset.C \
	../Element.d/Element.C \
	../Element.d/Elemset.C \
	../Element.d/State.C \
	../Element.d/SuperElement.C \
	../Element.d/Beam.d/beamRotation.f \
	../Element.d/Beam.d/e3dmas.f \
	../Element.d/Beam.d/EulerBeam.C \
	../Element.d/Beam.d/frame6.f \
	../Element.d/Beam.d/g_mass7.f \
	../Element.d/Beam.d/g_modmstif7.f \
	../Element.d/Beam.d/g_sands7.f \
	../Element.d/Beam.d/mass6.f \
	../Element.d/Beam.d/mass7.f \
	../Element.d/Beam.d/modmstif6.f \
	../Element.d/Beam.d/modmstif7.f \
	../Element.d/Beam.d/sands6.f \
	../Element.d/Beam.d/sands7.f \
	../Element.d/Beam.d/TimoshenkoBeam.C \
	../Element.d/Beam.d/transf.f \
	../Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.C \
	../Element.d/BelytschkoTsayShell.d/fint3dshl_hypo.f90 \
	../Element.d/BelytschkoTsayShell.d/base_mech_transf.f90 \
	../Element.d/BelytschkoTsayShell.d/base_mhd_fem.f90 \
	../Element.d/BelytschkoTsayShell.d/base.f90 \
	../Element.d/BelytschkoTsayShell.d/base_mech_tens.f90 \
	../Element.d/BelytschkoTsayShell.d/hypoelas.f90 \
	../Element.d/BelytschkoTsayShell.d/damage.f90 \
	../Element.d/BelytschkoTsayShell.d/base_mhd_hgc.f90 \
	../Element.d/BelytschkoTsayShell.d/fint3dshl_hgc.f90 \
	../Element.d/BelytschkoTsayShell.d/massl3dshl.f90 \
	../Element.d/BelytschkoTsayShell.d/parseinp.f90 \
	../Element.d/BelytschkoTsayShell.d/cri_effstr.f90 \
	../Element.d/BelytschkoTsayShell.d/cri_mpts.f90 \
	../Element.d/BelytschkoTsayShell.d/cri_xw_j2pstrs.f90 \
	../Element.d/BelytschkoTsayShell.d/bc2_3d.f90 \
	../Element.d/BelytschkoTsayShell.d/base_mhd_fem_cpp.cpp \
	../Element.d/Brick20.d/br20mas.f \
	../Element.d/Brick20.d/br20vmint.f \
	../Element.d/Brick20.d/Brick20.C \
	../Element.d/Brick20.d/brik20v.f \
	../Element.d/Brick20.d/g_br20mas.f \
	../Element.d/Brick20.d/g_br20vmint.f \
	../Element.d/Brick20.d/g_brik20v.f \
	../Element.d/Brick20.d/g_h20shpe.f \
	../Element.d/Brick20.d/g_vol20.f \
	../Element.d/Brick20.d/h20shpe.f \
	../Element.d/Brick20.d/hexa20.f \
	../Element.d/Brick20.d/Hexa20ShapeFct.C \
	../Element.d/Brick20.d/hxgaus20.f \
	../Element.d/Brick20.d/sands20.f \
	../Element.d/Brick20.d/vol20.f \
	../Element.d/Brick32.d/Brick32.C \
	../Element.d/Brick32.d/brk32.f \
	../Element.d/Brick32.d/Hexa32ShapeFct.C \
	../Element.d/Brick.d/br8mas.f \
	../Element.d/Brick.d/brik8v.f \
	../Element.d/Brick.d/brkcmt.f \
	../Element.d/Brick.d/EightNodeBrick.C \
	../Element.d/Brick.d/h8shpe.f \
	../Element.d/Brick.d/Hexa8ShapeFct.C \
	../Element.d/Brick.d/hxgaus.f \
	../Element.d/Brick.d/sands17c.f \
	../Element.d/Brick.d/sands17.f \
	../Element.d/Brick.d/thermbr8mas.f \
	../Element.d/Brick.d/ThermBrick.C \
	../Element.d/Brick.d/thermbrik8v.f \
	../Element.d/Brick.d/vol17.f \
	../Element.d/BulkFluid.d/PentaBulk.C \
	../Element.d/BulkFluid.d/TetraBulk.C \
	../Element.d/BulkFluid.d/TriangleBulk.C \
	../Element.d/CompShell.d/compatt.f \
	../Element.d/CompShell.d/compbBB.f \
	../Element.d/CompShell.d/compbBM.f \
	../Element.d/CompShell.d/compbMB.f \
	../Element.d/CompShell.d/compbMM.f \
	../Element.d/CompShell.d/compchk2.f \
	../Element.d/CompShell.d/compchk.f \
	../Element.d/CompShell.d/compcrd2.f \
	../Element.d/CompShell.d/compcrd.f \
	../Element.d/CompShell.d/compcst1.f \
	../Element.d/CompShell.d/compcst.f \
	../Element.d/CompShell.d/compfrot.f \
	../Element.d/CompShell.d/comphBB.f \
	../Element.d/CompShell.d/comphBM.f \
	../Element.d/CompShell.d/comphMB.f \
	../Element.d/CompShell.d/comphMM.f \
	../Element.d/CompShell.d/compjac.f \
	../Element.d/CompShell.d/complay.f \
	../Element.d/CompShell.d/compmrot.f \
	../Element.d/CompShell.d/compms.f \
	../Element.d/CompShell.d/compthmfr.f \
	../Element.d/CompShell.d/Compo3NodeShell.C \
	../Element.d/CompShell.d/Compo4NodeShell.C \
	../Element.d/CompShell.d/compst.f \
	../Element.d/CompShell.d/compvms.f \
	../Element.d/CompShell.d/g_compbBB.f \
	../Element.d/CompShell.d/g_compbBM.f \
	../Element.d/CompShell.d/g_compbMB.f \
	../Element.d/CompShell.d/g_compbMM.f \
	../Element.d/CompShell.d/g_compcrd2.f \
	../Element.d/CompShell.d/g_compcrd.f \
	../Element.d/CompShell.d/g_compcst.f \
	../Element.d/CompShell.d/g_compfrot.f \
	../Element.d/CompShell.d/g_comphBB.f \
	../Element.d/CompShell.d/g_comphBM.f \
	../Element.d/CompShell.d/g_comphMB.f \
	../Element.d/CompShell.d/g_comphMM.f \
	../Element.d/CompShell.d/g_complay.f \
	../Element.d/CompShell.d/g_compms1.f \
	../Element.d/CompShell.d/g_compms2.f \
	../Element.d/CompShell.d/g_compst.f \
	../Element.d/CompShell.d/g_compvms.f \
	../Element.d/CompShell.d/g_straineq.f \
	../Element.d/CompShell.d/g_transform.f \
	../Element.d/CompShell.d/g_trirotation.f \
	../Element.d/ContactResistance.d/BrickContact.C \
	../Element.d/ContactResistance.d/QuadContact.C \
	../Element.d/ContactResistance.d/PentaContact.C \
	../Element.d/Convection.d/BarConvec.C \
	../Element.d/Convection.d/convecquad.f \
	../Element.d/Convection.d/QuadConvec.C \
	../Element.d/Convection.d/trianarea.f \
	../Element.d/Convection.d/TriangleConvec.C \
	../Element.d/CtcVirtualElt.d/CtcVirtualElt.C \
	../Element.d/FelippaShell.d/EffMembraneTriangleInstance.cpp \
	../Element.d/FelippaShell.d/AndesBendingTriangleInstance.cpp \
	../Element.d/FelippaShell.d/ShellElementInstance.cpp \
	../Element.d/FelippaShell.d/FelippaShell.C \
	../Element.d/FelippaShell.d/FelippaShellX2.C \
	../Element.d/FelippaShell.d/ShellMaterial.cpp \
	../Element.d/FelippaShell.d/ShellMaterialType0.cpp \
	../Element.d/FelippaShell.d/ShellMaterialType1.cpp \
	../Element.d/FelippaShell.d/ShellMaterialTypes2And3.cpp \
	../Element.d/FelippaShell.d/ShellMaterialType4.cpp \
	../Element.d/FluidQuad.d/BarSloshFS.C \
	../Element.d/FluidTetra.d/HEVibTetra.C \
	../Element.d/FluidQuad.d/HEVibQuadGal.C \
	../Element.d/FluidTetra.d/SloshTetra.C \
	../Element.d/FluidQuad.d/SloshQuadGal.C \
	../Element.d/FluidTriangle3.d/SloshTriangleFS.C \
	../Element.d/FluidQuad.d/barsloshfs.f \
	../Element.d/FluidQuad.d/slsas2.f \
	../Element.d/FluidQuad.d/thermquad3b.f \
	../Element.d/Helm.d/GaussRule.C \
	../Element.d/Helm.d/GaussRules.C \
	../Element.d/Helm.d/genTetra10.C \
	../Element.d/Helm.d/helmbr8mas.f \
	../Element.d/Helm.d/HelmBrick32.C \
	../Element.d/Helm.d/HelmBrick.C \
	../Element.d/Helm.d/HelmBrickGLS.C \
	../Element.d/Helm.d/helmbrik8v.f \
	../Element.d/Helm.d/HelmElement.C \
	../Element.d/Helm.d/HelmIsoParamHexa.C \
	../Element.d/Helm.d/ThermIsoParamHexa.C \
	../Element.d/Helm.d/HelmIsoParamQuad.C \
	../Element.d/Helm.d/HelmIsoParamTetra.C \
	../Element.d/Helm.d/HelmIsoParamTri.C \
	../Element.d/Helm.d/HelmLagQuadGal.C \
	../Element.d/Helm.d/HelmPenta26.C \
	../Element.d/Helm.d/HelmPenta.C \
	../Element.d/Helm.d/HelmQuad8Gal.C \
	../Element.d/Helm.d/HelmQuadGal.C \
	../Element.d/Helm.d/HelmQuadGls.C \
	../Element.d/Helm.d/HelmSpectralIsoParamHexa.C \
	../Element.d/Helm.d/HelmSpectralIsoParamQuad.C \
	../Element.d/Helm.d/HelmTri3Gal.C \
	../Element.d/Helm.d/HelmTri3Gls.C \
	../Element.d/Helm.d/HelmTri6Gal.C \
	../Element.d/Helm.d/IsoParamUtils2d.C \
	../Element.d/Helm.d/IsoParamUtils.C \
	../Element.d/Helm.d/LEIsoParamHexa.C \
	../Element.d/Helm.d/LEIsoParamQuad.C \
	../Element.d/Helm.d/LEIsoParamTetra.C \
	../Element.d/Helm.d/LEIsoParamTri.C \
	../Element.d/Helm.d/linegaussq.f \
	../Element.d/Helm.d/lub3.C \
	../Element.d/Helm.d/lud3.C \
	../Element.d/Helm.d/q4d1dofmas.f \
	../Element.d/Helm.d/Q4shape2.C \
	../Element.d/Helm.d/quad1dofm.f \
	../Element.d/Helm.d/quad8mass1.f \
	../Element.d/Helm.d/quad8shape.f \
	../Element.d/Helm.d/quad8stif1.f \
	../Element.d/Helm.d/quadgaussq.f \
	../Element.d/Helm.d/Tetra10HelmGal.C \
	../Element.d/Helm.d/TetraHelmGal.C \
	../Element.d/Helm.d/TetraHelmGLS.C \
	../Element.d/Helm.d/trig6mass1.f \
	../Element.d/Helm.d/trig6shape.f \
	../Element.d/Helm.d/trig6stif1.f \
	../Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.C \
	../Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.C \
	../Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.C \
	../Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.C \
	../Element.d/Joint.d/BuildingBlocks.d/StraightLinePointFollowerConstraint.C \
	../Element.d/Joint.d/CylindricalJoint.C \
	../Element.d/Joint.d/PrismaticJoint.C \
	../Element.d/Joint.d/RevoluteJoint.C \
	../Element.d/Joint.d/WeldedJoint.C \
	../Element.d/Joint.d/SphericalJoint.C \
	../Element.d/Joint.d/TranslationalJoint.C \
	../Element.d/Joint.d/UniversalJoint.C \
	../Element.d/Joint.d/DotConstraintType1a.C \
	../Element.d/Joint.d/RevoluteActuator.C \
	../Element.d/Joint.d/PinInSlotJoint.C \
	../Element.d/Joint.d/PlanarJoint.C \
	../Element.d/Joint.d/LinearTranslationalSpring.C \
	../Element.d/Joint.d/NonlinearTorsionalSpring.C \
	../Element.d/Joint.d/NonlinearTranslationalSpring.C \
	../Element.d/Joint.d/NonlinearTranslationalSpring2.C \
	../Element.d/Joint.d/SphericalJointSpringCombo.C \
	../Element.d/Joint.d/TranslationalJointSpringCombo.C \
	../Element.d/Joint.d/UniversalJointSpringCombo.C \
	../Element.d/Joint.d/RevoluteJointSpringCombo.C \
	../Element.d/Joint.d/CylindricalJointSpringCombo.C \
	../Element.d/Joint.d/PrismaticJointSpringCombo.C \
	../Element.d/Joint.d/PinInSlotJointSpringCombo.C \
	../Element.d/Membrane.d/mass8.f \
	../Element.d/Membrane.d/Membrane.C \
	../Element.d/Membrane.d/sands19.f \
	../Element.d/Membrane.d/trimem.f \
	../Element.d/Membrane.d/trirotationx.f \
	../Element.d/Membrane.d/FourNodeMembrane.C \
	../Element.d/MpcElement.d/FsiElement.C \
	../Element.d/MpcElement.d/MpcElement.C \
	../Element.d/MpcElement.d/DistanceConstraintElement.C \
	../Element.d/MpcElement.d/DotType1ConstraintElement.C \
	../Element.d/MpcElement.d/DotType2ConstraintElement.C \
	../Element.d/MpcElement.d/DotType2v2ConstraintElement.C \
	../Element.d/MpcElement.d/AngleType1ConstraintElement.C \
	../Element.d/MpcElement.d/LineLineDistanceConstraintElement.C \
	../Element.d/MpcElement.d/LineVariLineDistanceConstraintElement.C \
	../Element.d/MpcElement.d/SegmentSegmentDistanceConstraintElement.C \
	../Element.d/MpcElement.d/SegVariSegDistanceConstraintElement.C \
	../Element.d/MpcElement.d/PointPointDistanceConstraintElement.C \
	../Element.d/MpcElement.d/PointLineDistanceConstraintElement.C \
	../Element.d/MpcElement.d/PointPlaneDistanceConstraintElement.C \
	../Element.d/MpcElement.d/PointVariPointDistanceConstraintElement.C \
	../Element.d/MpcElement.d/PointVariLineDistanceConstraintElement.C \
	../Element.d/MpcElement.d/PointVariPlaneDistanceConstraintElement.C \
        ../Element.d/MpcElement.d/PointVariPlaneSegmentDistanceConstraintElement.C \
	../Element.d/Function.d/Constraint.d/DistanceConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/DotType1ConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/DotType2ConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/DotType2v2ConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/AngleType1ConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/PointPointDistanceConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/PointLineDistanceConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/PointPlaneDistanceConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/PointVariPointDistanceConstraintFunction.C \
	../Element.d/Function.d/Constraint.d/exp-map.C \
	../Element.d/NonLinearity.d/2DMat.C \
	../Element.d/NonLinearity.d/ShapeFunction.C \
	../Element.d/NonLinearity.d/ElaLinIsoMat.C \
	../Element.d/NonLinearity.d/GaussIntgElem.C \
	../Element.d/NonLinearity.d/NLHexahedral.C \
	../Element.d/NonLinearity.d/NLPentahedral.C \
	../Element.d/NonLinearity.d/NLTetrahedral.C \
	../Element.d/NonLinearity.d/NLMembrane.C \
	../Element.d/NonLinearity.d/StrainEvaluator.C \
	../Element.d/Penta15.d/Penta15.C \
	../Element.d/Penta15.d/Penta15ShapeFct.C \
	../Element.d/Penta26.d/Penta26.C \
	../Element.d/Penta26.d/Penta26ShapeFct.C \
	../Element.d/Penta.d/eg3c2c.f \
	../Element.d/Penta.d/grav24.f \
	../Element.d/Penta.d/mass24.f \
	../Element.d/Penta.d/mstf24.f \
	../Element.d/Penta.d/Penta6ShapeFct.C \
	../Element.d/Penta.d/Pentahedral.C \
	../Element.d/Penta.d/sands24.f \
	../Element.d/Quad4.d/FourNodeQuad.C \
	../Element.d/Quad4.d/getcmt.f \
	../Element.d/Quad4.d/lgauss.f \
	../Element.d/Quad4.d/q4dmas.f \
	../Element.d/Quad4.d/q4shpe.f \
	../Element.d/Quad4.d/qgauss.f \
	../Element.d/Quad4.d/quad2d.f \
	../Element.d/Quad4.d/quad4m.f \
	../Element.d/Quad4.d/Quad.C \
	../Element.d/Quad4.d/sands2.f \
	../Element.d/Quad4.d/strainvmc.f \
	../Element.d/Quad4.d/strainvm.f \
	../Element.d/Radiation.d/BarRadiation.C \
	../Element.d/Radiation.d/TriangleRadiation.C \
	../Element.d/Radiation.d/QuadRadiation.C \
	../Element.d/Rigid.d/RigidBeam.C \
	../Element.d/Rigid.d/RigidEightNodeBrick.C \
	../Element.d/Rigid.d/RigidRotnSprlink.C \
	../Element.d/Rigid.d/RigidSolid6Dof.C \
	../Element.d/Rigid.d/RigidSolid.C \
	../Element.d/Rigid.d/RigidSpring.C \
	../Element.d/Rigid.d/RigidThreeNodeShell.C \
	../Element.d/Rigid.d/RigidTransSprlink.C \
	../Element.d/Rigid.d/RigidTwoNodeTruss.C \
	../Element.d/Rigid.d/RigidFourNodeShell.C \
	../Element.d/Rigid.d/FlexibleBeam.C \
	../Element.d/Rigid.d/FlexibleTwoNodeTruss.C \
	../Element.d/Shear.d/b3dstf.f \
	../Element.d/Shear.d/shearmass.f \
	../Element.d/Shear.d/ShearPanel.C \
	../Element.d/Shear.d/shearpanel.f \
	../Element.d/Shear.d/shearpanelF.f \
	../Element.d/Shear.d/spstress.f \
	../Element.d/Shear.d/vec.f \
	../Element.d/Shell.d/basico.f \
	../Element.d/Shell.d/ConnectedTri.C \
	../Element.d/Shell.d/FourNodeShell.C \
	../Element.d/Shell.d/ExpFourNodeShell.C \
	../Element.d/Shell.d/g_basico.f \
	../Element.d/Shell.d/g_dummy.f \
	../Element.d/Shell.d/g_mass8.f \
	../Element.d/Shell.d/g_membra.f \
	../Element.d/Shell.d/g_momen.f \
	../Element.d/Shell.d/g_rotation.f \
	../Element.d/Shell.d/g_sands8.f \
	../Element.d/Shell.d/g_sm3mb.f \
	../Element.d/Shell.d/g_sm3mhe.f \
	../Element.d/Shell.d/g_smcbh.f \
	../Element.d/Shell.d/g_straineq.f \
	../Element.d/Shell.d/g_transform.f \
	../Element.d/Shell.d/g_tria3d.f \
	../Element.d/Shell.d/g_trirotation.f \
	../Element.d/Shell.d/g_vonmis.f \
	../Element.d/Shell.d/localXY.f \
	../Element.d/Shell.d/mass8.f \
	../Element.d/Shell.d/membra.f \
	../Element.d/Shell.d/momen.f \
	../Element.d/Shell.d/rotati.f \
	../Element.d/Shell.d/rotation.f \
	../Element.d/Shell.d/sands8.f \
	../Element.d/Shell.d/sm3mb.f \
	../Element.d/Shell.d/sm3mhe.f \
	../Element.d/Shell.d/smcbh.f \
	../Element.d/Shell.d/straineq.f \
	../Element.d/Shell.d/Therm3NoShell.C \
	../Element.d/Shell.d/Therm4NoShell.C \
	../Element.d/Shell.d/ThreeNodeShell.C \
	../Element.d/Shell.d/transform.f \
	../Element.d/Shell.d/tria3d.f \
	../Element.d/Shell.d/trirotation.f \
	../Element.d/Shell.d/trirot.f \
	../Element.d/Shell.d/trithmfr.f \
	../Element.d/Shell.d/vonmis.f \
	../Element.d/Sommerfeld.d/conv.C \
	../Element.d/Sommerfeld.d/CurvedLine2SommerBC.C \
	../Element.d/Sommerfeld.d/IsoParamLineSommer.C \
	../Element.d/Sommerfeld.d/IsoParamQuadSommer.C \
	../Element.d/Sommerfeld.d/IsoParamTriLineSommer.C \
	../Element.d/Sommerfeld.d/IsoParamTriSommer.C \
	../Element.d/Sommerfeld.d/LagLineSommer.C \
	../Element.d/Sommerfeld.d/Line2SommerBC.C \
	../Element.d/Sommerfeld.d/LineSommerBC.C \
	../Element.d/Sommerfeld.d/q4d1sommas.f \
	../Element.d/Sommerfeld.d/QuadSommerBC.C \
	../Element.d/Sommerfeld.d/SommerElement.C \
	../Element.d/Sommerfeld.d/SpectralIsoParamQuadSommer.C \
	../Element.d/Sommerfeld.d/Tri6Gauss.C \
	../Element.d/Sommerfeld.d/Triangle6SommerBC.C \
	../Element.d/Sommerfeld.d/TriangleSommerBC.C \
	../Element.d/Sommerfeld.d/QuadPressureBC.C \
	../Element.d/Sommerfeld.d/TrianglePressureBC.C \
	../Element.d/Sommerfeld.d/Quad8PressureBC.C \
        ../Element.d/Sommerfeld.d/Quad9PressureBC.C \
	../Element.d/Sommerfeld.d/Triangle6PressureBC.C \
	../Element.d/Sommerfeld.d/Quad12PressureBC.C \
	../Element.d/Sommerfeld.d/Triangle10PressureBC.C \
	../Element.d/Spring.d/mstf11.f \
	../Element.d/Spring.d/mstf21.f \
	../Element.d/Spring.d/mstf22.f \
	../Element.d/Spring.d/outprd.f \
	../Element.d/Spring.d/RotnSprlink.C \
	../Element.d/Spring.d/TorSpring.C \
	../Element.d/Spring.d/TransSprlink.C \
	../Element.d/Tetra10.d/dpe25.C \
	../Element.d/Tetra10.d/mass25.f \
	../Element.d/Tetra10.d/mstf25.f \
	../Element.d/Tetra10.d/p2thcoupl.C \
	../Element.d/Tetra10.d/p2thcoupl_dat.C \
	../Element.d/Tetra10.d/sands25.f \
	../Element.d/Tetra10.d/TenNodeTetrahedral.C \
	../Element.d/Tetra10.d/Tetra10ShapeFct.C \
	../Element.d/Tetra10.d/vp1.C \
	../Element.d/Tetra.d/AddBtCBtoK3DSolid.C \
	../Element.d/Tetra.d/ec3c2c.f \
	../Element.d/Tetra.d/em3c2c.f \
	../Element.d/Tetra.d/er3c2c.f \
	../Element.d/Tetra.d/fobase.f \
	../Element.d/Tetra.d/mass23.f \
	../Element.d/Tetra.d/mstf23.f \
	../Element.d/Tetra.d/p1thcoupl.C \
	../Element.d/Tetra.d/p1thcoupl_dat.C \
	../Element.d/Tetra.d/sands23.f \
	../Element.d/Tetra.d/Tetra4ShapeFct.C \
	../Element.d/Tetra.d/Tetrahedral.C \
	../Element.d/Tetra.d/vmelmvc.f \
	../Element.d/Tetra.d/vmelmv.f \
	../Element.d/Tetra.d/ThermIsoParamTetra.C \
	../Element.d/ThermQuad.d/htsas2.f \
	../Element.d/ThermQuad.d/q4maslumpheat.f \
	../Element.d/ThermQuad.d/Therm3DQuad.C \
	../Element.d/ThermQuad.d/thermquad3a.f \
	../Element.d/ThermQuad.d/thermquad3b.f \
	../Element.d/ThermQuad.d/ThermQuadGal.C \
	../Element.d/ThermQuad.d/veccrs.f \
	../Element.d/Triangle3.d/ThermTriangle.C \
	../Element.d/Triangle3.d/Triangle3.C \
	../Element.d/Truss.d/Therm2NodeBar.C \
	../Element.d/Truss.d/TwoNodeTruss.C \
	../Element.d/Truss.d/TwoNodeTrussF.C \
	../Element.d/Utils.d/AddBtCBtoK3DSolid.C \
	../Element.d/Utils.d/CheckJacobian.C \
	../Element.d/Utils.d/RotateConstitutiveMatrix.C \
	../Element.d/Utils.d/StressAndStrain3DSolid.C \
	../Element.d/DEM.d/DEMElement.C \
	../Element.d/DEM.d/DEMHelm2d.C \
	../Element.d/DEM.d/DEMHelm3d.C \
	../Element.d/DEM.d/DEMLE2d.C \
	../Element.d/DEM.d/DEMLE3d.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := FetiLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Feti.d/CoarseSetCore.C \
	../Feti.d/DistrVectorCore.C \
	../Feti.d/FetiCore.C \
	../Feti.d/FetiOpCore.C \
	../Feti.d/NLFetiCore.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := GNU-getoptLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../GNU-getopt.d/my_getopt.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := HeteroLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Hetero.d/FlExchange.C \
	../Hetero.d/FilteredFile.C \
	../Hetero.d/DistFlExchange.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := LinpackLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Linpack.d/dsifa.f \
	../Linpack.d/dsisl.f \
	../Linpack.d/dsvdc.f \
	../Linpack.d/dpofa.f \
	../Linpack.d/dtrsl.f
#       ../Linpack.d/drot.f \
        ../Linpack.d/drotg.f

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := MaterialLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Material.d/Material.cpp \
	../Material.d/LinearElasticMaterial.cpp \
	../Material.d/NeoHookean.cpp \
	../Material.d/MooneyRivlin.cpp \
	../Material.d/IsotropicLinearElasticJ2PlasticMaterial.cpp \
	../Material.d/IsotropicLinearElasticJ2PlasticPlaneStressMaterial.cpp \
	../Material.d/KorkolisKyriakidesPlaneStressMaterial.cpp \
	../Material.d/KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding.cpp \
	../Material.d/KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding2.cpp

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := MathLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Math.d/spsmvp.f \
	../Math.d/sptmv.f \
	../Math.d/cspsmvp.f \
	../Math.d/cdspsmvp.f \
	../Math.d/IntFullM.C \
	../Math.d/SparseData.C \
	../Math.d/mathUtility.C \
	../Math.d/matrixFull.C \
	../Math.d/BigMatrixCore.C \
	../Math.d/BLKSparseMatrixCore.C \
	../Math.d/CuCSparseCore.C \
	../Math.d/DBSparseMatrixCore.C \
	../Math.d/FullRectMatrix.C \
	../Math.d/NBSparseMatrixCore.C \
	../Math.d/TTensor.C \
	../Math.d/EiSparseMatrixCore.C \
	../Math.d/Skyline.d/forbackr7.f \
	../Math.d/Skyline.d/BlockSkyCore.C \
	../Math.d/Skyline.d/forbackr7ns.f \
	../Math.d/Skyline.d/forbackr8.f \
	../Math.d/Skyline.d/forbackr8ns.f \
	../Math.d/Skyline.d/ffulpvm.f \
	../Math.d/Skyline.d/forbackc1.f \
	../Math.d/Skyline.d/forbackc2.f \
	../Math.d/Skyline.d/forbackc3.f \
	../Math.d/Skyline.d/forbackc4.f \
	../Math.d/Skyline.d/forbackr1.f \
	../Math.d/Skyline.d/forbackr1ns.f \
	../Math.d/Skyline.d/SkyData.C \
	../Math.d/Skyline.d/forbackr2.f \
	../Math.d/Skyline.d/forbackr2ns.f \
	../Math.d/Skyline.d/SkyMatrixCore.C \
	../Math.d/Skyline.d/forbackr3.f \
	../Math.d/Skyline.d/skymul.f \
	../Math.d/Skyline.d/forbackr3ns.f \
	../Math.d/Skyline.d/slacol.f \
	../Math.d/Skyline.d/forbackr4.f \
	../Math.d/Skyline.d/svbu4cb.f \
	../Math.d/Skyline.d/forbackr4ns.f \
	../Math.d/Skyline.d/svbu4gmb.f \
	../Math.d/Skyline.d/forbackr5.f \
	../Math.d/Skyline.d/svbu4rb.f \
	../Math.d/Skyline.d/forbackr5ns.f \
	../Math.d/Skyline.d/svbu4rb.save.f \
	../Math.d/Skyline.d/forbackr6.f \
	../Math.d/Skyline.d/utility.C \
	../Math.d/Skyline.d/forbackr6ns.f \
	../Math.d/ComplexFactorSparse.d/zassmb.f \
	../Math.d/ComplexFactorSparse.d/zblkslv.f \
	../Math.d/ComplexFactorSparse.d/zblkl2.f \
	../Math.d/ComplexFactorSparse.d/zblkslvp.f \
	../Math.d/ComplexFactorSparse.d/zblkldl.f \
	../Math.d/ComplexFactorSparse.d/zgecp.f \
	../Math.d/ComplexFactorSparse.d/zblkns.f \
	../Math.d/ComplexFactorSparse.d/zgers.f \
	../Math.d/ComplexFactorSparse.d/zblkslv2.f \
	../Math.d/ComplexFactorSparse.d/zmmpyi.f \
	../Math.d/ComplexFactorSparse.d/zblkslv3.f \
	../Math.d/ComplexFactorSparse.d/zpotf2_ldl.f \
	../Math.d/ComplexFactorSparse.d/zblkslv4.f \
	../Math.d/ComplexFactorSparse.d/zzeromat.f \
	../Math.d/FactorSparse.d/addMat.f \
	../Math.d/FactorSparse.d/epost2.f \
	../Math.d/FactorSparse.d/invinv.f \
	../Math.d/FactorSparse.d/addOne.f \
	../Math.d/FactorSparse.d/etordr.f \
	../Math.d/FactorSparse.d/ldindx.f \
	../Math.d/FactorSparse.d/assmb.f \
	../Math.d/FactorSparse.d/etpost.f \
	../Math.d/FactorSparse.d/lstats.f \
	../Math.d/FactorSparse.d/betree.f \
	../Math.d/FactorSparse.d/etree.f \
	../Math.d/FactorSparse.d/ltfrep.f \
	../Math.d/FactorSparse.d/bfinit.f \
	../Math.d/FactorSparse.d/fcnthn.f \
	../Math.d/FactorSparse.d/main.f \
	../Math.d/FactorSparse.d/blkl2.f \
	../Math.d/FactorSparse.d/fntsiz.f \
	../Math.d/FactorSparse.d/mmdelm.f \
	../Math.d/FactorSparse.d/blkldl.f \
	../Math.d/FactorSparse.d/fsup1.f \
	../Math.d/FactorSparse.d/mmdint.f \
	../Math.d/FactorSparse.d/blkns.f \
	../Math.d/FactorSparse.d/fsup2.f \
	../Math.d/FactorSparse.d/mmdnum.f \
	../Math.d/FactorSparse.d/genmmd.f \
	../Math.d/FactorSparse.d/mmdupd.f \
	../Math.d/FactorSparse.d/blkslv.f \
	../Math.d/FactorSparse.d/getadj.f \
	../Math.d/FactorSparse.d/mmpyi.f \
	../Math.d/FactorSparse.d/blkslvp.f \
	../Math.d/FactorSparse.d/getnrm.f \
	../Math.d/FactorSparse.d/ordmmd.f \
	../Math.d/FactorSparse.d/btree2.f \
	../Math.d/FactorSparse.d/getres.f \
	../Math.d/FactorSparse.d/ordnat.f \
	../Math.d/FactorSparse.d/chkns.f \
	../Math.d/FactorSparse.d/getrhs.f \
	../Math.d/FactorSparse.d/sfinit.f \
	../Math.d/FactorSparse.d/chordr.f \
	../Math.d/FactorSparse.d/gtimer.f \
	../Math.d/FactorSparse.d/symfc2.f \
	../Math.d/FactorSparse.d/dgecp.f \
	../Math.d/FactorSparse.d/igathr.f \
	../Math.d/FactorSparse.d/symfct.f \
	../Math.d/FactorSparse.d/dgers.f \
	../Math.d/FactorSparse.d/inpdof.f \
	../Math.d/FactorSparse.d/zero.f \
	../Math.d/FactorSparse.d/dpotf2_ldl.f \
	../Math.d/FactorSparse.d/inpmat.f \
	../Math.d/FactorSparse.d/dpotrf_ldl.f \
	../Math.d/FactorSparse.d/inpnv.f

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := MortarLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Mortar.d/FaceElement.d/FaceElemSet.C \
	../Mortar.d/FaceElement.d/AddFaceElemInFaceElemSet.C \
	../Mortar.d/FaceElement.d/SurfaceEntity.C \
	../Mortar.d/FaceElement.d/FaceElement.C \
	../Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.C \
	../Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.C \
	../Mortar.d/FaceElement.d/FaceQuad9.d/FaceQuad9.C \
	../Mortar.d/FaceElement.d/FaceQuad12.d/FaceQuad12.C \
	../Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.C \
	../Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.C \
	../Mortar.d/FaceElement.d/FaceTri10.d/FaceTri10.C \
	../Mortar.d/FaceElement.d/FacePoint1.d/FacePoint1.C \
	../Mortar.d/MortarDriver.d/MortarHandler.C \
	../Mortar.d/MortarElement.d/MortarElement.C \
	../Mortar.d/MortarElement.d/MortarQuad4.d/StdMortarQuad4.C \
	../Mortar.d/MortarElement.d/MortarQuad4.d/DualMortarQuad4.C \
	../Mortar.d/MortarElement.d/MortarQuad8.d/StdMortarQuad8.C \
	../Mortar.d/MortarElement.d/MortarQuad12.d/StdMortarQuad12.C \
	../Mortar.d/MortarElement.d/MortarTri3.d/StdMortarTri3.C \
	../Mortar.d/MortarElement.d/MortarTri3.d/DualMortarTri3.C \
	../Mortar.d/MortarElement.d/MortarTri6.d/StdMortarTri6.C \
	../Mortar.d/MortarElement.d/MortarTri10.d/StdMortarTri10.C \
	../Mortar.d/MortarElement.d/CreateMortarElement.C \
	../Mortar.d/Divers.d/Divers.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := ParalLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Paral.d/MDDynam.C \
	../Paral.d/MDNLStatic.C \
	../Paral.d/MDNLDynam.C 

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := ParserLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Parser.d/AuxDefs.C \
	../Parser.d/parser.cpp \
	../Parser.d/lexer.cpp

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := ProblemsLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Problems.d/EigenDescr.C \
	../Problems.d/DynamDescr.C \
	../Problems.d/NonLinStatic.C \
	../Problems.d/CondDescr.C \
	../Problems.d/NonLinDynam.C \
	../Problems.d/TempDescr.C \
	../Problems.d/ModalBase.C \
	../Problems.d/ModalGeomState.C \
	../Problems.d/DEMProblem.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := RomLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Rom.d/BasisBinaryFile.C \
	../Rom.d/BasisId.C \
	../Rom.d/BasisFileStream.C \
	../Rom.d/BasisInputFile.C \
	../Rom.d/BasisOutputFile.C \
	../Rom.d/BasisOrthoDriver.C \
	../Rom.d/CheckNonLinDynamic.C \
	../Rom.d/CholeskyUtils.C \
	../Rom.d/DistrBasisFile.C \
	../Rom.d/DistrBasisOrthoDriver.C \
	../Rom.d/DistrElementSamplingDriver.C \
	../Rom.d/DistrSnapshotNonLinDynamic.C \
	../Rom.d/DistrExplicitSnapshotNonLinDynamic.C \
	../Rom.d/DistrExplicitPodProjectionNonLinDynamicBase.C \
	../Rom.d/DistrExplicitPodProjectionNonLinDynamic.C \
	../Rom.d/DistrExplicitLumpedPodProjectionNonLinDynamic.C \
	../Rom.d/DistrROMPostProcessingDriver.C \
	../Rom.d/DistrNodeDof6Buffer.C \
	../Rom.d/DistrVecBasis.C \
	../Rom.d/DistrSvdOrthogonalization.C \
	../Rom.d/DistrVecNodeDof6Conversion.C \
	../Rom.d/ElementSamplingDriverInstance.C \
	../Rom.d/SubElementSamplingDriver.C \
	../Rom.d/FileNameInfo.C \
	../Rom.d/LumpedPodProjectionNonLinDynamic.C \
	../Rom.d/SparseNonNegativeLeastSquaresSolverInstance.C \
	../Rom.d/NodalRestrictionMapping.C \
	../Rom.d/MappedNodeDof6Buffer.C \
	../Rom.d/MeshOutput.C \
	../Rom.d/MeshDesc.C \
	../Rom.d/PodProjectionNonLinDynamic.C \
	../Rom.d/RenumberingUtils.C \
	../Rom.d/RestrictedVecNodeDof6Conversion.C \
	../Rom.d/SnapshotNonLinDynamic.C \
        ../Rom.d/SnapshotProjectionDriver.C \
	../Rom.d/SvdOrthogonalization.C \
	../Rom.d/VecBasis.C \
	../Rom.d/VecBasisFile.C \
	../Rom.d/VecBasisUtils.C \
	../Rom.d/VecNodeDof6Conversion.C \
	../Rom.d/VecNodeDof6Map.C \
	../Rom.d/DofSetUtils.C \
	../Rom.d/LawsonHanson.d/diff.cpp \
	../Rom.d/LawsonHanson.d/g1.cpp \
	../Rom.d/LawsonHanson.d/h12.cpp \
	../Rom.d/LawsonHanson.d/spnnls.cpp

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := SfemLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Sfem.d/Sfem.C \
	../Sfem.d/cijk.C \
	../Sfem.d/chaos.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := SolversLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Solvers.d/DSCsolver.C \
	../Solvers.d/MumpsCore.C \
	../Solvers.d/Rbm.C \
	../Solvers.d/MappedAssembledSolver.C \
	../Solvers.d/SolverCore.C \
	../Solvers.d/SpoolesCore.C \
	../Solvers.d/GoldfarbIdnani.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := ThreadsLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Threads.d/Paral.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := TimersLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Timers.d/Timing.C \
        ../Timers.d/StaticTimers.C \
        ../Timers.d/DistTimer.C \
        ../Timers.d/GetTime.C \
        ../Timers.d/NLTimers.C \
        ../Timers.d/FetiDPtimers.C

include $(BUILD_STATIC_LIBRARY)

include $(CLEAR_VARS)

LOCAL_MODULE     := UtilsLib
LOCAL_CPP_FEATURES := rtti exceptions
LOCAL_CPPFLAGS   := -D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS -DMAP_MIN_MEMORY -DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0 -DUSE_EIGEN3 -DUSE_EIGEN3_AUTODIFF -DEIGEN_QUATERNIONBASE_PLUGIN=\"Element.d/MpcElement.d/QuaternionBasePlugin.h\" -D__LINUX -fpermissive -fPIC
LOCAL_C_INCLUDES := $(LOCAL_PATH)/.. $(LOCAL_PATH)/../Acme.d/search $(LOCAL_PATH)/../Acme.d/enforcement $(LOCAL_PATH)/../../eigen
LOCAL_SRC_FILES  := ../Utils.d/Connectivity.C \
	../Utils.d/dofset.C \
	../Utils.d/BlockAlloc.C \
	../Utils.d/DistHelper.C \
	../Utils.d/MFTT.C \
	../Utils.d/BinaryOutputFile.C \
	../Utils.d/CompositeInfo.C \
	../Utils.d/pstress.C \
	../Utils.d/b2Rotation.C \
	../Utils.d/list.C \
	../Utils.d/BinaryResultFile.C \
	../Utils.d/NodeSpaceArray.C \
	../Utils.d/MathUtils.C \
	../Utils.d/dbg_alloca.C \
	../Utils.d/Conwep.d/BlastLoading.C

include $(BUILD_STATIC_LIBRARY)
