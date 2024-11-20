// ----------------------------------------------------------------
// HB - 03/02/04
//     Some compiler flags used in MortarHandler.[hC]
//     A priori ONLY for developpement purpose (i.e. should NOT be
//     accessed by the users)
// ----------------------------------------------------------------

// for using thread level parallelization in the computation of 
// the FFI M & N contributions 
// see MortarHandler::CreateFFIPolygons(...)
#define HB_THREAD_FFI_M_N

// for timing purpose
//#define HB_MORTAR_TIMER
//#define PJSA_ACME_TIMER

// for using thread level parallelization when assembling the nodal 
// mortar LMPC from the FFI M & N contributions
// see MortarHandler::CreateFFIPolygons(...)
// -> DO NOT SEEM TO BE EFFICIENT ...
//#define HB_THREAD_NODALMORTAR

