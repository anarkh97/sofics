#ifndef _BEAMELEMENTTEMPLATE_HPP_
#define _BEAMELEMENTTEMPLATE_HPP_

template<typename doublereal>
class BeamElementTemplate {
public:
	BeamElementTemplate() {}

	~BeamElementTemplate() {}

	void gForce(doublereal *x, doublereal *y, doublereal *z, doublereal *gravityAcceleration,
	            doublereal *_eframe, doublereal rho, doublereal A, doublereal *gravityForce, int gravflg);

	void modmstif7(doublereal *_estif, doublereal A, doublereal E,
	               doublereal *_eframe, doublereal Ix, doublereal Iy, doublereal Iz, doublereal alphay,
	               doublereal alphaz, doublereal C1, doublereal nu, doublereal *_x, doublereal *_y, doublereal *_z,
	               int flag);

	void sands6(doublereal area, doublereal e, int elm,
	            doublereal *_stress, int maxsze, int maxgus, int maxstr,
	            doublereal *_eframe,
	            doublereal ix, doublereal iy, doublereal iz,
	            doublereal nu, doublereal *_x, doublereal *_y, doublereal *_z,
	            doublereal *_ug, doublereal alpha, doublereal tref, doublereal *_temp);

	void vms6WRTdisp(doublereal area, doublereal e, int elm,
	                 doublereal *_vmsWRTdisp, int maxsze, int maxgus, int maxstr,
	                 doublereal *_eframe,
	                 doublereal ix, doublereal iy, doublereal iz,
	                 doublereal nu, doublereal *_x, doublereal *_y, doublereal *_z,
	                 doublereal *_ug, doublereal alpha, doublereal tref, doublereal *_temp);

	void sands7(int elm, doublereal A, doublereal E,
	            doublereal *eframe, doublereal Ix, doublereal Iy, doublereal Iz,
	            doublereal alphay, doublereal alphaz, doublereal C1, doublereal nu,
	            doublereal *x, doublereal *y, doublereal *z, doublereal *ug,
	            doublereal *stress, int numel, int maxgus, int maxstr,
	            int msize, doublereal alpha, doublereal tref, doublereal *temp);

	void vms7WRTdisp(int elm, doublereal A, doublereal E,
	                 doublereal *eframe, doublereal Ix, doublereal Iy, doublereal Iz,
	                 doublereal alphay, doublereal alphaz, doublereal C1, doublereal nu,
	                 doublereal *x, doublereal *y, doublereal *z, doublereal *ug,
	                 doublereal *vmsWRTdisp,
	                 doublereal alpha, doublereal tref, doublereal *_temp);


	void transform(doublereal *_l, doublereal *_g, doublereal *_str);

	void buildFrameInTemplate(doublereal *, doublereal *, doublereal *, doublereal *);

	doublereal getMassInTemplate(doublereal *x, doublereal *y, doublereal *z, doublereal rho, doublereal A);

private:
	void frame6(doublereal *frame, doublereal *u, doublereal *v, doublereal *w);

	void zeroLengthError();

	void wrongGeometricPropertyError();

	void rotationMatrixError();

	void missingRotationMatrix();

	void error100();

	void error200();

	void error300();

	void error400();

};

#endif
