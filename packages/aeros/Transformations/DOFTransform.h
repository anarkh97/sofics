/*
 * DOFTransform.h
 *
 *  Created on: Mar 16, 2009
 *      Author: michel
 */

#ifndef DOFTRANSFORM_H_
#define DOFTRANSFORM_H_

struct DOFID {
  int node, num;
};

struct Coef {
  DOFID dof;
  double coef;
};

class DOFTransform {
    DOFID target;
    Vec<Coef> components;
  public:
    DOFTransform();
    virtual ~DOFTransform();
};

#endif /* DOFTRANSFORM_H_ */
