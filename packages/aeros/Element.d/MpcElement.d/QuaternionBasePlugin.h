  /** Sets \c *this to be a quaternion corresponding to a rotation 
  * represented by a 3 DOF rotation vector \a v.
  *
  * \returns a reference to \c *this.
  */
  template<typename Derived1>
  Derived& setFromOneVector(const MatrixBase<Derived1>& v)
  {
    using std::sqrt;
    using std::cos;
    Vector3 e = v;

    Scalar thetaSquared = e.squaredNorm();
    Scalar theta;
    if(thetaSquared > Scalar(0)) theta = sqrt(thetaSquared); else theta = Scalar(0);
    double MIN_ANGLE = 1.0e-7; // crossover point to Taylor Series approximation
    if(theta > Scalar(MIN_ANGLE)) {
      e *= (1/theta);
      AngleAxisType aa(theta, e);
      *this = aa;
    }
    else { // Taylor series of 1/2*sinc(theta/2)
      e *= (Scalar(0.5) - thetaSquared/48 + thetaSquared*thetaSquared/3840);
      this->w() = cos(Scalar(0.5*theta));
      this->vec() = e;
    }

    //this->normalize();
    this->coeffs() *= (1/this->norm()); // workaround

    return derived();
  }
