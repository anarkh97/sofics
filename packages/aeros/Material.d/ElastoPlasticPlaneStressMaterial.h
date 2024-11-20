#ifndef ELASTOPLASTICPLANESTRESSMATERIAL
#define ELASTOPLASTICPLANESTRESSMATERIAL

#include <vector>
#include <limits>

class ElastoPlasticPlaneStressMaterial
{
 public:
  virtual ~ElastoPlasticPlaneStressMaterial() {};
  
  //! Compute the elastoplastic constitutive response.
  //! Returns true if calculations went well and false otherwise.
  //! \param Fnp1 Input. Deformation gradient at new state of material. Size 9x1.
  //! \param CauchyStress Output. Has size 9x1.
  //! \param Cep Output. Algorithmic elastoplastic tangent. If requested, has size 81x1.
  //! \param UpdateFlag Input. Material state updated if true. Note that by default, material state is updated.
  //! \param dt Input. Time increment.
  virtual bool ComputeElastoPlasticConstitutiveResponse(const std::vector<double> &Fnp1, 
                                                        std::vector<double> * CauchyStress, 
                                                        std::vector<double> * Cep = 0, 
                                                        const bool UpdateFlag = true,
                                                        const double dt = 0.) = 0;
  
  //! Returns the plastic strain in material (3x1 vector)
  virtual const std::vector<double> & GetMaterialPlasticStrain() const = 0;

  //! Returns equivalent plastic strain in material
  virtual double GetMaterialEquivalentPlasticStrain() const = 0;

  //! Returns back stress in material (3x1 vector)
  virtual const std::vector<double> & GetMaterialBackStress() const = 0;

  //! Returns dissipated energy in material
  virtual double GetDissipatedEnergy() const = 0;

  //! Returns the equivalent plastic strain at failure
  virtual double GetEquivalentPlasticStrainAtFailure() const { return std::numeric_limits<double>::max(); }

  //! Set the plastic strain in the material
  virtual void SetMaterialPlasticStrain(const std::vector<double> &EPSplastic) = 0;

  //! Set the equivalent plastic strain in the material
  virtual void SetMaterialEquivalentPlasticStrain(double equivEPSplastic) = 0;

  //! Set the back stress in the material
  virtual void SetMaterialBackStress(const std::vector<double> &BackStress) = 0;
};

#endif
