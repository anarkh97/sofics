// BlastLoading.C
// This file implements the CONWEP blast loading model.
#include "BlastLoading.h"
#include <cmath>
#include <iostream>
#include <fstream>
// ====================================================================================================
// This is the main function.

void BlastLoading::Conwep::SetUnitConversionsAndMassCubeRoot(BlastLoading::BlastData& P)
{
  switch (P.UnitConversionId)
  {
    case 1: {
      P.ScaleLength = 1.0; //per foot
      P.ScaleMass = 1.0; //per pound-mass
      P.ScaleTime = 1000; //per s
      P.ScalePressure = 1.0; //per Psi
      break;
    }
    case 2: {
      P.ScaleLength = 3.28084; //per m
      P.ScaleMass = 2.20462; //per kg
      P.ScaleTime = 1000; //per s
      P.ScalePressure = 0.000145038; //per Pa
      break;
    }
    case 3: {
      P.ScaleLength = 0.0328084; //per cm
      P.ScaleMass = 0.00220462; //per gm
      P.ScaleTime = 0.001; //per microseconds
      P.ScalePressure = 14503773.800722; //per megabar
      break;
    }
    case 4: {
      P.ScaleLength = 0.00328084; //per mm
      P.ScaleMass = 2.20462; //per kg
      P.ScaleTime = 1.0; //per ms
      P.ScalePressure = 145038; //per GPa
      break;
    }
    case 5: {
      P.ScaleLength = 0.00328084; //per mm
      P.ScaleMass = 0.00220462; //per gm
      P.ScaleTime = 1000; //per s
      P.ScalePressure = 0.000145038; //per Pa
      break;
    }
    default: {
      std::cerr << " ***WARNING: Unknown units conversion defined under the CONWEP card in the\n"
                << "             input file. Setting units to Kg-m-s-Pa.\n";
      P.ScaleLength = 3.28084; //per m
      P.ScaleMass = 2.20462; //per kg
      P.ScaleTime = 1000; //per s
      P.ScalePressure = 0.000145038; //per Pa
      break;
    }
  }

  // Compute mass cube root
  P.ExplosiveWeightCubeRoot = std::pow(P.ExplosiveWeight*P.ScaleMass, 1.0/3.0);

}

double BlastLoading::Conwep::Blast(const BlastLoading::BlastData& P,
                                   const double CurrentElementFacePosition[3],
                                   const double CurrentElementFaceNormalDirection[3],
                                   double CurrentTime) {
  // Calculate the distance between the current element face and the explosive:
  double DirectionFromElementFaceToExplosive[3] = {
    P.ExplosivePosition[0]-CurrentElementFacePosition[0],
    P.ExplosivePosition[1]-CurrentElementFacePosition[1],
    P.ExplosivePosition[2]-CurrentElementFacePosition[2]
  };
  double DistanceFromElementFaceCentroidToExplosive = sqrt( DirectionFromElementFaceToExplosive[0]*DirectionFromElementFaceToExplosive[0]
                  +DirectionFromElementFaceToExplosive[1]*DirectionFromElementFaceToExplosive[1]
                  +DirectionFromElementFaceToExplosive[2]*DirectionFromElementFaceToExplosive[2] );
  if(DistanceFromElementFaceCentroidToExplosive == 0 && !WarnedZeroDist) {
    std::cerr << " *** WARNING: Conwep blast distance is identically zero\n";
    WarnedZeroDist = true;
  }
  // Normalize the distance between the current element face and the explosive:
  DirectionFromElementFaceToExplosive[0] /= DistanceFromElementFaceCentroidToExplosive;
  DirectionFromElementFaceToExplosive[1] /= DistanceFromElementFaceCentroidToExplosive;
  DirectionFromElementFaceToExplosive[2] /= DistanceFromElementFaceCentroidToExplosive;
  // Convert this distance to feet:
  DistanceFromElementFaceCentroidToExplosive = DistanceFromElementFaceCentroidToExplosive*P.ScaleLength;
  // Calculate the current element's position cosine:
  double CurrentElementPositionCosine =  CurrentElementFaceNormalDirection[0]*DirectionFromElementFaceToExplosive[0]
                                        +CurrentElementFaceNormalDirection[1]*DirectionFromElementFaceToExplosive[1]
                                        +CurrentElementFaceNormalDirection[2]*DirectionFromElementFaceToExplosive[2];
  // Declare the blast wave parameters:
  double IncidentWaveArrivalTime;
  double PositivePhaseDuration;
  double IncidentWaveImpulse;
  double ReflectedWaveImpulse;
  double IncidentWavePressure;
  double ReflectedWavePressure;
  double IncidentWaveDecayExponent;
  double ReflectedWaveDecayExponent;
  // Call the parameters function to calculate all the parameters:
  Conwep::Parameters(P,
                 DistanceFromElementFaceCentroidToExplosive,
                 IncidentWaveArrivalTime,
                 PositivePhaseDuration,
                 IncidentWaveImpulse,
                 ReflectedWaveImpulse,
                 IncidentWavePressure,
                 ReflectedWavePressure,
                 IncidentWaveDecayExponent,
                 ReflectedWaveDecayExponent);
  // Calculate the time since the detonation and convert it to milliseconds:
  double CurrentTimeSinceDetonationTime = (CurrentTime- P.ExplosiveDetonationTime)*P.ScaleTime;
  // Call the pressure function to calculate the pressure on the current element:
  double Pressure = Conwep::Pressure(CurrentTimeSinceDetonationTime,
                                     IncidentWaveArrivalTime,
                                     PositivePhaseDuration,
                                     IncidentWavePressure,
                                     ReflectedWavePressure,
                                     CurrentElementPositionCosine,
                                     IncidentWaveDecayExponent,
                                     ReflectedWaveDecayExponent);
  return Pressure;
}
// ====================================================================================================
// Calculate the decay rate constants:
double BlastLoading::Conwep::Decay(double CurrentPressure,
                                   double CurrentImpulse,
                                   double PositivePhaseDuration) {
  // Calculate the pressure-to-impulse ratio:
    // CurrentPressure refers to either IncidentWavePressure or ReflectedWavePressure, depending on which wave called the function.
    // CurrentImpulse refers to either IncidentWaveImpulse or ReflectedWaveImpulse, depending on which wave called the function.
  double PressureToImpulseRatio = CurrentPressure*PositivePhaseDuration/CurrentImpulse;
  // Initialize the decay exponent using the pressure-to-impulse ratio:
  double DecayExponent = PressureToImpulseRatio-1.0;
  // Declare the errors:
  double ErrorUpper,ErrorLower;
  // Calculate the decay exponent by iteratively minimizing the errors:
  int numiter = 0, maxiter = 100;
  do {
    ErrorUpper = DecayExponent*DecayExponent-PressureToImpulseRatio*(-1.0+DecayExponent+exp(-DecayExponent));
    ErrorLower =           2.0*DecayExponent-PressureToImpulseRatio*( 1.0              -exp(-DecayExponent));
    DecayExponent = DecayExponent-ErrorUpper/ErrorLower;
    if(++numiter >= maxiter) {
      if(!WarnedDecayExp) {
        std::cerr << " *** WARNING: Conwep decay exponent calculation did not converge in " 
                  << maxiter << " iterations, error = " << fabs(ErrorUpper) << std::endl;
        WarnedDecayExp = true;
      }
      break;
    }
  } while (fabs(ErrorUpper) > 1.0e-6);
  // Return either the IncidentWaveDecayExponent or the ReflectedWaveDecayExponent, depending on which wave called the function:
  return DecayExponent;
}
// ====================================================================================================
// Calculate the indicent blast wave's pressure:
double BlastLoading::Conwep::IncidentPressure(const BlastLoading::BlastData& P,
                                              double ScaledStandoffDistanceLog10) {
  const static double Coefficients[2][12] = {  1.9422502013,
                                      -1.6958988741,
                                      -0.154159376846,
                                       0.514060730593,
                                       0.0988534365274,
                                      -0.293912623038,
                                      -0.0268112345019,
                                       0.109097496421,
                                       0.00162846756311,
                                      -0.0214631030242,
                                       0.0001456723382,
                                       0.00167847752266,
                                       1.77284970457,
                                      -1.69012801396,
                                       0.00804973591951,
                                       0.336743114941,
                                      -0.00516226351334,
                                      -0.0809228619888,
                                      -0.00478507266747,
                                       0.00793030472242,
                                       0.0007684469735,
                                       0.0,
                                       0.0,
                                       0.0 };
  double DistanceFactor = -0.756579301809 + 1.35034249993*ScaledStandoffDistanceLog10;
  int BlastType = (int)P.BlastType;
  double IncidentWavePressure = Coefficients[BlastType][11];
  for (int Counter = 10; Counter >= 0; --Counter) {
    IncidentWavePressure = IncidentWavePressure*DistanceFactor+Coefficients[BlastType][Counter];
  }
  return pow(10.0, IncidentWavePressure);
}
// ====================================================================================================
// Calculate the reflected blast wave's pressure:
double BlastLoading::Conwep::ReflectedPressure(const BlastLoading::BlastData& P,
                                               double ScaledStandoffDistanceLog10) {
  const static double CoefficientsSurfaceBlast[12] = {
     2.56431321138,
    -2.21030870597,
    -0.218536586295,
     0.895319589372,
     0.24989009775,
    -0.569249436807,
    -0.11791682383,
     0.224131161411,
     0.0245620259375,
    -0.0455116002694,
    -0.00190930738887,
     0.00361471193389 };
  const static double CoefficientsAirBlast[10] = {
     2.39106134946,
    -2.21400538997,
     0.035119031446,
     0.657599992109,
     0.0141818951887,
    -0.243076636231,
    -0.0158699803158,
     0.0492741184234,
     0.00227639644004,
    -0.00397126276058 };
  double DistanceFactor;
  double ReflectedWavePressure;
  int Counter;
  switch (P.BlastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    DistanceFactor = -0.789312405513+1.36637719229*ScaledStandoffDistanceLog10;
    ReflectedWavePressure = CoefficientsSurfaceBlast[11];
    for (Counter = 10; Counter >= 0; --Counter)
      ReflectedWavePressure = ReflectedWavePressure*DistanceFactor+CoefficientsSurfaceBlast[Counter];
    break;
  case BlastLoading::BlastData::AirBurst:
    DistanceFactor = -0.756579301809 + 1.35034249993*ScaledStandoffDistanceLog10;
    ReflectedWavePressure = CoefficientsAirBlast[9];
    for (Counter = 8; Counter >= 0; --Counter)
      ReflectedWavePressure = ReflectedWavePressure*DistanceFactor+CoefficientsAirBlast[Counter];
    break;
  }
  return pow(10.0, ReflectedWavePressure);
}
// ====================================================================================================
// Calculate the arrival time of the incident blast wave:
double BlastLoading::Conwep::ArrivalTime(const BlastLoading::BlastData& P,
                                         double ScaledStandoffDistanceLog10) {
  const static double CoefficientsSurfaceBlast[10] = { -0.173607601251,
                                     1.35706496258,
                                     0.052492798645,
                                    -0.196563954086,
                                    -0.0601770052288,
                                     0.0696360270891,
                                     0.0215297490092,
                                    -0.0161658930785,
                                    -0.00232531970294,
                                     0.00147752067524 };
  const static double CoefficientsAirBlast[8] = { -0.0423733936826,
                                    1.36456871214, 
                                   -0.0570035692784,
                                   -0.182832224796,
                                    0.0118851436014,
                                    0.0432648687627,
                                   -0.0007997367834,
                                   -0.00436073555033 };
  double DistanceFactor;
  double IncidentWaveArrivalTime;
  int Counter;
  switch (P.BlastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    DistanceFactor = -0.755684472698 + 1.37784223635*ScaledStandoffDistanceLog10;
    IncidentWaveArrivalTime = CoefficientsSurfaceBlast[9];
    for (Counter = 8; Counter >= 0; --Counter)
      IncidentWaveArrivalTime = IncidentWaveArrivalTime*DistanceFactor+CoefficientsSurfaceBlast[Counter];
    break;
  case BlastLoading::BlastData::AirBurst:
    DistanceFactor = -0.80501734056 + 1.37407043777*ScaledStandoffDistanceLog10;
    IncidentWaveArrivalTime = CoefficientsAirBlast[7];
    for (Counter = 6; Counter >= 0; --Counter)
      IncidentWaveArrivalTime = IncidentWaveArrivalTime*DistanceFactor+CoefficientsAirBlast[Counter];
    break;  
  }
  return pow(10.0, IncidentWaveArrivalTime);
}
// ====================================================================================================
// Calculate the time duration of the positive pressure after the incident blast wave hits:
double BlastLoading::Conwep::PositivePhaseDuration(const BlastLoading::BlastData& P,
                                                   double ScaledStandoffDistanceLog10) {
  const static double CoefficientsSurfaceBlastCase1[6] = { -0.728671776005,
                                                            0.130143717675,
                                                            0.134872511954,
                                                            0.0391574276906,
                                                           -0.00475933664702,
                                                           -0.00428144598008 };
  const static double CoefficientsSurfaceBlastCase2[9] = {  0.20096507334,
                                                           -0.0297944268976,
                                                            0.030632954288,
                                                            0.0183405574086,
                                                           -0.0173964666211,
                                                           -0.00106321963633,
                                                            0.00562060030977,
                                                            0.0001618217499,
                                                           -0.0006860188944 };
  const static double CoefficientsSurfaceBlastCase3[6] = {  0.572462469964,
                                                            0.0933035304009,
                                                           -0.0005849420883,
                                                           -0.00226884995013,
                                                           -0.00295908591505,
                                                            0.00148029868929 };
  const static double CoefficientsAirBlastCase1[9] = { -0.801052722864,
                                                        0.164953518069,
                                                        0.127788499497,
                                                        0.00291430135946,
                                                        0.00187957449227,
                                                        0.0173413962543,
                                                        0.00269739758043,
                                                       -0.00361976502798,
                                                       -0.00100926577934 };
  const static double CoefficientsAirBlastCase2[9] = {  0.115874238335,
                                                       -0.0297944268969,
                                                        0.0306329542941,
                                                        0.018340557407,
                                                       -0.0173964666286,
                                                       -0.00106321963576,
                                                        0.0056206003128,
                                                        0.0001618217499,
                                                       -0.0006860188944 };
  const static double CoefficientsAirBlastCase3[8] = {  0.50659210403,
                                                        0.0967031995552,
                                                       -0.00801302059667,
                                                        0.00482705779732,
                                                        0.00187587272287,
                                                       -0.00246738509321,
                                                       -0.000841116668,
                                                        0.0006193291052 };
  double DistanceFactor;
  double PositivePhaseDuration;
  int Counter;
  switch (P.BlastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    if (ScaledStandoffDistanceLog10 <= -0.34) {
      PositivePhaseDuration = -0.725;
    } else if (ScaledStandoffDistanceLog10 <= 0.4048337) {
      DistanceFactor = -0.1790217052 + 5.25099193925*ScaledStandoffDistanceLog10;
      PositivePhaseDuration = CoefficientsSurfaceBlastCase1[5];
      for (Counter = 4; Counter >= 0; --Counter)
        PositivePhaseDuration = PositivePhaseDuration*DistanceFactor+CoefficientsSurfaceBlastCase1[Counter];
    } else if (ScaledStandoffDistanceLog10 <= 0.845098) {
      DistanceFactor = -5.85909812338 + 9.2996288611*ScaledStandoffDistanceLog10;
      PositivePhaseDuration = CoefficientsSurfaceBlastCase2[8];
      for (Counter = 7; Counter >= 0; --Counter)
        PositivePhaseDuration = PositivePhaseDuration*DistanceFactor+CoefficientsSurfaceBlastCase2[Counter];
    } else  {
      DistanceFactor = -4.92699491141 + 3.46349745571*ScaledStandoffDistanceLog10;
      PositivePhaseDuration = CoefficientsSurfaceBlastCase3[5];
      for (Counter = 4; Counter >= 0; --Counter)
        PositivePhaseDuration = PositivePhaseDuration*DistanceFactor+CoefficientsSurfaceBlastCase3[Counter];
    }
    break;
  case BlastLoading::BlastData::AirBurst:  
    if (ScaledStandoffDistanceLog10 <= -0.34) {
      PositivePhaseDuration = -0.824;
    } else if (ScaledStandoffDistanceLog10 <= 0.350248) {
      DistanceFactor = 0.209440059933 + 5.11588554305*ScaledStandoffDistanceLog10;
      PositivePhaseDuration = CoefficientsAirBlastCase1[8];
      for (Counter = 7; Counter >= 0; --Counter)
        PositivePhaseDuration = PositivePhaseDuration*DistanceFactor+CoefficientsAirBlastCase1[Counter];
    } else if (ScaledStandoffDistanceLog10 <= 0.7596678) {
      DistanceFactor = -5.06778493835 + 9.2996288611*ScaledStandoffDistanceLog10;
      PositivePhaseDuration = CoefficientsAirBlastCase2[8];
      for (Counter = 7; Counter >= 0; --Counter)
        PositivePhaseDuration = PositivePhaseDuration*DistanceFactor+CoefficientsAirBlastCase2[Counter];
    } else  {
      DistanceFactor = -4.39590184126 + 3.1524725264*ScaledStandoffDistanceLog10;
      PositivePhaseDuration = CoefficientsAirBlastCase3[7];
      for (Counter = 6; Counter >= 0; --Counter)
        PositivePhaseDuration = PositivePhaseDuration*DistanceFactor+CoefficientsAirBlastCase3[Counter];
    }
    break;
  }
  return pow(10.0, PositivePhaseDuration);
}
// ====================================================================================================
// Calculate the impulse of the reflected blast wave:
double BlastLoading::Conwep::ReflectedImpulse(const BlastLoading::BlastData& P,
                                              double ScaledStandoffDistanceLog10) {
  const static double CoefficientsReflectedWaveImpulseSurfaceBlast[4] = { 1.75291677799,
                                                                         -0.949516092853,
                                                                          0.112136118689,
                                                                         -0.0250659183287 };
  const static double CoefficientsReflectedWaveImpulseAirBlast[4] = { 1.60579280091,
                                                                     -0.903118886091,
                                                                      0.101771877942,
                                                                     -0.0242139751146 };
  double AdjustmentFactor;
  double ReflectedWaveImpulse;
  int Counter;
  switch (P.BlastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    AdjustmentFactor = -0.781951689212 + 1.33422049854*ScaledStandoffDistanceLog10;
    ReflectedWaveImpulse = CoefficientsReflectedWaveImpulseSurfaceBlast[3];
    for (Counter = 2; Counter >= 0; --Counter)
      ReflectedWaveImpulse = ReflectedWaveImpulse*AdjustmentFactor+CoefficientsReflectedWaveImpulseSurfaceBlast[Counter];
    break;
  case BlastLoading::BlastData::AirBurst:
    AdjustmentFactor = -0.757659920369 + 1.37882996018*ScaledStandoffDistanceLog10;
    ReflectedWaveImpulse = CoefficientsReflectedWaveImpulseAirBlast[3];
    for (Counter = 2; Counter >= 0; --Counter)
      ReflectedWaveImpulse = ReflectedWaveImpulse*AdjustmentFactor+CoefficientsReflectedWaveImpulseAirBlast[Counter];
    break;
  }
  return pow(10.0, ReflectedWaveImpulse);
}
// ====================================================================================================
// Calculate the impulse of the incident blast wave:
double BlastLoading::Conwep::IncidentImpulse(const BlastLoading::BlastData& P,
                                             double ScaledStandoffDistanceLog10) {
  const static double CoefficientsIncidentWaveImpulseSurfaceBlast1[5] = { 1.57159240621,
                                   -0.502992763686,
                                    0.171335645235,
                                    0.0450176963051,
				   -0.0118964626402 };
  const static double CoefficientsIncidentWaveImpulseSurfaceBlast2[8] = { 0.719852655584,
                                   -0.384519026965,
                                   -0.0280816706301,
                                    0.00595798753822,
                                    0.014544526107,
                                   -0.00663289334734,
                                   -0.00284189327204,
                                    0.0013644816227 };
  const static double CoefficientsIncidentWaveImpulseAirBlast1[5] = { 1.43534136453,
                                   -0.443749377691,
                                    0.168825414684,
                                    0.0348138030308,
                                   -0.010435192824 };
  const static double CoefficientsIncidentWaveImpulseAirBlast2[9] = { 0.599008468099,
                                   -0.40463292088,
                                   -0.0142721946082,
                                    0.00912366316617,
                                   -0.0006750681404,
                                   -0.00800863718901,
                                    0.00314819515931,
                                    0.00152044783382,
                                   -0.0007470265899};
  double AdjustmentFactor;
  double IncidentWaveImpulse;
  int Counter;
  switch (P.BlastType) {
  case BlastLoading::BlastData::SurfaceBurst:
    if (ScaledStandoffDistanceLog10 <= 0.382017) {
      AdjustmentFactor = 0.832468843425 + 3.0760329666*ScaledStandoffDistanceLog10;
      IncidentWaveImpulse = CoefficientsIncidentWaveImpulseSurfaceBlast1[4];
      for (Counter = 3; Counter >= 0; --Counter)
        IncidentWaveImpulse = IncidentWaveImpulse*AdjustmentFactor+CoefficientsIncidentWaveImpulseSurfaceBlast1[Counter];
    } else  {
      AdjustmentFactor = -2.91358616806 + 2.40697745406*ScaledStandoffDistanceLog10;
      IncidentWaveImpulse = CoefficientsIncidentWaveImpulseSurfaceBlast2[7];
      for (Counter = 6; Counter >= 0; --Counter)
        IncidentWaveImpulse = IncidentWaveImpulse*AdjustmentFactor+CoefficientsIncidentWaveImpulseSurfaceBlast2[Counter];
    }
    break;
  case BlastLoading::BlastData::AirBurst:  
    if (ScaledStandoffDistanceLog10 <= 0.30103) {
      AdjustmentFactor = 1.04504877747 + 3.24299066475*ScaledStandoffDistanceLog10;
      IncidentWaveImpulse = CoefficientsIncidentWaveImpulseAirBlast1[4];
      for (Counter = 3; Counter >= 0; --Counter)
        IncidentWaveImpulse = IncidentWaveImpulse*AdjustmentFactor+CoefficientsIncidentWaveImpulseAirBlast1[Counter];
    } else  {
      AdjustmentFactor = -2.67912519532 + 2.30629231803*ScaledStandoffDistanceLog10;
      IncidentWaveImpulse = CoefficientsIncidentWaveImpulseAirBlast2[8];
      for (Counter = 7; Counter >= 0; --Counter)
        IncidentWaveImpulse = IncidentWaveImpulse*AdjustmentFactor+CoefficientsIncidentWaveImpulseAirBlast2[Counter];
    }
    break;
  }
  return pow(10.0, IncidentWaveImpulse);
}
// ====================================================================================================
// Calculate the Conwep parameters:
void BlastLoading::Conwep::Parameters(
  const BlastLoading::BlastData& P,
  double DistanceFromElementFaceCentroidToExplosive,
  double& IncidentWaveArrivalTime,
  double& PositivePhaseDuration,
  double& IncidentWaveImpulse,
  double& ReflectedWaveImpulse,
  double& IncidentWavePressure,
  double& ReflectedWavePressure,
  double& IncidentWaveDecayExponent,
  double& ReflectedWaveDecayExponent)
{
  double ScaledStandoffDistance = DistanceFromElementFaceCentroidToExplosive / P.ExplosiveWeightCubeRoot;
  double ScaledStandoffDistanceLog10 = log10(ScaledStandoffDistance);
  double LimitScaledStandoffDistance = (P.BlastType == BlastLoading::BlastData::SurfaceBurst ? 0.45 : 0.37);
  IncidentWaveArrivalTime = Conwep::ArrivalTime(P,ScaledStandoffDistanceLog10) * P.ExplosiveWeightCubeRoot;
  PositivePhaseDuration = Conwep::PositivePhaseDuration(P,ScaledStandoffDistanceLog10) * P.ExplosiveWeightCubeRoot;
  IncidentWaveImpulse = Conwep::IncidentImpulse(P, ScaledStandoffDistanceLog10) * P.ExplosiveWeightCubeRoot;
  ReflectedWaveImpulse = Conwep::ReflectedImpulse(P, ScaledStandoffDistanceLog10) * P.ExplosiveWeightCubeRoot;
  IncidentWavePressure = Conwep::IncidentPressure(P, ScaledStandoffDistanceLog10);
  ReflectedWavePressure = Conwep::ReflectedPressure(P, ScaledStandoffDistanceLog10);
  /*
  // For debugging:
  static int Counter = 0;
  static std::ofstream det("DetonationProperties.txt");
  if ((++Counter) < 122) {
    det << "DistanceFromElementFaceCentroidToExplosive = " << DistanceFromElementFaceCentroidToExplosive << "\n"
        << "IncidentWaveArrivalTime = " << IncidentWaveArrivalTime << "\n"
        << "PositivePhaseDuration = " << PositivePhaseDuration << "\n"
        << "IncidentWaveImpulse = " << IncidentWaveImpulse << "\n"
        << "ReflectedWaveImpulse = " << ReflectedWaveImpulse << "\n"
        << "IncidentWavePressure = " << IncidentWavePressure << "\n"
        << "ReflectedWavePressure = " << ReflectedWavePressure << std::endl;
  }
  */
  if (ScaledStandoffDistance >= ScaledStandoffDistance) {
    IncidentWaveDecayExponent  = Conwep::Decay(IncidentWavePressure, IncidentWaveImpulse, PositivePhaseDuration);
    ReflectedWaveDecayExponent = Conwep::Decay(ReflectedWavePressure, ReflectedWaveImpulse, PositivePhaseDuration); 
  } else {
    IncidentWaveDecayExponent = 0.0;
    ReflectedWaveDecayExponent = 0.0;
  }
}
// ====================================================================================================
// Calculate the pressure on the current element's face:
double BlastLoading::Conwep::Pressure(double CurrentTimeSinceExplosionTime,
                                      double IncidentWaveArrivalTime,
                                      double PositivePhaseDuration,
                                      double IncidentWavePressure,
                                      double ReflectedWavePressure,
                                      double CurrentElementPositionCosine,
                                      double IncidentWaveDecayExponent,
                                      double ReflectedWaveDecayExponent) {
  // Check if the current time exceeds the blast wave arrival time:
  if (CurrentTimeSinceExplosionTime >= IncidentWaveArrivalTime) {
    // If the blast wave has arrived, return the Conwep model's pressure:
    double IncidentWaveTime  = exp( -IncidentWaveDecayExponent*(CurrentTimeSinceExplosionTime-IncidentWaveArrivalTime)/PositivePhaseDuration);
    double ReflectedWaveTime = exp(-ReflectedWaveDecayExponent*(CurrentTimeSinceExplosionTime-IncidentWaveArrivalTime)/PositivePhaseDuration);
    // Check if the current element's position cosine is positive or negative:
    // If it is positive, then we are on the side of the element facing away from the explosive, so no pressure should be applied.
    double PositiveCurrentElementPositionCosine = (CurrentElementPositionCosine > 0.0 ? CurrentElementPositionCosine:0.0);
    // Compute the current element face's pressure:
    double Pressure = ( IncidentWavePressure*IncidentWaveTime*(1.0+PositiveCurrentElementPositionCosine
                       -2.0*PositiveCurrentElementPositionCosine*PositiveCurrentElementPositionCosine)
                       +ReflectedWavePressure*ReflectedWaveTime*PositiveCurrentElementPositionCosine*PositiveCurrentElementPositionCosine)
                       *(1.0-(CurrentTimeSinceExplosionTime-IncidentWaveArrivalTime)/PositivePhaseDuration);
    // Return the gage pressure:
    return (Pressure > -14.7 ? Pressure : -14.7);
  } else {
    // If the blast wave hasn't arrived yet, return zero pressure:
    return 0.0;
  } 
}
// ====================================================================================================
// Calculate the pressure at the centroid of the current 4 node quad (or degenerate) element:
double BlastLoading::ComputeShellPressureLoad(const double* CurrentElementNodePositions,
                                              double CurrentTime,
                                              const BlastLoading::BlastData& P) {
  // (AN) set unit convertors
  Conwep::SetUnitConversionsAndMassCubeRoot(const_cast<BlastLoading::BlastData&>(P));

  // Note: CurrentElementNodePositions contains the positions of the nodes of the current element.
  // CurrentElementNodePositions[0,1,2] = X,Y,Z position of node 1, CurrentElementNodePositions[3,4,5] = X,Y,Z position of node 2, etc.
  // Calculate 2 of the current element's edge directions:
  double CurrentElementEdge1Direction[3] = {
    CurrentElementNodePositions[6]-CurrentElementNodePositions[0],
    CurrentElementNodePositions[7]-CurrentElementNodePositions[1],
    CurrentElementNodePositions[8]-CurrentElementNodePositions[2]
  };
  double CurrentElementEdge2Direction[3] = {
    CurrentElementNodePositions[9] -CurrentElementNodePositions[3],
    CurrentElementNodePositions[10]-CurrentElementNodePositions[4],
    CurrentElementNodePositions[11]-CurrentElementNodePositions[5]
  };
  // Calculate the current element's normal vector by doing the cross-product of 2 of its edge directions:
  double CurrentElementNormalVector[3] = {
    CurrentElementEdge1Direction[1]*CurrentElementEdge2Direction[2]-CurrentElementEdge1Direction[2]*CurrentElementEdge2Direction[1],
    CurrentElementEdge1Direction[2]*CurrentElementEdge2Direction[0]-CurrentElementEdge1Direction[0]*CurrentElementEdge2Direction[2],
    CurrentElementEdge1Direction[0]*CurrentElementEdge2Direction[1]-CurrentElementEdge1Direction[1]*CurrentElementEdge2Direction[0]
  };
  // Calculate the magnitude of the current element's normal vector:
  double CurrentElementNormalVectorMagnitude =
    sqrt( CurrentElementNormalVector[0]*CurrentElementNormalVector[0]
         +CurrentElementNormalVector[1]*CurrentElementNormalVector[1]
         +CurrentElementNormalVector[2]*CurrentElementNormalVector[2]);
  // Normalize the current element's normal vector using the magnitude:
  CurrentElementNormalVector[0] /= CurrentElementNormalVectorMagnitude;
  CurrentElementNormalVector[1] /= CurrentElementNormalVectorMagnitude;
  CurrentElementNormalVector[2] /= CurrentElementNormalVectorMagnitude;
  // Check if the current element is degenerate (is a triangle instead of a quad):
  bool CurrentElementIsDegenerate = (   CurrentElementNodePositions[6]==CurrentElementNodePositions[9]
                                     && CurrentElementNodePositions[7]==CurrentElementNodePositions[10]
                                     && CurrentElementNodePositions[8]==CurrentElementNodePositions[11] );
  // Calculate the current element's number of nodes:
  int CurrentElementNumberOfNodes = CurrentElementIsDegenerate ? 3 : 4;
  // Calculate the current element's centroid coordinates:
  double CurrentElementCentroidCoordinates[3] = {0,0,0};
  for (int CurrentNodeOfCurrentElement = 0; CurrentNodeOfCurrentElement < CurrentElementNumberOfNodes; ++CurrentNodeOfCurrentElement) { 
    for (int CurrentDimension = 0; CurrentDimension < 3; ++CurrentDimension)
      CurrentElementCentroidCoordinates[CurrentDimension] += CurrentElementNodePositions[CurrentNodeOfCurrentElement*3+CurrentDimension];
  }
  for (int CurrentDimension = 0; CurrentDimension < 3; ++CurrentDimension) { 
    CurrentElementCentroidCoordinates[CurrentDimension] *= 1.0/CurrentElementNumberOfNodes;
  }
  // Calculate the pressure at the centroid of the current element:
  double CurrentElementPressure = Conwep::Blast(P,CurrentElementCentroidCoordinates,CurrentElementNormalVector,CurrentTime);
  // Return the pressure at the centroid of the current element:
  // Note that the pressure is in psi: convert to appropriate units
  return -CurrentElementPressure/P.ScalePressure;
}
// ====================================================================================================
// Calculate the pressure at a gauss point of the current element:
double BlastLoading::ComputeGaussPointPressure(const double CurrentElementGaussPointCoordinates[3],
                                              const double CurrentElementGaussPointNormalVector[3],
                                              double CurrentTime,
                                              const BlastLoading::BlastData& P) {
  // (AN) set unit convertors
  Conwep::SetUnitConversionsAndMassCubeRoot(const_cast<BlastLoading::BlastData&>(P));

  // Calculate the pressure at a gauss point of the current element:
  double CurrentElementPressure = Conwep::Blast(P,CurrentElementGaussPointCoordinates,CurrentElementGaussPointNormalVector,CurrentTime);
  // Return the pressure at a gauss point of the current element:
  // Note that the pressure is in psi: convert it to Pa (6.89e3 factor), then use ScaleLength,
  // ScaleTime and ScaleMass to convert it to the correct pressure units.
  return -CurrentElementPressure/P.ScalePressure;
}
// ====================================================================================================
// Print the BlastLoading::BlastData member variables to screen:
void BlastLoading::BlastData::print() {
  std::cerr << "ExplosivePosition = " << ExplosivePosition[0] << " " << ExplosivePosition[1] << " " << ExplosivePosition[2] << std::endl;
  std::cerr << "ExplosiveDetonationTime = " << ExplosiveDetonationTime << std::endl;
  if(BlastType == SurfaceBurst) std::cerr << "BlastType = SurfaceBurst\n"; else std::cerr << "BlastType = AirBurst\n";
  std::cerr << "ExplosiveWeight = " << ExplosiveWeight << std::endl;
  std::cerr << "ExplosiveWeightCubeRoot = " << ExplosiveWeightCubeRoot << std::endl;
  std::cerr << "ScaleLength = " << ScaleLength << std::endl;
  std::cerr << "ScaleTime = " << ScaleTime << std::endl;
  std::cerr << "ScaleMass = " << ScaleMass << std::endl;
  std::cerr << "ScalePressure = " << ScalePressure << std::endl;
}
// ====================================================================================================
// Initialize the BlastLoading::InputFileData structure and other static member variables:
BlastLoading::BlastData BlastLoading::InputFileData = {{0.0,0.0,0.0},0.0,
                                                BlastLoading::BlastData::AirBurst,1.0,2,0.3048,1.0,1.0,1.0,0.0};
bool BlastLoading::WarnedZeroDist = false;
bool BlastLoading::WarnedDecayExp = false;
// ====================================================================================================
// End of file.
