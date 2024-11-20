// If the BlastLoading class is not defined, define it:
#ifndef _BLASTLOADING_H_
#define _BLASTLOADING_H_

class BlastLoading {
    public:
        struct BlastData {
            double ExplosivePosition[3];
            double ExplosiveDetonationTime;
            enum {SurfaceBurst, AirBurst} BlastType;
            double ExplosiveWeight;

            /*
            (AN) Now user specifies an integer ID for unit conversions.
            Locally, we need units to be in pound-mass, foot, ms, psi.
            */

            int UnitConversionId;
            double ScaleLength; // units specified by user to foot
            double ScaleMass; // units specified by user to pound-mass
            double ScaleTime; // units specified by user to ms
            double ScalePressure; // units specified by user to psi
            double ExplosiveWeightCubeRoot;
            void print();
        };
    private:
        class Conwep {
            public:
                static void SetUnitConversionsAndMassCubeRoot(
                    BlastLoading::BlastData& P
                );
                static double Blast(
                    const BlastLoading::BlastData& P,
                    const double CurrentElementFaceCentroidPosition[3],
                    const double CurrentElementFaceNormalDirection[3],
                    double CurrentTime);
                static double Decay(
                    double CurrentPressure,
                    double CurrentImpulse,
                    double PositivePhaseDuration);
                static double IncidentPressure(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double ReflectedPressure(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double ArrivalTime(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10) ;
                static double PositivePhaseDuration(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double ReflectedImpulse(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10);
                static double IncidentImpulse(
                    const BlastLoading::BlastData& P,
                    double ScaledStandoffDistanceLog10) ;
                static void Parameters(
                    const BlastLoading::BlastData& P,
                    double DistanceFromElementFaceCentroidToExplosive,
                    double& IncidentWaveArrivalTime,
                    double& PositivePhaseDuration,
                    double& IncidentWaveImpulse,
                    double& ReflectedWaveImpulse,
                    double& IncidentWavePressure,
                    double& ReflectedWavePressure,
                    double& IncidentWaveDecayExponent,
                    double& ReflectedWaveDecayExponent);
                static double Pressure(
                    double CurrentTimeSinceExplosionTime,
                    double IncidentWaveArrivalTime,
                    double PositivePhaseDuration,
                    double IncidentWavePressure,
                    double ReflectedWavePressure,
                    double CurrentElementPositionCosine,
                    double IncidentWaveDecayExponent,
                    double ReflectedWaveDecayExponent);
        };
    public:
        static double ComputeShellPressureLoad(
            const double* CurrentElementNodePositions,
            double CurrentTime,
            const BlastLoading::BlastData& P);
        static double ComputeGaussPointPressure(const double CurrentElementGaussPointCoordinates[3],
            const double CurrentElementGaussPointNormalVector[3],
            double CurrentTime,
            const BlastLoading::BlastData& P);
        static BlastData InputFileData;
        static bool WarnedZeroDist;
        static bool WarnedDecayExp;

};
#endif
