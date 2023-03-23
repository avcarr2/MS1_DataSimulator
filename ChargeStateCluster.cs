using Easy.Common.Extensions;

namespace MS1_DataSimulator
{
    internal class ChargeStateIsotopeCluster
    {
        public readonly double firstMzValue;
        public readonly double mzSpacing;
        public readonly double mostAbundantMz;
        public readonly int peakCount;
        public readonly int chargeState;
        public readonly (double[], double[]) unitSpectrum;

        public ChargeStateIsotopeCluster((double[], double[]) spectrum, int charge)
        {
            this.firstMzValue = spectrum.Item1[0];
            this.chargeState = charge;
            if(spectrum.Item1.Length== 1 ) 
            {
                this.mzSpacing = 1.0 / (double)charge;
            } 
                else
            {
                this.mzSpacing = spectrum.Item1[1] - spectrum.Item1[0];
            }
            int maxIntensityIndex = spectrum.Item2.IndexOf(spectrum.Item2.Max());
            this.mostAbundantMz = spectrum.Item1[maxIntensityIndex];
            this.peakCount = spectrum.Item1.Length;
            this.unitSpectrum = spectrum;
        }
    }
}
