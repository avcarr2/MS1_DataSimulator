

using System.Collections.Generic;
using System;
using Chemistry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using Proteomics.ProteolyticDigestion;

namespace MS1_DataSimulator
{
    internal class ChargeStatesAndEnvelopeAbundances
    {
        public int[] ChargeStates;
        public double[] EnvelopeAbundances;

        public ChargeStatesAndEnvelopeAbundances(int[] chargeStates, double[] envelopeAbundances)
        {
            ChargeStates = chargeStates;
            EnvelopeAbundances = envelopeAbundances;
        }

        public void AssignEnvelopeAbundances(int[] chargeStates, int medianChargeState)
        {
            ChargeStates = chargeStates;
            Random seed = new();

            int chargeStateWidth = chargeStates[0] - chargeStates[^1]; 

            List<double> envelopeAbundances = new List<double>();
            for (double i = 0; i < chargeStates.Length; i++)
            {
                envelopeAbundances.Add(Normal.PDF((double)medianChargeState, seed.NextDouble() * chargeStateWidth + 2.5,
                    i)); 
            }

            EnvelopeAbundances = envelopeAbundances.ToArray(); 
        }

        public ChargeStatesAndEnvelopeAbundances()
        {

        }
    }

    public static class Extensions
    {
        /// <summary>
        /// Returns an array of mzVals and the median charge state present. 
        /// </summary>
        /// <param name="lowMz"></param>
        /// <param name="highMz"></param>
        /// <param name="lowZ"></param>
        /// <param name="highZ"></param>
        /// <param name="mzVals"></param>
        /// <returns></returns>
        public static int CreateChargeStates(this PeptideWithSetModifications peptide, double lowMz, double highMz, int lowZ, int highZ, out int[] charges)
        {
            List<double> mzVals = new List<double>();
            List<int> chargesPresent = new List<int>();
            for (int i = highZ; i >= lowZ; i--)
            {
                double mzVal = peptide.MonoisotopicMass.ToMz(i);
                if (mzVal >= lowMz && mzVal <= highMz)
                {
                    mzVals.Add(mzVal);
                    chargesPresent.Add(i);
                }
            }

            charges = chargesPresent.ToArray();
            return (int)charges.Select(i => (double)i).Median();
        }
    }
}
