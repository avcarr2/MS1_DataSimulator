using Easy.Common.Extensions;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.Metrics;
using System.Linq;
using System.Text;
using System.Threading.Tasks; 

namespace MS1_DataSimulator
{
    internal class PeptideSpectrum
    {
        public readonly PeptideWithSetModifications peptideWithSetModifications;
        public readonly ChargeStatesAndEnvelopeAbundances chargeStatesAndEnvelopeAbundances;
        public readonly double[] mzValues;
        public readonly double[] intensityValues;
        public double TotalSpectrumIntensity { get; private set; }
        public readonly string peptideSpectrumLabel;
        public readonly List<ChargeStateIsotopeCluster> ChargeStateClusters;

        public PeptideSpectrum(PeptideWithSetModifications peptideWithSetModifications, int maxNumberChargeStates, int minChargeState, int maxChargeState, double minEnvelopeAbundance = 0.1, double totalSpectrumIntensity = 1)
        {
            this.peptideWithSetModifications = peptideWithSetModifications;
            this.chargeStatesAndEnvelopeAbundances = new ChargeStatesAndEnvelopeAbundances(maxNumberChargeStates, minChargeState, maxChargeState, minEnvelopeAbundance);
            this.TotalSpectrumIntensity = totalSpectrumIntensity;
            this.ChargeStateClusters = PopulateClusters();
            var mzAndIntensityValues = PopulateSpectrum();
            this.mzValues = mzAndIntensityValues.Item1.ToArray();
            this.intensityValues = mzAndIntensityValues.Item2.ToArray();
            this.peptideSpectrumLabel = GetLabel();
        }

        public void UpdateSpectrumWithNewTotalIntensity(double newTotalIntensity) 
        {
            double intensityMutiplier = newTotalIntensity / TotalSpectrumIntensity;

            for (int i = 0; i < intensityValues.Length; i++)
            {
                intensityValues[i] *= intensityMutiplier;
            }
            TotalSpectrumIntensity = intensityValues.Sum();
        }

        private List<ChargeStateIsotopeCluster> PopulateClusters() 
        {
            List<ChargeStateIsotopeCluster> chargeStateClusters = new();

            foreach (int chargeState in chargeStatesAndEnvelopeAbundances.ChargeStates)
            {
                IsotopicMassesAndNormalizedAbundances unchargedParentCluster = new IsotopicMassesAndNormalizedAbundances(peptideWithSetModifications);
                (double[], double[]) spectrum = unchargedParentCluster.ComputeMzAndIntensity(chargeState);
                chargeStateClusters.Add(new ChargeStateIsotopeCluster(spectrum,chargeState));
            }
            return chargeStateClusters;
        }
        private (List<double>,List<double>) PopulateSpectrum() 
        {
            List<double> mzs = new();
            List<double> intValues = new();

            for (int i = 0; i < ChargeStateClusters.Count; i++)
            {
                foreach (double mz in ChargeStateClusters[i].unitSpectrum.Item1)
                {
                    mzs.Add(mz);
                }
                foreach (double intensity in ChargeStateClusters[i].unitSpectrum.Item2)
                {
                    intValues.Add(intensity * chargeStatesAndEnvelopeAbundances.EnvelopeAbundances[i]);
                }
            }

            double intensitySum = intValues.Sum();
            double intensityMultipler = TotalSpectrumIntensity / intensitySum;

            for (int i = 0; i < intValues.Count; i++)
            {
                intValues[i] *= intensityMultipler;
            }

            return (mzs, intValues);
        }

        private string? GetLabel()
        {
            return string.Empty;
        }

    }
}
