using Proteomics.ProteolyticDigestion;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common;
using Easy.Common.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.Statistics;

namespace MS1_DataSimulator
{
    internal class PeptideSpectrum
    {
        public readonly PeptideWithSetModifications peptideWithSetModifications;
        public readonly ChargeStatesAndEnvelopeAbundances chargeStatesAndEnvelopeAbundances;
        public readonly double[] mzValues;
        public readonly double[] intensityValues;
        public double TotalSpectrumIntensity { get; private set; }
        public readonly List<ChargeStateIsotopeCluster> ChargeStateClusters;

        public PeptideSpectrum(PeptideWithSetModifications peptideWithSetModifications, ChargeStatesAndEnvelopeAbundances csea, double totalSpectrumIntensity = 1)
        {
            this.peptideWithSetModifications = peptideWithSetModifications;
            this.chargeStatesAndEnvelopeAbundances = csea;
            this.TotalSpectrumIntensity = totalSpectrumIntensity;
            this.ChargeStateClusters = PopulateClusters();
            var mzAndIntensityValues = PopulateSpectrum();
            this.mzValues = mzAndIntensityValues.Item1.ToArray();
            this.intensityValues = mzAndIntensityValues.Item2.ToArray();
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
        /// <summary>
        /// Returns an array of mzVals and the median charge state present. 
        /// </summary>
        /// <param name="lowMz"></param>
        /// <param name="highMz"></param>
        /// <param name="lowZ"></param>
        /// <param name="highZ"></param>
        /// <param name="mzVals"></param>
        /// <returns></returns>
        public int CreateChargeStates(double lowMz, double highMz, int lowZ, int highZ, out int[] charges)
        {
            List<double> mzVals = new List<double>();
            List<int> chargesPresent = new List<int>();
            for (int i = highZ; i >= lowMz; i--)
            {
                double mzVal = peptideWithSetModifications.MonoisotopicMass.ToMz(i);
                if (mzVal >= lowMz && mzVal <= highMz)
                {
                    mzVals.Add(mzVal); 
                    chargesPresent.Add(i);
                }
            }

            charges = chargesPresent.ToArray(); 
            return (int)charges.Select(i => (double)i).Median();
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
        private List<ChargeStateIsotopeCluster> PopulateClusters(int[] chargeStates)
        {
            List<ChargeStateIsotopeCluster> chargeStateClusters = new();

            foreach (int chargeState in chargeStates)
            {
                IsotopicMassesAndNormalizedAbundances unchargedParentCluster = new IsotopicMassesAndNormalizedAbundances(peptideWithSetModifications);
                (double[], double[]) spectrum = unchargedParentCluster.ComputeMzAndIntensity(chargeState);
                chargeStateClusters.Add(new ChargeStateIsotopeCluster(spectrum, chargeState));
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
        
        public string Label()
        {
            StringBuilder sb = new();

            //base sequence
            sb.Append(this.peptideWithSetModifications.BaseSequence + "\t");

            //full sequence
            sb.Append(this.peptideWithSetModifications.FullSequence + "\t");

            //monoisotopic mass
            sb.Append(this.peptideWithSetModifications.MonoisotopicMass + "\t");

            //number of displayed charge states
            sb.Append(String.Join(",",this.ChargeStateClusters.Select(c => c.chargeState).ToList()) + "\t");

            //list of spacing for each charge state
            sb.Append(String.Join(",", this.ChargeStateClusters.Select(c => c.mzSpacing).ToList()) + "\t");

            //list of most abundant mz value for each charge state
            sb.Append(String.Join(",", this.ChargeStateClusters.Select(c => c.mostAbundantMz).ToList()));

            return sb.ToString();
        }

        public string LabelHeader()
        {
            StringBuilder sb = new();

            //scan number
            sb.Append("Scan#" + "\t");

            //base sequence
            sb.Append("Base Sequence" + "\t");

            //full sequence
            sb.Append("Full Sequence" + "\t");

            //monoisotopic mass
            sb.Append("Monoisotopic Mass" + "\t");

            //number of displayed charge states
            sb.Append("Displayed Charge States" + "\t");

            //list of spacing for each charge state
            sb.Append("mz Spacing for each Charge State" + "\t");

            //list of most abundant mz value for each charge state
            sb.Append("Most Abundant mz for each Charge State");

            return sb.ToString();
        }
    }
}
