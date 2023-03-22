using Easy.Common.Extensions;
using Proteomics.ProteolyticDigestion;

namespace MS1_DataSimulator
{
    internal class Scan
    {
        public readonly List<PeptideWithSetModifications> peptides;
        public readonly List<PeptideSpectrum> peptideSpectra;

        public Scan(List<PeptideWithSetModifications> peptides, List<PeptideSpectrum> peptideSpectra)
        {
            this.peptides = peptides;
            this.peptideSpectra = peptideSpectra;
        }

        public (double[], double[]) Spectrum(double minimumIntensity = 0)
        {
            List<double> mzs = new();
            List<double> intensities = new();

            for (int i = 0; i < peptideSpectra.Count; i++)
            {
                for (int j = 0; j < peptideSpectra[i].intensityValues.Length; j++)
                {
                    if (peptideSpectra[i].intensityValues[j] > minimumIntensity)
                    {
                        mzs.Add(peptideSpectra[i].mzValues[j]);
                        intensities.Add(peptideSpectra[i].intensityValues[j]);
                    }
                }
            }
            intensities = intensities.SortLike(mzs.ToArray()).ToList();
            mzs.Sort();
            return (mzs.ToArray(), intensities.ToArray());
        }
        public void AddPeptideToScan(PeptideWithSetModifications peptide, PeptideSpectrum spectrum)
        {
            peptides.Add(peptide);
            peptideSpectra.Add(spectrum);
        }
        private static string GetLabel()
        {
            return string.Empty;
        }        
    }
}
