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

        public (double[], double[]) Spectrum()
        {
            List<double> mzs = new();
            List<double> intensities = new();

            for (int i = 0; i < peptideSpectra.Count; i++)
            {
                mzs.AddRange(peptideSpectra[i].mzValues);
                intensities.AddRange(peptideSpectra[i].intensityValues);
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
