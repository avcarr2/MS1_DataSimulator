using Chemistry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace MS1_DataSimulator
{
    internal class DigestProteinsIntoPeptides
    {
        public static string? DataDir { get; private set; }
        public readonly PeptideWithSetModifications[] Peptides;

        public DigestProteinsIntoPeptides(List<string> fullPathsToFastaFiles, int minPeptideLength = 7, int maxMissedCleavages = 0)
        {
            
            Loaders.LoadElements();
            DataDir = SetUpDataDirectory();
            var protDic = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(DataDir, @"ProteolyticDigestion", @"proteases.tsv"));

            List<PeptideWithSetModifications> peptides = new();
            foreach (string path in fullPathsToFastaFiles)
            {
                if (File.Exists(path))
                {
                    List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out List<string> errors);

                    foreach (Protein protein in proteinList)
                    {
                        peptides.AddRange(protein.Digest(new DigestionParams(protease: "trypsin", minPeptideLength: minPeptideLength, maxMissedCleavages: maxMissedCleavages), new List<Modification>(), new List<Modification>()).ToList());
                    }
                } 
            }

            Peptides = peptides.ToArray();
            peptides.Clear();
        }
        private static string SetUpDataDirectory()
        {
            // get data directory
            var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles)
                && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins"))
            {
                return Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
            }
            else
            {
                return AppDomain.CurrentDomain.BaseDirectory;
            }
        }
    }
}
