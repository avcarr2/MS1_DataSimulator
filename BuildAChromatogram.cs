using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MS1_DataSimulator
{
    internal class BuildAnMzmlChromatogram
    {
        List<string> FullPathsToFastas;
        double TrainingFraction;
        /// <summary>
        /// The ultimate goal here is to create an mzML file of simultated mass spec data
        /// Each scan will be a theoretical MS1 spectrum of intact peptides or proteoforms
        /// There will be from 1 to n different peptides/proteforms in each scan
        /// Each peptide/proteoform can have multiple charge states in a single scan
        /// The intensities or peptides and proteoforms will vary over up to four orders of magnitude
        /// Consecutive scans will have the same peptides/proteoforms coming in and going out as if it was a real chromatography peak
        /// Peptides/proteoforms will begin eluting at random times.
        /// Training and test data will appear as totally separate mzML files to ease analysis.
        /// </summary>
        /// <param name="fullPathsToFastas"></param>
        /// <param name="trainingFraction"></param>
        public BuildAnMzmlChromatogram(List<string> fullPathsToFastas, double trainingFraction)
        {
            this.FullPathsToFastas = fullPathsToFastas;
            this.TrainingFraction = trainingFraction;
            DigestProteinsIntoPeptides proteinDigest = new(fullPathsToFastas);
            SplitPeptides split = new((0, proteinDigest.Peptides.Length), trainingFraction);
            int maxNumScans = 100;
            Scan[] scans = new Scan[maxNumScans];

            Random rnd = new Random(42);

            int centerOfPeptide = 50;
            int[] chargeStates = new int[2] { 2, 3 };

            foreach (int trainingIndex in split.IndicesOfTrainingSet)
            {
                double peakArea = 10000 * rnd.NextDouble();
                GeneratePeak peak = new GeneratePeak(peakArea, centerOfPeptide);

                for (int i = 0; i < peak.RelativeScanPositions.Count; i++)
                {
                    if (peak.PeakHeights[i] > 0)
                    {
                        if(!((centerOfPeptide + peak.RelativeScanPositions[i]) > maxNumScans))
                        {
                            scans[centerOfPeptide + peak.RelativeScanPositions[i]].AddPeptideToScan(proteinDigest.Peptides[trainingIndex], peak.PeakHeights[i], chargeStates);
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                centerOfPeptide++;
            }
        }

        
    }
}
