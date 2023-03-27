using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;

namespace MS1_DataSimulator
{
    internal class BuildAnMzmlChromatogram
    {
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

        public BuildAnMzmlChromatogram(List<string> fullPathsToFastas, string fullPathToOutput, int maxNumScans, double trainingFraction, double minimumIntensity, double lowestObservedMz, double highestObservedMz)
        {
            DigestProteinsIntoPeptides proteinDigest = new(fullPathsToFastas);
            SplitPeptides split = new((0, proteinDigest.Peptides.Length), trainingFraction);
            
            Random rnd = new Random(42);

            

            List<List<int>> splits = new() { split.IndicesOfTrainingSet, split.IndicesOfTestSet };

            for (int j = 0; j < splits.Count; j++)
            {
                int centerOfPeptide = 50; 
                Scan[] scans = new Scan[maxNumScans];
                foreach (int trainingIndex in splits[j])
                {
                    double peakArea = 10000 * rnd.NextDouble() + 10.0; //having a minimum of 10 means that we'll have at least 1 charge state
                    List<int> peakShifts = Enumerable.Range(0, 7).ToList();
                    centerOfPeptide += peakShifts[rnd.Next(peakShifts.Count)];
                    GeneratePeak peak = new GeneratePeak(peakArea, centerOfPeptide);
                    int maxChargeStateCount = (int)Math.Round(Math.Log10(peakArea), 0);
                    (int, int) minMax = MinAndMaxChargeStatesForPeptide(proteinDigest.Peptides[trainingIndex], lowestObservedMz,highestObservedMz);
                    ChargeStatesAndEnvelopeAbundances csea = new ChargeStatesAndEnvelopeAbundances(maxChargeStateCount, minMax.Item1, minMax.Item2);
                    PeptideSpectrum genericSpectrum = new(proteinDigest.Peptides[trainingIndex], csea);

                    for (int i = 0; i < peak.RelativeScanPositions.Count; i++)
                    {
                        genericSpectrum.UpdateSpectrumWithNewTotalIntensity(peak.PeakHeights[i]);
                        int projectedScanNumber = centerOfPeptide + peak.RelativeScanPositions[i];
                        if (!(projectedScanNumber > (maxNumScans - 1)))
                        {
                            if (scans[centerOfPeptide + peak.RelativeScanPositions[i]] == null)
                            {
                                scans[centerOfPeptide + peak.RelativeScanPositions[i]] = new Scan(new List<Proteomics.ProteolyticDigestion.PeptideWithSetModifications>() { proteinDigest.Peptides[trainingIndex] }, new List<PeptideSpectrum>() { genericSpectrum });
                            }
                            else
                            {
                                scans[centerOfPeptide + peak.RelativeScanPositions[i]].AddPeptideToScan(proteinDigest.Peptides[trainingIndex], genericSpectrum);
                            }
                        }
                        else
                        {
                            break;
                        }

                    }
                    centerOfPeptide++;
                }

                if (j == 0)
                {
                    CreateMzMlAndLabelsFile(Path.Combine(fullPathToOutput, "train.mzML"), scans, minimumIntensity);
                }
                else
                {
                    CreateMzMlAndLabelsFile(Path.Combine(fullPathToOutput, "test.mzML"), scans, minimumIntensity);
                }      
            }
        }

        public static void CreateMzMlAndLabelsFile(string fullPathToOutput, Scan[] scans, double minimumIntensity)
        {
            int scanNumber = 0;
            List<MsDataScan> dataScans = new();
            List<string > labels = new();
            labels.Add(scans.Where(s => s != null).First().peptideSpectra[0].LabelHeader());
            for (int i = 0; i < scans.Length; i++)
            {
                if (!(scans[i] == null) && scans[i].Spectrum(minimumIntensity).Item1.Length > 0)
                {
                    scanNumber++;
                    MzSpectrum mzSpectrum = new(scans[i].Spectrum(minimumIntensity).Item1, scans[i].Spectrum(minimumIntensity).Item2, true);
                    int oneBasedScanNumber = scanNumber;
                    int msnOrder = 1;
                    bool isCentroid = true;
                    Polarity polarity = Polarity.Positive;
                    double retentionTime = (double)i;
                    MzRange scanWindowRange = new MzRange(scans[i].Spectrum(minimumIntensity).Item1.Min(), scans[i].Spectrum(minimumIntensity).Item1.Max());
                    string scanFilter = "";
                    MZAnalyzerType mzAnalyzer = MZAnalyzerType.Orbitrap;
                    double totalIonCurrent = scans[i].Spectrum(minimumIntensity).Item2.Sum();
                    double? injectionTime = null;
                    double[,]? noiseData = null;
                    string nativeId = "controllerType=0 controllerNumber=1 scan=" + oneBasedScanNumber.ToString();

                    dataScans.Add(new MsDataScan(mzSpectrum, oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent, injectionTime, noiseData, nativeId));

                    foreach (PeptideSpectrum spectrum in scans[i].peptideSpectra)
                    {
                        labels.Add(scanNumber + "\t" + spectrum.Label());
                    }
                }

            }
            SourceFile sourceFile = new("no nativeID format", "mzML format", null, null, null);

            MsDataFile dataFile = new(dataScans.ToArray(), sourceFile);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(dataFile, fullPathToOutput, false);

            File.WriteAllLines(fullPathToOutput.Replace("mzML","txt"), labels);
        }
        
        private (int,int) MinAndMaxChargeStatesForPeptide(PeptideWithSetModifications peptide, double minObservedMz, double maxObservedMz)
        {
            double approximateMinChargeStateProtonsMass = Math.Max(1.0, Math.Round(peptide.MonoisotopicMass / maxObservedMz - 0.5, 0)) * 1.003;
            double approximateMaxChargeStateProtonsMass = Math.Round(peptide.FullChemicalFormula.AverageMass / minObservedMz + 0.5, 0) * 1.003;

            int min = (int)Math.Max(1, Math.Round((peptide.MonoisotopicMass + approximateMinChargeStateProtonsMass) / maxObservedMz - 0.5, 0));
            int max = (int)Math.Round((peptide.FullChemicalFormula.AverageMass + approximateMaxChargeStateProtonsMass) / minObservedMz + 0.5,0);
            
            return (min, max);
        }
    }
}
