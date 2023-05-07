// See https://aka.ms/new-console-template for more information

using MS1_DataSimulator;

List<string> fullPathToFastas = new List<string>() { @"C:\Users\Austin\Desktop\uniprot-download_true_format_fasta_query__28human_29_20AND_20_28mode-2023.05.04-16.21.30.18.fasta" };
string fullPathToOutput = @"C:\Users\Austin\Documents\Projects\MsDataSimulatorOutput";
int maxNumScans = 1000;
double trainingFraction = 0.01;
double minimumIntensity = 1E-6;
double lowestObservedMz = 800.0;
double highestObservedMz = 2000.0;

BuildAnMzmlChromatogram chromatogram = new(fullPathToFastas, fullPathToOutput, maxNumScans, trainingFraction, minimumIntensity,
    lowestObservedMz, highestObservedMz);

Console.WriteLine("All done, buddy.");
Console.Read(); 