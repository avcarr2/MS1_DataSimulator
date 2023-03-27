// See https://aka.ms/new-console-template for more information

using MS1_DataSimulator;

List<string> fullPathToFastas = new List<string>() { @"C:\Users\mrsho\Documents\Projects\MSDataSimulator\yeast.fasta" };
string fullPathToOutput = @"C:\Users\mrsho\Documents\Projects\MSDataSimulator";
int maxNumScans = 3000;
double trainingFraction = 0.5;
double minimumIntensity = 0.1;
double lowestObservedMz = 200.0;
double highestObservedMz = 2200.0;

BuildAnMzmlChromatogram chromatogram = new(fullPathToFastas, fullPathToOutput, maxNumScans, trainingFraction, minimumIntensity, lowestObservedMz, highestObservedMz);

Console.WriteLine("Hello, World!");
