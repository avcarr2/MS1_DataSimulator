// See https://aka.ms/new-console-template for more information

using MS1_DataSimulator;

List<string> fullPathToFastas = new List<string>() { @"C:\Users\mrsho\Documents\Projects\MSDataSimulator\tenProtein.fasta" };
string fullPathToOutput = @"C:\Users\mrsho\Documents\Projects\MSDataSimulator";
int maxNumScans = 6000;
double trainingFraction = 0.5;
double minimumIntensity = 0.1;

BuildAnMzmlChromatogram chromatogram = new(fullPathToFastas, fullPathToOutput, maxNumScans, trainingFraction, minimumIntensity);

Console.WriteLine("Hello, World!");
