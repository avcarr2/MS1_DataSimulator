// See https://aka.ms/new-console-template for more information

using MS1_DataSimulator;

List<string> fullPathToFastas = new List<string>() { @"C:\Users\mrsho\Documents\Projects\MSDataSimulator\tenProtein.fasta" };
string fullPathToOutput = @"C:\Users\mrsho\Documents\Projects\MSDataSimulator";
BuildAnMzmlChromatogram chromatogram = new(fullPathToFastas, fullPathToOutput, 0.5);

Console.WriteLine("Hello, World!");
