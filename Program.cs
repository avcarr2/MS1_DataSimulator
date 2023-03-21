// See https://aka.ms/new-console-template for more information

using MS1_DataSimulator;

List<string> fullPathToFastas = new List<string>() { @"C:\Users\mrsho\Documents\Projects\MSDataSimulator\tenProtein.fasta" };
BuildAnMzmlChromatogram chromatogram = new (fullPathToFastas, 0.5);

Console.WriteLine("Hello, World!");
