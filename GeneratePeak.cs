using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MS1_DataSimulator
{
    internal class GeneratePeak
    {
        public readonly double PeakArea;
        public readonly double PeakHeightAtCenter;
        public readonly double StandardDeviation;
        public readonly int CenterScanNumber;
        public readonly double MinimumFraction;
        public readonly List<double> PeakHeights = new();
        public readonly List<int> RelativeScanPositions = new();
        
        public GeneratePeak(double peakArea, int centerScanNumber = 0, double minimumFraction = 0.1)
        {
            this.MinimumFraction= minimumFraction;
            this.CenterScanNumber = centerScanNumber;
            this.PeakArea = peakArea;
            this.PeakHeightAtCenter = peakArea / 2.0;
            this.StandardDeviation = 0.3989 * this.PeakArea / this.PeakHeightAtCenter;

            //The parameter
            //  a: is the height of the curve's peak,
            //  b: is the position of the center of the peak, and
            //  c: (the standard deviation, sometimes called the Gaussian RMS width) controls the width of the "bell". 

            List<(int, double)> scanAndHeight = new();
            for (int i = -100; i <= 100; i++)
            {
                //this is the gaussia formula
                double x = Convert.ToDouble(i);
                double h = this.PeakHeightAtCenter * Math.Exp(-Math.Pow(x,2)/(2.0 * Math.Pow(this.StandardDeviation,2)));

                if (!double.IsNaN(h) && h > minimumFraction)
                {
                    scanAndHeight.Add((i, h));
                }
            }

            double peakSum = scanAndHeight.Select(s => s.Item2).Sum();
            foreach ((int,double) scan in scanAndHeight)
            {
                 RelativeScanPositions.Add(scan.Item1);
                 PeakHeights.Add(scan.Item2);  
            }
        }
    }
}
