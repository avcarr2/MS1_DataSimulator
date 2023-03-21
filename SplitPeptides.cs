using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MS1_DataSimulator
{
    internal class SplitPeptides
    {
        public readonly List<int> IndicesOfTrainingSet = new();
        public readonly List<int> IndicesOfTestSet = new();

        public SplitPeptides((int, int) contiguousRange, double trainingFraction)
        {
            Random rnd = new(42);

            for (int i = contiguousRange.Item1; i < contiguousRange.Item2; i++)
            {
                if (rnd.NextDouble() < trainingFraction)
                {
                    IndicesOfTrainingSet.Add(i);
                }
                else
                {
                    IndicesOfTestSet.Add(i);
                }
            }
        }
    }
}
