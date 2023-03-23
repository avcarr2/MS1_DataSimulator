

using System.Collections.Generic;
using System;

namespace MS1_DataSimulator
{
    internal class ChargeStatesAndEnvelopeAbundances
    {
        public readonly int[] ChargeStates;
        public readonly double[] EnvelopeAbundances;

        public ChargeStatesAndEnvelopeAbundances(int[] chargeStates, double[] envelopeAbundances)
        {
            ChargeStates = chargeStates;
            EnvelopeAbundances = envelopeAbundances;
        }

        public ChargeStatesAndEnvelopeAbundances(int maxNumberChargeStates, int minChargeState, int maxChargeState, double minEnvelopeAbundance = 0.1) 
        { 
            List<int> chargeStates = new();
            
            if(maxNumberChargeStates > (maxChargeState - minChargeState))
            {
                maxNumberChargeStates = maxChargeState - minChargeState;   
            }

            Random rnd = new Random();

            int numChargeStates = (int)Math.Round(maxNumberChargeStates * rnd.NextDouble() + 0.5,0);//add 0.5 makes it so we get at least 1 charge state.
            List<int> possibleChargeStateList = Enumerable.Range(minChargeState, maxChargeState - numChargeStates).ToList();
            int index = rnd.Next(possibleChargeStateList.Count);
            int firstChargeState = possibleChargeStateList[index];

            for (int i = 0; i < numChargeStates; i++)
            {
                chargeStates.Add(firstChargeState + i);
            }
            ChargeStates = chargeStates.ToArray();

            //abundance of each envelope is chosen randomly with the min abundance definced in the constructor. 
            //the sum of all envelope abundances equals 1
            double[] envelopeAbundances= new double[numChargeStates];
            for (int i = 0; i < envelopeAbundances.Length; i++)
            {
                envelopeAbundances[i] = minEnvelopeAbundance + rnd.Next();
            }
            double arraySum = envelopeAbundances.Sum();

            for (int i = 0; i < envelopeAbundances.Length; i++)
            {
                envelopeAbundances[i] /= arraySum;
            }
            EnvelopeAbundances= envelopeAbundances;
        }
    }
}
