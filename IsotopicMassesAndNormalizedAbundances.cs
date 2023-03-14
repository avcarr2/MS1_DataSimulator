using Chemistry;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text; 
using System.Threading.Tasks;

namespace MS1_DataSimulator
{
    internal class IsotopicMassesAndNormalizedAbundances
    {
        public readonly PeptideWithSetModifications peptideWithSetModifications;
        public readonly double intensityFractionRequired;
        public readonly double[] massShiftsFromMonoisotopicMass;
        public readonly double[] isotopeAbundances;
        public IsotopicMassesAndNormalizedAbundances(PeptideWithSetModifications pwsm, double minIntensityFraction = 0.01) 
        {
            peptideWithSetModifications = pwsm;
            intensityFractionRequired = minIntensityFraction;
            (massShiftsFromMonoisotopicMass, isotopeAbundances) = ComputeMassesAndAbundances();
        }

        /// <summary>
        /// Yields one isotope envelope and the abundances. you must add the monoisotopic mass and then compute m/z
        /// </summary>
        /// <returns></returns>
        private (double[], double[]) ComputeMassesAndAbundances()
        {
            ChemicalFormula formula = peptideWithSetModifications.FullChemicalFormula;

            var isotopicMassesAndNormalizedAbundances = new List<(double massShift, double abundance)>();

            var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);

            double[] masses = isotopicDistribution.Masses.ToArray();
            double[] abundances = isotopicDistribution.Intensities.ToArray();

            for (int i = 0; i < masses.Length; i++)
            {
                masses[i] += (peptideWithSetModifications.MonoisotopicMass - formula.MonoisotopicMass);
            }

            double highestAbundance = abundances.Max();

            for (int i = 0; i < masses.Length; i++)
            {
                // expected isotopic mass shifts for this peptide
                masses[i] -= peptideWithSetModifications.MonoisotopicMass;

                // normalized abundance of each isotope
                abundances[i] /= highestAbundance;

                // look for these isotopes
                if (abundances[i] > intensityFractionRequired)
                {
                    isotopicMassesAndNormalizedAbundances.Add((masses[i], abundances[i]));
                }
            }
            return (masses,abundances);
        }

        public (double[], double[]) ComputeMzAndIntensity(int charageState)
        {
            return (massShiftsFromMonoisotopicMass.Select(i => (i + this.peptideWithSetModifications.MonoisotopicMass).ToMz(charageState)).ToArray(), isotopeAbundances);
        }
    }
}
