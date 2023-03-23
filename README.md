# MS1_DataSimulator

        /// The ultimate goal here is to create an mzML file of simultated mass spec data
        /// Each scan will be a theoretical MS1 spectrum of intact peptides or proteoforms
        /// There will be from 1 to n different peptides/proteforms in each scan
        /// Each peptide/proteoform can have multiple charge states in a single scan
        /// The intensities or peptides and proteoforms will vary over up to four orders of magnitude
        /// Consecutive scans will have the same peptides/proteoforms coming in and going out as if it was a real chromatography peak
        /// Peptides/proteoforms will begin eluting at random times.
        /// Training and test data will appear as totally separate mzML files to ease analysis.
