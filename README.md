# Time-Homogeneity_of_Infection
Checks constant probability of infection

README

Code is written in Python.
Installation of the appropriate Python packages is assumed.

This current working directory contains:  
    1. This file README  
    1. The directory 5_llr_test/ (abbreviating "Log Likelihood-Ratio Test")  

5_llr_test/ contains:  
    1. Data/ contains:  
        1. Bib/bib.csv  
            contains metadata for SIV-macaque database, stored as CSV files.  
        1. CSV files  
            named [PubMed ID]_[Source of data, e.g., figure or table].csv  
    1. Executable/ contains:  
        1. jls_animal_format.py, 1. jls_animal_io.py, 1. jls_animal_model.py  
            auxiliary subroutines  
        1. to_summary.log  
            output log file produced as to_summary.py runs, to identify computational problems 
        1. to_summary.py  
            the program file that inputs Data/ and computes Output/  
        1. to_summary_controls_make.py  
            drives the program in to_summary.py by providing parameters  
    1. Output/Backup/summary.csv is the output of the program, appearing as Output/summary.csv.  
        Note that the maximization of the frailty models is not perfectly reproducible, 
        because their likelihood surfaces are multimodal. Maximization therefore requires a 
        randomized algorithm, which may not be perfectly reproducible across hardware. The
        Supplementary Information of the article associated with this repository explains the 
        columns of summary.csv.
        
