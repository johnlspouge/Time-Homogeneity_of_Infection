# Time-Homogeneity_of_Infection
Checks constant probability of infection

README

Code is written in Python.
Installation of the appropriate Python packages is assumed.

This current working directory contains:
    <ol>
        <li>This file README</li>
        <li>The directory 5_llr_test/ (abbreviating "Log Likelihood-Ratio Test")</li>
    </ol>
5_llr_test/ contains:  
    <ol>
    <li>Data/ contains:  
        <ul>
            <li>Bib/bib.csv  
            contains metadata for SIV-macaque database, stored as CSV files.</li>
            <li>CSV files
            named [PubMed ID]_[Source of data, e.g., figure or table].csv</li> 
        </ul>
    <li>Executable/ contains:  
        <ul>
            <li>jls_animal_format.py, jls_animal_io.py, jls_animal_model.py  
                auxiliary subroutines</li>
            <li>to_summary.log  
            output log file produced as to_summary.py runs, to identify computational problems</li>
           <li>to_summary.py  
            the program file that inputs Data/ and computes Output/</li>
            <li>to_summary_controls_make.py  
            drives the program in to_summary.py by providing parameters</li>
        </ul>
    </li> 
    <li>Output/Backup/summary.csv is the output of the program, appearing as Output/summary.csv.  
        Note that the maximization of the frailty models is not perfectly reproducible, 
        because their likelihood surfaces are multimodal. Maximization therefore requires a 
        randomized algorithm, which may not be perfectly reproducible across hardware. The
        Supplementary Information of the article associated with this repository explains the 
        columns of summary.csv.</li>
    </ol>
        
