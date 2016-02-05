----------code/-------------
Directory containing some scripts for computing hitran curves and fitting them to data.



---------Standalone Scripts-------------

masterFittingScript.py
    ->	Wrapper script to automate fitting of data.
    ->	Prompts for data file, prompts for molecules and isotopologues to fit
	extracts proper range from HITRAN .par files (using extract.py), and 
	finally initiates the fitting routine (using optimizeConc.py).
    ->	For usage instructions run: python masterFittingScript.py --help

masterHitranCurve.py
    ->	Script for generating absorption coefficient curves from HITRAN database.
    ->	Prompts for molecules and isotopologues and their respective concentrations,
	temperature and pressure.
    ->	For usage instructions run: python masterHitranCurve.py --help


to be completed........