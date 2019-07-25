# 5Beam-Turbulence-Methods

A set of codes to process turbulence data from a 5-beam Accoustic Doppler Current Profiler (specifically the Nortek Signature 1000 AD2CP).

This work is published on the Journal of Atmospheric and Oceanic Technology as Guerra, M. and Thomson, J. "Turbulence measurements from 5-beam acoustic Doppler profilers", and is available for open access at https://doi.org/10.1175/JTECH-D-16-0148.1 

The following turbulence parameters can be estimated using this scripts and functions:

1. Turbulent kinetic energy dissipation rate from TKE frequency spectra and structure function
2. Reynolds stresses using all 5 beam information and tilt corrections 
3. Turbulent kinetic energy vertical shear production 

Plotting functions are also included.

The 5beamMethods_Document.pdf file contains a detailed explanation of each of the scripts and functions.

The following working directory is recomended to use this scripts:

    5Beam-Turbulence-Methods/
    	README.md
    	5BeamMethods_Document.pdf
    	Guerra_Thomson_FiveBeamMethods_Submitted_2016.pdf
    	RS_5beam.m
    	RS_VT.m
    	Signature_E_P_Profile_Plots.m
    	Signature_Production_5beam.m
    	Signature_QC.m
    	Signature_RawData_ENU.m
    	Signature_Reynolds_Stress_5Beam.m
    	Signature_Spectra_DissipationRate.m
    	Signature_Spectra_Plots.m
    	Signature_Spectra_w.m
    	Signature_StructureFunction_Plots.m
    	Signature_StructureFunction_dissipation.m
    	Signature_Usigned.m
    	Signature_join_files.m
    	dissipationMG_SF.m
    	sign_speed.m
    	signatureAD2CP_beam2xyz_enu.m
    	spectralestimate.m
    	structureFunction.m
    	BinDataSignature/ 
    	RawData/ %put raw data from MIDAS software here
    	SDF_Signtaure/
    	Signature_Dissipation/
    	Signature_Production/

Please use these codes carefully and check if functions and processing steps are suitable for your data.

For any questions, comments and suggestions please contact me at mguerrap@uw.edu or at mguerra@dal.ca
