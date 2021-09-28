# Membrane-Fluctuation
Grid-based and Least-squares analysis of surface and direction fields of fluctuating flat membrane

Following text may not be accurate

## Fourier_GB_2D
Inputs:

        Head, tail, and surface positions of each lipid molecule in a bilayer
Outputs:

        Time-series data of spectra (Gridded surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Related wave numbers [nm]

## Fourier_LS
Inputs:

        Head, tail, and surface positions of each lipid molecule in a bilayer
Outputs:

        Time-series data of spectra (Least-squares fitted surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Related wave numbers [nm]

## Blocking Spectra
Inputs:

        Time-series data of spectra (Gridded surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Related wave numbers [nm]
Outputs:

        Spectra (Surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Error bars of spectra
Carries:

        Related wave numbers [nm]

## Get_Params
Inputs:

        Spectra (Surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Error bars of spectra
        Related wave numbers [nm]
Outputs:

        Nboot=2000 bootstrapped fitting parameters
Carries:

        Spectra (Surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Error bars of spectra
        Related wave numbers [nm]

## Save Results
Inputs:

        Nboot=2000 bootstrapped fitting parameters
        Spectra (Surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Error bars of spectra
        Related wave numbers [nm]
Outputs:

        Best fitting parameters
        Error bars of best fitting parameters
Carries:

        Nboot
        Number of used wave numbers
        Nboot=2000 bootstrapped fitting parameters
        Spectra (Surface, parallel, and perpendicular direction fluctuations) [$kT, nm$]
        Error bars of spectra
        Related wave numbers [nm]
        

