"""
AU Mic Project
05/25/16
Author: Cail Daley
"""


This repository contains my work on AU_Microscopii, including data analysis, modeling, and fitting.


Directories in AU Mic folder:
    Imaging
    Modeling




Visibilities:

    Contains visibilities and image files for AU Mic. Observations were taken on three different nights:
        26mar2014
        18aug2014
        24jun2015

    For each night, there is a <date>_reduced directory containing the original data files that were the result of the data reduction performed by Meredith Hughes and I at the NRAO headquarters in October of 2015. There is also a <date>_main directory for each night containing images and visibilities created in the process of modeling. The _reduced and _main directories contain a variety of filetypes and naming conventions, which I will explain below:


        spw:    An abbreviation of spectral window; there are four spectral windows in this dataset.

        _flare/_noflare:    A flare occured during observation on the night of 24jun2015; part of the reason we performed data reduction at NRAO was to correct for this. We divided this data set into two categories, the part of the observation before the flare occured (_noflare) and the part of the observation during which the flare occured (_flare). We split the _flare subset into even smaller time windows, correcting for the flare in each window; this is the reason for the _3, _4, _5, etc. naming convention. NOTE: the 24jun2015_main directory contains concatenations of the _flare and _noflare time windows, and thus contains the files used for modeling.

        uid__*:     Visibilities of this type contain the raw data that was given to us by ALMA/NRAO.

        .fixvis:    The phase center of visibilities with this extension have been changed so as to account for the proper motion of the star.

        .uvsub:     Visibilities with this extension have had a model point source subtracted from the visibility data to remove any star flux. For the 26mar2014 and 18aug2014 dates, files with this extension have been succesfully reduced and can be used for modeling.

        .concat:   Visibilities with this extension are concatenations of the different time windows for the flare date, 24jun2015. Files with this extension have been sucessfully reduced and can be used for modeling.

        .image, .mask, .model, .psf, .residual, .sumwt: Files with these extensions were created by the CASA tclean task; I have little experience with any of these extensions except for .image, an image file,  which can be converted into a FITS file or opened with CASA's viewer.

        .vis:   Miriad's visibility extension; visibilities of this type can be processed using miriad.

        .bm, .cl, .cm, .mp:     Files with these extensions by miriad's cleaning procedure; .cm, an abbrevation of clean map, is an image file and can be imaged using miriad's cgdisp command.

        .uvf:   A U-V FITS file; visibilities in this format can be manipulated in python.

        aumic_ALL/: The visibilities/images in this directory are concatenations of all spectral windows and dates.



    NOTE: The 0th and 3rd spectral windows for the flare date (24jun2015) deviate from the model much more than any of the other spectral windows/dates (at the moment, all datasets have a reduced chi^2 of roughly 2, except for the 0th spw, which has a value of ~5, and the 3th spw, which has a value of almost 5000. I have yet to figure out where this is coming from.




Modeling:

    This directory contains everything needed for creating a model of AU Mic; running the script "Modeling_Code.py" generates models titled 'aumicmodel#' with numbers ranging from 0~11 (one for each of the four spectral windows on each of the three dates), creates both model visibilities and cleaned images, and returns the chi squared of the model. The model is created using the radiative transfer code radmc3dpy (http://www.ast.cam.ac.uk/~juhasz/radmc3dPyDoc/index.html), the code for which can be found in the directory "radmc3dPy-0.28.1".

    The model makes use of ramc3dpy's 'ppdisk' parameter set; the parameters I have changed are listed at the top of "Modeling_Code.py", but a full list of the parameters can be found in the file problem_params.inp. Many of the default parameters have been replaced with Macgregor et al. 2013's best fit model parameters or have been modified based on suggestions my research advisor, Professor Meredith Hughes; star mass and luminosity were taken from Plavchan et al. 2009.

    I have copied below a partial list of the input files from documentation in 'setup.py' in the radmc3dPy-0.28.1/radmc3dPy directory; not all are necessarily present/created in this specific model:

    Files written by problemSetupDust() for RADMC-3D:

        * dustopac.inp             : Dust opacity master file.

        * wavelength_micron.inp    : Wavelength grid.

        * amr_grid.inp             : Spatial grid.

        * stars.inp                : Input radiation field (discrete stellar sources).

        * stellarsrc_density.inp   : Input radiation field (continuous stellar sources).

        * stellarsrc_templates.inp : Input radiation field (continuous stellar sources).

        * dust_density.inp         : Dust density distribution.

        * radmc3d.inp              : Parameters for RADMC-3D (e.g. Nr of photons to be used, scattering type, etc).

        The dustkappa_silicate file is generated when the model is created, and is determined by the dustkappa_ext parameter under dust opacity.

    The code ouputs at least three files, radmc3d.out, spectrum.out, and image.out; however, I have not needed to use any of these files thus far in the modeling process.

    The names of the data files used in the modeling code are included in the modeling script under the section 'Variables and input files'.




    The script 'var_vis.py' was written by Kevin Flaherty, and I have used it to generate weights to be used in the modeling process.
