# chemokine_buffering_paper

This repo hosts the code for the python-based quantitative image and data analysis used in the paper entitled "Dynamic buffering of extracellular chemokine enables robust adaptation during directed tissue migration" by Wong and colleagues, as described in the `METHOD DETAILS` chapter in the sections `High resolution imaging and processing` and `3D image analysis and quantification`.  All code was written by Jonas Hartmann.


### Analyses Included

- Intensity analysis (as shown in figures `3J` and `3L`)
- Point cloud analysis (as shown in figures `2L`, `4B`, `6D`)
- Colocalization analysis (as shown in figures `5J` and `6N`)


### Code Structure

- The code in this repo is structured in three 'layers':
	- Modules (`quant`, `coloc`, `util`) contain refactored functions, mainly for the `RUN` notebooks
	- `RUN` notebooks are pipelines built on the modules; they take in images and output quantitative measures
	- `ANALYSIS` notebooks take in quantitative measures and produce various plots and statistics


### Data Availability

This repository only hosts a couple of example images (in `data_ex`) to help follow through the `RUN` notebooks. These examples have been zipped to comply with GitHub's 100MB maximum file size policy; they need to be unzipped before use.

The `ANALYSIS` notebooks are set up to reproduce the analyses presented in the paper and expect the full data to be present in a directory named `data_full`. This large dataset is available only on request via the Lead Contact, Darren Gilmour (![Email DG](_images/email_DG.png)).


### Workflow

1. Convert raw images to 8bit `.tif` images using `macro_8bit.ijm` in [Fiji](https://fiji.sc/) with a suitable intensity range (logged in `metadata.xlsx`).
	- The example data provided in `data_ex` has already been converted.
	- These example images must be unzipped before use and the resulting `.tif` files must be located directly in the `data_ex` directory.
2. Manually identify the coordinates of the tissue's apical focus point, note them in `metadata.xlsx` and save the file as a tab-separated text file called `metadata.txt`.
	- For the example data provided and the study's full dataset, coordinates and `metadata.txt` are already provided.
3. For both point cloud and intensity analyses, run `RUN_preprocessing.ipynb` to mask the overall tissue.
4. For intensity analysis, run `RUN_intensity_quantification.ipynb`.
5. For point cloud analysis, run `RUN_landmark_extraction.ipynb`.
6. For colocalization analysis, run `RUN_colocalization.ipynb`.
7. Use the the various `ANALYSIS` notebooks to visualize and analzye the respective results.


### Dependencies

- Python 2.7.13 (we recommend the [Anaconda distribution](https://www.anaconda.com/distribution/))
- Scientific python stack including numpy, scipy, scikit-image, matplotlib
	- The versions used by us were numpy 1.11.3, scipy 2.0, scikit-image 0.13.0, and matplotlib 1.5.1
- The tifffile module (can be installed from [Anaconda cloud](https://anaconda.org/conda-forge/tifffile))


### Contact and Support

- The study's corresponding authors are Mie Wong (![Email MW](_images/email_MW.png)) and Darren Gilmour (![Email DG](_images/email_DG.png)).
- For questions regarding the code in this repository, contact Jonas Hartmann (![Email JH](_images/email_JH.png)) or open an issue on GitHub. Note that we cannot promise support for any use cases other than direct reproduction of the study's results.
- To request the full data or any other materials, resources or reagents, please get in touch with the study's Lead Contact, Darren Gilmour (![Email DG](_images/email_DG.png)).