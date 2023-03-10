# Intro
Minerva is a software tool that processes raw GC-MS data from multiple samples and then aligns peaks into a database, facilitating multi-sample comparison of multiple peaks at once, and enabling downstream analysis with machine learning approaches. The main Minerva features are:

- Alignment of peaks from multiple GC-MS samples into a single database.
- Detection of peak overlaps and deconvolution.
- Visualization of peak alignments for quality control and further analysis.
- Graphic User Interface.
- Configurable to different sample matrices and data analysis approaches.

Minerva is an open source software released under the MIT license. It was developed by Ryan Snowdon at the University of Calgary, Department of Geoscience, under the supervision of Renzo Silva and funding from Dr. Steve Larter’s group (PRG, via Project Rip van Winkle), and a GRI grant. Visit the [PRG group](https://ucalgary.ca/prg) at the University of Calgary webpage for more information on the academic context Minerva was created. 


# Instructions

## Download and installation
Download the .zip file under */windows_build* and execute **main.exe**.
Once Minerva is running you should see two windows:

**status**: generally the status window can be ignored or minimized unless there is a problem, or you would like more details about what the software is currently doing.  

![Minerva empty terminal window](docs/Minerva_w_terminal.png?raw=true)

**main interface**: the main interface is grouped into three sections of configuration options, and buttons to start processing and display results, detailed below.  

![Minerva main window](docs/Minerva_w_main_menu.png?raw=true)

## Experiment
The first section is the most commonly changed values. Once you have chosen values in the other sections, processing new batches of samples will only require changes here.

**Experiment Name**: A unique name that will be added to all output files and used to keep track of the results from this set of samples. If data is later uploaded to a database, this name will be used as a key to keep track of these results. Note that because alignment is specific to an experiment, comparing data from a different run (even with the same method setup) cannot be assured.

**Source Directory**: The directory with all the CDF files (e.g. produced by an AIA export) to be aligned. All files in subdirectories are also added. Any other files will be ignored, with a message in the status window.

**Destination Directory**: Directory for output files, which must have be created before Minerva begins processing. The directory can be empty or have prior results. If the same <Experiment Name> is used, any existing files will be overwritten! Results include .csv which can be opened in Excel and an sqlite .db file for import to other software:  

|experiment name|_aligned_peaks.bin  
|experiment name|_area.csv  
|experiment name|_area_common_ion.csv  
|experiment name|_grouped.csv  
|experiment name|_grouping.csv  
|experiment name|_Output.db  
|experiment name|_results.csv  
|experiment name|_rt.csv  

The |experiment name|_results.csv file is the last step and it includes identification information, aligned areas, and possible compound identifications.

## Method
Method items determine how samples are processed in Minerva, which peaks are considered and how samples are aligned.

**Method Name**: A description of these settings which is stored in the results.  
**Minimum RT, Maximum RT**: Time, ending in s or m (for seconds or minutes), to begin and end processing peaks.  
**Minimum Area**: Minimum peak size of significance (noise baseline) to reduce the number of results.  
**RT sensitivity (sec)**: Maximum variance of RT between samples. If peaks are not aligning, check the RT of a peak in all samples to make sure this is large enough, make sure to convert to seconds.  
**Gap Penalty**: Value from 0-1. Higher values consider 'more' different mass spectra to be the same for alignment, lower values will split 'more' different mass spectra to separate groups with a smaller gap in RT.  
**Min Peak %**: Peaks with intensity below this percent of the most intense peak will be excluded (in testing settings the minimum area has already removed all of these)  
**Minimum Ion count**: Peaks with fewer fragments can be excluded, note the minimum fragment area is not the same as the minimum peak area.

## File Locations
Files and directories that are not typically changed when a new project is run are listed in this tab. These should be set up on first install and can then be basically ignored.

**Excel Library**: Excel file with a user library of compounds and their identifying RT and largest MS fragments. Should list CAS number to reduce ambiguity. Can also contain other compound data like formula, molecular weight, etc.  
**Temporary Directory**: Location to save checkpoint files that are used to speed the processing of files or experiments multiple times. Any files here can be deleted and they will be recreated as required, but alignment of large sets can take hours, which will be skipped if the cached results file exists (keyed to the exact file and experiment name). These can be large files, make sure there is ~300MB per file and ~1GB per alignment experiment.  
**Recalculate Alignment Cache**: When checked, processing will be done regardless of the cached file. This is required for method changes to be used when an experiment with the current name already exists but you require results with different settings.

## Results
When a sample set is completed a spreadsheet of summary results is created with the name:
<experiment name>_summary.csv

This includes 3 sections:
- Run summary of the project and method names, date and total aligned groups for the entire set.  
- Sample summary of each sample in the run, giving the number of peaks and some data of the top 50* peaks.  
- Compound summary of the largest aligned groups (by total area of all samples).  
**User setting*

## Other Tools
**Utilities/minerva_xl_compare.py**  
Result comparison tool, including a comparison of [peak_ID, retention time, corrected area, and others].

## Troubleshooting
We always appreciate feedback and suggestions. If you encounter any problems or have any questions, please feel free to submit an issue to our GitHub repo.

# Acknowledgments
This research was undertaken thanks in part to funding from the Canada First Research Excellence.   
Minerva partially uses the [PyMassSpec · PyPI](https://pypi.org/project/PyMassSpec/) library, which was released under the GNU General Public License version 2.

We would like to acknowledge the support of the following people and organizations:
- Lloyd Snowdon, Kim Nightingale, and Mengsha Yin as early adopters.  
- Dr. Steve Larter and the PRG group, including the Project Rip van Winkle team and sponsors.  
- The University of Calgary: GRI and the Department of Geoscience. 
