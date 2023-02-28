# gcms-minerva

Once Minerva is running (execute main.exe in the .zip file under windows_build) you should see two windows, status:
![Minerva empty terminal window](https://github.com/snowdonr/gcms-minerva/tree/main/docs/Minerva_w_terminal.png?raw=true)
<img src="https://github.com/snowdonr/gcms-minerva/tree/main/docs/Minerva_w_terminal.png?raw=true">

And the main interface:

![Minerva main window](https://github.com/snowdonr/gcms-minerva/tree/main/docs/Minerva_w_main_menu.png?raw=true)
<img src="https://github.com/snowdonr/gcms-minerva/tree/main/docs/Minerva_w_main_menu.png?raw=true">

Generally the status window can be ignored or minimized unless there is a problem, or you would like more details about what the software is currently doing.

The main interface is grouped into three sections of configuration options, and buttons to start processing and display results, detailed below.


# Experiment
The first section is the most commonly changed values. Once you have chosen values in the other sections, processing new batches of samples will only require changes here.

Experiment Name: A unique name that will be added to all output files and used to keep track of this method and set of samples. If data is later uploaded to a database, this name will be used as a key to keep track of these results. Note that because alignment is specific to an experiment, comparing data from a different run (even with the same method setup) cannot be assured.

Source Directory: The directory with all the CDF files (e.g. produced by an AIA export) to be aligned. All files in subdirectories are also added. Any other files will be ignored, with a message in the status window.

Sample Metadata: Spreadsheet file with sample data unrelated to GCMS.

Destination Directory: Directory for output files, this must have been created before Minerva begins processing. The directory can be empty or have prior results. If the same Experiment Name is used, files will be overwritten. Results include .csv which can be opened in Excel and an sqlite .db file for import to other software:
LL_merged_97_aligned_peaks.bin
LL_merged_97_area.csv
LL_merged_97_area_common_ion.csv
LL_merged_97_grouped.csv
LL_merged_97_grouping.csv
LL_merged_97_Output.db
LL_merged_97_results.csv
LL_merged_97_rt.csv

The _results.csv file is the last step and it includes identification information, aligned areas, and possible compound identifications.

# Method
Method items determine how samples are processed in Minerva, which peaks are considered and how samples are aligned.

Method Name: A description of these settings which is stored in the results.

Minimum RT, Maximum RT: Time, ending in s or m (for seconds or minutes), to begin and end processing peaks.

Minimum Area: Minimum peak size of significance (noise baseline) to reduce the number of results.

RT sensitivity (sec): Maximum variance of RT between samples. If peaks are not aligning, check the RT of a peak in all samples to make sure this is large enough, make sure to convert to seconds.

Gap Penalty: Value from 0-1. Higher values consider 'more' different mass spectra to be the same for alignment, lower values will split 'more' different mass spectra to separate groups with a smaller gap in RT.

Min Peak %: Peaks with intensity below this percent of the most intense peak will be excluded (in testing settings the minimum area has already removed all of these)

Minimum Ion count: Peaks with fewer fragments can be excluded, note the minimum fragment area is not the same as the minimum peak area.

# File Locations

Files and directories that are not typically changed when a new project is run are listed here. These should be set up on first install and can then be basically ignored.


Nist Library Directory: Directory containing the nist.db and related files, e.g. MSSEARCH/mainlib

Nist Working Dir: Temporary storage directory for NIST related data

Excel Library: Excel file with a user library of compounds and their identifying RT and largest MS fragments. Should list CAS number to reduce ambiguity. Can also contain other compound data like formula, MW, odor etc.

Temporary Directory: Location to save checkpoint files that are used to speed the processing of files or experiments multiple times. Any files here can be deleted and they will be recreated as required, but alignment of large sets can take hours, which will be skipped if the cached results file exists (keyed to the exact file and experiment name). These can be large files, make sure there is ~300MB per file and ~1GB per alignment experiment.

Recalculate Alignment Cache: When checked, processing will be done regardless of the cached file. This is required for method changes to be used when an experiment with the current name already exists but you require results with different settings.

# Results

When a sample set is completed a spreadsheet of summary results is created with the name:
<experiment name>_summary.csv

This includes 3 sections:
Run summary of the project and method names, date and total aligned groups for the entire set
Sample summary of each sample in the run, giving the number of peaks and some data of the top 50* peaks
Compound summary of the largest aligned groups (by total area of all samples)
*User setting

Other Tools
Utilities/minerva_xl_compare.py
Result comparison tool, including a comparison of [peak_ID, retention time, corrected area, and others].




# Acknowledgments
