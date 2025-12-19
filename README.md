# Isotope ratio analysis of Phytane fragments via Orbitrap mass spectrometry

This repository contains an implementation of the isotopologue calculations and data culling methodology presented in "Measuring carbon isotopic signatures of mass spectrometric fragment ions of phytane using ultra high resolution orbitrap mass spectrometry". It contains a .R scrpt which contains all functions used and can be run to reproduce the major results of the paper. 

# files

- isotopologs.tsv -- text file used to extract .isox files from RAW aquisition files using IsoX. File contains exact masses of target isotopologues of intesest
  
- isox_files -- folder that containes all input data files of individual aquisition runs for all analytes. RAW aquisition files were extracted using IsoX.
  
- qual.figs -- partial output of Phytane_GC_Orbi_data.R. Qualitiy check figures with signal braceting, TICxIT trace, 13C/12C Nio ratio, SNL figure. see figure 1
  
- output.tables -- partial output of Phytane_GC_Orbi_data.R. Contains quality.all.csv file (results of filtering steps), Table1.csv adn Table2.csv
