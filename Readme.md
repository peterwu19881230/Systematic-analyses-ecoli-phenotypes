# Information on running the code

* To generate the complete results in this study: 
	1. Download the complete folder
    2. Run the script Nichols_preload.R. This will install all required packages and load important files. 
    3. Change the directories in each .R script based on your the place where the folder was downloaded. Run the scripts to reproduce the results
    Run phenotypeData.R to generate various versions of the phenotypic profile data.
    4. In the id_mapping folder, run these files to generate necessary data before running scripts in the map folder:
        + generate_id_operon.R, get_uniProtID.R, parse_CGSC_Strains.R, parse_regulonDB.R, parse.genes.dat.R, parse.pwy-genes.dat.R, parse.transunits.dat.R, proteinComplex.R, scrapeKEGG.R.
    5. In the map folder, run map1, map2 map3 in sequence to generate the table that maps strains to their annotations.
    6. In the pairwise_obj file, run the following files in sequence to obtain pairwise co-annotation information and various similarity metrics: generate_strain1_strain2_GO_separated.R, generate_strain1strain2_allAnnotations.R, generate_strain1strain2_allDistances.R, generate_strain1strain2_allAnnotations_allDistances.R.

* Notes:
    + Results are summarized in the Jupyter Notebook: paper_figs.ipynb    
    + Change the designated directories for the .R files based on the downloaded path of the local folder, if needed.
    + All data (files with .RData, excel, .txt...etc., extensions ) required to run the scripts as well as intermediate data files are all provided in the Data  folder
    + Self-defined functions are in the functions folder
    + .gitignore and .gitattributes are configuration files for git, which should be ignored from users' end
    
  
