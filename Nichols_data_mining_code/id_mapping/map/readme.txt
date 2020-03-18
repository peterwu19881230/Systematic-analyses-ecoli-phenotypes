#run map1.R, map2.R, map3.R to generate id_allAttributes


* ids: from row no. of Nichols' original data

* Original Name: original row names from Nichols' data

* ECKs: painfully, manually corrected ECKs based on "original ECKs"

* associated_gene_names: mainly based on "original ECKs" with some manual correction

* other.synonyms: based on EcoCyc's genes.dat

* JW: from 2006 Baba;s paper: inline-supplementary-material-4(Baba et al., 2006).xls, mapped by ECK

* JW2: from scraping CGSC (http://cgsc2.biology.yale.edu)

* position: from scraping CGSC (http://cgsc2.biology.yale.edu)

* EcoCyc ID: based on genes.dat joined by ECKs

* bNumber: based on genes.dat joined by EcoCyc ID

* strain_availability: from scraping CGSC (http://cgsc2.biology.yale.edu), joined by JW numbers

* GO: map from the .gaf file Suzi made, joined by EcoCyc ID

* Pwy: from EcoCyc's pathwayscol, joined by EcoCyc ID

* pcomplex: from EcoCyc's protcplxs.col, joined by EcoCyc ID

* operon: from EcoCyc's transunits.dat, joined by EcoCyc ID

* regulator_name:   from regulonDB's object_synonym.txt get ECK12 numbers, joined by bNumbers from regulonDB's gene.txt get ECK12 numbers for genes, joined by ECK12 numbers. From regulonDB's generegulation_tmp.txt get regulator names, joined by ECK12 numbers

* kegg_modules: from scrapping genome.jp/kegg, joined by bNumber (11/11/2019 The result of KEGG modules have drastically changed compared to my previous download. I used the newest version)

* UniProtID: retrieve uniprotID from a package by EcoCyc product ID. EcoCyc product ID was mapped using genes.dat, joined by EcoCyc ID 























