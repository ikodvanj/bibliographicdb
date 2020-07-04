## World Health Organization (WHO) COVID-19 database: WHO needs it?

Accompanying git repository with code and data used in the article "World Health Organization (WHO) COVID-19 database: WHO needs it?" by Kodvanj I, Homolak J, Virag D, Trkulja V. 

-----

#### R script dependencies:
    ggplot2
    dplyr
    ggsci
    stringr
    devEMF
    pubmedR
    VennDiagram
    bibliometrix
    naniar
    
User has to provide API key to run `pmApiRequest()` function in R script. API key can be obtained on PubMed.

-----

#### Raw data

All raw data is contined in the data folder. 

Data retrieved from PubMed with pubmedR (2020-06-26):

    * pubmed_at_26062020.RDS

Data exported from Scopus (2020-06-27):

    * scopus_top2000.bib 
    * scopus_bottom2000.bib 

WHO database (2020-06-26):

    * WHOdb.RDS
    
Search results on specific topic (PubMed database, 2020-06-26):

    * Azithromycin.RDS
    * Chloroquine.RDS
    * Depression.RDS
    * Diabetes.RDS
    * Hypertension.RDS
    * Quarantine.RDS
    * Shock.RDS
    * Tocilizumab.RDS

Search results on specific topic (Scopus database, 2020-06-26):

    * s_azithromycin.csv
    * s_chloroquine.csv
    * s_depression.csv
    * s_diabetes.csv
    * s_hypertension.csv
    * s_quarantine.csv
    * s_shock.csv
    * s_tocilizumab.csv

Search results on specific topic (WHO database, 2020-06-26):

    * w_azithromycin.csv
    * w_chloroquine.csv
    * w_depression.csv
    * w_diabetes.csv
    * w_hypertension.csv
    * w_quarantine.csv
    * w_shock.csv
    * w_tocilizumab.csv

Older version of WHO database (2020-05-19):

    * CSV as at 19 May 2020-Full database.xlsx

-----

#### Processed data

`WHO - final db - duplicated based on URL.xlsx`, in the Results folder, contains duplicate entries from the WHO database.


