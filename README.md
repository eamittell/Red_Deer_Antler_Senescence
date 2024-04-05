## Antler size in red deer: declining selection and increasing genetic variance with age, but little senescence

The code in this repository is to run the models presented in the manuscript above. In this study we assess quantitative genetic parameters to gain an evolutionary explanation as to why antler size traits show little sensescnce within this wild population of red deer.

Elizabeth A. Mittell^1^*$\dag$, Priyam Mandaliya^1^$\dag$, Josephine M. Pemberton^1^, Alison Morris^1^, Sean Morris^1^, Susan E. Johnston^1^ & Loeske E. B. Kruuk^1^

$\dag$Joint first authors

^1^Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, UK

*Corresponding author: Elizabeth A. Mittell, e.mittell@ed.ac.uk, e.mittell@gmail.com; Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, EH9 2LD, UK

J.M.P, A.M., S.M, L.E.B.K, S.E.J. and E.A.M. collected the data. E.A.M., L.E.B.K. and P.M. analysed the data. E.A.M. is responsible for the code in this repository.

### Data and Scripts

DataFormRandomIds.csv contains the phenotypic data, and PedigreeRandomIds.csv contains the pedigree for the red deer population. The IDs within both these files have been randomised. This means that they should link to each other, but will not link to any other data published from the Rum Red Deer Project.

Models_Antler_Senescence.R contains the models that were used in the manuscript. The packages used are shown within the script. These models were run in various versions of R. All run and are installable in R version 4.3.2 as of 5th April 2024 on macOS Monterey version 12.7.1.

### RUM RED DEER PROJECT DATA REUSE STATEMENT:
		
The attached files contain data derived from the long term field project monitoring individual red deer on Rum and their environment. This is a request to please let us know if you use them. Several people have spent the best part of their careers collecting the data. If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact the authors or Professor Josephine Pemberton (J.Pemberton@ed.ac.uk) before doing so.
		
[NB. If you are interested in analysing the detailed project data in any depth you may find it helpful to have our full relational database rather than the file(s) available here. If so, then we have a simple process for bringing you onto the project as a collaborator.]
		
1) The data can be subject to change due to updates in the pedigree, merging of	records, occasional errors and so on.
		
2) The data are complex and workers who do not know the study system may benefit from advice when interpreting it.
		
3) At any one time a number of people within the existing project collaboration	are analysing data from this project. Someone else may already be conducting the analysis you have in mind and it is desirable to prevent duplication of effort. Frequently these projects will be the only research of an early-career researcher.
		
4) In order to maintain funding for the project, every few years we have to write proposals for original analyses to funding agencies. It is therefore very helpful for those running the project to know what data analyses are in progress.
		
5) Individual identifiers may vary relative to other data archives from papers using	the individual-level data. Therefore, contacting us is advised.
		
6) Data are here for the explicit purpose of reproducible research.
