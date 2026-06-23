Reference Information
=====================
Provenance for this README
--------------------------
* File name: README_Dataset-AphidThermalToleranceLandscape.md
* Authors: Yuan-Jie Li
* Other contributors:Si-Yang Chen, Lisa Bjerregaard Jørgensen, Johannes Overgaard, David Renault, Hervé Colinet, Chun-Sen Ma 
* Date created: 2023-03-08
* Date modified: 2023-05-23

Dataset Attribution and Usage
-----------------------------

* Dataset Title: Data for the article "Interspecific differences in thermal tolerance landscape explain aphid community abundance under climate change"

* Persistent Identifier: https://doi.org/10.5061/dryad.mcvdnck4j

* Dataset Contributors: Yuan-Jie Li,Si-Yang Chen, Lisa Bjerregaard Jørgensen, Johannes Overgaard, David Renault, Hervé Colinet, Chun-Sen Ma 

* License: Use of these data is covered by the following license:
  * Title: CC0 1.0 Universal (CC0 1.0)
  * Specification: https://creativecommons.org/publicdomain/zero/1.0/; the authors respectfully request to be contacted by researchers interested in the re-use of these data so that the possibility of collaboration can be discussed. 

* Suggested Citations:
* Dataset citation:
> Li, Yuan-Jie et al. (2023), Interspecific differences in thermal tolerance landscape explain aphid community abundance under climate change, Dryad, Dataset, https://doi.org/10.5061/dryad.mcvdnck4j

* Corresponding publication:
> Li, Y.-J., S.-Y. Chen, L. B. Jørgensen, J. Overgaard, D. Renault, H. Colinet, and C.-S. Ma. 2023. Interspecific differences in thermal tolerance landscape explain aphid community abundance under climate change. Journal of Thermal Biology:103583.

Contact Information
-------------------

  * Name: Yuan-Jie Li
  * Affiliations: University of Rennes
  * Email: li.yuanjie@aliyun.com
  * Alternate Email: yuanjie.li2021@gmail.com


* Alternative Contact: PI
  * Name: Chun-Sen Ma
  * Affiliations: Hebei university
  * Email: machunsen@caas.cn

- - -

Additional Dataset Metadata
===========================

Acknowledgements
----------------
* Funding sources: This work was supported by National Natural Science Foundation of China 31620103914, the Fundamental Research Funds of CAAS (Y2017LM10), and Innovation Program of CAAS (CAAS-ZDRW202012). The study was funded by the International Research Project (IRP) “Phenomic Responses of Invertebrates to Changing Environments and Multiple Stress” (PRICES) supported by InEE-CNRS. This project was supported by ARED PhD grant, provided by Region Bretagne (project climPuce COH20023).

Dates and Locations
-------------------

* Dates of data collection: Wet lab work performed between November 2020 and November 2021.
* Geographic locations of data collection: Wet lab work performed in Institute of Plant Protection, Chinese Academy of Agricultural Sciences, Beijing, China.

- - -

Methodological Information
==========================
* Methods of data collection/generation: see manuscript for details

- - -

Data and File Overview
======================

Summary Metrics
---------------
* File count: 10
* Total file size: 417 KB
* Range of individual file sizes: 2KB - 367 KB
* File formats: .csv, .txt, .Rmd

Table of Contents
-----------------

*surv.txt 
*TDT analysis.Rmd
*Thermal injury.Rmd
*tdt.csv
*wuhan 2016-05-01 to 2016-08-31.csv
*xinxiang 2016-05-01 to 2016-08-31.csv
*beijing 2016-05-01 to 2016-08-31.csv
*wuhan 2016-11-01 to 2017-02-28.csv
*xinxiang 2016-11-01 to 2017-02-28.csv
*beijing 2016-11-01 to 2017-02-28.csv

Setup
-----
* Recommended software/tools:  RStudio 2021.09.2
* Relationships between files/folders: The first folder"Thermal death time analysis" includes the Rscript and the input data file (surv.txt). the second folder “Thermal injury calculation” contains the Rscript(Thermal injury.Rmd) and the input data files(tdt.csv;wuhan 2016-05-01 to 2016-08-31.csv; xinxiang 2016-05-01 to 2016-08-31.csv; beijing 2016-05-01 to 2016-08-31.csv; wuhan 2016-11-01 to 2017-02-28.csv;xinxiang 2016-11-01 to 2017-02-28.csv; beijing 2016-11-01 to 2017-02-28.csv)
- - -

File/Folder Details
===================

Details for: surv.txt 
---------------------------------------
* Description:This dataset describes the survival proportions of three aphid species with three ages across a broad range of stressful high (34~40°C) and low (-3~-11°C) temperatures 
* Format(s): .txt
* Size(s): 79 KB
* Dimensions: 3042 rows x 8 columns
* Variables:
 *spec: species name, MD = Metopolophium dirhodum, SD = Sitobion avenae, RP = Rhopalosiphum padi
 * age: age of the focal species, 2 = 2-day-old; 6 = 6-day-old; 12 = 12-day-old
 * temp: the tested temperature (degrees Celsius)
 * dur: the exposure duration (min)
 * rep: replication
 * surv: the number of survived aphid individuals after treatment and recover
 * total: the total number of treated aphid individuals
 * survival_rate: the survival rate of tested aphids ranges from 0-1.

Details for: TDT analysis.Rmd
---------------------------------------
* Description: This is the r script to fit the survival curve using the survival dataset and derive the TDT parameters from the thermal death time (TDT) curve. The input data file is surv.txt for our manuscript.
* Format(s): .Rmd

* Size(s):  20 KB


Details for: Thermal injury.Rmd
---------------------------------------
* Description: This code package includes the R script to calculate the thermal injury under fluctuate temperature regimes. The the input data files are temperature files and TDT parameters， which could be find in the same folder.
The mathematical framework and the original r script are referred to Jørgensen, L. B., H. Malte, M. Orsted, N. A. Klahn, and J. Overgaard. 2021. A unifying model to estimate thermal tolerance limits in ectotherms across static, dynamic and fluctuating exposures to thermal stress. Scientific Reports 11:12840

* Format(s): .Rmd

* Size(s):  44 KB


Details for: tdt.csv
---------------------------------------
* Description: This file contains all the parameters of thermal death time curves for each species and each age, which are calculated using r script (TDT analysis.Rmd)

* Format(s): .csv
* Size(s):  1 KB
* Dimensions: 19 rows x 10 columns
* Variables:
 *species: species name, MD = Metopolophium dirhodum, SD = Sitobion avenae, RP = Rhopalosiphum padi
 * age: age of the focal species, 2 = 2-day-old; 6 = 6-day-old; 12 = 12-day-old
 * slope: the slope of thermal death time curve
 * z: thermal sensitivity, z = -1/slope at high temperatures, z = 1/slope at low temperatures.
 * ctm: Ctm is the extrapolated average temperature that would result in a median lethal time = 1 min 
 * ctm1h: Ctm is the extrapolated average temperature that would result in a median lethal time = 1 h
 * rsquare: r square shows how well the data fit the linear regression model 
 * pvalue: p value shows if the line regression relationships are statistically significant
 * intercept: the intercept of the thermal death time curve (the linear regression model)
 * type: "warm" indicates high temperatures. "cold" indicates low temperatures



Details for: wuhan 2016-05-01 to 2016-08-31.csv
---------------------------------------
* Description: the hourly temperature data in Wuhan, N 30.78°, E 114.21° between the period from 2016-05-01 to 2016-08-31, which was downloaded from a global weather API (https://www.visualcrossing.com/weather-data) in May 2022.

* Format(s): .csv
* Size(s):  90 KB
* Dimensions: 2954 rows x 3 columns
* Variables:
 * name: the name of the site (city)
 * datetime: the date and time when the temperature is recorded.
 * temp: temperature (degrees Celsius)
 
 Details for: xinxiang 2016-05-01 to 2016-08-31.csv
---------------------------------------
* Description: the hourly temperature data in Xinxiang, N 35.30°, E 113.92° between the period from 2016-05-01 to 2016-08-31, which was downloaded from a global weather API (https://www.visualcrossing.com/weather-data) in May 2022.

* Format(s): .csv
* Size(s):  103 KB
* Dimensions: 2954 rows x 3 columns
* Variables:
 * name: the name of the site (city)
 * datetime: the date and time when the temperature is recorded.
 * temp: temperature (degrees Celsius)


 Details for: beijing 2016-05-01 to 2016-08-31.csv
---------------------------------------
* Description: the hourly temperature data in Beijing, N 40.07°, E 116.58° between the period from 2016-05-01 to 2016-08-31, which was downloaded from a global weather API (https://www.visualcrossing.com/weather-data) in May 2022.

* Format(s): .csv
* Size(s):  190 KB
* Dimensions: 2954 rows x 3 columns
* Variables:
 * name: the name of the site (city)
 * datetime: the date and time when the temperature is recorded.
 * temp: temperature (degrees Celsius)


Details for: wuhan 2016-11-01 to 2017-02-28.csv
---------------------------------------
* Description: the hourly temperature data in Wuhan, N 30.78°, E 114.21° between the period from 2016-11-01 to 2017-02-28, which was downloaded from a global weather API (https://www.visualcrossing.com/weather-data) in May 2022.

* Format(s): .csv
* Size(s):  87 KB
* Dimensions: 2905 rows x 3 columns
* Variables:
 * name: the name of the site (city)
 * datetime: the date and time when the temperature is recorded.
 * temp: temperature (degrees Celsius)
 
 Details for: xinxiang 2016-11-01 to 2017-02-28.csv
---------------------------------------
* Description: the hourly temperature data in Xinxiang, N 35.30°, E 113.92° between the period from 2016-11-01 to 2017-02-28, which was downloaded from a global weather API (https://www.visualcrossing.com/weather-data) in May 2022.

* Format(s): .csv
* Size(s):  99 KB
* Dimensions: 2905 rows x 3 columns
* Variables:
 * name: the name of the site (city)
 * datetime: the date and time when the temperature is recorded.
 * temp: temperature (degrees Celsius)


 Details for: beijing 2016-11-01 to 2017-02-28.csv
---------------------------------------
* Description: the hourly temperature data in Beijing, N 40.07°, E 116.58° between the period from 2016-11-01 to 2017-02-28, which was downloaded from a global weather API (https://www.visualcrossing.com/weather-data) in May 2022.

* Format(s): .csv
* Size(s):  93 KB
* Dimensions: 2905 rows x 3 columns
* Variables:
 * name: the name of the site (city)
 * datetime: the date and time when the temperature is recorded.
 * temp: temperature (degrees Celsius)

- - -
END OF README


