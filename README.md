# dnaMethyAge
A user friendly **R package** to predict epigenetic age and calculate age acceleration from DNA methylation data. Currently, supported age clocks are listed below:

| Name | Published Name | First Author (Published Year) | Trained Phenotype | Num. of CpGs | Tissues Derived |
| ------------- | ------------- |------------- | ------------- |------------- |------------- |
| HannumG2013 | | [Gregory Hannum (2013)](https://www.sciencedirect.com/science/article/pii/S1097276512008933) | Chronological age | 71 | Whole blood|
| HorvathS2013 | Multi-tissue age estimator| [Steve Horvath (2013)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115) | Chronological age | 353 | 27 tissues/cells |
| YangZ2016 | epiTOC | [Zhen Yang (2016)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1064-3)| Mitotic divisions | 385 | Whole blood |
| ZhangY2017 | | [Yan Zhang (2017)](https://www.nature.com/articles/ncomms14617) | Mortality risk | 10 | Whole blood |
| HorvathS2018 | Skin & Blood Clock | [Steve Horvath (2018)](https://www.aging-us.com/article/101508/text) | Chronological age | 391 | Skin, blood, buccal cells and 5 other tissues
| LevineM2018 | PhenoAge | [Morgan E. Levine (2018)](https://www.aging-us.com/article/101414/text) | Mortality risk | 513 | Whole blood |
| McEwenL2019 | PedBE | [Lisa M. McEwen (2019)](https://www.pnas.org/doi/10.1073/pnas.1820843116) | Chronological age | 94 | Buccal epithelial cells |
| ZhangQ2019 | | [Qian Zhang (2019)](https://link.springer.com/article/10.1186/s13073-019-0667-1) | Chronological age | 514 | Whole blood |
| LuA2019 | DNAmTL | [Ake T. Lu (2019)](https://www.aging-us.com/article/102173/text) | Leukocyte telomere length | 140 | Whole blood |
| epiTOC2 | epiTOC2 | [Andrew E. Teschendorff (2020)](https://doi.org/10.1186/s13073-020-00752-3) | Mitotic divisions | 163 | Whole blood |
| ShirebyG2020 | Cortical clock | [Gemma L Shireby (2020)](https://academic.oup.com/brain/article/143/12/3763/5942151?login=true) | Chronological age | 347 | Brain cortex |
| DunedinPACE | DunedinPACE | [Daniel W Belsky (2022)](https://elifesciences.org/articles/73420) | Pace of ageing | 173 | Blood |
| PCGrimAge | PCGrimAge | [Albert T. Higgins-Chen (2022)](https://doi.org/10.1038/s43587-022-00248-2) | [GrimAge](https://doi.org/10.18632/aging.101684) estimated from DNAm data| 78464 | Blood |
| BernabeuE2023c | cAge | [Elena Bernabeu (2023)](https://doi.org/10.1186/s13073-023-01161-y) | Chronological age |3225 | Blood |
 
## 1. Usage

### 1.1 Installation

**Install from Github**
```R
## Make sure 'devetools' is installed in your R
# install.packages("devtools")
devtools::install_github("yiluyucheng/dnaMethyAge")
```

### 2.1 How to use

Start a R work environment

#### (1) Predict epigenetic age from DNA methylation data

```R
library('dnaMethyAge')

## prepare betas dataframe
data('subGSE174422') ## load example betas

print(dim(betas)) ## probes in row and samples in column
# 485577 8

availableClock() ## List all supported clocks
# "HannumG2013"  "HorvathS2013" "LevineM2018"  "ZhangQ2019"   "ShirebyG2020"  "YangZ2016"    "ZhangY2017"

clock_name <- 'HorvathS2013'  # Select one of the supported clocks.
## Use Horvath's clock with adjusted-BMIQ normalisation (same as Horvath's paper)
horvath_age <- methyAge(betas, clock=clock_name)

print(horvath_age)
#                         Sample     mAge
# 1 GSM5310260_3999979009_R02C02 74.88139
# 2 GSM5310261_3999979017_R05C01 62.36400
# 3 GSM5310262_3999979018_R02C02 68.04759
# 4 GSM5310263_3999979022_R02C01 61.62691
# 5 GSM5310264_3999979027_R02C01 59.65161
# 6 GSM5310265_3999979028_R01C01 60.95991
# 7 GSM5310266_3999979029_R04C02 52.48954
# 8 GSM5310267_3999979031_R06C02 64.29711

```

More age models will be added in the future, please get contact if you would like me to add a new clock.


#### (2)  Predict epigenetic age and calculate age acceleration

```R
library('dnaMethyAge')

## prepare betas dataframe
data('subGSE174422') ## load example betas and info

print(dim(betas)) ## probes in row and samples in column
# 485577 8
print(info) ##  info should be a dataframe which includes at least two columns: Sample, Age.
#                         Sample  Age    Sex
# 1 GSM5310260_3999979009_R02C02 68.8 Female
# 2 GSM5310261_3999979017_R05C01 45.6 Female
# 3 GSM5310262_3999979018_R02C02 67.4 Female
# 4 GSM5310263_3999979022_R02C01 45.6 Female
# 5 GSM5310264_3999979027_R02C01 62.5 Female
# 6 GSM5310265_3999979028_R01C01 45.1 Female
# 7 GSM5310266_3999979029_R04C02 53.2 Female
# 8 GSM5310267_3999979031_R06C02 63.8 Female


clock_name <- 'HorvathS2013'  # Select one of the supported clocks, try: availableClock()
## Apply Horvath's clock and calculate age acceleration
## Use Horvath's clock with adjusted-BMIQ normalisation (same as Horvath's paper)
horvath_age <- methyAge(betas, clock=clock_name, age_info=info, fit_method='Linear', do_plot=TRUE)

print(horvath_age)
#                         Sample  Age    Sex     mAge Age_Acceleration
# 1 GSM5310260_3999979009_R02C02 68.8 Female 74.88139         7.334461
# 2 GSM5310261_3999979017_R05C01 45.6 Female 62.36400         3.318402
# 3 GSM5310262_3999979018_R02C02 67.4 Female 68.04759         1.013670
# 4 GSM5310263_3999979022_R02C01 45.6 Female 61.62691         2.581311
# 5 GSM5310264_3999979027_R02C01 62.5 Female 59.65161        -5.586763
# 6 GSM5310265_3999979028_R01C01 45.1 Female 60.95991         2.097534
# 7 GSM5310266_3999979029_R04C02 53.2 Female 52.48954        -9.340977
# 8 GSM5310267_3999979031_R06C02 63.8 Female 64.29711        -1.417638
```

By default, methyAge would plot the age prediction results and the distribution of age acceleration, to save the plot:
```R
pdf('savename.pdf', width=4.3, height=6)
horvath_age <- methyAge(betas, clock=clock_name, age_info=info, fit_method='Linear', do_plot=TRUE)
dev.off()
```
Here is the result plot(very nice!):

<img width="320" src="https://github.com/yiluyucheng/dnaMethyAge/blob/main/test_res/test_result1.png">

Now you can try other clocks by simply redefine the 'clock_name' and keep other codes unedited. Normally, the clock of Horvath2013 costs the highest amount of time due to its unefficient normalisation steps, for other clocks such as Hannum2013, the overall running time is very short. Please refer to the below code to learn more about how to use the method.
```R
library('dnaMethyAge')

help(methyAge)
```


Lastly, I provide the epiginetic age prediction results in four clocks for [GSE147221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147221) in below(click on the image to have a more clear look). 

<img width="960" src="https://github.com/yiluyucheng/dnaMethyAge/blob/main/test_res/test_result2.png">

## 2. Attention
The original model of Hannum2013 not only uses the 72 CpG sites, but also includes covariates gender, BMI, diabetes status, ethnicity and batch, howevever, the authors did not provide the coefficients for those covariates. In this method, I only used the 72 CpG sites to calculate the Hannum2013 age, and this is the commom practice in this field.

The four clocks' prediciton performance may vary in different datasets, and the Levine2018 also known as PhenoAge was not directly trained on chronological age.

## 3. Citation
If you used this package in your research, please cite us:
[Wang et al., 2023](https://doi.org/10.1007/s11357-023-00871-w)
```
@article{Wang2023,
  title={Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock},
  author={Wang, Yucheng and Grant, Olivia A and Zhai, Xiaojun and Mcdonald-Maier, Klaus D and Schalkwyk, Leonardo C},
  journal={GeroScience},
  pages={1--18},
  year={2023},
  publisher={Springer}
}
```


## 4. Contact me

yw19282@essex.ac.uk



