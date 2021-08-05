# dnaMethAge
Predict epigenetic age from DNA methylation data, for example, Horvath age, Hannum age, PhenoAge(Levine), Zhang2019.

## Usage

#### Installation
```
git clone git@github.com:yiluyucheng/dnaMethAge.git
```

or download the ZIP file(dnaMethAge-main.zip), after unzip the file.

#### How to use

1. Change the work directory into /dnaMethAge/scripts/

2. Open R environment

* 2.1 Predict Horvath age:

```R
source('MethAge.R')

horvath_age <- methyAge(betas, clock='Horvath2013')
hannum_age <- methyAge(betas, clock='Hannum2013')
pheno_age <- methyAge(betas, clock='Levine2018')
zhang_age <- methyAge(betas, clock='Zhang2019')

```

The $betas$ should be a dataframe which samples in the columns and probe in the rows.

Currently, support age clocks are 'Hannum2013', 'Horvath2013', 'Levine2018', 'Zhang2019'. 

More age models will be included in the future.


* 2.2 Predict Horvath age and calculate age acceleration:

```
source('MethAge.R')

horvath_age <- methyAge(betas, clock='Horvath2013', age_info=info, fit_method='Linear')
hannum_age <- methyAge(betas, clock='Hannum2013', age_info=info, fit_method='Linear')
pheno_age <- methyAge(betas, clock='Levine2018', age_info=info, fit_method='Linear')
zhang_age <- methyAge(betas, clock='Zhang2019', age_info=info, fit_method='Linear')
```

$info$ should be a dataframe which contains sample ID and age information, like:

Sample | Age
-------- | -----
name_1 | 31
name_2 | 60
name_3 | 42
... | ...
name_n | 53


Please also check '/scripts/predict_on_betas.R' to get better known as how to use methyAge.

#### Reference

Hannum2013: Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. Hannum et al., 2013

Horvath2013: DNA methylation age of human tissues and cell types. Steve Horvath, 2013

Levine2018: An epigenetic biomarker of aging for lifespan and healthspan. Levine et al., 2018

Zhang2019: Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing, Zhang et al., 2019

## Contact me

yw19282@essex.ac.uk



