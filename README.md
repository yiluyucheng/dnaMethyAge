# dnaMethAge
Predict epigenetic age from DNA methylation data, for example, Horvath age, Hannum age, PhenoAge(Levine).

## Usage

#### Installation
git clone git@github.com:yiluyucheng/dnaMethAge.git

#### How to use

1. Change the work directory into /dnaMethAge/scripts

2. Open R environment

* 2.1 Predict Horvath age:

```R
source('MethAge.R')

res <- methyAge(betas, clock='Horvath2013')
```

The *betas* should be a dataframe which samples in the columns and probe in the rows.

Currently, support age clocks are 'Hannum2013', 'Horvath2013', 'Levine2018'. More age models will be included in the future.


* 2.2 Predict Horvath age and calculate age acceleration:

```
source('MethAge.R')

res <- methyAge(betas, clock='Horvath2013', age_info=info, fit_method='Linear')
```

*info* should be a dataframe which contains sample ID and age information, like:

Sample | Age
-------- | -----
name_1 | 31
name_2 | 60
name_3 | 42
... | ...
name_n | 53


Please also check '/scripts/predict_on_betas.R' to get better known as how to use methyAge.

## Contact me

yw19282@essex.ac.uk



