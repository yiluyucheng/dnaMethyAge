# dnaMethyAge
A user friendly **R package** to predict epigenetic age and calculate age acceleration from DNA methylation data. Supported age clocks are Horvath's clock, Hannum's clock, Levine's PhenoAge, and Zhang's clock.

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
print(dim(betas))
# 485577 85

hannum_age <- methyAge(betas, clock='Hannum2013')

## Directly use Horvath's clock without adjusted-BMIQ normalisation 
horvath_age <- methyAge(betas, clock='Horvath2013', fast_mode=TRUE)
## Use Horvath's clock with adjusted-BMIQ normalisation (same as Horvath's paper)

horvath_age <- methyAge(betas, clock='Horvath2013')

pheno_age <- methyAge(betas, clock='Levine2018')

zhang_age <- methyAge(betas, clock='Zhang2019')

```

The above variable **betas** should be a dataframe which samples in the columns and probe in the rows.

Currently, supported age clocks are 'Hannum2013', 'Horvath2013', 'Levine2018', 'Zhang2019'. More age models will be added in the future.


#### (2)  Predict epigenetic age and calculate age acceleration

```R
library('dnaMethyAge')

## prepare betas dataframe
print(dim(betas))
# 485577 85

## prepare age_info
print(dim(info))
# 85 2

# Apply Hannum's clock and calculate age acceleration
hannum_age <- methyAge(betas, clock='Hannum2013', age_info=info)

# Apply Horvath's clock and calculate age acceleration
horvath_age <- methyAge(betas, clock='Horvath2013', age_info=info, fit_method='Linear')

# Apply Levine's PhenoAge and calculate age acceleration
pheno_age <- methyAge(betas, clock='Levine2018', age_info=info)

# Apply Zhang's clock and calculate age acceleration
zhang_age <- methyAge(betas, clock='Zhang2019', age_info=info, fit_method='Linear')
```

The above variable **info** should be a dataframe which contains sample ID and age information, like:

Sample | Age
-------- | -----
name_1 | 31
name_2 | 60
name_3 | 42
... | ...
name_n | 53


Please refer to the below code to learn more about how to use the method.
```R
library('dnaMethyAge')

help(methyAge)
```

## 2. Reference

Hannum2013: Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. Hannum et al., 2013

Horvath2013: DNA methylation age of human tissues and cell types. Steve Horvath, 2013

Levine2018: An epigenetic biomarker of aging for lifespan and healthspan. Levine et al., 2018

Zhang2019: Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing, Zhang et al., 2019

## 3. Contact me

yw19282@essex.ac.uk



