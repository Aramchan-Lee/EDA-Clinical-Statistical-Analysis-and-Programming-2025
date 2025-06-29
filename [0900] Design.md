# [0900] Kaplan-Meier Log-rank test and Cox-PH: Design 

---

- [Input Dataset](#Input-Dataset)
    - [Content](#Content)
    - [Key Variables](#Key-Variables)
    - [Rationale for Use](#Rationale-for-Use)
- [Numeric Agreement Criteria](#Numeric-Agreement-Criteria)
    - [Kaplan-Meier Estimator](#Kaplan-Meier-Estimator-Crit)
    - [Log-rank Test](#Log-rank-Test-Crit)
    - [Cox Proportional Hazards Model](#Cox-Proportional-Hazards-Model-Crit)
- [Expectations for Supported Statistics](#Expectations-for-Supported-Statistics)
    - [SAS (LIFETEST)](#SAS-LIFETEST-Ex)
    - [SAS (PHREG)](#SAS-PHREG-Ex)
    - [R (SURVIVAL)](#R-SURVIVAL-Ex)
    - [Python (LIFELINES)](#Python-LIFELINES-Ex)
- [Critical Arguments](#Arguments)
    - [SAS (LIFETEST)](#SAS-LIFETEST-Arg)
    - [SAS (PHREG)](#SAS-PHREG-Arg)
    - [R (SURVIVAL)](#R-SURVIVAL-Arg)
    - [Python (LIFELINES)](#Python-LIFELINES-Arg)
- [Known Incompatabilities](#Incompatabilities)
    - [SAS (LIFETEST)](#SAS-LIFETEST-Incomp)
    - [SAS (PHREG)](#SAS-PHREG-Incomp)
    - [R (SURVIVAL)](#R-SURVIVAL-Incomp)
    - [Python (LIFELINES)](#Python-LIFELINES-Incomp)
- [Comparison Protocol](#Comparison-Protocol)
    - [Kaplan-Meier Estimator](#Kaplan-Meier-Estimator-Prot)
    - [Log-rank Test](#Log-rank-Test-Prot)
    - [Cox Proportional Hazards Model](#Cox-Proportional-Hazards-Model-Prot)


---

# Input Dataset

The input dataset for this validation is the **lung** dataset, which is included as a built-in example in the R survival package.

- This dataset originates from the North Central Cancer Treatment Group (NCCTG) study of advanced lung cancer patients and has been widely used in textbooks, tutorials, and comparative studies of survival analysis methods.

- Notably, the CAMIS Survival Analysis Comparison project used this dataset to validate and demonstrate consistency between R and SAS implementations of Kaplan-Meier estimation, log-rank testing, and Cox proportional hazards regression.

- By using the same dataset, this project ensures that the Python implementation (lifelines) can be directly compared against both R and SAS outputs with known reference results.

## Content

The lung dataset contains clinical trial data recording time to death or censoring, along with patient characteristics such as age, sex, and performance status.
It includes censored observations (patients who were still alive at the end of the observation period) and allows testing of all major survival methods:

- Kaplan-Meier survival function estimation

- Group comparisons via log-rank test

- Cox proportional hazards regression modeling

## Key Variables

| Variable      | Description                                                    |
|---------------|----------------------------------------------------------------|
| `time`        | Survival time in days                                          |
| `status`      | Event indicator (1=censored, 2=dead). For modeling, recoded as 0/1. |
| `age`         | Patient age in years                                           |
| `sex`         | Sex of patient (1=male, 2=female)                              |
| `ph.ecog`     | ECOG performance status                                        |
| `ph.karno`    | Karnofsky score (physician-rated)                              |
| `pat.karno`   | Karnofsky score (patient-rated)                                |
| `meal.cal`    | Average daily calorie intake                                   |
| `wt.loss`     | Weight loss in pounds over the last 6 months                   |


## Rationale for Use

- The dataset is public, reproducible, and commonly accepted in examples.

- It is small enough for rapid testing across platforms.

- It has clear expected outputs published in CAMIS, supporting numeric agreement validation.

# Numeric Agreement Criteria

Numeric agreement criteria define the acceptable range of difference between results produced by the new Python implementation and reference implementations (SAS and R). They are set to account for expected minor variation due to computational precision and method defaults, while ensuring statistical and clinical equivalence of the outputs.

## Kaplan-Meier Estimator <span style="display:none">(Crit)</span>

- Survival probability estimates should agree within 0.001 to account for minor differences in step function computation and rounding.

- Median survival times should agree within 0.01.
 
## Log-rank Test <span style="display:none">(Crit)</span>

Test statistic and p-values should agree within 0.001, considering SAS’s default rounding to four decimal places and ensuring consistency across tools.
  
## Cox Proportional Hazards Model <span style="display:none">(Crit)</span>

- Estimated coefficients (β) should agree within 0.01, allowing for algorithmic variation while maintaining interpretation consistency. Minor deviations in coefficients typically do not alter interpretation.

- Hazard Ratios (HRs) and their confidence intervals should agree within 0.01.
 
- p-values should agree within 0.001 to ensure reproducible significance levels.

# Expectations for Supported Statistics

## SAS LIFETEST <span style="display:none">(Ex)</span>

- Outputs: Survival estimates, Standard Errors, Confidence Intervals, Number at Risk, Log-Rank Test Statistic, Degrees of Freedom, P-value (Formatted as Tables)

- Uses: PROC LIFETEST for Kaplan-Meier estimation and Log-Rank test for comparing survival curves

## SAS PHREG <span style="display:none">(Ex)</span>

- Outputs: Parameter Estimates (Coefficients), Standard Errors, Hazard Ratios, Confidence Intervals, Likelihood Ratio / Wald / Score Test Statistics, P-values (Formatted as Tables)

- Uses: PROC PHREG for Cox Proportional Hazards Regression Model

## R SURVIVAL <span style="display:none">(Ex)</span>

- Outputs: Fitted Model Object (with Coefficients, Standard Errors, Hazard Ratios, p-values), Summary Tables for Kaplan-Meier curves or Cox Models (Formatted as List or Data Frames)

- Uses:
    - survfit() for Kaplan-Meier estimates
    - coxph() for Cox Proportional Hazards Model

## Python LIFELINES <span style="display:none">(Ex)</span>

- Outputs: Summary DataFrame with Coefficients, Standard Errors, Hazard Ratios, Confidence Intervals, P-values (Formatted as pandas DataFrame)

- Uses: 
    - KaplanMeierFitter() for survival curves
    - CoxPHFitter() for Cox Proportional Hazards Regression

# Critical Arguments

## SAS LIFETEST <span style="display:none">(Arg)</span>

- data: dataset containing survival times, event indicator, and optionally grouping variable

- time: time *event(censoring code): specifies survival time and censoring status (e.g., time*status(0) where 0 = censored)

- strata: grouping variable to compare survival curves (e.g., treatment group)

## SAS PHREG <span style="display:none">(Arg)</span>

- data: dataset containing survival times, event indicator, predictors

- model: time *event(censoring code) = predictors: specifies survival time, censoring, and covariates (e.g., time*status(0) = age sex treatment)

- strata: optional variable to stratify baseline hazards

## R SURVIVAL <span style="display:none">(Arg)</span>

- formula: Surv(time, status) ~ predictors: specifies survival time, event indicator, and model predictors

    - note: for Kaplan-Meier estimation, predictors is typically a grouping variable, e.g., treatment group

- data: data frame containing variables used in the formula

- subset: optional logical vector to specify a subset of rows

## Python LIFELINES <span style="display:none">(Arg)</span>

- KaplanMeierFitter().fit()
    - durations: array-like of survival times

    - event_observed: array-like indicator of event occurrence (1=event, 0=censored)

    - label: optional label for the curve

- CoxRegressionFitter().fit()
    - df: pandas DataFrame with all columns (time, event, covariates)

    - duration_col: column name for survival time

    - event_col: column name for event indicator

# Known Incompatibilities

## SAS LIFETEST <span style="display:none">(Incomp)</span>

### Output Format

- The output is primarily formatted as SAS tables and printed output, which is less convenient for direct programmatic extraction compared to R or Python.

### Model Complexity

- Only supports univariate survival curve estimation (Kaplan-Meier) and stratified comparisons—cannot fit multivariate Cox models (requires PROC PHREG).

## SAS PHREG <span style="display:none">(Incomp)</span>

### Tied Events Handling

- By default, ties are handled using the Breslow approximation; handling ties exactly (Efron method) requires explicit options, which may differ from R defaults.

### Output Format

- Model results are in printed output and ODS tables; extra steps are needed to export results to data sets or CSV.

## R SURVIVAL <span style="display:none">(Incomp)</span>

### Output Complexity

- The fitted objects (survfit, coxph) are S3 objects that may require additional steps (e.g., summary(), as.data.frame()) to extract human-readable summaries.

### Defaults

- Different default tie handling in Cox models (efron in R vs. breslow in SAS), which can lead to small discrepancies.

## Python LIFELINES <span style="display:none">(Incomp)</span>

### Input Format

- Requires Pandas DataFrames or arrays; less forgiving of raw lists or other data types compared to R formula interfaces.

### Limited Multistate Modeling

- Does not natively support multi-state models or more advanced survival modeling options (e.g., time-varying covariates require extra work).

### Output

- While outputs are DataFrames, formatting and exporting often need additional manual steps compared to R’s tidyverse integration.

# Comparison Protocol

## Kaplan-Meier Estimator <span style="display:none">(Prot)</span>

- Survival probabilities at key time points  
- Median survival time  
- Number at risk  
- Confidence interval bounds (e.g., lower and upper CI)

## Log-rank Test <span style="display:none">(Prot)</span>

- Chi-square test statistic  
- p-value  
- Degrees of freedom (df)

## Cox Proportional Hazards Model <span style="display:none">(Prot)</span>

- Coefficient estimates (β)  
- Wald test statistic and p-values for each coefficient
- Hazard ratios (HR = exp(β))  
- Standard errors of β  
- 95% confidence intervals for coefficients or HRs
- Log-likelihood and AIC  
- p-values from Likelihood Ratio tests
