---
title: '[0900] Kaplan-Meier Log-rank test and Cox-PH - Requirements Analysis'

---


# [0900] Kaplan-Meier Log-rank test and Cox-PH: Requirements Analysis

---

- [Background](#Background)
  - [Hypotheses](#Hypotheses)
  - [Assumptions](#Assumptions)
  - [Test Statistic & P-value](#Test-Statistic)
- [Package Implementations](#Package-Implementations)
  - [SAS (LIFETEST)](#SAS-LIFETEST)
    - [Function/Procedure](#Function-Lifetest)
    - [Inputs](#Inputs-Lifetest)
    - [Outputs](#Outputs-Lifetest)
    - [Sample Code](#Sample-Code-Lifetest)
    - [Limitations](#Limitations-Lifetest)
  - [SAS (PHREG)](#SAS-PHREG)
    - [Function/Procedure](#Function-PHREG)
    - [Inputs](#Inputs-PHREG)
    - [Outputs](#Outputs-PHREG)
    - [Sample Code](#Sample-Code-PHREG)
    - [Limitations](#Limitations-PHREG)
  - [R (SURVIVAL)](#R-SURVIVAL)
    - [Function/Procedure](#Function-Survival)
    - [Inputs](#Inputs-Survival)
    - [Outputs](#Outputs-Survival)
    - [Sample Code](#Sample-Code-Survival)
    - [Limitations](#Limitations-Survival)
  - [Python (LIFELINES)](#Python-LIFELINES)
    - [Function/Procedure](#Function-Lifelines)
    - [Inputs](#Inputs-Lifelines)
    - [Outputs](#Outputs-Lifelines)
    - [Sample Code](#Sample-Code-Lifelines)
    - [Limitations](#Limitations-Lifelines)
- [Summary](#Summary)
- [Extra Resources](#Extra-Resourses)


---


# Background

Survival analysis encompasses a set of statistical methods designed to analyze time-to-event data—where the outcome is the time until the occurrence of an event of interest, such as death, relapse, equipment failure, or recovery. These methods are uniquely capable of handling censored data, which occurs when the event has not been observed for all subjects during the observation period.

Three fundamental tools in survival analysis are:

- Kaplan-Meier Estimator: A non-parametric method that estimates the survival function from lifetime data. It calculates the probability of surviving past certain time points and accounts for right-censored observations. It is particularly useful for visualizing survival curves and comparing them descriptively.

- Log-rank Test: A hypothesis test used to formally compare two or more Kaplan-Meier survival curves. It evaluates whether there is a statistically significant difference in the survival distributions between groups. The test is based on the observed versus expected number of events over time and assumes proportional hazards.

- Cox Proportional Hazards Model (Cox-PH): A semi-parametric regression model that estimates the hazard (or risk) of an event occurring, adjusting for covariates. It does not assume a specific baseline hazard function but assumes that hazard ratios between groups remain constant over time (proportional hazards). This model is powerful for examining the impact of continuous and categorical covariates on survival.

These methods are widely used in clinical trials, epidemiology, engineering, and other fields where time-to-event data is critical. Their ability to accommodate censoring, handle covariates, and provide interpretable results make them essential tools in both exploratory and confirmatory data analysis to estimate the probability of survival past a certain point in time.

## Hypotheses

<ins>**Log-rank test:**</ins>

  - H₀: Survival functions are equal across groups.

  - H₁: At least one group differs.

<ins>**Cox-PH:**</ins>

  - H₀: Hazard ratio for each covariate = 1 (no effect).

  - H₁: At least one covariate has HR ≠ 1.
  
## Assumptions

<ins>**Kaplan-Meier (KM)**</ins>

  - Independent censoring: Censored subjects are assumed to have the same survival prospects as those who remain under observation.

  - Non-informative censoring: Censoring is unrelated to the event risk.

  - Exact event times: Event times are known and measured without error.

  - Random sampling: Subjects are representative of the population of interest.

<ins>**Log-Rank Test**</ins>

  - Proportional hazards: Hazard ratios between groups are constant over time.

  - Independent censoring: Same as in KM.
  
  - Independent survival times: Between individuals and across groups.

  - Correct group assignment: Every individual belongs to one (and only one) group.

<ins>**Cox Proportional Hazards (Cox-PH)**</ins>

  - Proportional hazards: The hazard ratio between individuals is constant over time.

  - Linearity of covariates in log-hazard: The effect of continuous variables is linear on the log hazard scale.

  - Independent survival times: Survival times are independent between individuals.

  - Independent censoring: Censoring is unrelated to survival probability.

  - No omitted confounders: All relevant predictors are included in the model (not formally testable but crucial for valid inference).

## Test Statistic

<ins>**Kaplan-Meier (KM):**</ins>

  The Kaplan-Meier survival estimate at time t is:

    S(t) = product over all tᵢ ≤ t of [1 - dᵢ / nᵢ]

    where:

      dᵢ = number of events at time tᵢ

      nᵢ = number at risk just before tᵢ

    Variance of the survival estimate (Greenwood's formula):

      Var[S(t)] ≈ S(t)² × sum over tᵢ ≤ t of [dᵢ / (nᵢ × (nᵢ - dᵢ))]

<ins>**Log-rank Test:**</ins>

  The test statistic is:

    Q = (O - E)² / Var(O)

    where:

      O = observed number of events

      E = expected number of events under the null hypothesis

    Var(O) = variance of observed events under the null

    Q follows a chi-squared distribution with degrees of freedom = (number of groups - 1)

<ins>**Cox Proportional Hazards (Cox-PH):**</ins>

  The partial likelihood for estimating regression coefficients:

    L(β) = product over all events i of [exp(βXᵢ) / sum over j in R(tᵢ) of exp(βXⱼ)]

    where R(tᵢ) is the risk set at time tᵢ

  The log partial likelihood is:

    log L(β) = sum over i of [βXᵢ - log sum over j in R(tᵢ) of exp(βXⱼ)]

  Common test statistics:

    Wald test: Z = β̂ / SE(β̂)

    Likelihood Ratio Test (LRT): Compare log-likelihoods of nested models

    Score test: Uses derivative of the log-likelihood at β = 0

    Hazard Ratio (HR): HR = exp(β̂)

# Package Implementations

Dependencies:
All methods rely only on standard, widely available components in their respective environments:

- SAS procedures require no additional modules beyond SAS/STAT.

- R functions use only base packages (e.g., survival).

- Python (lifelines) dependencies (numpy, pandas, scipy, matplotlib, patsy) are installed automatically via pip.

## SAS LIFETEST

  ### Function <span style="display:none">(Lifetest)</span>
  
  - Use PROC LIFETEST for Kaplan-Meier survival curves and log-rank test.
    
  ### Inputs <span style="display:none">(Lifetest)</span>
  <ins>**Required:**</ins>
  
  - DATA= (input dataset), wide format
  
  - TIME statement for survival time and censoring
  
  <ins>**Optional**:</ins>
  
  - STRATA= (for group comparison)
  
  - PLOTS= (e.g. survival)
  
  - TEST (equality test)

  ### Outputs <span style="display:none">(Lifetest)</span>
  
  - Kaplan-Meier curves
  
  - Log-rank test statistic & p-value
  
  - Median survival & CI

  ### Sample Code <span style="display:none">(Lifetest)</span>
  ```sas
  proc lifetest data=study plots=survival;
    time time*status(0);
    strata treatment;
  run;
  ```

  ### Limitations <span style="display:none">(Lifetest)</span>
  - No built-in test for proportional hazards; requires user-driven diagnostic procedures.

  - Default method for handling ties is Efron's, which differs from some R/Python defaults.

  - Cannot handle time-varying covariates without macro extensions.

  - Limited support for custom baseline hazard specification.

## SAS PHREG    

  ### Function <span style="display:none">(PHREG)</span>
  
  - Use PROC PHREG from the SAS STAT module for Cox Proportional Hazards modeling.

  ### Inputs <span style="display:none">(PHREG)</span>
    
  <ins>**Required:**</ins>
  
  - DATA= (input dataset), wide format
  
  - MODEL (time and censoring indicator, covariates)
  
  <ins>**Optional:**</ins>
  
  - CLASS= (categorical predictors)
  
  - STRATA= (stratify baseline hazards)
  
  - HAZARDRATIO= (requests hazard ratios)
  
  - ID= (subject identifier)
  
  ### Outputs <span style="display:none">(PHREG)</span>
  
  - Hazard ratios (HR)
  
  - Standard errors and confidence intervals
  
  - Wald test, LRT statistics
  
  - Diagnostic outputs (residuals, influence stats)
    
### Sample Code <span style="display:none">(PHREG)</span>
  ```sas
    proc phreg data=study;
      model time*status(0) = age sex treatment;
    run;
  ```

### Limitations <span style="display:none">(PHREG)</span>  
  
  - Cannot incorporate covariates—limited to stratified or grouped comparisons.
  
  - Output can be verbose and require manual post-processing for reporting.
  
  - Graphical outputs depend on enabling ODS graphics manually.
  
## R SURVIVAL

  ### Function <span style="display:none">(Survival)</span>
  - Use survfit(), survdiff(), and coxph() from the survival package.
  
  ### Inputs <span style="display:none">(Survival)</span>
  <ins>**Required:**</ins>
  
  - Surv(time, status) object
  
  - Formula input for groups/covariates
  
  <ins>**Optional:**</ins>
  
  - cluster, weights, ties, na.action
  
  ### Outputs <span style="display:none">(Survival)</span>
  
  - survfit: KM estimates
  
  - survdiff: log-rank p-value
  
  - coxph: HRs, CIs, test stats
    
### Sample Code <span style="display:none">(Survival)</span>
  ```r
  library(survival)
  fit_km <- survfit(Surv(time, status) ~ group, data = df)
  survdiff(Surv(time, status) ~ group, data = df)
  fit_cox <- coxph(Surv(time, status) ~ age + treatment, data = df)
  summary(fit_cox)
  ```
### Limitations <span style="display:none">(Survival)</span>

  - Diagnostics (e.g., proportional hazards via cox.zph) must be explicitly run and interpreted.
  
  - Default handling of ties is Efron's method; Breslow available, but not default.
  
  - Formulas can become complex when many interactions or strata are involved.
  
## Python LIFELINES

  ### Function <span style="display:none">(Lifelines)</span>
  - Use KaplanMeierFitter, logrank_test, and CoxPHFitter from the lifelines package.

  ### Inputs <span style="display:none">(Lifelines)</span>
  <ins>Required:**</ins>
  
  - duration_col, event_col
  
  - Pandas DataFrame with covariates
  
  <ins>**Optional:**</ins>
  
  - robust, step_size, weights, formula

  ### Outputs <span style="display:none">(Lifelines)</span>
  
  - KM survival estimates
  
  - Log-rank p-values
  
  - HRs, standard errors, confidence intervals

  ### Sample Code <span style="display:none">(Lifelines)</span>
  ```python
  from lifelines import KaplanMeierFitter, CoxPHFitter
  from lifelines.statistics import logrank_test
  
  kmf = KaplanMeierFitter()
  kmf.fit(df['time'], df['event'])
  kmf.plot()
  
  logrank_test(df_A['time'], df_B['time'], df_A['event'], df_B['event']).print_summary()
  
  cph = CoxPHFitter()
  cph.fit(df, duration_col='time', event_col='event')
  cph.print_summary()
  ```
  
### Limitations <span style="display:none">(Lifelines)</span>

  - Assumption checking (e.g., proportional hazards) must be performed via .check_assumptions(); not automatic.
  
  - Results may vary depending on tie-breaking method used (Efron, Breslow, or exact).
  
  - Does not support time-dependent covariates as robustly as R.
  
  - Some plotting features require knowledge of Matplotlib.
  
# Summary

- Kaplan-Meier, Log-rank, and Cox Proportional Hazards (Cox-PH) are foundational tools in survival analysis, each suited to different analytic needs. Kaplan-Meier provides non-parametric survival estimates, Log-rank tests detect differences between survival curves, and Cox-PH models estimate the impact of covariates on time-to-event outcomes under the proportional hazards assumption. These methods are supported in validated statistical packages (SAS, R, Python) and widely accepted in regulatory settings, with varying strengths and limitations in implementation, flexibility, and diagnostics.

# Extra Resourses

**Regulatory and Industry References**

  - FDA Statistical Guidance for Clinical Trials: Discusses statistical principles in clinical trials; includes KM and Cox-PH.
    - https://www.fda.gov/media/71145/download

  - ICH E9 Guideline: Standardizes survival analysis approaches internationally.
    - https://www.ema.europa.eu/en/documents/scientific-guideline/ich-e-9-statistical-principles-clinical-trials-step-5_en.pdf

  - SAS Documentation: Procedures PROC LIFETEST and PROC PHREG used in FDA-reviewed submissions. SAS Lifetest, SAS PHREG
    - https://documentation.sas.com/doc/en/pgmsascdc/v_063/statug/statug_lifetest_overview.htm

  - R Survival Package: Official reference from CRAN, widely used in academic and industry publications. CRAN Survival
    - https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf

  - Python Lifelines Documentation: Supports accepted statistical methods for open-science usage. Lifelines Docs
    - https://lifelines.readthedocs.io/en/latest/

  - EMA Biostatistics Guidance: Endorses Cox and stratified log-rank for regulatory use.
    - https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-adjustment-covariates-randomised-clinical-trials_en.pdf

**Notes and SOPs**

  - UCLA SAS Survival Seminar
    - https://stats.oarc.ucla.edu/sas/seminars/sas-survival/

  - UC Davis Cox Model Presentation
    - https://health.ucdavis.edu/media-resources/ctsc/documents/pdfs/cph-model-presentation.pdf

  - FDA ICH E9 Statistical Principles
    - https://www.fda.gov/regulatory-information/search-fda-guidance-documents/e9-statistical-principles-clinical-trials

# [0900] Kaplan-Meier Log-rank test and Cox-PH: Requirements Analysis

---

- [Background](#Background)
  - [Hypotheses](#Hypotheses)
  - [Assumptions](#Assumptions)
  - [Test Statistic & P-value](#Test-Statistic)
- [Package Implementations](#Package-Implementations)
  - [SAS (LIFETEST)](#SAS-LIFETEST)
    - [Function/Procedure](#Function-Lifetest)
    - [Inputs](#Inputs-Lifetest)
    - [Outputs](#Outputs-Lifetest)
    - [Sample Code](#Sample-Code-Lifetest)
    - [Limitations](#Limitations-Lifetest)
  - [SAS (PHREG)](#SAS-PHREG)
    - [Function/Procedure](#Function-PHREG)
    - [Inputs](#Inputs-PHREG)
    - [Outputs](#Outputs-PHREG)
    - [Sample Code](#Sample-Code-PHREG)
    - [Limitations](#Limitations-PHREG)
  - [R (SURVIVAL)](#R-SURVIVAL)
    - [Function/Procedure](#Function-Survival)
    - [Inputs](#Inputs-Survival)
    - [Outputs](#Outputs-Survival)
    - [Sample Code](#Sample-Code-Survival)
    - [Limitations](#Limitations-Survival)
  - [Python (LIFELINES)](#Python-LIFELINES)
    - [Function/Procedure](#Function-Lifelines)
    - [Inputs](#Inputs-Lifelines)
    - [Outputs](#Outputs-Lifelines)
    - [Sample Code](#Sample-Code-Lifelines)
    - [Limitations](#Limitations-Lifelines)
- [Summary](#Summary)
- [Extra Resources](#Extra-Resourses)


---


# Background

Survival analysis encompasses a set of statistical methods designed to analyze time-to-event data—where the outcome is the time until the occurrence of an event of interest, such as death, relapse, equipment failure, or recovery. These methods are uniquely capable of handling censored data, which occurs when the event has not been observed for all subjects during the observation period.

Three fundamental tools in survival analysis are:

- Kaplan-Meier Estimator: A non-parametric method that estimates the survival function from lifetime data. It calculates the probability of surviving past certain time points and accounts for right-censored observations. It is particularly useful for visualizing survival curves and comparing them descriptively.

- Log-rank Test: A hypothesis test used to formally compare two or more Kaplan-Meier survival curves. It evaluates whether there is a statistically significant difference in the survival distributions between groups. The test is based on the observed versus expected number of events over time and assumes proportional hazards.

- Cox Proportional Hazards Model (Cox-PH): A semi-parametric regression model that estimates the hazard (or risk) of an event occurring, adjusting for covariates. It does not assume a specific baseline hazard function but assumes that hazard ratios between groups remain constant over time (proportional hazards). This model is powerful for examining the impact of continuous and categorical covariates on survival.

These methods are widely used in clinical trials, epidemiology, engineering, and other fields where time-to-event data is critical. Their ability to accommodate censoring, handle covariates, and provide interpretable results make them essential tools in both exploratory and confirmatory data analysis to estimate the probability of survival past a certain point in time.

## Hypotheses

<ins>**Log-rank test:**</ins>

  - H₀: Survival functions are equal across groups.

  - H₁: At least one group differs.

<ins>**Cox-PH:**</ins>

  - H₀: Hazard ratio for each covariate = 1 (no effect).

  - H₁: At least one covariate has HR ≠ 1.
  
## Assumptions

<ins>**Kaplan-Meier (KM)**</ins>

  - Independent censoring: Censored subjects are assumed to have the same survival prospects as those who remain under observation.

  - Non-informative censoring: Censoring is unrelated to the event risk.

  - Exact event times: Event times are known and measured without error.

  - Random sampling: Subjects are representative of the population of interest.

<ins>**Log-Rank Test**</ins>

  - Proportional hazards: Hazard ratios between groups are constant over time.

  - Independent censoring: Same as in KM.
  
  - Independent survival times: Between individuals and across groups.

  - Correct group assignment: Every individual belongs to one (and only one) group.

<ins>**Cox Proportional Hazards (Cox-PH)**</ins>

  - Proportional hazards: The hazard ratio between individuals is constant over time.

  - Linearity of covariates in log-hazard: The effect of continuous variables is linear on the log hazard scale.

  - Independent survival times: Survival times are independent between individuals.

  - Independent censoring: Censoring is unrelated to survival probability.

  - No omitted confounders: All relevant predictors are included in the model (not formally testable but crucial for valid inference).

## Test Statistic

<ins>**Kaplan-Meier (KM):**</ins>

  The Kaplan-Meier survival estimate at time t is:

    S(t) = product over all tᵢ ≤ t of [1 - dᵢ / nᵢ]

    where:

      dᵢ = number of events at time tᵢ

      nᵢ = number at risk just before tᵢ

    Variance of the survival estimate (Greenwood's formula):

      Var[S(t)] ≈ S(t)² × sum over tᵢ ≤ t of [dᵢ / (nᵢ × (nᵢ - dᵢ))]

<ins>**Log-rank Test:**</ins>

  The test statistic is:

    Q = (O - E)² / Var(O)

    where:

      O = observed number of events

      E = expected number of events under the null hypothesis

    Var(O) = variance of observed events under the null

    Q follows a chi-squared distribution with degrees of freedom = (number of groups - 1)

<ins>**Cox Proportional Hazards (Cox-PH):**</ins>

  The partial likelihood for estimating regression coefficients:

    L(β) = product over all events i of [exp(βXᵢ) / sum over j in R(tᵢ) of exp(βXⱼ)]

    where R(tᵢ) is the risk set at time tᵢ

  The log partial likelihood is:

    log L(β) = sum over i of [βXᵢ - log sum over j in R(tᵢ) of exp(βXⱼ)]

  Common test statistics:

    Wald test: Z = β̂ / SE(β̂)

    Likelihood Ratio Test (LRT): Compare log-likelihoods of nested models

    Score test: Uses derivative of the log-likelihood at β = 0

    Hazard Ratio (HR): HR = exp(β̂)

# Package Implementations

Dependencies:
All methods rely only on standard, widely available components in their respective environments:

- SAS procedures require no additional modules beyond SAS/STAT.

- R functions use only base packages (e.g., survival).

- Python (lifelines) dependencies (numpy, pandas, scipy, matplotlib, patsy) are installed automatically via pip.

## SAS LIFETEST

  ### Function <span style="display:none">(Lifetest)</span>
  
  - Use PROC LIFETEST for Kaplan-Meier survival curves and log-rank test.
    
  ### Inputs <span style="display:none">(Lifetest)</span>
  <ins>**Required:**</ins>
  
  - DATA= (input dataset), wide format
  
  - TIME statement for survival time and censoring
  
  <ins>**Optional**:</ins>
  
  - STRATA= (for group comparison)
  
  - PLOTS= (e.g. survival)
  
  - TEST (equality test)

  ### Outputs <span style="display:none">(Lifetest)</span>
  
  - Kaplan-Meier curves
  
  - Log-rank test statistic & p-value
  
  - Median survival & CI

  ### Sample Code <span style="display:none">(Lifetest)</span>
  ```sas
  proc lifetest data=study plots=survival;
    time time*status(0);
    strata treatment;
  run;
  ```

  ### Limitations <span style="display:none">(Lifetest)</span>
  - No built-in test for proportional hazards; requires user-driven diagnostic procedures.

  - Default method for handling ties is Efron's, which differs from some R/Python defaults.

  - Cannot handle time-varying covariates without macro extensions.

  - Limited support for custom baseline hazard specification.

## SAS PHREG    

  ### Function <span style="display:none">(PHREG)</span>
  
  - Use PROC PHREG from the SAS STAT module for Cox Proportional Hazards modeling.

  ### Inputs <span style="display:none">(PHREG)</span>
    
  <ins>**Required:**</ins>
  
  - DATA= (input dataset), wide format
  
  - MODEL (time and censoring indicator, covariates)
  
  <ins>**Optional:**</ins>
  
  - CLASS= (categorical predictors)
  
  - STRATA= (stratify baseline hazards)
  
  - HAZARDRATIO= (requests hazard ratios)
  
  - ID= (subject identifier)
  
  ### Outputs <span style="display:none">(PHREG)</span>
  
  - Hazard ratios (HR)
  
  - Standard errors and confidence intervals
  
  - Wald test, LRT statistics
  
  - Diagnostic outputs (residuals, influence stats)
    
### Sample Code <span style="display:none">(PHREG)</span>
  ```sas
    proc phreg data=study;
      model time*status(0) = age sex treatment;
    run;
  ```

### Limitations <span style="display:none">(PHREG)</span>  
  
  - Cannot incorporate covariates—limited to stratified or grouped comparisons.
  
  - Output can be verbose and require manual post-processing for reporting.
  
  - Graphical outputs depend on enabling ODS graphics manually.
  
## R SURVIVAL

  ### Function <span style="display:none">(Survival)</span>
  - Use survfit(), survdiff(), and coxph() from the survival package.
  
  ### Inputs <span style="display:none">(Survival)</span>
  <ins>**Required:**</ins>
  
  - Surv(time, status) object
  
  - Formula input for groups/covariates
  
  <ins>**Optional:**</ins>
  
  - cluster, weights, ties, na.action
  
  ### Outputs <span style="display:none">(Survival)</span>
  
  - survfit: KM estimates
  
  - survdiff: log-rank p-value
  
  - coxph: HRs, CIs, test stats
    
### Sample Code <span style="display:none">(Survival)</span>
  ```r
  library(survival)
  fit_km <- survfit(Surv(time, status) ~ group, data = df)
  survdiff(Surv(time, status) ~ group, data = df)
  fit_cox <- coxph(Surv(time, status) ~ age + treatment, data = df)
  summary(fit_cox)
  ```
### Limitations <span style="display:none">(Survival)</span>

  - Diagnostics (e.g., proportional hazards via cox.zph) must be explicitly run and interpreted.
  
  - Default handling of ties is Efron's method; Breslow available, but not default.
  
  - Formulas can become complex when many interactions or strata are involved.
  
## Python LIFELINES

  ### Function <span style="display:none">(Lifelines)</span>
  - Use KaplanMeierFitter, logrank_test, and CoxPHFitter from the lifelines package.

  ### Inputs <span style="display:none">(Lifelines)</span>
  <ins>Required:**</ins>
  
  - duration_col, event_col
  
  - Pandas DataFrame with covariates
  
  <ins>**Optional:**</ins>
  
  - robust, step_size, weights, formula

  ### Outputs <span style="display:none">(Lifelines)</span>
  
  - KM survival estimates
  
  - Log-rank p-values
  
  - HRs, standard errors, confidence intervals

  ### Sample Code <span style="display:none">(Lifelines)</span>
  ```python
  from lifelines import KaplanMeierFitter, CoxPHFitter
  from lifelines.statistics import logrank_test
  
  kmf = KaplanMeierFitter()
  kmf.fit(df['time'], df['event'])
  kmf.plot()
  
  logrank_test(df_A['time'], df_B['time'], df_A['event'], df_B['event']).print_summary()
  
  cph = CoxPHFitter()
  cph.fit(df, duration_col='time', event_col='event')
  cph.print_summary()
  ```
  
### Limitations <span style="display:none">(Lifelines)</span>

  - Assumption checking (e.g., proportional hazards) must be performed via .check_assumptions(); not automatic.
  
  - Results may vary depending on tie-breaking method used (Efron, Breslow, or exact).
  
  - Does not support time-dependent covariates as robustly as R.
  
  - Some plotting features require knowledge of Matplotlib.
  
# Summary

- Kaplan-Meier, Log-rank, and Cox Proportional Hazards (Cox-PH) are foundational tools in survival analysis, each suited to different analytic needs. Kaplan-Meier provides non-parametric survival estimates, Log-rank tests detect differences between survival curves, and Cox-PH models estimate the impact of covariates on time-to-event outcomes under the proportional hazards assumption. These methods are supported in validated statistical packages (SAS, R, Python) and widely accepted in regulatory settings, with varying strengths and limitations in implementation, flexibility, and diagnostics.

# Extra Resourses

**Regulatory and Industry References**

  - FDA Statistical Guidance for Clinical Trials: Discusses statistical principles in clinical trials; includes KM and Cox-PH.
    - https://www.fda.gov/media/71145/download

  - ICH E9 Guideline: Standardizes survival analysis approaches internationally.
    - https://www.ema.europa.eu/en/documents/scientific-guideline/ich-e-9-statistical-principles-clinical-trials-step-5_en.pdf

  - SAS Documentation: Procedures PROC LIFETEST and PROC PHREG used in FDA-reviewed submissions. SAS Lifetest, SAS PHREG
    - https://documentation.sas.com/doc/en/pgmsascdc/v_063/statug/statug_lifetest_overview.htm

  - R Survival Package: Official reference from CRAN, widely used in academic and industry publications. CRAN Survival
    - https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf

  - Python Lifelines Documentation: Supports accepted statistical methods for open-science usage. Lifelines Docs
    - https://lifelines.readthedocs.io/en/latest/

  - EMA Biostatistics Guidance: Endorses Cox and stratified log-rank for regulatory use.
    - https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-adjustment-covariates-randomised-clinical-trials_en.pdf

**Notes and SOPs**

  - UCLA SAS Survival Seminar
    - https://stats.oarc.ucla.edu/sas/seminars/sas-survival/

  - UC Davis Cox Model Presentation
    - https://health.ucdavis.edu/media-resources/ctsc/documents/pdfs/cph-model-presentation.pdf

  - FDA ICH E9 Statistical Principles
    - https://www.fda.gov/regulatory-information/search-fda-guidance-documents/e9-statistical-principles-clinical-trials
