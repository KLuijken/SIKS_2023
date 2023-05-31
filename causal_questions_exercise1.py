import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyreadr
from sklearn.linear_model import LogisticRegression
from statsmodels.api import Logit

# Load data
vaccine_data = pyreadr.read_r('vaccine_data.rds')[None]

# Describe data
print(vaccine_data.describe())

# ---------------------------------------------------------------------------- #
# Example analysis 1
# ---------------------------------------------------------------------------- #
# Fit propensity score model
ps_1 = LogisticRegression()
ps_1.fit(vaccine_data[['sex', 'age', 'cvd', 'pulm', 'dm']], vaccine_data['vacc'])

# Create density plot (unweighted sample)
plt.figure()
vaccine_data.groupby('vacc')['propensity_score_1'].plot(kind='density', legend=True)
plt.title('Unweighted Sample')
plt.xlabel('Propensity Score')
plt.ylabel('Density')
plt.legend(title='Vaccination')
plt.ylim(0, 5.5)
plt.show()

# Estimate weights
vaccine_data['weight_1'] = np.where(vaccine_data['vacc'] == 1, 1 / ps_1.predict_proba(vaccine_data[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm']])[:, 1], 1 / (1 - ps_1.predict_proba(vaccine_data[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm']])[:, 1]))

# Create density plot (weighted sample)
plt.figure()
vaccine_data.groupby('vacc')['propensity_score_1'].plot(kind='density', legend=True, weights=vaccine_data['weight_1'])
plt.title('Weighted Sample')
plt.xlabel('Propensity Score')
plt.ylabel('Density')
plt.legend(title='Vaccination')
plt.ylim(0, 5.5)
plt.show()

# Estimate effect of vaccination on mortality
mod_1 = Logit(vaccine_data['mort'], vaccine_data['vacc'], weights=vaccine_data['weight_1']).fit()
print(mod_1.summary())

# Compute marginal risk difference
vacc_yes_data = vaccine_data.copy()
vacc_yes_data['vacc'] = 1
mort_vacc_yes_1 = mod_1.predict(vacc_yes_data)

vacc_no_data = vaccine_data.copy()
vacc_no_data['vacc'] = 0
mort_vacc_no_1 = mod_1.predict(vacc_no_data)

rd_1 = np.mean(mort_vacc_yes_1 - mort_vacc_no_1)
print('Marginal Risk Difference:', rd_1)

# Bootstrap for standard error estimation
b_rep = 1000
effect_estimates_1 = np.zeros((b_rep, 1))

for i in range(b_rep):
    bootstrap_sample = vaccine_data.sample(n=len(vaccine_data), replace=True)
    ps_bs_1 = LogisticRegression()
    ps_bs_1.fit(bootstrap_sample[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm']], bootstrap_sample['vacc'])
    bootstrap_sample['weight_1'] = np.where(bootstrap_sample['vacc'] == 1, 1 / ps_bs_1.predict_proba(bootstrap_sample[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm']])[:, 1], 1 / (1 - ps_bs_1.predict_proba(bootstrap_sample[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm']])[:, 1]))
    mod_bs_1 = Logit(bootstrap_sample['mort'], bootstrap_sample['vacc'], weights=bootstrap_sample['weight_1']).fit(disp=0)
    vacc_yes_bs = bootstrap_sample.copy()
    vacc_yes_bs['vacc'] = 1
    mort_vacc_yes_bs = mod_bs_1.predict(vacc_yes_bs)
    vacc_no_bs = bootstrap_sample.copy()
    vacc_no_bs['vacc'] = 0
    mort_vacc_no_bs = mod_bs_1.predict(vacc_no_bs)
    effect_estimates_1[i] = np.mean(mort_vacc_yes_bs - mort_vacc_no_bs)

se_1 = np.std(effect_estimates_1)
print('Standard Error:', se_1)

# ---------------------------------------------------------------------------- #
# Example analysis 2
# ---------------------------------------------------------------------------- #
# Fit propensity score model
ps_2 = LogisticRegression()
ps_2.fit(vaccine_data[['sex', 'age', 'cvd', 'pulm', 'dm', 'wt']], vaccine_data['vacc'])

# Estimate weights
vaccine_data['weight_2'] = np.where(vaccine_data['vacc'] == 1, 1 / ps_2.predict_proba(vaccine_data[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm', 'wt']])[:, 1], 1 / (1 - ps_2.predict_proba(vaccine_data[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm', 'wt']])[:, 1]))

# Create density plot (weighted sample)
plt.figure()
vaccine_data.groupby('vacc')['propensity_score_2'].plot(kind='density', legend=True, weights=vaccine_data['weight_2'])
plt.title('Weighted Sample')
plt.xlabel('Propensity Score')
plt.ylabel('Density')
plt.legend(title='Vaccination')
plt.ylim(0, 5.5)
plt.show()

# Estimate effect of vaccination on mortality
mod_2 = Logit(vaccine_data['mort'], vaccine_data['vacc'], weights=vaccine_data['weight_2']).fit()
print(mod_2.summary())

# Compute marginal risk difference
vacc_yes_data = vaccine_data.copy()
vacc_yes_data['vacc'] = 1
mort_vacc_yes_2 = mod_2.predict(vacc_yes_data)

vacc_no_data = vaccine_data.copy()
vacc_no_data['vacc'] = 0
mort_vacc_no_2 = mod_2.predict(vacc_no_data)

rd_2 = np.mean(mort_vacc_yes_2 - mort_vacc_no_2)
print('Marginal Risk Difference:', rd_2)

# Bootstrap for standard error estimation
b_rep = 1000
effect_estimates_2 = np.zeros((b_rep, 1))

for i in range(b_rep):
    bootstrap_sample = vaccine_data.sample(n=len(vaccine_data), replace=True)
    ps_bs_2 = LogisticRegression()
    ps_bs_2.fit(bootstrap_sample[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm', 'wt']], bootstrap_sample['vacc'])
    bootstrap_sample['weight_2'] = np.where(bootstrap_sample['vacc'] == 1, 1 / ps_bs_2.predict_proba(bootstrap_sample[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm', 'wt']])[:, 1], 1 / (1 - ps_bs_2.predict_proba(bootstrap_sample[['sex', 'age', 'age_sq', 'cvd', 'pulm', 'dm', 'wt']])[:, 1]))
    mod_bs_2 = Logit(bootstrap_sample['mort'], bootstrap_sample['vacc'], weights=bootstrap_sample['weight_2']).fit(disp=0)
    vacc_yes_bs = bootstrap_sample.copy()
    vacc_yes_bs['vacc'] = 1
    mort_vacc_yes_bs = mod_bs_2.predict(vacc_yes_bs)
    vacc_no_bs = bootstrap_sample.copy()
    vacc_no_bs['vacc'] = 0
    mort_vacc_no_bs = mod_bs_2.predict(vacc_no_bs)
    effect_estimates_2[i] = np.mean(mort_vacc_yes_bs - mort_vacc_no_bs)

se_2 = np.std(effect_estimates_2)
print('Standard Error:', se_2)


