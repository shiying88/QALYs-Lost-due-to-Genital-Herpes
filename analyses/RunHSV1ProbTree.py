import pickle
import numpy as np
import pandas as pd
from classes import DecisionTreeClass as dt
from supports.RunProbTreeSupport import get_summary_stats
from classes.ParameterClass import ParametersTypeOne
from classes.ProbTreeClasses import buildHSV1Tree

# Parameter initiation
DISCOUNT = 0.03      # (float) discounting rate
NUM_PSA = 1000       # (int) number of probability sensitivity analysis iterations
virus_type = 'HSV-1'
sex_list = ['male', 'female']
age_list = [21, 27, 32, 42]
age_text_list = ['18-24', '25-29', '30-34', '35-49']
hsv1_params = ParametersTypeOne()
mf_id_lists = []                # male and female incidence list
mf_utility_lists = []           # male and female QALYs list
mf_rate_of_recur_lists = []     # male and female rate of recurrence list

# QALYs breakdowns
# psychological loss, recurrent outbreaks, urinary retention, primary outbreak, aseptic meningitis
colnames = ['psych', 'recur', 'ur', 'primary', 'am', 'sex', 'age', 'iter_num']
component_df = pd.DataFrame(columns=colnames)

# Main Analysis and Probability Sensitivity Analysis
for i in range(NUM_PSA):
    print('current number of iterations:', i)
    # random number generator
    rng = np.random.RandomState(seed=i)
    # resample parameter independent by age and sex
    hsv1_params.resample_hsv1_non_age_sex_params(rng=rng)
    f_id_distr = hsv1_params.f_id_distr  # incidence - female
    m_id_distr = hsv1_params.m_id_distr  # incidence - male
    mf_id_lists.append(m_id_distr + f_id_distr)
    utility_list = []
    rate_of_recur_each_sex_list = []
    for sex in sex_list:
        # resample parameters dependent on sex
        hsv1_params.resample_hsv1_sex_params(sex=sex, rng=rng)
        rate_of_recur_age_list = []
        for age in age_list:
            # adjust disutility based on background utility (by age)
            hsv1_params.update_hsv1_age_params(age_of_infection=age)
            # construct HSV-1 probabilistic tree
            # rate_of_recur: record the rate of recurrences per year until no future recurrence
            dictDecisions, dictChances, dictTerminals, rate_of_recur = \
                buildHSV1Tree(params=hsv1_params, sex=sex, discount=DISCOUNT, age_of_infection=age)
            # rate of recurrences for each sex&age subpopulation
            rate_of_recur_age_list.append(rate_of_recur)
            # calculate expected QALYs loss
            myDT = dt.DecisionNode('d1', dictDecisions, dictChances, dictTerminals)
            myDT.evaluate()
            cost_utility = myDT.get_cost_utility()
            utility_list.append(cost_utility['c0'][-1])     # only get utility, do not want cost
            # age- and sex- specific component QALYs lost
            breakdown_uti_dic = myDT.get_component_loss()['c0']
            df = pd.DataFrame([[breakdown_uti_dic[colnames[0]], breakdown_uti_dic[colnames[1]],
                                breakdown_uti_dic[colnames[2]], breakdown_uti_dic[colnames[3]],
                                breakdown_uti_dic[colnames[4]], sex, age, i]],
                              columns=colnames)
            component_df = component_df.append(df)
        rate_of_recur_each_sex_list.append(rate_of_recur_age_list[-1])    # rate of recur is independent of age
    mf_utility_lists.append(utility_list)
    mf_rate_of_recur_lists.append(rate_of_recur_each_sex_list)


age_sex_specific_qaly, non_age_sex_qaly_avg, t_qaly_f_avg, t_qaly_m_avg, t_qaly_hsv_avg, qaly_id_non_age_sex_list = \
    get_summary_stats(mf_utility_lists=mf_utility_lists,
                      mf_id_lists=mf_id_lists,
                      age_list=age_list,
                      sex_list=sex_list,
                      num_psa=NUM_PSA,
                      virus_type='HSV-1')

# output dictionary
hsv1_Dic = {'age_sex_specific_qaly': age_sex_specific_qaly,  # age- and sex-specific QALYs lost per case
            'non_age_sex_qaly_avg': non_age_sex_qaly_avg,   # avg. QALYs lost per case
            't_qaly_f_avg': t_qaly_f_avg,   # total QALYs lost for female
            't_qaly_m_avg': t_qaly_m_avg,   # total QALYs lost for male
            't_qaly_hsv_avg': t_qaly_hsv_avg,   # total QALYs lost not by sex
            'rate_of_recur': mf_rate_of_recur_lists,    # rage of recurrences by sex
            'qaly_id_non_age_sex_list': qaly_id_non_age_sex_list}  # incidence

# save QALYs lost outputs
output = open('tree_outputs/dics/DicHSV1.pkl', 'wb')
pickle.dump(hsv1_Dic, output)
output.close()
# save component QALYs lost
component_df.to_csv('tree_outputs/component_utl/hsv1.csv')
