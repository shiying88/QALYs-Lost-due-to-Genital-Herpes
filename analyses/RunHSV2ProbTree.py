import pickle
import numpy as np
import pandas as pd
from classes import DecisionTreeClass as dt
from supports.RunProbTreeSupport import get_summary_stats
from classes.ParameterClass import ParametersTypeTwo
from classes.ProbTreeClasses import buildHSV2Tree

# Parameter initiation
DISCOUNT = 0.06      # (float) discounting rate
NUM_PSA = 1000       # (int) number of probability sensitivity analysis iterations
virus_type = 'HSV-2'
sex_list = ['male', 'female']
age_list = [21, 27, 32, 42]
age_text_list = ['18-24', '25-29', '30-34', '35-49']

hsv2_params = ParametersTypeTwo()
mf_id_lists = []        # male amd female incidence list
mf_utility_lists = []   # male and female QALYs list
mf_infreq_lists = []    # male and female # of infrequent recurrence list
mf_freq_lists = []      # male and female # of frequent recurrence list
mf_freq_cst_lists = []  # male and female # of frequent (with CST) recurrence list

# QALYs breakdowns
# primary outbreak, aseptic meningitis, urinary retention, psychological,
# infrequent recurrence, frequent recurrence (no CST), frequent recurrence (with CST)
colnames = ['primary', 'am', 'ur', 'psych', 'infreq', 'freq_no_cst', 'freq_cst', 'rm', 'sex', 'age', 'iter_num']
component_df = pd.DataFrame(columns=colnames)

# Main Analysis and Probability Sensitivity Analysis
for i in range(NUM_PSA):
    print('current number of iterations:', i)
    # random number generator
    rng = np.random.RandomState(seed=i)
    # resample age and sex independent parameters
    hsv2_params.resample_hsv2_non_age_sex_params(rng=rng)
    f_id_distr = hsv2_params.f_id_distr     # incidence - female
    m_id_distr = hsv2_params.m_id_distr     # incidence - male
    mf_id_lists.append(m_id_distr + f_id_distr)
    utility_list = []
    infreq_r_sex_list = []
    freq_r_sex_list = []
    freq_cst_r_sex_list = []
    # resample sex dependent parameters
    for sex in sex_list:
        utility_list_each_sex = []
        infreq_r_age_list = []
        freq_r_age_list = []
        freq_cst_r_age_list = []
        hsv2_params.resample_hsv2_sex_params(sex=sex, rng=rng)
        for age in age_list:
            # adjust disutility based on background utility (by age)
            hsv2_params.update_hsv2_age_params(age_of_infection=age)
            # construct HSV-2 probabilistic tree
            dictDecisions, dictChances, dictTerminals, infreq_r, freq_r, freq_cst_r = \
                buildHSV2Tree(params=hsv2_params, sex=sex, discount=DISCOUNT, age_of_infection=age)
            infreq_r_age_list.append(infreq_r)
            freq_r_age_list.append(freq_r)
            freq_cst_r_age_list.append(freq_cst_r)
            # calculate expected QALYs loss
            myDT = dt.DecisionNode('d1', dictDecisions, dictChances, dictTerminals)
            myDT.evaluate()
            cost_utility = myDT.get_cost_utility()
            age_sex_specific_utility = cost_utility['c0'][-1]
            utility_list_each_sex.append(age_sex_specific_utility)
            # age- and sex- specific component QALYs lost
            breakdown_uti_dic = myDT.get_component_loss()['c0']
            df = pd.DataFrame([[breakdown_uti_dic[colnames[0]], breakdown_uti_dic[colnames[1]],
                                breakdown_uti_dic[colnames[2]], breakdown_uti_dic[colnames[3]],
                                breakdown_uti_dic[colnames[4]], breakdown_uti_dic[colnames[5]],
                                breakdown_uti_dic[colnames[6]], breakdown_uti_dic[colnames[7]],
                                sex, age, i]],
                              columns=colnames)
            component_df = component_df.append(df)
        utility_list.append(utility_list_each_sex)
        infreq_r_sex_list.append(infreq_r_age_list[-1])
        freq_r_sex_list.append(freq_r_age_list[-1])
        freq_cst_r_sex_list.append(freq_cst_r_age_list[-1])

    mf_utility_lists.append(utility_list[0] + utility_list[1])
    mf_infreq_lists.append(infreq_r_sex_list)
    mf_freq_lists.append(freq_r_sex_list)
    mf_freq_cst_lists.append(freq_cst_r_sex_list)

age_sex_specific_qaly, non_age_sex_qaly_avg, t_qaly_f_avg, t_qaly_m_avg, t_qaly_hsv_avg, qaly_id_non_age_sex_list = \
    get_summary_stats(mf_utility_lists=mf_utility_lists,
                      mf_id_lists=mf_id_lists,
                      age_list=age_list,
                      sex_list=sex_list,
                      num_psa=NUM_PSA,
                      virus_type='HSV-2')
# output dictionary
hsv2_Dic = {'age_sex_specific_qaly': age_sex_specific_qaly,         # QALYs lost per case by age and sex
            'non_age_sex_qaly_avg': non_age_sex_qaly_avg,           # QALYs lost per case - general
            't_qaly_f_avg': t_qaly_f_avg,                           # total QALYs lost - female
            't_qaly_m_avg': t_qaly_m_avg,                           # total QALYs lost - male
            't_qaly_hsv_avg': t_qaly_hsv_avg,                       # total QALYs lost - general
            'infreq_recur': mf_infreq_lists,                        # num of infrequent recurrences by sex
            'freq_no_cst_recur': mf_freq_lists,                     # num of frequent recurrences by sex
            'freq_cst_recur': mf_freq_cst_lists,                    # num of frequent recurrences (CST) by sex
            'qaly_id_non_age_sex_list': qaly_id_non_age_sex_list}   # incidence

# Save outputs
output = open('tree_outputs/dics/DicHSV2.pkl', 'wb')
pickle.dump(hsv2_Dic, output)
output.close()
# save component QALYs lost
component_df.to_csv('tree_outputs/component_utl/hsv2.csv')
