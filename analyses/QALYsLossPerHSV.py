import pickle
from SimPy.Statistics import SummaryStat

# read data
read_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/trees/dics/'
pkl_file = open('{}DicHSV1_dr3_psycho1_noGBD_20221010_dN100.pkl'.format(read_path), 'rb')
hsv1_dic = pickle.load(pkl_file)
pkl_file = open('{}DicHSV2_dr3_psycho1_noGBD_20221010_dN100.pkl'.format(read_path), 'rb')
hsv2_dic = pickle.load(pkl_file)
pkl_file.close()

hsv1_outcome_list = hsv1_dic['qaly_id_non_age_sex_list']
hsv2_outcome_list = hsv2_dic['qaly_id_non_age_sex_list']

# within each PSA iteration
avg_qaly_list = []
total_qaly_list = []
avg_hsv1_id_list = []
avg_hsv2_id_list = []
for i in range(len(hsv1_outcome_list)):
    hsv1_qaly_per_psa = hsv1_outcome_list[i][0]
    hsv1_id_per_psa = hsv1_outcome_list[i][1]
    hsv2_qaly_per_psa = hsv2_outcome_list[i][0]
    hsv2_id_per_psa = hsv2_outcome_list[i][1]
    # total QALYs loss
    qaly_per_psa = sum(hsv1_qaly_per_psa) + sum(hsv2_qaly_per_psa)
    # total HSV incidence
    id_per_psa = sum(hsv1_id_per_psa) + sum(hsv2_id_per_psa)
    # append
    avg_hsv1_id_list.append(sum(hsv1_id_per_psa))
    avg_hsv2_id_list.append(sum(hsv2_id_per_psa))
    avg_qaly_list.append(qaly_per_psa / id_per_psa)
    total_qaly_list.append(qaly_per_psa)

######################
# SUMMARY STATISTICS #
######################
DECIMAL = 2
# mean & 95% UI: avg. HSV-1 incidence
avg_hsv1_id_stat = SummaryStat(name='average genital HSV-1 incidence', data=avg_hsv1_id_list)
avg_hsv1_id = avg_hsv1_id_stat.get_formatted_mean_and_interval(deci=DECIMAL, interval_type='p')
print('avg. HSV-1 incidence', avg_hsv1_id)

# mean & 95% UI: avg. HSV-2 incidence
avg_hsv2_id_stat = SummaryStat(name='average HSV-2 incidence', data=avg_hsv2_id_list)
avg_hsv2_id = avg_hsv2_id_stat.get_formatted_mean_and_interval(deci=DECIMAL, interval_type='p')
print('avg. HSV-2 incidence', avg_hsv2_id)

# mean & 95% UI: avg. QALYs loss per infection
avg_qaly_stat = SummaryStat(name='(overall) QALYs lost per HSV infection in the US', data=avg_qaly_list)
avg_avg_qaly = avg_qaly_stat.get_formatted_mean_and_interval(deci=DECIMAL, interval_type='p')
print('(Overall) avg. QALYs lost per HSV infection:', avg_avg_qaly)

# mean & 95% UI: total QALYs loss per infection
total_qaly_stat = SummaryStat(name='total QALYs lost in the US', data=total_qaly_list)
avg_total_qaly = total_qaly_stat.get_formatted_mean_and_interval(deci=DECIMAL, interval_type='p')
print(' avg. total QALYs lost:', avg_total_qaly)
