import pandas as pd
from SimPy.Statistics import SummaryStat


def convert_list_to_dic(list_data, age_list, sex_list):
    """ convert [[m1, m2, m3, m4, f1, f2, f3, f4], [...], ..., [...]] into dictionary """
    colnames = []
    for sex in sex_list:
        for age in age_list:
            name = '{}_{}'.format(sex, age)
            colnames.append(name)
    df = pd.DataFrame(list_data, columns=colnames)
    uListDic = {'male': {}, 'female': {}}
    uStatDic = {'male': {}, 'female': {}}
    for sex in sex_list:
        for age in age_list:
            uListDic[sex][age] = df['{}_{}'.format(sex, age)].tolist()
            uStatDic[sex][age] = SummaryStat(name='QALYs loss per infection for {} at age of {}'.format(sex, age),
                                             data=df['{}_{}'.format(sex, age)].tolist())
    return uListDic, uStatDic


def get_pair_product(list1, list2):
    inner_product = []
    for a, b in zip(list1, list2):
        inner_product.append(a * b)
    return inner_product


def get_summary_stats(mf_utility_lists, mf_id_lists, age_list, sex_list, num_psa, virus_type):
    """ calculate summary stats, including:
    1) QALYs loss per infection within each age- and sex-subgroups
    2) QALYs loss per infection for all people with HSV-1 or HSV-2 infection
    3) total QALYs loss per gender (e.g., total QALYs loss due to female HSV-1 infection in 2018)
    4) total QALYs loss per virus type (e.g., total QALYs loss due to HSV-1 in 2018)
    :param mf_utility_lists: sex and age specific QALYs loss per infection
    :param mf_id_lists: sex and age specific incidence
    :param age_list: list of age groups
    :param sex_list: list of sex groups
    :param num_psa: number of probability sensitivity analysis
    :param virus_type: 'HSV-1' or 'HSV-2'
    :return: several outcomes
    """

    # 1) QALYs loss per infection within each age and sex subgroup
    qaly_list_dic, qaly_stat_dic = convert_list_to_dic(list_data=mf_utility_lists, age_list=age_list, sex_list=sex_list)
    # Incidence per infection within each age and sex subgroup
    id_list_dic, id_stat_dic = convert_list_to_dic(list_data=mf_id_lists, age_list=age_list, sex_list=sex_list)
    # for each (age- and sex-specific) id_stat_dic
    age_sex_specific_qaly = []  # [m1, m2, m3, m4, f1, f2, f3, f4]
    for sex in sex_list:
        for age in age_list:
            avg = qaly_stat_dic[sex][age].get_formatted_mean_and_interval(deci=5, interval_type='p')
            age_sex_specific_qaly.append(avg)
            # print('QALYs loss per {} HSV infection with an avg age of {}:'.format(sex, age), avg)

    # 2) QALYs loss per infection among the general population
    non_age_sex_qaly = []
    qaly_id_non_age_sex_list = []
    for i in range(num_psa):
        qaly_non_age_sex_per_psa = []
        id_non_age_sex_per_psa = []
        for sex in sex_list:
            for age in age_list:
                qaly_this_group = qaly_list_dic[sex][age][i]
                id_this_group = id_list_dic[sex][age][i]
                total_qaly_this_group = qaly_this_group * id_this_group
                id_non_age_sex_per_psa.append(id_this_group)
                qaly_non_age_sex_per_psa.append(total_qaly_this_group)
        non_age_sex_qaly.append(sum(qaly_non_age_sex_per_psa) / sum(id_non_age_sex_per_psa))
        qaly_id_non_age_sex_list.append([qaly_non_age_sex_per_psa, id_non_age_sex_per_psa])
    non_age_sex_qaly_stat = SummaryStat(name='QALYs loss per {} infection among the general population'.format(virus_type),
                                        data=non_age_sex_qaly)
    non_age_sex_qaly_avg = non_age_sex_qaly_stat.get_formatted_mean_and_interval(deci=5, interval_type='p')
    print()
    print('QALYs loss per {} infection among the general population:'.format(virus_type), non_age_sex_qaly_avg)

    # 3) Total QALYs loss each gender (summing up age groups)
    # 4) Total QALYs loss per virus type (summing up sex groups)
    t_qaly_dic = {'female': [], 'male': []}
    for age in age_list:
        t_qaly_dic['female'].append(get_pair_product(list1=qaly_list_dic['female'][age], list2=id_list_dic['female'][age]))
        t_qaly_dic['male'].append(get_pair_product(list1=qaly_list_dic['male'][age], list2=id_list_dic['male'][age]))
    t_qaly_loss_for_m = []
    t_qaly_loss_for_f = []
    t_qaly_loss_for_hsv = []
    for i in range(num_psa):
        t_qaly_loss_m = 0
        t_qaly_loss_f = 0
        t_qaly_loss = 0
        for f_qaly, m_qaly in zip(t_qaly_dic['female'], t_qaly_dic['male']):
            t_qaly_loss_m += m_qaly[i]
            t_qaly_loss_f += f_qaly[i]
            t_qaly_loss += m_qaly[i] + f_qaly[i]
        t_qaly_loss_for_m.append(t_qaly_loss_m)
        t_qaly_loss_for_f.append(t_qaly_loss_f)
        t_qaly_loss_for_hsv.append(t_qaly_loss)
    print()
    t_qaly_m_stat = SummaryStat(name='total QALYs loss among males due to {} infection in 2018'.format(virus_type),
                                data=t_qaly_loss_for_m)
    t_qaly_m_avg = t_qaly_m_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
    print('Total QALYs loss among males due to {} infection in 2018:'.format(virus_type), t_qaly_m_avg)
    t_qaly_f_stat = SummaryStat(name='total QALYs loss among females due to {} infection in 2018'.format(virus_type),
                                data=t_qaly_loss_for_f)
    t_qaly_f_avg = t_qaly_f_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
    print('Total QALYs loss among females due to {} infection in 2018:'.format(virus_type), t_qaly_f_avg)

    print()
    t_qaly_hsv_stat = SummaryStat(name='total QALYs loss for {} infection in 2018'.format(virus_type),
                                  data=t_qaly_loss_for_hsv)
    t_qaly_hsv_avg = t_qaly_hsv_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
    print('Total QALYs loss due to {} infection in 2018:'.format(virus_type), t_qaly_hsv_avg)
    return age_sex_specific_qaly, non_age_sex_qaly_avg, t_qaly_f_avg, t_qaly_m_avg, t_qaly_hsv_avg, qaly_id_non_age_sex_list
