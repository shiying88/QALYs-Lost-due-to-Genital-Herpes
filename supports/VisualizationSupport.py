import re
import pickle
import numpy as np
import math
import matplotlib.pyplot as plt

sex_list = ['male', 'female']
age_list = [21, 27, 32, 42]


def get_mean_min_max(string):
    split_list = re.split('[(,)]', string)
    mean = float(split_list[0])
    min_value = float(split_list[1])
    max_value = float(split_list[2])
    return mean, min_value, max_value


def get_mean_min_max_list(list_of_strings):
    mean_list = []
    min_list = []
    max_list = []
    for item in list_of_strings:
        mean_v, min_v, max_v = get_mean_min_max(item)
        mean_list.append(mean_v)
        min_list.append(min_v)
        max_list.append(max_v)
    return mean_list, min_list, max_list


def get_age_sex_qaly_loss(dic):
    """ get age and sex specific QALYs loss"""
    # age- and sex- specific QALYs loss
    age_sex_outcome_list = dic['age_sex_specific_qaly']
    age_sex_Dic = {'male': {}, 'female': {}}
    i = 0
    for sex in sex_list:
        for age in age_list:
            age_sex_Dic[sex][age] = age_sex_outcome_list[i]
            i += 1
    age_sex_male_list = []
    age_sex_female_list = []

    for age in age_list:
        age_sex_male_list.append(age_sex_Dic['male'][age])
        age_sex_female_list.append(age_sex_Dic['female'][age])

    age_sex_mean_male, age_sex_min_male, age_sex_max_male = get_mean_min_max_list(age_sex_male_list)
    age_sex_mean_female, age_sex_min_female, age_sex_max_female = get_mean_min_max_list(age_sex_female_list)

    ci_male = []
    ci_female = []

    for i in range(len(age_sex_min_male)):
        ci_male.append(age_sex_max_male[i] - age_sex_min_male[i])
        ci_female.append(age_sex_max_female[i] - age_sex_min_female[i])

    return age_sex_mean_male, ci_male, age_sex_mean_female, ci_female


def get_total_qaly_loss(dic, scalar=1000):
    mean_list = []
    # min_list = []
    # max_list = []
    ci_list = []        # confidence interval list of male, female, total
    # 3) Total QALYs loss per gender
    mean_m_t_qaly, min_m_t_qaly, max_m_t_qaly = get_mean_min_max(dic['t_qaly_m_avg'])       # men
    mean_list.append(mean_m_t_qaly / scalar)
    # min_list.append(min_m_t_qaly / scalar)
    # max_list.append(max_m_t_qaly / scalar)
    ci_list.append((max_m_t_qaly - min_m_t_qaly) / scalar)
    mean_f_t_qaly, min_f_t_qaly, max_f_t_qaly = get_mean_min_max(dic['t_qaly_f_avg'])       # women
    mean_list.append(mean_f_t_qaly / scalar)
    # min_list.append(min_f_t_qaly / scalar)
    # max_list.append(max_f_t_qaly / scalar)
    ci_list.append((max_f_t_qaly - min_f_t_qaly) / scalar)
    # 4) Total QALYs loss per virus type
    mean_t_qaly, min_t_qaly, max_t_qaly = get_mean_min_max(dic['t_qaly_hsv_avg'])
    mean_list.append(mean_t_qaly / scalar)
    # min_list.append(min_t_qaly / scalar)
    # max_list.append(max_t_qaly / scalar)
    ci_list.append((max_t_qaly - min_t_qaly) / scalar)
    return mean_list, ci_list


def get_neonate_stats(dic_5, dic_10, dic_15, dic_lifetime):
    """
    :param dic_5: dictionary outcome with sim_length = 5 years
    :param dic_10: dictionary outcome with sim_length = 10 years
    :param dic_15: dictionary outcome with sim_length = 15 years
    :param dic_lifetime: dictionary outcome with sim_length = lifetime
    :return: (4 lists of strings) per neonatal case loss, total loss, maternal_per case loss, maternal total loss
    """
    outcome_list = [dic_5, dic_10, dic_15, dic_lifetime]
    loss_per_case = [dic['loss_per_infection'] for dic in outcome_list]
    loss_total = [dic['loss_total'] for dic in outcome_list]
    m_loss_per_case = [dic['m_loss_per_infection'] for dic in outcome_list]
    m_loss_total = [dic['m_loss_total'] for dic in outcome_list]

    mean_n_per_case_list = []
    ci_n_per_case_list = []
    mean_m_per_case_list = []
    ci_m_per_case_list = []
    mean_total_list = []
    ci_total_list = []
    mean_m_total_list = []
    ci_m_total_list = []
    for i in range(len(loss_per_case)):
        # loss per case
        mean_v, min_v, max_v = get_mean_min_max(loss_per_case[i])
        mean_n_per_case_list.append(mean_v)
        ci_n_per_case_list.append(max_v - min_v)
        # maternal loss per case
        mean_v, min_v, max_v = get_mean_min_max(m_loss_per_case[i])
        mean_m_per_case_list.append(mean_v)
        ci_m_per_case_list.append(max_v - min_v)
        # total loss
        mean_v, min_v, max_v = get_mean_min_max(loss_total[i])
        mean_total_list.append(mean_v)
        ci_total_list.append(max_v - min_v)
        # maternal total loss
        mean_v, min_v, max_v = get_mean_min_max(m_loss_total[i])
        mean_m_total_list.append(mean_v)
        ci_m_total_list.append(max_v - min_v)

    return mean_n_per_case_list, ci_n_per_case_list, mean_m_per_case_list, ci_m_per_case_list,\
           mean_total_list, ci_total_list, mean_m_total_list, ci_m_total_list


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def get_recur_traj_by_sex(list_of_recur_outcomes):
    """ the input list could be directly accessed from the saved dictionary for HSV1/HSV2 """
    m_traj_list = []
    m_x_list = []
    m_max_x = 0
    m_max_rate = 0
    f_traj_list = []
    f_x_list = []
    f_max_x = 0
    f_max_rate = 0
    for iteration_outcome in list_of_recur_outcomes:
        # male
        male_traj = iteration_outcome[0]
        male_traj.append(0)
        m_traj_list.append(male_traj)
        m_x_list.append(np.arange(len(male_traj)))
        m_max_x = max(m_max_x, len(male_traj))
        m_max_rate = max(m_max_rate, max(male_traj))
        # female
        female_traj = iteration_outcome[1]
        female_traj.append(0)
        f_traj_list.append(female_traj)
        f_x_list.append(np.arange(len(female_traj)))
        f_max_x = max(f_max_x, len(female_traj))
        f_max_rate = max(f_max_rate, max(female_traj))
    # find xlim and ylim values
    x_max = max(m_max_x, f_max_x)
    y_max = math.ceil(max(m_max_rate, f_max_rate))
    return m_traj_list, m_x_list, f_traj_list, f_x_list, x_max, y_max


def get_hsv_sa_stats(hsv_dic_list):
    """ return: six list, [QALYs loss / infection under each scenario] and [total QALYs loss under each scenario] """
    # order: base case, psycho0, psycho0.02, (genital5, genital25, genital75, genital95), GBD, discount0, discount0.05
    qaly_loss_per_case_list = [dic['non_age_sex_qaly_avg'] for dic in hsv_dic_list]
    total_loss_list = [dic['t_qaly_hsv_avg'] for dic in hsv_dic_list]

    # QALYs loss per case
    qaly_loss_per_case_mean_list = []
    qaly_loss_per_case_ci_list = []
    for element in qaly_loss_per_case_list:
        mean_var, min_var, max_var = get_mean_min_max(element)
        qaly_loss_per_case_mean_list.append(mean_var)
        qaly_loss_per_case_ci_list.append(max_var - min_var)

    # total QALYs loss
    total_loss_mean_list = []
    total_loss_ci_list = []
    for element in total_loss_list:
        mean_var, min_var, max_var = get_mean_min_max(element)
        total_loss_mean_list.append(mean_var)
        total_loss_ci_list.append(max_var - min_var)

    return qaly_loss_per_case_mean_list, qaly_loss_per_case_ci_list, total_loss_mean_list, total_loss_ci_list


# set figure titles BOLD
def setBold(txt): return r"$\bf{" + str(txt) + "}$"
