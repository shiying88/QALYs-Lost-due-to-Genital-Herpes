from source.DemographicValues import get_adjusted_disu
import pandas as pd


def get_recur_initiate_rate_helper(recur_type, parameter, long_term_therapy):
    """ generate a list of recurrence rates for the first and second year after primary outbreak
    :param recur_type: infrequent or frequent. For HSV-1, we only have infrequent type.
    :param parameter: the Parameter object.
    :param long_term_therapy: whether the person is on CST (chronic suppressive therapy). Only applicable for HSV-2.
    """
    if recur_type == 'infrequent':
        first_year_recur_rate = parameter.infreq_first_year_recur_rate_sample
    elif recur_type == 'frequent':
        if long_term_therapy:
            first_year_recur_rate = parameter.freq_cst_first_year_recur_rate_sample
        else:
            first_year_recur_rate = parameter.freq_nocst_first_year_recur_rate_sample
    else:
        raise ValueError('wrong type of recurrences')
    second_year_recur_rate = first_year_recur_rate - parameter.avg_yearly_recur_reduction_sample
    return first_year_recur_rate, second_year_recur_rate


def get_long_term_qaly_loss(start_age, length_of_condition, disease_disu, year_idx, discount, disease_dura=1):
    """ estimate long-term QALY loss, adjusted by background disutility
    :param start_age: initial age of certain condition
    :param length_of_condition: length of years people in this condition
    :param disease_disu: disutility of disease, unadjusted
    :param disease_dura: duration of disease within one year, default = 1 (i.e., exist for the full year)
    :param year_idx: number of years after primary HSV infection. in the same year of primary infection, year_idx = 0
    :param discount: discounting rate
    """
    idx = 0  # the number of years after start_age
    qaly_list = []  # list of QALY loss from start age to end age
    for age in range(start_age, round(length_of_condition) + start_age):
        # adjust disutility based on age in the current year
        current_year_adj_disu = get_adjusted_disu(current_age=age, disu=disease_disu)
        # calculate QALYs lost in the current year (discounting)
        current_year_qaly = current_year_adj_disu * disease_dura / (1 + discount)**(year_idx + idx)
        # attach QALYs lost into the list
        qaly_list.append(current_year_qaly)
        idx += 1
    # return total QALY loss due to this long-term condition
    return sum(qaly_list)


def separate_cst_helper(np_array, cst_prob):
    """ list_recur_type_after_sympt_primary include [normal, infrequent, frequent]
    separate frequent into frequent_no_cst and frequent_cst based on cst_prob
    :returns (list) [normal, infrequent, frequent_no_cst, frequent_cst]
    """
    list_of_array = list(np_array)
    # get % having frequent recurrence
    frequent = list_of_array[-1]
    # remove last element in the list
    list_of_array.pop()
    # calculate % having cst and % not having cst, append it to the list
    frequent_cst = frequent * cst_prob
    frequent_no_cst = frequent - frequent_cst
    list_of_array.append(frequent_no_cst)
    list_of_array.append(frequent_cst)
    return list_of_array


def life_table_get_long_term_loss(sex, age_onset, seq_disu, discount, seq_dura=0,
                                  acute_sympt_disu=0, emort=0, emort_dura=None, tmort=0, if_utility=False):
    """ calculate total QALYs loss for people with no/mild/moderate/severe sequelae due to encephalitis
    :param sex: 'female', or 'male', or 'general'
    :param age_onset: (float) age of disease onset
    :param acute_sympt_disu: disutility of having acute symptoms
    :param discount: discounting rate
    :param seq_disu: disutility of sequelae
    :param seq_dura: duration of sequelae. (default 0=lifetime), otherwise, specify length of years (1, 2, etc.)
    :param emort_dura: duration of having excess mortality (start from age onset)
    :param emort: (float) excess mortality due to living with sequelae
    :param tmort: (float) total mortality (= nqx + ad_nqx)
    :param if_utility: (bool) whether we calculate QALYs loss (using utility=false) or total QALYs (using utility=True)
    """
    if sex == 'female':
        df = pd.read_csv('/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/data/life_table_for_qaly_female.csv')
    elif sex == 'male':
        df = pd.read_csv('/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/data/life_table_for_qaly_male.csv')
    elif sex == 'general':
        df = pd.read_csv('/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/data/life_table_for_qaly.csv')
    else:
        raise ValueError('wrong type of sex input')
    # round age onset
    age_onset = round(age_onset)
    # excess mortality during disease sequelae period
    age_list = df['x']
    ad_nqx_list = []
    nqx_list = df['nqx']
    for pair in zip(age_list, nqx_list):
        if emort_dura is not None:
            if age_onset <= pair[0] < (age_onset + emort_dura):
                if emort > 0:
                    ad_nqx_list.append(emort)
                elif tmort > 0:
                    ad_nqx_list.append(max(tmort - pair[1], 0))
                else:
                    raise ValueError('emort_dura should be NoneType')
            else:
                ad_nqx_list.append(0)
        else:
            if pair[0] >= age_onset:
                if emort > 0:
                    ad_nqx_list.append(emort)
                elif tmort > 0:
                    ad_nqx_list.append(max(tmort - pair[1], 0))
                else:
                    ad_nqx_list.append(0)
            else:
                ad_nqx_list.append(0)
    df['ad_nqx'] = ad_nqx_list
    # total mortality/survival rate within time period x and x+n
    df['aggre_nqx'] = df.apply(lambda row: row['nqx'] + row['ad_nqx'], axis=1)  # mortality = npx + excess mortality
    df['npx'] = df.apply(lambda row: 1-row['aggre_nqx'], axis=1)                # survival rate = 1 - mortality
    # calculate % alive starting at the age of onset
    aggre_nqx_list = df['aggre_nqx']
    nlx_list = []
    seq_disu_list = []
    prev_row = None
    for pair in zip(age_list, aggre_nqx_list):
        # age of sequelae onset
        if pair[0] == age_onset:
            nlx_list.append(1)
            seq_disu_list.append(seq_disu + acute_sympt_disu)
        # age after sequelae onset
        elif pair[0] > age_onset:
            nlx_list.append(nlx_list[-1] * (1 - prev_row[1]))
            seq_disu_list.append(seq_disu)
        # age before sequelae onset
        else:
            nlx_list.append(0)
            seq_disu_list.append(seq_disu)
        prev_row = pair
    df['nlx'] = nlx_list
    # disutility & utility
    df['seq_disu'] = seq_disu_list
    df['seq_utl'] = df.apply(lambda row: 1 - row['seq_disu'], axis=1)
    # avg. life years people live within each interval
    df['nLx_ud'] = df.apply(lambda row: (row['nax']*row['aggre_nqx'] + row['nbx']*row['npx'])*row['nlx'], axis=1)
    # cumulative life years lived at age x/x+n
    nLx_ud_list = df['nLx_ud']
    xn_T_ud_list = []
    x_T_ud_list = []
    ud_T_total = 0
    for value_ud in nLx_ud_list:
        x_T_ud_list.append(ud_T_total)
        ud_T_total += value_ud
        xn_T_ud_list.append(ud_T_total)
    df['(x+n)T_ud'] = xn_T_ud_list
    df['(x)T_ud'] = x_T_ud_list
    # discounted life years people live within each interval = discounted_(x+n)T_ud - discounted_(x)T_ud
    if discount == 0:
        # df['nLx_d'] = df.apply(lambda row: row['nLx_ud'])
        df['nLx_d'] = df['nLx_ud']
    else:
        df['nLx_d'] = df.apply(lambda row: ((1-(1-discount)**row['(x+n)T_ud'])/discount - (1-(1-discount)**row['(x)T_ud'])/discount), axis=1)
    # life years lived with sequelae within each interval
    if seq_dura > 0:
        nLx_d_list = df['nLx_d']
        nlx_d_with_sequelae_list = []
        count_n_years = 0
        for pair in zip(age_list, nLx_d_list):
            if pair[0] >= age_onset:
                count_n_years += 1
                if count_n_years <= seq_dura:
                    nlx_d_with_sequelae_list.append(pair[1])
                else:
                    diff = count_n_years - seq_dura
                    nlx_d_with_sequelae_list.append(max(0, 1-diff))
            else:
                nlx_d_with_sequelae_list.append(0)
        df['nlx_d_with_sequelae'] = nlx_d_with_sequelae_list
    else:
        df['nlx_d_with_sequelae'] = df['nLx_d']
    # calculate QALYs loss with sequelae
    df['nQALYx_loss'] = df.apply(lambda row: row['nlx_d_with_sequelae'] * row['bg_utl'] * row['seq_disu'], axis=1)
    df['nQALYs'] = df.apply(lambda row: row['nlx_d_with_sequelae'] * row['bg_utl'] * row['seq_utl'], axis=1)
    if if_utility:
        return sum(df['nQALYs'])
    else:
        return sum(df['nQALYx_loss'])


def life_table_get_long_term_loss_mixed_sex(age_onset, seq_disu, discount,
                                            seq_dura=0, emort=0, acute_sympt_disu=0, if_utility=False):
    female_loss = life_table_get_long_term_loss(sex='female', age_onset=age_onset, acute_sympt_disu=acute_sympt_disu,
                                                seq_disu=seq_disu, discount=discount, seq_dura=seq_dura,
                                                emort=emort, if_utility=if_utility)
    male_loss = life_table_get_long_term_loss(sex='male', age_onset=age_onset, acute_sympt_disu=acute_sympt_disu,
                                              seq_disu=seq_disu, discount=discount, seq_dura=seq_dura,
                                              emort=emort, if_utility=if_utility)
    # sex ratio at birth
    male_ratio = 0.511480215
    female_ratio = 0.488519785

    return male_ratio * male_loss + female_ratio * female_loss

