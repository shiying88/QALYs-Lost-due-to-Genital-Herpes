from classes.RecurrentPeriodClass import RecurrentPeriodTypeTwo, RecurrentPeriodTypeOne


def initiate_recur_object_helper_type_two(params, rm, age_of_infection, sex, discount):
    """ for HSV-2 probabilistic tree, return three objects for no, infreq, freq, freq + cst terminal nodes """
    no = RecurrentPeriodTypeTwo(parameters=params, recur_type='no', rm=rm,
                                age_of_infection=age_of_infection, sex=sex, discount=discount)
    infreq = RecurrentPeriodTypeTwo(parameters=params, recur_type='infrequent', rm=rm,
                                    age_of_infection=age_of_infection, sex=sex, discount=discount)
    freq = RecurrentPeriodTypeTwo(parameters=params, recur_type='frequent', rm=rm,
                                  age_of_infection=age_of_infection, sex=sex, discount=discount)
    freq_cst = RecurrentPeriodTypeTwo(parameters=params, recur_type='frequent', cst=True, rm=rm,
                                      age_of_infection=age_of_infection, sex=sex, discount=discount)
    return no, infreq, freq, freq_cst


def buildHSV1Tree(params, age_of_infection, sex, discount):
    """ construct a probabilistic tree for HSV-1 infection """
    # disutility values during primary infection
    p_untreat = params.primary_sympt_undia_notreat_qaly
    p_treat = params.primary_sympt_diag_treat_qaly
    am = params.am_qaly
    ur = params.ur_qaly

    # cost, disutility, component_disutility_dics, [future nodes], probabilities
    dictChances = {'c0': [0, 0, {}, ['c1', 'c2'], params.list_primary_type],
                   'c1': [0, 0, {}, ['c3', 'c6', 'c7'], params.nor_am_ur_after_primary_prob_sample, {}],
                   'c2': [0, 0, {}, ['cn5', 'cr5'], params.list_recur_type_after_asympt_primary, {}],
                   'c3': [0, 0, {}, ['c4', 'c5'], params.list_diagnose, {}],
                   'c4': [0, p_untreat, {'primary': p_untreat}, ['cn1', 'cr1'], params.list_recur_type_after_sympt_primary],
                   'c5': [0, p_treat, {'primary': p_treat}, ['cn2', 'cr2'], params.list_recur_type_after_sympt_primary],
                   'c6': [0, p_treat, {'primary': p_treat}, ['c8'], [1]],
                   'c7': [0, p_treat, {'primary': p_treat}, ['c9'], [1]],
                   'c8': [0, am, {'am': am}, ['cn3', 'cr3'], params.list_recur_type_after_sympt_primary],
                   'c9': [0, ur, {'ur': ur}, ['cn4', 'cr4'], params.list_recur_type_after_sympt_primary]
                   }

    # no/undetected recurrences (ence=False is newly added during revision)
    recur_no = RecurrentPeriodTypeOne(parameters=params, recur_type='no', ence=False,
                                      age_of_infection=age_of_infection, sex=sex, discount=discount)
    # infrequent recurrences (ence=False is newly added during revision)
    recur_infreq = RecurrentPeriodTypeOne(parameters=params, recur_type='infrequent', ence=False,
                                          age_of_infection=age_of_infection, sex=sex, discount=discount)

    # recurrent disutility component dictionary
    # for the first version (before revision), include QALYs lost associated with ence: recur_no.ence_qaly
    # during the first revision, exclude encephalitis, assuming it has a 0 loss.
    # the dictionary is used to generate breakdown estimates.
    total_loss_hsv1_components_no_recur = {'ence': recur_no.ence_qaly}
    total_loss_hsv1_components_infreq = {'psych': recur_infreq.psych_qaly,
                                         'recur': recur_infreq.recur_only_qaly,
                                         'ur': recur_infreq.ur_qaly,
                                         'ence': recur_no.ence_qaly}

    # cost, disutility, component_disutility_dics
    dictTerminals = \
        {'cn1': [0, recur_no.total_loss_hsv1, total_loss_hsv1_components_no_recur],  # ence
         'cn2': [0, recur_no.total_loss_hsv1, total_loss_hsv1_components_no_recur],  # ence
         'cn3': [0, recur_no.total_loss_hsv1, total_loss_hsv1_components_no_recur],  # ence
         'cn4': [0, recur_no.total_loss_hsv1, total_loss_hsv1_components_no_recur],  # ence
         'cn5': [0, recur_no.total_loss_hsv1, total_loss_hsv1_components_no_recur],  # ence
         # recur type: infrequent recurrences
         'cr1': [0, recur_infreq.total_loss_hsv1, total_loss_hsv1_components_infreq],   # psycho + recur + ence
         'cr2': [0, recur_infreq.total_loss_hsv1, total_loss_hsv1_components_infreq],   # psycho + recur + ence
         'cr3': [0, recur_infreq.total_loss_hsv1, total_loss_hsv1_components_infreq],   # psycho + recur + ence
         'cr4': [0, recur_infreq.total_loss_hsv1, total_loss_hsv1_components_infreq],   # psycho + recur + ence
         'cr5': [0, recur_infreq.total_loss_hsv1, total_loss_hsv1_components_infreq]    # psycho + recur + ence
         }

    # cost, utility, component_utility_dic, future_node
    dictDecisions = {'d1': [0, 0, {}, ['c0']]}

    # rate of recurrences for people with infrequent recurrent outbreaks
    rate_of_recur = recur_infreq.rate_of_recurrence

    return dictDecisions, dictChances, dictTerminals, rate_of_recur


def buildHSV2Tree(params, age_of_infection, sex, discount):
    """ construct a probabilistic tree for HSV-2 infection """
    # disutility values during primary infection
    p_untreat = params.primary_sympt_undia_notreat_qaly
    p_treat = params.primary_sympt_diag_treat_qaly
    ur = params.ur_qaly
    am = params.am_qaly

    # cost, disutility, [future nodes], probabilities
    dictChances = {'c0': [0, 0, {}, ['c1', 'c2'], params.list_primary_type],
                   'c1': [0, 0, {}, ['c3', 'c4', 'c5'], params.nor_am_ur_after_primary_prob_sample],
                   'c2': [0, 0, {}, ['t6', 'ri6'], params.list_recur_type_after_asympt_primary],
                   'c3': [0, 0, {}, ['c6', 'c7'], params.list_diagnose],
                   'c4': [0, p_treat, {'primary': p_treat}, ['c10'], [1]],
                   'c5': [0, p_treat, {'primary': p_treat}, ['c11'], [1]],
                   'c6': [0, p_untreat, {'primary': p_untreat}, ['t1', 'ri1', 'rf1', 'rfc1'], params.list_recur_type_after_sympt_primary],
                   'c7': [0, p_treat, {'primary': p_treat}, ['t2', 'ri2', 'rf2', 'rfc2'], params.list_recur_type_after_sympt_primary],
                   'c8': [0, 0, {}, ['t3', 'ri3', 'rf3', 'rfc3'], params.list_recur_type_after_sympt_primary],
                   'c9': [0, 0, {}, ['t4', 'ri4', 'rf4', 'rfc4'], params.list_recur_type_after_sympt_primary],
                   'c10': [0, ur, {'ur': ur}, ['c8', 'c9'], params.list_recur_meningitis],
                   'c11': [0, am, {'am': am}, ['t5', 'ri5', 'rf5', 'rfc5'], params.list_recur_type_after_sympt_primary]
                   }

    # define recurrent objects
    # after undiagnosed symptomatic primary HSV
    t1, ri1, rf1, rfc1 = initiate_recur_object_helper_type_two(params=params, rm=False, sex=sex, discount=discount,
                                                               age_of_infection=age_of_infection)
    # after diagnosed symptomatic primary HSV
    t2, ri2, rf2, rfc2 = initiate_recur_object_helper_type_two(params=params, rm=False, sex=sex, discount=discount,
                                                               age_of_infection=age_of_infection)
    # after symptomatic primary HSV (with aseptic meningitis) + no recurrent meningitis
    t3, ri3, rf3, rfc3 = initiate_recur_object_helper_type_two(params=params, rm=False, sex=sex, discount=discount,
                                                               age_of_infection=age_of_infection)
    # after symptomatic primary HSV (with aseptic meningitis) + recurrent meningitis
    t4, ri4, rf4, rfc4 = initiate_recur_object_helper_type_two(params=params, rm=True, sex=sex, discount=discount,
                                                               age_of_infection=age_of_infection)
    # after symptomatic primary HSV (with urinary retention)
    t5, ri5, rf5, rfc5 = initiate_recur_object_helper_type_two(params=params, rm=False, sex=sex, discount=discount,
                                                               age_of_infection=age_of_infection)
    # after asymptomatic primary HSV
    ri6 = RecurrentPeriodTypeTwo(parameters=params, recur_type='infrequent', age_of_infection=age_of_infection,
                                 sex=sex, discount=discount)

    dictTerminals = \
        {'t1': [0, t1.total_loss_hsv2, {}],  # no sympt recur -> no psycho loss
         't2': [0, t2.total_loss_hsv2, {}],  # no sympt recur -> no psycho loss
         't3': [0, t3.total_loss_hsv2, {}],  # no sympt recur -> no psycho loss
         't4': [0, t4.total_loss_hsv2, {'rm': t4.rm_qaly}],     # risk of recurrent meningitis
         't5': [0, t5.total_loss_hsv2, {}],  # no sympt recur -> no psycho loss
         't6': [0, 0, {}],                   # not diagnosed & no need for treatment
         # symptomatic primary, undiagnosed, no complication associated with primary HSV
         'ri1': [0, ri1.total_loss_hsv2, {'infreq': ri1.recur_only_qaly, 'ur': ri1.ur_qaly, 'psych': ri1.psych_qaly}],
         'rf1': [0, rf1.total_loss_hsv2, {'freq_no_cst': rf1.recur_only_qaly, 'ur': rf1.ur_qaly, 'psych': rf1.psych_qaly}],
         'rfc1': [0, rfc1.total_loss_hsv2, {'freq_cst': rfc1.recur_only_qaly, 'ur': rfc1.ur_qaly, 'psych': rfc1.psych_qaly}],
         # symptomatic primary, diagnosed, no complication associated with primary HSV
         'ri2': [0, ri2.total_loss_hsv2, {'infreq': ri2.recur_only_qaly, 'ur': ri2.ur_qaly, 'psych': ri2.psych_qaly}],
         'rf2': [0, rf2.total_loss_hsv2, {'freq_no_cst': rf2.recur_only_qaly, 'ur': rf2.ur_qaly, 'psych': rf2.psych_qaly}],
         'rfc2': [0, rfc2.total_loss_hsv2, {'freq_cst': rfc2.recur_only_qaly, 'ur': rfc2.ur_qaly, 'psych': rfc2.psych_qaly}],
         # symptomatic primary, diagnosed, meningitis_primary + no meningitis_recurrent
         'ri3': [0, ri3.total_loss_hsv2, {'infreq': ri3.recur_only_qaly, 'ur': ri3.ur_qaly, 'psych': ri3.psych_qaly}],
         'rf3': [0, rf3.total_loss_hsv2, {'freq_no_cst': rf3.recur_only_qaly, 'ur': rf3.ur_qaly, 'psych': rf3.psych_qaly}],
         'rfc3': [0, rfc3.total_loss_hsv2, {'freq_cst': rfc3.recur_only_qaly, 'ur': rfc3.ur_qaly, 'psych': rfc3.psych_qaly}],
         # symptomatic primary, diagnosed, meningitis_primary + meningitis_recurrent
         'ri4': [0, ri4.total_loss_hsv2, {'infreq': ri4.recur_only_qaly, 'ur': ri4.ur_qaly, 'psych': ri4.psych_qaly, 'rm': ri4.rm_qaly}],
         'rf4': [0, rf4.total_loss_hsv2, {'freq_no_cst': rf4.recur_only_qaly, 'ur': rf4.ur_qaly, 'psych': rf4.psych_qaly, 'rm': rf4.rm_qaly}],
         'rfc4': [0, rfc4.total_loss_hsv2, {'freq_cst': rfc4.recur_only_qaly, 'ur': rfc4.ur_qaly, 'psych': rfc4.psych_qaly, 'rm': rfc4.rm_qaly}],
         # symptomatic primary, diagnosed, urinary retention_primary
         'ri5': [0, ri5.total_loss_hsv2, {'infreq': ri5.recur_only_qaly, 'ur': ri5.ur_qaly, 'psych': ri5.psych_qaly}],
         'rf5': [0, rf5.total_loss_hsv2, {'freq_no_cst': rf5.recur_only_qaly, 'ur': rf5.ur_qaly, 'psych': rf5.psych_qaly}],
         'rfc5': [0, rfc5.total_loss_hsv2, {'freq_cst': rfc5.recur_only_qaly, 'ur': rfc5.ur_qaly, 'psych': rfc5.psych_qaly}],
         # asymptomatic primary, undiagnosed, no need to treat, no complication associated with primary HSV
         'ri6': [0, ri6.total_loss_hsv2, {'infreq': ri6.recur_only_qaly, 'ur': ri6.ur_qaly, 'psych': ri6.psych_qaly}]
         }

    dictDecisions = {'d1': [0, 0, {}, ['c0']]}

    infreq_recur = ri1.rate_of_recurrence
    freq_no_cst_recur = rf1.rate_of_recurrence
    freq_cst_recur = rfc1.rate_of_recurrence

    return dictDecisions, dictChances, dictTerminals, infreq_recur, freq_no_cst_recur, freq_cst_recur


def buildNeonatalTree(params):
    """ construct a structure for neonatal HSV
    :param params: (object) input parameters
    # :param sim_time: (int) simulation length, 0=life time (default, neonatal period), ≥1: number of years
    """

    # cost, disutility, [future nodes], probabilities
    dict_chances = {'c0': [0, 0, {}, ['c1', 'c2'], params.list_infection_types],
                    'c1': [0, 0, {}, ['t1', 't2'], params.list_intrauterine_demise_or_infected],
                    'c2': [0, 0, {}, ['c4', 'c5', 'c6'], params.list_intrapartum_short_term_outcomes],
                    # 'c3': [0, 0, ['t2', 't3', 't4'], params.list_intrauterine_long_term_outcomes],
                    'c4': [0, 0, {}, ['t5', 't6', 't7'], params.list_intrapartum_long_term_sem],
                    'c5': [0, 0, {}, ['t10', 't11', 't12', 't13', 't14'], params.list_intrapartum_long_term_cns],
                    'c6': [0, 0, {}, ['t15', 't16', 't17', 't18', 't19'], params.list_intrapartum_long_term_dis]
                    }

    # Maternal + neonatal
    dict_terminals = \
        {
         # intrauterine
         't1': [0, params.get_death_qaly_loss(), {}],
         't2': [0, params.get_intrauterine_qaly_loss(), {}],
         # intrapartum, SEM
         't5': [0, params.get_intrapartum_qaly_loss(acute='SEM', sequelae='nor'), {}],
         't6': [0, params.get_intrapartum_qaly_loss(acute='SEM', sequelae='mil'), {}],
         't7': [0, params.get_intrapartum_qaly_loss(acute='SEM', sequelae='mod'), {}],
         # intrapartum, CNS
         't10': [0, params.get_intrapartum_qaly_loss(acute='CNS', sequelae='nor'), {}],
         't11': [0, params.get_intrapartum_qaly_loss(acute='CNS', sequelae='mil'), {}],
         't12': [0, params.get_intrapartum_qaly_loss(acute='CNS', sequelae='mod'), {}],
         't13': [0, params.get_intrapartum_qaly_loss(acute='CNS', sequelae='sev'), {}],
         't14': [0, params.get_death_qaly_loss(), {}],
         # intrapartum, DIS
         't15': [0, params.get_intrapartum_qaly_loss(acute='DIS', sequelae='nor'), {}],
         't16': [0, params.get_intrapartum_qaly_loss(acute='DIS', sequelae='mil'), {}],
         't17': [0, params.get_intrapartum_qaly_loss(acute='DIS', sequelae='mod'), {}],
         't18': [0, params.get_intrapartum_qaly_loss(acute='DIS', sequelae='sev'), {}],
         't19': [0, params.get_death_qaly_loss(), {}]
         }

    dict_maternal_terminals = \
        {
         # intrauterine
         't1': [0, params.get_maternal_loss(neonate_outcome='dea'), {}],
         't2': [0, params.get_maternal_loss(neonate_outcome='neuro'), {}],
         # intrapartum, SEM
         't5': [0, params.get_maternal_loss(neonate_outcome='nor'), {}],
         't6': [0, params.get_maternal_loss(neonate_outcome='mil'), {}],
         't7': [0, params.get_maternal_loss(neonate_outcome='mod'), {}],
         # intrapartum, CNS
         't10': [0, params.get_maternal_loss(neonate_outcome='nor'), {}],
         't11': [0, params.get_maternal_loss(neonate_outcome='mil'), {}],
         't12': [0, params.get_maternal_loss(neonate_outcome='mod'), {}],
         't13': [0, params.get_maternal_loss(neonate_outcome='sev'), {}],
         't14': [0, params.get_maternal_loss(neonate_outcome='dea'), {}],
         # intrapartum, DIS
         't15': [0, params.get_maternal_loss(neonate_outcome='nor'), {}],
         't16': [0, params.get_maternal_loss(neonate_outcome='mil'), {}],
         't17': [0, params.get_maternal_loss(neonate_outcome='mod'), {}],
         't18': [0, params.get_maternal_loss(neonate_outcome='sev'), {}],
         't19': [0, params.get_maternal_loss(neonate_outcome='dea'), {}]
         }

    # print('maternal loss due to child death', params.get_maternal_loss(neonate_outcome='dea'))
    # print('maternal loss due to child with mid sequelae', params.get_maternal_loss(neonate_outcome='mil'))
    # print('maternal loss due to child with mod sequelae', params.get_maternal_loss(neonate_outcome='mod'))
    # print('maternal loss due to child with sev sequelae', params.get_maternal_loss(neonate_outcome='sev'))

    dict_decisions = {'d1': [0, 0, {}, ['c0']]}

    return dict_decisions, dict_chances, dict_terminals, dict_maternal_terminals
