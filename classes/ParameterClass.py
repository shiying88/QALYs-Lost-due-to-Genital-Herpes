from supports.UtilityOrderingSupport import *
from supports.DemographicValuesSupport import get_adjusted_disu
from supports.RandomValueGenerators import BetaValueGenerator, LogNormalValueGenerator, DirichletValueGenerator
from supports.ParameterAndRecurrentPeriodSupport import separate_cst_helper, \
    life_table_get_long_term_loss_mixed_sex, life_table_get_long_term_loss


class Incidence:
    def __init__(self):
        # Incidence / disease burden for HSV-1
        self.hsv1_genital_id = \
            LogNormalValueGenerator(log_mean=11.556, log_std=0.766,
                                    note='genital hsv-1 incidence assuming 20% of total infection to be genital')
        self.hsv1_female_prob = BetaValueGenerator(mean=0.49, stdv=0.05, note='proportion of female genital cases')
        self.hsv1_female_age_specific_incidence = \
            DirichletValueGenerator(N=1000, prob_list=[54.0, 21.9, 7.4, 16.7],
                                    note='age distribution for female subpopulation')
        self.hsv1_male_age_specific_incidence = \
            DirichletValueGenerator(N=1000, prob_list=[30.0, 23.2, 11.2, 35.6],
                                    note='age distribution for male subpopulation')

        # Incidence / disease burden for HSV-2
        self.hsv2_male_incidence = LogNormalValueGenerator(log_mean=12.60, log_std=0.37, note='male incident hsv2')
        self.hsv2_female_incidence = LogNormalValueGenerator(log_mean=12.46, log_std=0.35, note='female incident hsv2')
        self.hsv2_female_age_specific_incidence = \
            DirichletValueGenerator(N=1000, prob_list=[54.0, 21.9, 7.4, 16.7],
                                    note='age distribution for female subpopulation')
        self.hsv2_male_age_specific_incidence = \
            DirichletValueGenerator(N=1000, prob_list=[30.0, 23.2, 11.2, 35.6],
                                    note='age distribution for male subpopulation')

        # outcome initiation
        self.female_age_specific_incidence_sample = 0
        self.male_age_specific_incidence_sample = 0
        self.genital_id_sample = 0      # only for HSV-1 incidence, used for prob(encephalitis) calculation

    def resample_hsv1_incidence(self, rng):
        # total number of genital incident cases
        self.genital_id_sample = self.hsv1_genital_id.sample(rng=rng)
        # prob of cases to be female
        f_prob_sample = self.hsv1_female_prob.sample(rng=rng)
        # total incident cases for male and female, respectively
        total_id_f =self.genital_id_sample * f_prob_sample
        total_id_m = self.genital_id_sample * (1 - f_prob_sample)
        # case distribution among age groups [18-24, 25-29, 30-34, 35-49] for male and female
        prob_f_age_spcfc_id_sample = self.hsv1_female_age_specific_incidence.sample(rng=rng)
        prob_m_age_spcfc_id_sample = self.hsv1_male_age_specific_incidence.sample(rng=rng)
        # number of cases among age groups
        self.female_age_specific_incidence_sample = [prob * total_id_f for prob in prob_f_age_spcfc_id_sample]
        self.male_age_specific_incidence_sample = [prob * total_id_m for prob in prob_m_age_spcfc_id_sample]

    def resample_hsv2_incidence(self, rng):
        # incidence by gender
        total_id_f = self.hsv2_female_incidence.sample(rng=rng)
        total_id_m = self.hsv2_male_incidence.sample(rng=rng)
        # case distribution among age groups [18-24, 25-29, 30-34, 35-49] for male and female
        prob_f_age_spcfc_id_sample = self.hsv2_female_age_specific_incidence.sample(rng=rng)
        prob_m_age_spcfc_id_sample = self.hsv2_male_age_specific_incidence.sample(rng=rng)
        # number of cases among age groups
        self.female_age_specific_incidence_sample = [prob * total_id_f for prob in prob_f_age_spcfc_id_sample]
        self.male_age_specific_incidence_sample = [prob * total_id_m for prob in prob_m_age_spcfc_id_sample]


class Parameters:
    def __init__(self):
        """ 1) age and sex insensitive parameters & initiations; 2) age and sex sensitive initiations """
        # INCIDENCE
        self.incidence = Incidence()
        # PROBABILITY
        # diagnostic test sensitivity
        self.diagnose_prob = BetaValueGenerator(mean=0.89, stdv=0.045, note='diagnosis rate')  # 0.89 [0.9, 0.98]

        # DURATION
        # primary infection - HSV symptoms
        self.primary_treat_duration = \
            LogNormalValueGenerator(mean=12 / 365, stdv=1 / 365,
                                    note='duration of primary sympt outbreak (with treatment)')
        self.primary_notreat_duration = \
            LogNormalValueGenerator(mean=18 / 365, stdv=2.65 / 365,
                                    note='duration of primary sympt outbreak (no treatment)')

        # complications
        self.urinary_retention_duration = LogNormalValueGenerator(mean=5 / 365, stdv=2.5 / 365,
                                                                  note='duration of sacral radiculitis UR')
        self.aseptic_meningitis_duration = LogNormalValueGenerator(mean=8 / 365, stdv=2.5 / 365,
                                                                   note='duration of aseptic meningitis')

        # DISUTILITY
        # long-term psychosocial disutility
        self.diagnosis_disu = BetaValueGenerator(mean=0.01, stdv=0.0025, note='disu due to hsv diagnosis')
        self.diagnosis_disu_reduction = \
            BetaValueGenerator(mean=0.005, stdv=0.00265, note='annual reduction in diagnosis disu as age increase')
        # primary infection period
        # From literatures
        self.sympt_outbreak_disu = BetaValueGenerator(mean=0.15, stdv=0.0375, note='disu of sympt outbreak')
        # # From GBD
        # self.primary_sympt_notreat_disu = \
        #     BetaValueGenerator(mean=0.051, stdv=0.0105, note='disutility due to symptomatic primary infection (no treat)')
        # self.primary_sympt_treat_disu = \
        #     BetaValueGenerator(mean=0.051, stdv=0.0105, note='disutility due to symptomatic primary infection (treated)')
        # self.recurrence_sympt_treated_disu = \
        #     BetaValueGenerator(mean=0.006, stdv=0.0025, note='disutility per sympt recurrence (treated)')
        # complications
        self.aseptic_meningitis_disu = BetaValueGenerator(mean=0.05, stdv=0.0125,
                                                          note='disutility due to aseptic meningitis')
        self.urinary_retention_disu = BetaValueGenerator(mean=0.03, stdv=0.01,
                                                         note='disutility due to urinary retention (recurrence)')

        ##################
        # initialization #
        ##################
        # Non Age and Sex Parameters #
        # incidence
        self.f_id_distr = None
        self.m_id_distr = None
        # incidence / rate of recurrence initiation
        self.avg_yearly_recur_reduction_sample = 0

        # primary and diagnosis
        self.list_diagnose = []
        self.diagnosis_disu_sample = 0
        self.diagnosis_disu_reduction_sample = 0
        self.list_primary_type = []
        self.primary_sympt_treat_disu_sample = 0
        self.primary_sympt_notreat_disu_sample = 0
        self.primary_treat_duration_sample = 0
        self.primary_notreat_duration_sample = 0

        # complication related to primary infection
        self.urinary_retention_disu_sample = 0
        self.urinary_retention_duration_sample = 0
        self.aseptic_meningitis_disu_sample = 0
        self.aseptic_meningitis_duration_sample = 0

        # recurrent period
        self.list_recur_type_after_asympt_primary = []
        self.list_recur_type_after_sympt_primary = []
        self.recur_sympt_treated_disu_sample = 0
        self.recur_treat_duration_sample = 0

        ######################################
        # initiate QALYs / probability lists #
        ######################################
        """ age dependent parameters """
        # QALYs
        self.primary_sympt_undia_notreat_qaly = 0
        self.primary_sympt_diag_treat_qaly = 0
        self.am_qaly = 0
        self.ur_qaly = 0

    def resample_non_age_sex_params(self, rng):
        # primary infection diagnosis [undiagnosed vs. diagnosed]
        diagnose_prob_sample = self.diagnose_prob.sample(rng=rng)
        self.list_diagnose = [(1 - diagnose_prob_sample), diagnose_prob_sample]

        # psychosocial discomfort due to diagnosis
        self.diagnosis_disu_sample = self.diagnosis_disu.sample(rng=rng)
        # self.diagnosis_disu_sample = 0
        self.diagnosis_disu_reduction_sample = self.diagnosis_disu_reduction.sample(rng=rng)
        # self.diagnosis_disu_reduction_sample = 0

        # primary symptomatic infection
        self.primary_sympt_notreat_disu_sample = self.sympt_outbreak_disu.sample(rng=rng)
        self.primary_sympt_treat_disu_sample = self.sympt_outbreak_disu.sample(rng=rng)
        self.recur_sympt_treated_disu_sample = self.sympt_outbreak_disu.sample(rng=rng)
        # GBD
        # self.primary_sympt_notreat_disu_sample = self.primary_sympt_notreat_disu.sample(rng=rng)
        # self.primary_sympt_treat_disu_sample = self.primary_sympt_treat_disu.sample(rng=rng)
        # self.recur_sympt_treated_disu_sample = self.recurrence_sympt_treated_disu.sample(rng=rng)
        self.primary_treat_duration_sample = self.primary_treat_duration.sample(rng=rng)
        self.primary_notreat_duration_sample = self.primary_notreat_duration.sample(rng=rng)
        # complication
        self.urinary_retention_disu_sample = self.urinary_retention_disu.sample(rng=rng)
        self.urinary_retention_duration_sample = self.urinary_retention_duration.sample(rng=rng)
        self.aseptic_meningitis_disu_sample = self.aseptic_meningitis_disu.sample(rng=rng)
        self.aseptic_meningitis_duration_sample = self.aseptic_meningitis_duration.sample(rng=rng)

    def update_age_params(self, age_of_infection):
        """ resample parameter values sensitive to age and sex based on pre-defined distributions """
        # qaly loss of patients with symptomatic primary infection (+ diagnosis + treatment)
        self.primary_sympt_diag_treat_qaly = \
            get_adjusted_disu(current_age=age_of_infection,
                              disu=self.primary_sympt_treat_disu_sample) * self.primary_treat_duration_sample
        # qaly loss of patients with symptomatic primary infection (undiagnosed + no treatment)
        self.primary_sympt_undia_notreat_qaly = \
            get_adjusted_disu(current_age=age_of_infection,
                              disu=self.primary_sympt_notreat_disu_sample) * self.primary_notreat_duration_sample

        # qaly loss of complications associated with primary infection
        self.ur_qaly = \
            get_adjusted_disu(current_age=age_of_infection,
                              disu=self.urinary_retention_disu_sample) * self.urinary_retention_duration_sample
        self.am_qaly = \
            get_adjusted_disu(current_age=age_of_infection,
                              disu=self.aseptic_meningitis_disu_sample) * self.aseptic_meningitis_duration_sample


class ParametersTypeOne(Parameters):
    def __init__(self):
        super().__init__()
        ####################
        # INCIDENCE / RATE #
        ####################
        # recurrence rate and avg. reduction rate in number of recurrences
        self.m_infreq_first_year_recur_rate = LogNormalValueGenerator(mean=1.14, stdv=0.05,
                                                                      note='first-year rate of recurrences')
        self.m_infreq_second_year_recur_rate = LogNormalValueGenerator(mean=1, stdv=0.05,
                                                                       note='second year rate of recur')
        self.f_infreq_first_year_recur_rate = LogNormalValueGenerator(mean=2.62, stdv=0.05,
                                                                      note='first-year rate of recurrences')
        self.f_infreq_second_year_recur_rate = LogNormalValueGenerator(mean=2.29, stdv=0.05,
                                                                       note='second year rate of recur')
        self.avg_yearly_recur_reduction = LogNormalValueGenerator(mean=0.5, stdv=0.025, note='avg rate of reduction/y')

        ###############
        # PROBABILITY #
        ###############
        # primary infection period
        self.primary_sympt_prob = BetaValueGenerator(mean=0.25, stdv=.075, note='primary symptomatic infection')
        self.sympt_primary_outcome = DirichletValueGenerator(N=100,  # [normal, aseptic meningitis, urinary retention]
                                                             # N=1000,
                                                             prob_list=[0.8575, 0.05, 0.0925],
                                                             note='outcome after primary sympt infection')
        # recurrent period
        self.recur_after_asympt_primary_prob = \
            BetaValueGenerator(mean=0.3, stdv=0.05, note='prob of infrequent recurrence after primary asympt infection')
        self.recur_after_sympt_primary_prob = \
            BetaValueGenerator(mean=0.57, stdv=0.05, note='prob of infrequent recurrence after primary sympt infection')
        self.urinary_retention_recur_prob = \
            BetaValueGenerator(mean=0.045, stdv=0.0025, note='probability of urinary retention after recurrence')
        self.num_ence_case = LogNormalValueGenerator(mean=3, stdv=0.5, note='number of ence cases among 1m population')
        self.hsv1_ence = BetaValueGenerator(mean=0.9, stdv=0.002, note='proportion of ence attributed to type1 virus')
        self.adult_ence = BetaValueGenerator(mean=0.915, stdv=0.002, note='proportion of ence happen among adults')
        self.ence_outcomes = \
            DirichletValueGenerator(N=182,  # [normal, mild, moderate, severe, death]
                                    prob_list=[0.187, 0.275, 0.231, 0.192, 0.115],
                                    note='long term outcome of encephalitis')

        ##############
        # DISUTILITY #
        ##############
        # disutilities for primary, recurrent, and common complications are defined in parent class
        # encephalitis at later age
        # self.ence_disu = BetaValueGenerator(mean=0.13, stdv=0.0255, note='disutility due to acute encephalitis')
        # self.ence_mil_disu = BetaValueGenerator(mean=0.03, stdv=0.0075, note='mild long-term sequelae of encephalitis')
        # self.ence_mod_disu = BetaValueGenerator(mean=0.20, stdv=0.04,
        #                                         note='moderate long-term sequelae of encephalitis')
        # self.ence_sev_disu = BetaValueGenerator(mean=0.54, stdv=0.0825,
        #                                         note='severe long-term sequelae of encephalitis')
        # self.ence_dead_disu = 1

        ############
        # DURATION #
        ############
        # primary infection period - treated and untreated - (in parent class)
        # recurrent treated duration
        self.recur_treat_duration = \
            LogNormalValueGenerator(mean=7/365, stdv=2.15/365, note='duration of recurrent sympt outbreak (treated)')
        # complications
        # urinary retention & aseptic meningitis (in parent class)
        # # encephalitis
        # self.ence_duration = LogNormalValueGenerator(mean=25 / 365, stdv=2.5 / 365, note='duration of encephalitis')
        # self.ence_sev_life_expectancy = \
        #     LogNormalValueGenerator(mean=3.8, stdv=1.25,
        #                             note='life expectancy for people with severe neurological sequelae after ence')

        ##################
        # INITIALIZATION #
        ##################
        # incidence/rate
        self.infreq_first_year_recur_rate_sample = 0
        self.infreq_second_year_recur_rate_sample = 0
        # primary infection
        self.nor_am_ur_after_primary_prob_sample = None
        # urinary retention, recurrent period
        self.urinary_retention_recur_prob_sample = None
        # # encephalitis
        # self.ence_elderly_prob_sample = None
        # self.ence_sequelae_prob_sample = None
        # self.ence_disu_sample = 0
        # self.ence_mil_disu_sample = 0
        # self.ence_mod_disu_sample = 0
        # self.ence_sev_disu_sample = 0
        # self.ence_duration_sample = 0
        # self.ence_sev_mort_dur_sample = 0
        # self.ence_sev_life_expectancy_sample = 0

    def resample_hsv1_non_age_sex_params(self, rng):
        self.resample_non_age_sex_params(rng=rng)
        ##############################
        # Incidence / disease burden #
        ##############################
        self.incidence.resample_hsv1_incidence(rng=rng)
        self.f_id_distr = self.incidence.female_age_specific_incidence_sample
        self.m_id_distr = self.incidence.male_age_specific_incidence_sample

        #############################################
        # Recurrence rate and annual reduction rate #
        #############################################
        self.avg_yearly_recur_reduction_sample = self.avg_yearly_recur_reduction.sample(rng=rng)

        #################
        # Probabilities #
        #################
        # Primary infection period
        # diagnosis prob (parent class)
        primary_sympt_prob_sample = self.primary_sympt_prob.sample(rng=rng)
        self.nor_am_ur_after_primary_prob_sample = self.sympt_primary_outcome.sample(rng=rng)
        infrequent_recur_after_asympt_primary_sample = self.recur_after_asympt_primary_prob.sample(rng=rng)
        infrequent_recur_after_sympt_primary_sample = self.recur_after_sympt_primary_prob.sample(rng=rng)

        # Recurrent period
        # urinary retention during recurrent period
        self.urinary_retention_recur_prob_sample = self.urinary_retention_recur_prob.sample(rng=rng)
        # # encephalitis at age of 60
        # num_ence_case_sample = self.num_ence_case.sample(rng=rng) / 1000000
        # hsv1_id = self.incidence.genital_id_sample
        # hsv1_ence = self.hsv1_ence.sample(rng=rng)
        # adult_ence = self.adult_ence.sample(rng=rng)
        # pop_size = 327115193
        # self.ence_elderly_prob_sample = hsv1_ence * adult_ence * num_ence_case_sample * pop_size / hsv1_id
        # # self.ence_elderly_prob_sample = self.ence_elderly_prob.sample(rng=rng)
        # self.ence_sequelae_prob_sample = self.ence_outcomes.sample(rng=rng)

        # SUMMARY LIST
        # primary infection diagnosis [undiagnosed vs. diagnosed] (parent class)
        # primary infection [symptomatic vs. asymptomatic manifestation]
        self.list_primary_type = [primary_sympt_prob_sample, (1 - primary_sympt_prob_sample)]
        # after asymptomatic primary infection [no/undetected recurrence vs. infrequent recurrence]
        self.list_recur_type_after_asympt_primary = [(1 - infrequent_recur_after_asympt_primary_sample),
                                                     infrequent_recur_after_asympt_primary_sample]
        # after symptomatic primary infection [no/undetected recurrence vs. infrequent recurrence]
        self.list_recur_type_after_sympt_primary = [(1 - infrequent_recur_after_sympt_primary_sample),
                                                    infrequent_recur_after_sympt_primary_sample]
        # # encephalitis total mortality for people with severe sequelae
        # self.ence_sev_mort_sample = self.ence_sev_mort.sample(rng=rng)

        ##############
        # Disutility #
        ##############
        # # primary infection and associated complication, recurrent period and associated complication (parent class)
        # self.ence_disu_sample = self.ence_disu.sample(rng=rng)
        # self.ence_mil_disu_sample = self.ence_mil_disu.sample(rng=rng)
        # self.ence_mod_disu_sample = self.ence_mod_disu.sample(rng=rng)
        # self.ence_sev_disu_sample = self.ence_sev_disu.sample(rng=rng)

        ############
        # Duration #
        ############
        # psychosocial discomfort due to diagnosis (parent class)
        # primary infection period (parent class)
        # recurrent period
        self.recur_treat_duration_sample = self.recur_treat_duration.sample(rng=rng)
        # complication (parent class - urinary retention & aseptic meningitis)
        # self.ence_duration_sample = self.ence_duration.sample(rng=rng)
        # self.ence_sev_life_expectancy_sample = self.ence_sev_life_expectancy.sample(rng=rng)

    def resample_hsv1_sex_params(self, sex, rng):
        # self.resample_sex_params(sex=sex, rng=rng)
        if sex == 'male':
            self.infreq_first_year_recur_rate_sample = self.m_infreq_first_year_recur_rate.sample(rng=rng)
            self.infreq_second_year_recur_rate_sample = self.m_infreq_second_year_recur_rate.sample(rng=rng)
        elif sex == 'female':
            self.infreq_first_year_recur_rate_sample = self.f_infreq_first_year_recur_rate.sample(rng=rng)
            self.infreq_second_year_recur_rate_sample = self.f_infreq_second_year_recur_rate.sample(rng=rng)
        else:
            raise ValueError('wrong sex type')

    def update_hsv1_age_params(self, age_of_infection):
        self.update_age_params(age_of_infection=age_of_infection)


class ParametersTypeTwo(Parameters):
    def __init__(self):
        super().__init__()
        ############################
        # SEX AGE INSENSITIVE PART #
        ############################
        ####################
        # INCIDENCE / RATE #
        ####################
        self.avg_yearly_recur_reduction = LogNormalValueGenerator(mean=0.5, stdv=0.025, note='avg rate of reduction /y')

        # complication (recurrent meningitis)
        self.total_num_recur_meningitis = LogNormalValueGenerator(mean=4.6, stdv=1.5,
                                                                  note='avg. total number of recurrent meningitis')

        ###############
        # PROBABILITY #
        ###############
        # primary infection period (the second row is updated during paper revision)
        self.primary_sympt_prob = BetaValueGenerator(mean=0.17, stdv=0.04, note='primary symptomatic infection')
        # recurrent period
        self.recur_after_asympt_primary_prob = \
            BetaValueGenerator(mean=0.62, stdv=0.05,
                               note='prob of infrequent recurrence after primary asympt infection')
        self.recur_after_sympt_primary_prob = \
            DirichletValueGenerator(N=203, prob_list=[0.11, 0.51, 0.38],
                                    note='prob of no/infrequent/frequent recurrence after primary sympt infection')
        self.prob_cst = BetaValueGenerator(mean=0.133, stdv=0.004,
                                           note='proportion people with frequent recur take CST')
        # recurrent meningitis for people having aseptic meningitis associated with primary HSV-2 infection
        self.recur_meningitis_prob = BetaValueGenerator(mean=0.3, stdv=0.06,
                                                        note='prob of people having recurrent meningitis')

        ############
        # DURATION #
        ############
        # primary infection period - treated & untreated - (in parent class)
        # recurrent treated duration
        self.recur_treat_duration = \
            LogNormalValueGenerator(mean=8.5 / 365, stdv=1.75 / 365,
                                    note='duration of recurrent sympt outbreak (treated)')
        # complications
        # urinary retention & aseptic meningitis (in parent class)
        self.recur_meningitis_duration = \
            LogNormalValueGenerator(mean=6 / 365, stdv=2.75 / 365, note='duration of recurrent meningitis')
        self.total_year_with_meningitis = \
            LogNormalValueGenerator(mean=8.4, stdv=4.75, note='total years experiencing meningitis symptoms')

        ##############
        # DISUTILITY #
        ##############
        # disutilities for primary, recurrent, and common complications are defined in parent class
        # recurrent meningitis
        self.recur_meningitis_disu = BetaValueGenerator(mean=0.13, stdv=0.0255,
                                                        note='disutility due to 1 recurrent meningitis')

        ##############
        # INITIATION #
        ##############
        self.total_num_recur_meningitis_sample = 0
        self.list_recur_meningitis = None
        self.list_cst = None
        self.recur_meningitis_disu_sample = None
        self.recur_meningitis_duration_sample = None
        self.total_year_with_meningitis_sample = None

        ######################
        # SEX SENSITIVE PART #
        ######################
        ####################
        # INCIDENCE / RATE #
        ####################
        # recurrent rate and avg. reduction rate in number of recurrences
        self.m_infreq_first_year_recur_rate = \
            LogNormalValueGenerator(mean=3.09, stdv=0.05, note='male, first-year rate of infrequent recur')
        self.m_freq_nocst_first_year_recur_rate = \
            LogNormalValueGenerator(mean=9.76, stdv=0.05, note='male, first-year rate of frequent recur without CST')
        self.m_freq_cst_first_year_recur_rate = \
            LogNormalValueGenerator(mean=2.44, stdv=0.05, note='male, first-year rate of frequent recur with CST')
        self.f_infreq_first_year_recur_rate = \
            LogNormalValueGenerator(mean=2.29, stdv=0.05, note='female, first-year rate of infrequent recur')
        self.f_freq_nocst_first_year_recur_rate = \
            LogNormalValueGenerator(mean=9.87, stdv=0.05, note='female, first-year rate of frequent recur without CST')
        self.f_freq_cst_first_year_recur_rate = \
            LogNormalValueGenerator(mean=2.47, stdv=0.05, note='female, first-year rate of frequent recur with CST')

        ###############
        # PROBABILITY #
        ###############
        self.m_sympt_primary_outcome = DirichletValueGenerator(N=100,
                                                               # N=63,
                                                               # [normal, aseptic meningitis, urinary retention]
                                                               prob_list=[0.924, 0.016, 0.06],
                                                               note='male outcome after primary sympt infection')
        self.f_sympt_primary_outcome = DirichletValueGenerator(N=100,
                                                               # N=126,
                                                               # [normal, aseptic meningitis, urinary retention]
                                                               prob_list=[0.811, 0.064, 0.125],
                                                               note='female outcome after primary sympt infection')
        # complication for people without aseptic meningitis associated with primary HSV-2 infection
        self.m_urinary_retention_recur_prob = \
            BetaValueGenerator(mean=0.03, stdv=0.015, note='male prob of urinary retention after recurrence')
        self.f_urinary_retention_recur_prob = \
            BetaValueGenerator(mean=0.06, stdv=0.03, note='female prob of urinary retention after recurrence')

        ##############
        # INITIATION #
        ##############
        self.nor_am_ur_after_primary_prob_sample = None
        self.urinary_retention_recur_prob_sample = None
        self.infreq_first_year_recur_rate_sample = None
        self.freq_nocst_first_year_recur_rate_sample = None
        self.freq_cst_first_year_recur_rate_sample = None

    def resample_hsv2_non_age_sex_params(self, rng):
        self.resample_non_age_sex_params(rng=rng)
        ####################
        # INCIDENCE / RATE #
        ####################
        # Incidence / disease burden
        self.incidence.resample_hsv2_incidence(rng=rng)
        self.f_id_distr = self.incidence.female_age_specific_incidence_sample
        self.m_id_distr = self.incidence.male_age_specific_incidence_sample

        # primary infection
        ###############
        # PROBABILITY #
        ###############
        # diagnosis prob (parent class)
        primary_sympt_prob_sample = self.primary_sympt_prob.sample(rng=rng)
        infrequent_recur_after_asympt_primary_sample = self.recur_after_asympt_primary_prob.sample(rng=rng)

        recur_meningitis_prob_sample = self.recur_meningitis_prob.sample(rng=rng)
        use_cst_prob_sample = self.prob_cst.sample(rng=rng)

        # recurrent period
        self.recur_treat_duration_sample = self.recur_treat_duration.sample(rng=rng)
        self.list_recur_type_after_sympt_primary = \
            separate_cst_helper(np_array=self.recur_after_sympt_primary_prob.sample(rng=rng),
                                cst_prob=use_cst_prob_sample)
        self.total_num_recur_meningitis_sample = self.total_num_recur_meningitis.sample(rng=rng)

        # SUMMARY LIST
        # primary infection diagnosis [undiagnosed vs. diagnosed] (parent class)
        # primary infection [symptomatic vs. asymptomatic]
        self.list_primary_type = [primary_sympt_prob_sample, (1 - primary_sympt_prob_sample)]
        # after asymptomatic primary infection [no/undetected recurrence vs. infrequent recurrence]
        self.list_recur_type_after_asympt_primary = [(1 - infrequent_recur_after_asympt_primary_sample),
                                                     infrequent_recur_after_asympt_primary_sample]
        # recurrent meningitis [no recurrent vs. with recurrent attacks]
        self.list_recur_meningitis = [(1 - recur_meningitis_prob_sample), recur_meningitis_prob_sample]
        # people with frequent recurrence type choosing CST [not taking CST vs. taking CST]
        self.list_cst = [(1 - use_cst_prob_sample), use_cst_prob_sample]

        ##############
        # DISUTILITY #
        ##############
        # primary (parent class)
        # recurrent (parent class)
        self.recur_meningitis_disu_sample = self.recur_meningitis_disu.sample(rng=rng)

        ############
        # DURATION #
        ############
        # psychosocial discomfort due to diagnosis (parent class)
        # primary (parent class)
        # recurrent (parent class)
        # complication (parent class: urinary retention & aseptic meningitis)
        self.recur_meningitis_duration_sample = self.recur_meningitis_duration.sample(rng=rng)
        self.total_year_with_meningitis_sample = self.total_year_with_meningitis.sample(rng=rng)

        #############################################
        # recurrence rate and annual reduction rate #
        #############################################
        self.avg_yearly_recur_reduction_sample = self.avg_yearly_recur_reduction.sample(rng=rng)

    def resample_hsv2_sex_params(self, sex, rng):
        # self.resample_sex_params(sex=sex, rng=rng)
        if sex == 'male':
            # outcome after primary infection
            self.nor_am_ur_after_primary_prob_sample = self.m_sympt_primary_outcome.sample(rng=rng)

            # urinary retention during recurrent period
            self.urinary_retention_recur_prob_sample = self.m_urinary_retention_recur_prob.sample(rng=rng)

            #############################################
            # recurrence rate and annual reduction rate #
            #############################################
            self.infreq_first_year_recur_rate_sample = self.m_infreq_first_year_recur_rate.sample(rng=rng)
            self.freq_nocst_first_year_recur_rate_sample = self.m_freq_nocst_first_year_recur_rate.sample(rng=rng)
            self.freq_cst_first_year_recur_rate_sample = self.m_freq_cst_first_year_recur_rate.sample(rng=rng)
        elif sex == 'female':
            # outcome after primary infection
            self.nor_am_ur_after_primary_prob_sample = self.f_sympt_primary_outcome.sample(rng=rng)

            # urinary retention during recurrent period
            self.urinary_retention_recur_prob_sample = self.f_urinary_retention_recur_prob.sample(rng=rng)

            #############################################
            # recurrence rate and annual reduction rate #
            #############################################
            self.infreq_first_year_recur_rate_sample = self.f_infreq_first_year_recur_rate.sample(rng=rng)
            self.freq_nocst_first_year_recur_rate_sample = self.f_freq_nocst_first_year_recur_rate.sample(rng=rng)
            self.freq_cst_first_year_recur_rate_sample = self.f_freq_cst_first_year_recur_rate.sample(rng=rng)
        else:
            raise ValueError('wrong sex type')

    def update_hsv2_age_params(self, age_of_infection):
        self.update_age_params(age_of_infection=age_of_infection)


class ParametersNeonatal:
    def __init__(self, discount=0.03, sim_time=0, num_psa=1000):
        """
        :param discount: discounting rate
        :param sim_time: simulation length. 0: lifetime, ≥1: the actual number of years to simulation
        """
        self.num_psa = num_psa

        # discounting
        self.discount = discount

        # simulation length
        self.sim_time = sim_time

        ####################
        # INCIDENCE / RATE #
        ####################
        self.incidence = LogNormalValueGenerator(mean=11.5, stdv=1.5, note='annual neonatal incidence/100,000 births')

        ###############
        # PROBABILITY #
        ###############
        # type of transmission
        self.route_of_transmission = \
            DirichletValueGenerator(N=500,
                                    prob_list=[0.05, 0.85, 0.1],
                                    note='prob. of transmission for intrauterine vs. intrapartum vs. postpartum')
        # clinical manifestations and long-term outcomes of intrauterine infection
        self.abortion = BetaValueGenerator(mean=0.078, stdv=0.015,
                                           note='prob of abortion due to intrauterine infection')
        self.intrauterine_sem = BetaValueGenerator(mean=0.32, stdv=0.058,
                                                   note='prob of SEM for neonates infected with intrauterine infection')
        self.intrauterine_long_term_outcomes = \
            DirichletValueGenerator(N=51,  # [Normal, Neurological impair, Death]
                                    prob_list=[0.16, 0.25, 0.59], note='long-term outcome of intrauterine infection')
        # clinical manifestations and long-term outcomes of intrapartum infection
        self.intrapartum_sem_cns_dis = DirichletValueGenerator(N=200,  # [SEM, CNS, Disseminated]
                                                               prob_list=[.45, .3, .25],
                                                               note='clinical manifestations of intrapartum infection')
        self.intrapartum_long_term_after_sem = \
            DirichletValueGenerator(N=200,
                                    prob_list=[0.98, 0.01, 0.01],  # [normal, mild, moderate]
                                    note='long-term outcome of intrapartum SEM infection')
        self.intrapartum_long_term_after_cns = \
            DirichletValueGenerator(N=200,
                                    prob_list=[0.3, 0.13, 0.17, 0.36, 0.04],  # [no, mi, mod, severe, death]
                                    note='long-term outcome of intrapartum CNS infection')
        self.intrapartum_long_term_after_dis = \
            DirichletValueGenerator(N=200,
                                    prob_list=[0.59, 0.01, 0.03, 0.08, 0.29],  # [n, mi, mo, s, d]
                                    note='long-term outcome of intrapartum disseminated infection')

        ############
        # DURATION #
        ############
        self.neonate_sev_life_expectancy = \
            LogNormalValueGenerator(mean=20, stdv=10, note='life expectancy of severely impaired neonates')
        self.maternal_dea_duration = \
            LogNormalValueGenerator(mean=0.6, stdv=0.475, note='initial grief period for losing a child')
        # after revision: assuming mother having moderate impaired child have higher utility after 1st year
        self.maternal_mod_initial_period_duration = \
            LogNormalValueGenerator(mean=1, stdv=0.5, note='initial grief period for moderate impaired child')
        self.m_life_expectancy = LogNormalValueGenerator(mean=53.1, stdv=1, note='avg life expectancy of mothers')
        self.m_avg_age_pregnancy = LogNormalValueGenerator(mean=29.03, stdv=1, note='avg age of pregnancy')

        ##############
        # DISUTILITY #
        ##############
        # long-term outcome
        # maternal side
        self.maternal_dea_disu = \
            BetaValueGenerator(mean=0.08, stdv=0.035, note='disutility of losing a child due to neonatal HSV')
        self.maternal_mil_disu = \
            BetaValueGenerator(mean=0.06, stdv=0.025, note='disu due to child has mild neurological sequelae')
        self.maternal_mod_disu = \
            BetaValueGenerator(mean=0.13, stdv=0.025, note='disu due to child has moderate neurological sequelae')
        self.maternal_sev_disu = \
            BetaValueGenerator(mean=0.24, stdv=0.07, note='disu due to child has severe neurological sequelae')
        # neonatal side
        self.intrauterine_neuro_disu = \
            BetaValueGenerator(mean=0.38, stdv=0.1075, note='disutility due to neurological sequelae')
        self.intrapartum_mild_disu = \
            BetaValueGenerator(mean=0.18, stdv=0.05, note='disutility due to intrapartum mild neurological sequelae')
        self.intrapartum_moderate_disu = \
            BetaValueGenerator(mean=0.48, stdv=0.05,
                               note='disutility due to intrapartum moderate neurological sequelae')
        self.intrapartum_severe_disu = \
            BetaValueGenerator(mean=0.84, stdv=0.05, note='disutility due to intrapartum severe neurological sequelae')

        self.maternal_disu_matrix = self._sample_by_utility_order(mil=self.maternal_mil_disu,
                                                                  mod=self.maternal_mod_disu,
                                                                  sev=self.maternal_sev_disu)

        ##############
        # INITIATION #
        ##############
        # incidence/rate
        self.incidence_sample = 0
        self.first_year_recur_rate_sample = 0
        self.avg_yearly_recur_reduction_sample = 0
        # probabilities
        self.intrauterine_sem_sample = None
        self.list_infection_types = []
        self.list_intrauterine_demise_or_infected = []
        self.list_intrauterine_long_term_outcomes = None
        self.list_intrapartum_short_term_outcomes = None
        self.list_intrapartum_long_term_sem = None
        self.list_intrapartum_long_term_cns = None
        self.list_intrapartum_long_term_dis = None
        # maternal outcome
        self.m_avg_age_pregnancy_sample = 0  # average age of pregnancy
        # neonatal outcome
        self.DicAcute = {}
        self.DicSequelae = {}
        self.qaly_death = 0

    def resample_by_distr(self, seed):
        i = int(seed)
        rng = np.random.RandomState(seed=i)
        ###############
        # PROBABILITY #
        ###############
        # mother-to-child transmission types
        infection_types = self.route_of_transmission.sample(rng=rng)
        prob_utero = infection_types[0] / (infection_types[0] + infection_types[1])
        prob_prenatal = infection_types[1] / (infection_types[0] + infection_types[1])
        self.list_infection_types = [prob_utero, prob_prenatal]
        # incidence of (intrauterine infection + intrapartum infection)
        self.incidence_sample = self.incidence.sample(rng=rng) * (infection_types[0] + infection_types[1])
        # intrauterine infection outcomes
        abortion_sample = self.abortion.sample(rng=rng)
        self.list_intrauterine_demise_or_infected = [abortion_sample, (1 - abortion_sample)]
        self.intrauterine_sem_sample = self.intrauterine_sem.sample(rng=rng)
        self.list_intrauterine_long_term_outcomes = self.intrauterine_long_term_outcomes.sample(rng=rng)
        # intrapartum infection outcomes
        self.list_intrapartum_short_term_outcomes = self.intrapartum_sem_cns_dis.sample(rng=rng)
        self.list_intrapartum_long_term_sem = self.intrapartum_long_term_after_sem.sample(rng=rng)
        self.list_intrapartum_long_term_cns = self.intrapartum_long_term_after_cns.sample(rng=rng)
        self.list_intrapartum_long_term_dis = self.intrapartum_long_term_after_dis.sample(rng=rng)

        ############
        # DURATION #
        ############
        neonate_sev_life_expec_sample = self.neonate_sev_life_expectancy.sample(rng=rng)

        ##############
        # DISUTILITY #
        ##############
        # neonatal
        intrauterine_neuro_disu_sample = self.intrauterine_neuro_disu.sample(rng=rng)
        intrapartum_mild_disu_sample = self.intrapartum_mild_disu.sample(rng=rng)
        intrapartum_moderate_disu_sample = self.intrapartum_moderate_disu.sample(rng=rng)
        intrapartum_severe_disu_sample = self.intrapartum_severe_disu.sample(rng=rng)

        #####################
        # Maternal outcomes #
        #####################
        # age and life expectancy
        maternal_dea_duration_sample = self.maternal_dea_duration.sample(rng=rng)  # duration: impact of losing a child
        maternal_mod_duration_sample = self.maternal_mod_initial_period_duration.sample(rng=rng)
        self.m_avg_age_pregnancy_sample = self.m_avg_age_pregnancy.sample(rng=rng)
        maternal_life_expec_sample = self.m_life_expectancy.sample(rng=rng)
        # disutilities
        maternal_mil_disu_sample = self.maternal_disu_matrix[i, 0]
        maternal_mod_disu_sample = self.maternal_disu_matrix[i, 1]
        maternal_sev_disu_sample = self.maternal_disu_matrix[i, 2]
        maternal_dea_disu_sample = self.maternal_dea_disu.sample(rng=rng)

        # Quality adjusted life expectancy for babies die before birth
        self.qaly_death = \
            life_table_get_long_term_loss_mixed_sex(age_onset=0, seq_disu=0, discount=self.discount, if_utility=True)

        ###################################
        # Dictionary of neonatal outcomes #
        ###################################
        # calculate maternal loss based on simulation length
        # QALYs loss due to death is a one-time loss, while QALYs loss of having an impaired child is life-time

        self.DicSequelae = {'nor': {'disu': 0,
                                    'maternal_disu': 0,
                                    'maternal_dura': 0},
                            'mil': {'disu': intrapartum_mild_disu_sample,
                                    'maternal_disu': maternal_mil_disu_sample,
                                    'maternal_dura': maternal_life_expec_sample},
                            'mod': {'disu': intrapartum_moderate_disu_sample,
                                    'maternal_disu': maternal_mod_disu_sample,
                                    'maternal_dura': maternal_mod_duration_sample},
                            'sev': {'disu': intrapartum_severe_disu_sample,
                                    'dura': neonate_sev_life_expec_sample,
                                    'maternal_disu': maternal_sev_disu_sample,
                                    'maternal_dura': min(neonate_sev_life_expec_sample, maternal_life_expec_sample)},
                            'dea': {'disu': 0,
                                    'maternal_disu': maternal_dea_disu_sample,
                                    'maternal_dura': maternal_dea_duration_sample},
                            'neuro': {'disu': intrauterine_neuro_disu_sample,
                                      'maternal_disu': maternal_mod_disu_sample,
                                      'maternal_dura': maternal_life_expec_sample}
                            }

    def get_intrauterine_qaly_loss(self):
        """ get total QALYs lost due to intrauterine HSV infection within simulation duration (sim_d)"""
        # long-term QALYs loss (normal, sequelae, death)
        # QALYs loss of death = total QALYs for a child with normal outcome (discounted)
        prob_outcomes = self.list_intrauterine_long_term_outcomes
        neuro_loss = life_table_get_long_term_loss_mixed_sex(age_onset=0, acute_sympt_disu=0,
                                                             seq_disu=self.DicSequelae['neuro']['disu'],
                                                             discount=self.discount, seq_dura=self.sim_time,
                                                             emort=0, if_utility=False)
        sequelae_and_death_qaly_loss = neuro_loss * prob_outcomes[1] + self.qaly_death * prob_outcomes[2]

        # maternal QALYs loss (having an impaired child, losing a child), dependent on simulation length
        m_loss_neuro = self.get_maternal_loss(neonate_outcome='neuro')
        m_loss_death = self.get_maternal_loss(neonate_outcome='dea')
        maternal_intrauterine_qaly_loss = m_loss_neuro * prob_outcomes[1] + m_loss_death * prob_outcomes[2]

        # SUM
        total_qaly_loss = sequelae_and_death_qaly_loss + maternal_intrauterine_qaly_loss
        return total_qaly_loss

    def get_intrapartum_qaly_loss(self, sequelae, utility=False, excess_mort=0, acute='type'):
        """ get total QALYs loss due to intrapartum HSV infection + SEM acute symptoms
        :param acute: acute infection type: SEM/CNS/DIS
        :param sequelae: severity of sequelae: (string) nor, mil, mod, sev -> normal, mild, moderate, severe
        :param utility: whether we calculate QALYs (utility=True) or QALYs loss (utility=False)
        :param excess_mort: (float) excess mortality with sequelae
        # :param sim_d: simulation length: (float) simulation length, sim_d = 0 means life time
        """
        # calculate QALYs loss based on sequelae type
        sim_time = self.sim_time
        if sequelae == 'sev':
            if sim_time > 0:
                sim_time = min(self.sim_time, self.DicSequelae['sev']['dura'])
        # Neonatal perspective
        loss = life_table_get_long_term_loss_mixed_sex(age_onset=0,
                                                       acute_sympt_disu=0,
                                                       seq_disu=self.DicSequelae[sequelae]['disu'],
                                                       discount=self.discount,
                                                       seq_dura=sim_time,
                                                       emort=excess_mort,
                                                       if_utility=utility)
        # maternal perspective
        m_loss = self.get_maternal_loss(neonate_outcome=sequelae)
        total_loss = loss + m_loss
        return total_loss

    def get_death_qaly_loss(self):
        """ return loss of both neonatal and maternal loss due to neonatal death of HSV """
        m_loss = self.get_maternal_loss(neonate_outcome='dea')
        return self.qaly_death + m_loss

    def get_maternal_loss(self, neonate_outcome):
        """ return (un)discounted/age-adjusted maternal QALYs loss due to neonatal herpes
        :param neonate_outcome: 'nor', 'mid', 'mod', 'sev', or 'dea'
        # :param sim_time: simulation length (0: lifetime, >0: actual number of years)
        """
        if neonate_outcome not in ['nor', 'mil', 'mod', 'sev', 'dea', 'neuro']:
            raise ValueError('Wrong type of neonatal outcome')
        pregnant_age = self.m_avg_age_pregnancy_sample
        if neonate_outcome == 'dea':
            # 0.6 years as 'grief period'
            adjusted_disu = get_adjusted_disu(current_age=pregnant_age, disu=self.DicSequelae['dea']['maternal_disu'])
            dea_duration = self.DicSequelae['dea']['maternal_dura']
            loss_grief_period = adjusted_disu * dea_duration
            # long-term loss (same as mothers with children having mild impairment
            seq_dura = self.DicSequelae['mil']['maternal_dura']
            seq_disu = self.DicSequelae['mil']['maternal_disu']
            if self.sim_time != 0:
                seq_dura = self.sim_time
            loss = life_table_get_long_term_loss(sex='female', age_onset=pregnant_age + dea_duration,
                                                 seq_disu=seq_disu, seq_dura=seq_dura, discount=self.discount)
            return loss_grief_period + loss
        elif neonate_outcome == 'mod':
            # moderately impaired child -> mothers also experience an initial 'shock' period
            # same logic as neonate_outcome == dea
            # 1 year as 'grief period'
            adjusted_disu = get_adjusted_disu(current_age=pregnant_age, disu=self.DicSequelae['dea']['maternal_disu'])
            mod_duration = self.DicSequelae['mod']['maternal_dura']
            loss_grief_period = adjusted_disu * mod_duration
            # long-term loss (same as mothers with children having mild impairment)
            seq_dura = self.DicSequelae['mil']['maternal_dura']
            seq_disu = self.DicSequelae['mil']['maternal_disu']
            if self.sim_time != 0:
                seq_dura = self.sim_time
            loss = life_table_get_long_term_loss(sex='female', age_onset=pregnant_age + mod_duration,
                                                 seq_disu=seq_disu, seq_dura=seq_dura, discount=self.discount)
            return loss_grief_period + loss
        else:
            seq_dura = self.DicSequelae[neonate_outcome]['maternal_dura']
            seq_disu = self.DicSequelae[neonate_outcome]['maternal_disu']
            if self.sim_time != 0:
                seq_dura = self.sim_time
            loss = life_table_get_long_term_loss(sex='female', age_onset=pregnant_age, seq_disu=seq_disu,
                                                 seq_dura=seq_dura, discount=self.discount)
            return loss

    def _sample_by_utility_order(self, mil, mod, sev):
        # random sample num_psa times
        seed_list = np.linspace(0, self.num_psa - 1, self.num_psa)
        mil_list = []
        mod_list = []
        sev_list = []
        for i in range(self.num_psa):
            rng = np.random.RandomState(seed=int(seed_list[i]))
            mil_list.append(mil.sample(rng=rng))
            mod_list.append(mod.sample(rng=rng))
            sev_list.append(sev.sample(rng=rng))

        # construct matrix
        matrix_u = np.asarray([mil_list, mod_list, sev_list]).T   # col1 mil, col2 mod, col3 sev
        matrix_q = np.asarray([[0, 0, 0], [-1, 0, 0], [-1, -1, 0]])

        epsilon = 0.05  # tolerance threshold
        delta = 0.01  # decrement correlation (rho)

        matrix_u_star = CorrelateUtils(matrix_u=matrix_u, matrix_q=matrix_q, epsilon=epsilon, delta=delta)
        return matrix_u_star
