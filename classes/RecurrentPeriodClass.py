from source.DemographicValues import get_conditional_survival_rate, get_adjusted_disu
from supports.ParameterAndRecurrentPeriodSupport import get_recur_initiate_rate_helper
import math


class RecurrentPeriod:
    def __init__(self, parameters, recur_type, age_of_infection, sex, discount=0.03):
        """ record all events happening during recurrent period and calculate QALY loss
        :param parameters: the object stored all attribute parameters for a specific cohort
        :param recur_type: (string) 'no', 'infrequent', 'frequent'
        :param age_of_infection: age of HSV infection
        :param sex: 'male' or 'female'
        :param discount: discount if specify discounting rate
        """
        # the first year of recurrences happen in the same year of primary infection
        self.age_of_infection = age_of_infection
        self.sex = sex
        self.discount = discount
        self.params = parameters
        self.recur_type = recur_type    # type of recurrence, infrequent vs. frequent
        self.sympt_recur_qaly = 0       # total QALYs loss attributed to symptomatic recurrences (+ complication)
        self.recur_only_qaly = 0        # total QALYs loss attributed to symptomatic recurrences (no complication)
        self.ur_qaly = 0                # total QALYs loss attributed to urinary retention
        self.psych_qaly = 0             # total QALYs loss attributed to long-term psychosocial impacts of diagnosis
        self.rate_of_recurrence = []    # rate of recurrences each year

    def _get_sympt_recur_loss(self, first_year_recur_rate, second_year_recur_rate, annual_reduction_rate):
        """
        calculate QALY loss due to symptomatic recurrences, which decreases in a linear fashion over time
        :param first_year_recur_rate: rate of symptomatic recurrences in the first year in recur period
        :param second_year_recur_rate: rate of symptomatic recurrences in the second year in recur period
        :param annual_reduction_rate: annual reduction in the rate of recurrence
        """
        ##################################################
        # Calendar 1: # of years after primary infection #
        ##################################################
        # cal_1 related events: current recurrence rate; utility discounting
        # initiate calender 1, first- and second- year recurrence rate
        cal_1 = 0  # initiate calender 1 (same as the exponent number for discounting equation)

        # number of years with sympt recurrences (first year is year_0)
        num_years_with_recur = int(math.ceil(second_year_recur_rate / annual_reduction_rate))

        # recurrences rate
        self.rate_of_recurrence = [first_year_recur_rate, second_year_recur_rate]
        for i in range(1, num_years_with_recur - 1):
            self.rate_of_recurrence.append(second_year_recur_rate - i * annual_reduction_rate)

        # the last year with sympt recurrences
        cal_1_end = num_years_with_recur - 1

        ###################################
        # Calendar 2: age in current year #
        ###################################
        # cal_2 related events: survive rate; adjusted for background utility
        cal_2 = self.age_of_infection  # age of primary infection (same year having first-year sympt recur)
        # the age when the person does not experience sympt recurrences = cal_2 + num_year_with_recur
        # for each iteration during the interval with symptomatic recurrences
        psycho_disu_a = self.params.diagnosis_disu_sample
        # one-time reduction
        psycho_disu_b = psycho_disu_a - self.params.diagnosis_disu_reduction_sample
        for i in range(cal_1, cal_1_end + 1):
            # surviving rate & recurrence rate
            if cal_1 == 0:
                sur_rate = 1
                psycho_disu = psycho_disu_a
            else:
                sur_rate = get_conditional_survival_rate(younger_age=self.age_of_infection,
                                                         older_age=cal_2,
                                                         sex=self.sex)
                psycho_disu = psycho_disu_b

            # recurrent rate
            recur_rate = self.rate_of_recurrence[cal_1]

            # discounting
            discount_index = cal_1

            # adjusted disutility
            recur_adjusted_disu = get_adjusted_disu(current_age=cal_2, disu=self.params.recur_sympt_treated_disu_sample)
            ur_adjusted_disu = get_adjusted_disu(current_age=cal_2, disu=self.params.urinary_retention_disu_sample)
            psycho_adjusted_disu = get_adjusted_disu(current_age=cal_2, disu=psycho_disu)

            # total QALYs lost at age cal_2 = s * d * (# of recur * QALYs_recur + # of recur * prob_ur * QALYs_ur)
            recur_qaly = recur_rate * recur_adjusted_disu * self.params.recur_treat_duration_sample
            ur_prob = self.params.urinary_retention_recur_prob_sample
            ur_qaly = recur_rate * ur_prob * ur_adjusted_disu * self.params.urinary_retention_duration_sample
            psycho_qaly = psycho_adjusted_disu * 1
            annual_qaly = sur_rate * (recur_qaly + ur_qaly + psycho_qaly) / (1 + self.discount) ** discount_index

            # update calendar times
            cal_1 += 1
            cal_2 += 1
            self.sympt_recur_qaly += annual_qaly    # total QALYs lost
            self.recur_only_qaly += sur_rate * recur_qaly / (1 + self.discount) ** discount_index
            self.ur_qaly += sur_rate * ur_qaly / (1 + self.discount) ** discount_index
            self.psych_qaly += sur_rate * psycho_qaly / (1 + self.discount) ** discount_index


class RecurrentPeriodTypeOne(RecurrentPeriod):
    def __init__(self, parameters, recur_type, age_of_infection, sex, discount=0.03):
        super().__init__(parameters=parameters, recur_type=recur_type, age_of_infection=age_of_infection,
                         sex=sex, discount=discount)
        if self.recur_type == 'frequent':
            raise ValueError('Wrong type of recurrence for HSV-1 infection')

        self.total_loss_hsv1 = 0        # total QALYs loss due to HSV-1 recurrence + encephalitis
        self.ence_qaly = 0              # total QALYs loss attributed to encephalitis at age 60 years old

        # calculate total QALYs loss
        if self.recur_type == 'infrequent':
            # self._add_qaly_loss_of_psychosocial_impact()   # only people with sympt recur have long psychosocial loss
            self._add_qaly_loss_of_sympt_recur()
        # if ence:
        #     self._add_qaly_loss_of_ence()

    def _add_qaly_loss_of_sympt_recur(self):
        first_year_rate = self.params.infreq_first_year_recur_rate_sample  # recur rate in year 0 (after primary hsv)
        second_year_rate = self.params.infreq_second_year_recur_rate_sample  # recur rate in year 1
        annual_reduction_rate = self.params.avg_yearly_recur_reduction_sample
        self._get_sympt_recur_loss(first_year_recur_rate=first_year_rate, second_year_recur_rate=second_year_rate,
                                   annual_reduction_rate=annual_reduction_rate)
        # sum up to total QALYs loss
        self.total_loss_hsv1 += self.sympt_recur_qaly

    # def _add_qaly_loss_of_ence(self, ence_age=61):
    #     """ can probably happen at age of 61
    #     1. QALYs lost due to acute symptoms (at age of 61)
    #     2. QALYs lost due to long-term sequelae (with a shortened life expectancy for severe group)
    #     3. QALYs lost for those dead people (with discounted life expectancy of normal people)
    #     """
    #     nor_long_term_loss = life_table_get_long_term_loss(sex=self.sex, age_onset=ence_age, discount=self.discount,
    #                                                        acute_sympt_disu=self.params.ence_disu_sample, seq_disu=0)
    #     mil_long_term_loss = life_table_get_long_term_loss(sex=self.sex, age_onset=ence_age, discount=self.discount,
    #                                                        acute_sympt_disu=self.params.ence_disu_sample,
    #                                                        seq_disu=self.params.ence_mil_disu_sample)
    #     mod_long_term_loss = life_table_get_long_term_loss(sex=self.sex, age_onset=ence_age, discount=self.discount,
    #                                                        acute_sympt_disu=self.params.ence_disu_sample,
    #                                                        seq_disu=self.params.ence_mod_disu_sample)
    #     dea_long_term_loss = life_table_get_long_term_loss(sex=self.sex, age_onset=ence_age, discount=self.discount,
    #                                                        acute_sympt_disu=self.params.ence_disu_sample,
    #                                                        seq_disu=0, if_utility=True)
    #
    #     sev_life_e = self.params.ence_sev_life_expectancy_sample
    #     sev_disu = self.params.ence_sev_disu_sample
    #     num_years_after_infection = ence_age - self.age_of_infection
    #     sev_long_term_loss = get_long_term_qaly_loss(start_age=ence_age, length_of_condition=sev_life_e,
    #                                                  disease_disu=sev_disu, year_idx=num_years_after_infection,
    #                                                  discount=self.discount)
    #
    #     # parameters
    #     risk_of_ence = self.params.ence_elderly_prob_sample
    #     # for people developed with encephalitis (a small probabilistic tree)
    #     # // key: [cost, disutility, component_dic, [future nodes], [probabilities]
    #     dictChancesEnce = {'E0': [0, 0, {}, ['t0', 'E1'], [1 - risk_of_ence, risk_of_ence]],
    #                        'E1': [0, 0, {}, ['t1', 't2', 't3', 't4', 't5'], self.params.ence_sequelae_prob_sample]
    #                        }
    #     # // key: [cost, disutility, component_dic]
    #     dictTerminalsEnce = {'t0': [0, 0, {}],
    #                          't1': [0, nor_long_term_loss, {'ence': nor_long_term_loss}],
    #                          't2': [0, mil_long_term_loss, {'ence': mil_long_term_loss}],
    #                          't3': [0, mod_long_term_loss, {'ence': mod_long_term_loss}],
    #                          't4': [0, sev_long_term_loss, {'ence': sev_long_term_loss}],
    #                          't5': [0, dea_long_term_loss, {'ence': dea_long_term_loss}],
    #                          }
    #     # // key: [cost, disutility, component_dic, [future nodes]]
    #     dictDecisionsEnce = {'d1': [0, 0, {}, ['E0']]}
    #     # build the decision tree
    #     myDT = dt.DecisionNode('d1', dictDecisionsEnce, dictChancesEnce, dictTerminalsEnce)
    #     myDT.evaluate()
    #     cost_utility = myDT.get_cost_utility()
    #
    #     # get total QALY loss due to encephalitis (discounting)
    #     self.ence_qaly = cost_utility['E0'][-1]
    #     # sum up to total QALYs loss
    #     self.total_loss_hsv1 += self.ence_qaly


class RecurrentPeriodTypeTwo(RecurrentPeriod):
    def __init__(self, parameters, recur_type, age_of_infection, sex, rm=False, discount=0.03, cst=False):
        """
        :param parameters: (object) parameters
        :param recur_type: (string) 'no', 'infrequent', 'frequent'
        :param age_of_infection: (int) age of infection
        :param sex: (string) 'male' or 'female'
        :param rm: (bool) whether have recurrent meningitis
        :param discount: (float) discounting rate
        :param cst: (bool) whether receive CST, only useful when recur_type is 'frequent'
        """
        super().__init__(parameters=parameters, recur_type=recur_type, age_of_infection=age_of_infection,
                         sex=sex, discount=discount)
        self.long_term_therapy = cst
        # total QALYs loss due to infrequent HSV-2 recurrence (no recurrent meningitis)
        self.total_loss_hsv2 = 0
        self.rm_qaly = 0                            # total QALYs loss attributed to recurrent meningitis

        # calculate total QALYs loss
        if recur_type in ['infrequent', 'frequent']:
            self._add_qaly_loss_of_sympt_recur()
        # people with no symptomatic recurrences can also experience recurrent meningitis
        if rm:
            self._add_qaly_loss_of_recur_meningitis()

    def _add_qaly_loss_of_sympt_recur(self):
        first_year_rate, second_year_rate = get_recur_initiate_rate_helper(recur_type=self.recur_type,
                                                                           parameter=self.params,
                                                                           long_term_therapy=self.long_term_therapy)
        annual_reduction_rate = self.params.avg_yearly_recur_reduction_sample
        self._get_sympt_recur_loss(first_year_recur_rate=first_year_rate, second_year_recur_rate=second_year_rate,
                                   annual_reduction_rate=annual_reduction_rate)
        # sum up to total QALYs loss
        self.total_loss_hsv2 += self.sympt_recur_qaly

    def _get_recur_meningitis_loss(self):
        """ calculate QALYs lost due to recurrent meningitis, which is independent of the presence of HSV relapses """
        # (avg.) total number of years that people experience meningitis relapses
        total_year_with_meningitis = self.params.total_year_with_meningitis_sample
        # (avg.) total number of meningitis people will experience
        total_num_meningitis = round(self.params.total_num_recur_meningitis_sample)
        # duration of each recurrent meningitis
        rm_duration = self.params.recur_meningitis_duration_sample
        # unadjusted disutility for each recurrent meningitis
        rm_unadjusted_disu = self.params.recur_meningitis_disu_sample
        # assuming the intermittent time between two relapses are the same
        # initiate calendar time
        cal_1 = total_year_with_meningitis / total_num_meningitis  # (constant) time to the next recurrent meningitis
        cal_2 = self.age_of_infection  # age when having a meningitis

        # calculate QALYs lost for each iteration
        for iterate in range(0, total_num_meningitis):
            cal_2_prev_round = round(cal_2)  # (int) rounded age used for survival rate calculation
            cal_2 += cal_1  # (float) exact age having one recurrent meningitis
            cal_2_round = round(cal_2)  # (int) rounded age used for survival rate and discounting calculation
            sur_rate = get_conditional_survival_rate(younger_age=cal_2_prev_round, older_age=cal_2_round,
                                                     sex=self.sex)
            rm_adjusted_disu = get_adjusted_disu(current_age=cal_2_round, disu=rm_unadjusted_disu)
            # QALYs lost = s(t) * d(t) * disutility * duration
            per_rm_qaly_loss = sur_rate * rm_adjusted_disu * rm_duration / (1 + self.discount) ** (
                        cal_2_round - self.age_of_infection)
            self.rm_qaly += per_rm_qaly_loss

    def _add_qaly_loss_of_recur_meningitis(self):
        self._get_recur_meningitis_loss()
        # sum up to total QALYs loss
        self.total_loss_hsv2 += self.rm_qaly
