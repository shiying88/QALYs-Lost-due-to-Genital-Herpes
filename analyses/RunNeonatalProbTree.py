import numpy as np
import pickle
from classes import DecisionTreeClass as dt
from classes.ParameterClass import ParametersNeonatal
from classes.ProbTreeClasses import buildNeonatalTree
from SimPy.Statistics import SummaryStat

# Parameter initiation
DISCOUNT = 0.03     # (float) discounting rate
NUM_PSA = 1000      # (int) number of repetition for PSA

###################
# NEONATAL HERPES #
###################
SIM_TIME = 15        # (int) 0=life time, ≥1: exact years for simulation duration

# specify parameter distributions
params_neonate = ParametersNeonatal(discount=DISCOUNT, sim_time=SIM_TIME)
# initiating results
incidence_list = []                 # incidence
neonate_qaly_loss_list = []         # neonatal + maternal: per case QLAYs lost
neonate_total_qaly_loss_list = []   # neonatal + maternal: total QLAYs lost
maternal_qaly_loss_list = []        # maternal: per case QALYs lost
maternal_total_qaly_loss_list = []  # maternal: total QALYs lost

# Main Analysis and Probability Sensitivity Analysis
for i in range(NUM_PSA):
    print('current number of iterations:', i)
    # random number generator
    rng = np.random.RandomState(seed=i)
    # resample parameters
    params_neonate.resample_by_distr(seed=i)
    # incidence
    incidence = params_neonate.incidence_sample/100000 * 3791712
    incidence_list.append(incidence)
    # construct probabilistic tree
    dictDecisions, dictChances, dictTerminals, dictTerminals_maternal = buildNeonatalTree(params=params_neonate)
    # calculate expected QALYs loss for neonatal + maternal
    myDT = dt.DecisionNode('d1', dictDecisions, dictChances, dictTerminals)
    myDT.evaluate()
    cost_utility = myDT.get_cost_utility()
    neonate_qaly_loss_list.append(cost_utility['c0'][-1])
    neonate_total_qaly_loss_list.append(incidence * cost_utility['c0'][-1])
    # calculate expected QALYs loss of maternal
    myDT_maternal = dt.DecisionNode('d1', dictDecisions, dictChances, dictTerminals_maternal)
    myDT_maternal.evaluate()
    cost_utility_maternal = myDT_maternal.get_cost_utility()
    maternal_qaly_loss = cost_utility_maternal['c0'][-1]
    maternal_qaly_loss_list.append(maternal_qaly_loss)
    maternal_total_qaly_loss_list.append(incidence * maternal_qaly_loss)

print('simulation length', SIM_TIME)
print('discount rate', DISCOUNT)
# incidence
incidence_stat = SummaryStat(name='incidence of neonatal herpes in 2018', data=incidence_list)
incidence_avg = incidence_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
print('total number of births infected by HSV in 2018:', incidence_avg)
# neonatal + maternal
qaly_loss_per_infection_stat = SummaryStat(name='QALYs loss due to neonatal herpes', data=neonate_qaly_loss_list)
qaly_loss_total_stat = SummaryStat(name='Total QALYs loss due to neonatal herpes', data=neonate_total_qaly_loss_list)
qaly_loss_per_infection_avg = qaly_loss_per_infection_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
qaly_loss_total_avg = qaly_loss_total_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
print('The average and 95% CI of QALYs loss due to per neonatal herpes:', qaly_loss_per_infection_avg)
print('The average and 95% CI of QALYs loss due to neonatal herpes in 2018:', qaly_loss_total_avg)
print('quality-adjusted life expectancy for stillbirth', round(params_neonate.qaly_death, 2))
# maternal
maternal_qaly_loss_per_infection_stat = SummaryStat(name='maternal QALYs loss due to per neonatal herpes',
                                                    data=maternal_qaly_loss_list)
maternal_qaly_loss_total_stat = SummaryStat(name='maternal QALYs loss due to neonatal herpes in 2018',
                                            data=maternal_total_qaly_loss_list)
maternal_qaly_loss_per_infection_avg = \
    maternal_qaly_loss_per_infection_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
maternal_qaly_loss_total_avg = \
    maternal_qaly_loss_total_stat.get_formatted_mean_and_interval(deci=2, interval_type='p')
print('The average and 95% CI maternal QALYs loss due to per neonatal herpes:', maternal_qaly_loss_per_infection_avg)
print('The average and 95% CI maternal QALYs loss due to neonatal herpes in 2018:', maternal_qaly_loss_total_avg)

# output dictionary
N_dic = {'simulation length': SIM_TIME,                                 # simulation length
         'num_case': incidence_avg,                                     # incidence
         'loss_per_infection': qaly_loss_per_infection_avg,             # QALYs lost per infection: neonatal+maternal
         'loss_total': qaly_loss_total_avg,                             # QALYs lost total: neonatal+maternal
         'm_loss_per_infection': maternal_qaly_loss_per_infection_avg,  # QALYs lost per infection: maternal
         'm_loss_total': maternal_qaly_loss_total_avg}                  # QALYs lost total: maternal
# save output
output = open('tree_outputs/dics/neonatalGH_{}.pkl'.format(SIM_TIME), 'wb')
pickle.dump(N_dic, output)
output.close()
