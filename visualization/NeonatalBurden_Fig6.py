import matplotlib.pyplot as plt
import numpy as np
from supports.VisualizationSupport import *

neonatal_color = '#5D62D5'
maternal_color = '#2596be'

# outcomes: age_sex_specific_qaly, non_age_sex_qaly_avg, t_qaly_f_avg, t_qaly_m_avg, t_qaly_hsv_avg
output_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/visualization/graphs/'

sim_length = ['5 years', '10 years', '15 years', 'lifetime']

# Neonatal + Maternal
per_neonatal_loss_mean = [4.31, 4.99, 5.51, 7.93]
per_neonatal_loss_min = [3.54, 4.14, 4.56, 6.63]
per_neonatal_loss_max = [5.13, 5.9, 6.52, 9.19]
t_neonatal_loss_mean = [1707.62, 1976.08, 2182.26, 3136.39]
t_neonatal_loss_min = [1211.32, 1400.02, 1567.27, 2261.55]
t_neonatal_loss_max = [2294.82, 2666.8, 2921.75, 4143.99]
# Maternal
per_maternal_loss_mean = [0.21, 0.37, 0.51, 0.77]
per_maternal_loss_min = [0.12, 0.22, 0.3, 0.42]
per_maternal_loss_max = [0.31, 0.57, 0.8, 1.22]
t_maternal_loss_mean = [81.62, 148.17, 203.58, 305.03]
t_maternal_loss_min = [43.88, 78.45, 109.58, 153.36]
t_maternal_loss_max = [132.95, 243.98, 333.67, 510.98]


ci_per_neonatal = []
ci_t_neonatal = []
ci_per_maternal = []
ci_t_maternal = []
for i in range(len(per_maternal_loss_max)):
    ci_per_neonatal.append(per_neonatal_loss_max[i] - per_neonatal_loss_min[i])
    ci_t_neonatal.append(t_neonatal_loss_max[i] - t_neonatal_loss_min[i])
    ci_per_maternal.append(per_maternal_loss_max[i] - per_maternal_loss_min[i])
    ci_t_maternal.append(t_maternal_loss_max[i] - t_maternal_loss_min[i])

x_pos = np.arange(4)
width = 0.6

err_c_scalar = 1.4
capsize = 3

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 4.5))

ax[0].bar(x_pos, per_neonatal_loss_mean, width=width, color=neonatal_color, label='Neonatal',
          yerr=[ci/4 for ci in ci_per_neonatal], ecolor=lighten_color(neonatal_color, err_c_scalar), capsize=capsize)
ax[0].bar(x_pos, per_maternal_loss_mean, width=width, color=maternal_color, label='Maternal',
          yerr=[ci/4 for ci in ci_per_maternal], ecolor=lighten_color(maternal_color, err_c_scalar), capsize=capsize)

ax[1].bar(x_pos, t_neonatal_loss_mean, width=width, color=neonatal_color, label='Neonatal',
          yerr=[ci/4 for ci in ci_t_neonatal], ecolor=lighten_color(neonatal_color, err_c_scalar), capsize=capsize)
ax[1].bar(x_pos, t_maternal_loss_mean, width=width, color=maternal_color, label='Maternal',
          yerr=[ci/4 for ci in ci_t_maternal], ecolor=lighten_color(maternal_color, err_c_scalar), capsize=capsize)

plt.rcParams["font.family"] = "Times New Roman"  # set legend font

for i in range(2):
    ax[i].set_xticks(x_pos)
    ax[i].set_xticklabels(sim_length, fontname='Times New Roman', fontsize=10)
    ax[i].legend(loc='upper left')
    ax[i].set_xlabel('Simulation Length', **{'fontname': 'Times New Roman'}, fontsize=10)  # set x-label font
    ax[i].set_ylabel('QALYs lost', **{'fontname': 'Times New Roman'}, fontsize=10)  # set x-label font
    ax[i].grid(axis='y', alpha=0.5, linestyle='--')
ax[0].set_yticklabels(np.arange(0, 10, 1), fontname='Times New Roman', fontsize=10)
ax[1].set_yticklabels(np.arange(0, 4500, 500), fontname='Times New Roman', fontsize=10)
ax[0].set_title('QALYs lost of one neonatal HSV infection (2018)',
                fontname='Times New Roman', fontsize=10, fontweight='bold')
ax[1].set_title('Total QALYs lost of neonatal HSV in the US (2018)',
                fontname='Times New Roman', fontsize=10, fontweight='bold')

# save file
plt.tight_layout()
plt.savefig('{}neonatal_burden_20221010.png'.format(output_path), dpi=350, bbox_inches='tight')
plt.show()

