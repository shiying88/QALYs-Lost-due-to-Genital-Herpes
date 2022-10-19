import pandas as pd
from SimPy.Statistics import SummaryStat
from supports.VisualizationSupport import *

# read dataset
read_path = '/tree_outputs/component_utl/'
output_path = 'visualization/graphs/'
df_hsv1 = pd.read_csv('{}hsv1_dr3_psycho1_noGBD_20221010_dN100.csv'.format(read_path))
df_hsv2 = pd.read_csv('{}hsv2_dr3_psycho1_noGBD_20221010_dN100.csv'.format(read_path))

# enumerate component names, sex list, and age lists
sex_list = ['male', 'female']
age_list = ['21', '27', '32', '42']
components_hsv1 = ['psych', 'recur', 'ur', 'primary', 'am']
components_hsv2 = ['primary', 'am', 'ur', 'psych', 'infreq', 'freq_no_cst', 'freq_cst', 'rm']

summary_stat_dic_hsv1 = {'psych': {'male': [], 'female': []},
                         'recur': {'male': [], 'female': []},
                         'ur': {'male': [], 'female': []},
                         'primary': {'male': [], 'female': []},
                         'am': {'male': [], 'female': []}}
summary_stat_dic_hsv2 = {'rm': {'male': [], 'female': []},
                         'psych': {'male': [], 'female': []},
                         'infreq': {'male': [], 'female': []},
                         'freq_no_cst': {'male': [], 'female': []},
                         'freq_cst': {'male': [], 'female': []},
                         'ur': {'male': [], 'female': []},
                         'primary': {'male': [], 'female': []},
                         'am': {'male': [], 'female': []}}

# save statistics into summary dictionary
# HSV-1
for name in components_hsv1:
    print('component:', name)
    for sex in sex_list:
        print('sex:', sex)
        df_sex_com = df_hsv1[df_hsv1["sex"] == sex][[name, 'age']]
        for age in age_list:
            print('age:', age)
            age_sex_com_list = df_sex_com[df_sex_com['age'] == int(age)][name].tolist()
            summary_stat = SummaryStat(name='{}_{}_{}_QALYs'.format(name, sex, age), data=age_sex_com_list)
            summary_stat_dic_hsv1[name][sex].append(summary_stat.get_mean())
# HSV-2
for name in components_hsv2:
    print('component:', name)
    for sex in sex_list:
        print('sex:', sex)
        df_sex_com = df_hsv2[df_hsv2["sex"] == sex][[name, 'age']]
        for age in age_list:
            print('age:', age)
            age_sex_com_list = df_sex_com[df_sex_com['age'] == int(age)][name].tolist()
            summary_stat = SummaryStat(name='{}_{}_{}_QALYs'.format(name, sex, age), data=age_sex_com_list)
            summary_stat_dic_hsv2[name][sex].append(summary_stat.get_mean())


# plot stacked bar charts
# create data
x = ['18-24', '25-29', '30-34', '35-49']
# HSV-1
primary_male_hsv1 = np.array(summary_stat_dic_hsv1['primary']['male'])
psych_male_hsv1 = np.array(summary_stat_dic_hsv1['psych']['male'])
recur_male_hsv1 = np.array(summary_stat_dic_hsv1['recur']['male'])
ur_male_hsv1 = np.array(summary_stat_dic_hsv1['ur']['male'])
am_male_hsv1 = np.array(summary_stat_dic_hsv1['am']['male'])
primary_female_hsv1 = np.array(summary_stat_dic_hsv1['primary']['female'])
psych_female_hsv1 = np.array(summary_stat_dic_hsv1['psych']['female'])
recur_female_hsv1 = np.array(summary_stat_dic_hsv1['recur']['female'])
ur_female_hsv1 = np.array(summary_stat_dic_hsv1['ur']['female'])
am_female_hsv1 = np.array(summary_stat_dic_hsv1['am']['female'])
# HSV-2
primary_male_hsv2 = np.array(summary_stat_dic_hsv2['primary']['male'])
psych_male_hsv2 = np.array(summary_stat_dic_hsv2['psych']['male'])
infreq_male_hsv2 = np.array(summary_stat_dic_hsv2['infreq']['male'])
freq_no_cst_male_hsv2 = np.array(summary_stat_dic_hsv2['freq_no_cst']['male'])
freq_cst_male_hsv2 = np.array(summary_stat_dic_hsv2['freq_cst']['male'])
ur_male_hsv2 = np.array(summary_stat_dic_hsv2['ur']['male'])
am_male_hsv2 = np.array(summary_stat_dic_hsv2['am']['male'])
rm_male_hsv2 = np.array(summary_stat_dic_hsv2['rm']['male'])
primary_female_hsv2 = np.array(summary_stat_dic_hsv2['primary']['female'])
psych_female_hsv2 = np.array(summary_stat_dic_hsv2['psych']['female'])
infreq_female_hsv2 = np.array(summary_stat_dic_hsv2['infreq']['female'])
freq_no_cst_female_hsv2 = np.array(summary_stat_dic_hsv2['freq_no_cst']['female'])
freq_cst_female_hsv2 = np.array(summary_stat_dic_hsv2['freq_cst']['female'])
ur_female_hsv2 = np.array(summary_stat_dic_hsv2['ur']['female'])
am_female_hsv2 = np.array(summary_stat_dic_hsv2['am']['female'])
rm_female_hsv2 = np.array(summary_stat_dic_hsv2['rm']['female'])

# plot statistics
ALPHA = 0.8
WIDTH = 0.6
# primary_color = 'burlywood'
primary_color = 'orange'
psych_color = 'olivedrab'
recur_color = 'lightblue'  # recurrence for HSV-1 and infre recur for HSV-2
freq_no_cst_color = 'royalblue'

# four subplots
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))

# plot bars in stack manner
# HSV-1
# male
ax[0][0].bar(x, psych_male_hsv1, color=psych_color, alpha=ALPHA, width=WIDTH, label='Psycho')
ax[0][0].bar(x, recur_male_hsv1, bottom=psych_male_hsv1, color=recur_color, alpha=ALPHA, width=WIDTH, label='Recur')
ax[0][0].bar(x, primary_male_hsv1 + ur_male_hsv1 + am_male_hsv1, bottom=psych_male_hsv1 + recur_male_hsv1,
             color=primary_color, alpha=ALPHA, width=WIDTH, label='All other')
# female
ax[0][1].bar(x, psych_female_hsv1, color=psych_color, alpha=ALPHA, width=WIDTH, label='Psycho')
ax[0][1].bar(x, recur_female_hsv1, bottom=psych_female_hsv1, color=recur_color, alpha=ALPHA, width=WIDTH, label='Recur')
ax[0][1].bar(x, primary_female_hsv1 + ur_female_hsv1 + am_female_hsv1, bottom=psych_female_hsv1 + recur_female_hsv1,
             color=primary_color, alpha=ALPHA, width=WIDTH, label='All other')
# HSV-2
# male
ax[1][0].bar(x, psych_male_hsv2, color=psych_color, alpha=ALPHA, width=WIDTH, label='Psycho')
ax[1][0].bar(x, infreq_male_hsv2, bottom=psych_male_hsv2, color=recur_color, alpha=ALPHA, width=WIDTH, label='Infreq')
ax[1][0].bar(x, freq_no_cst_male_hsv2, bottom=psych_male_hsv2+infreq_male_hsv2,
             color=freq_no_cst_color, alpha=ALPHA, width=WIDTH, label='Freq (episodic T)')
ax[1][0].bar(x, primary_male_hsv2+freq_cst_male_hsv2+ur_male_hsv2+am_male_hsv2+rm_male_hsv2,
             bottom=psych_male_hsv2+infreq_male_hsv2+freq_no_cst_male_hsv2,
             color=primary_color, alpha=ALPHA, width=WIDTH, label='All other')
# female
ax[1][1].bar(x, psych_female_hsv2, color=psych_color, alpha=ALPHA, width=WIDTH, label='Psycho')
ax[1][1].bar(x, infreq_female_hsv2, bottom=psych_female_hsv2, color=recur_color, alpha=ALPHA, width=WIDTH, label='Infreq')
ax[1][1].bar(x, freq_no_cst_female_hsv2, bottom=psych_female_hsv2+infreq_female_hsv2,
             color=freq_no_cst_color, alpha=ALPHA, width=WIDTH, label='Freq (episodic T)')
ax[1][1].bar(x, primary_female_hsv2+freq_cst_female_hsv2+ur_female_hsv2+am_female_hsv2+rm_female_hsv2,
             bottom=psych_female_hsv2+infreq_female_hsv2+freq_no_cst_female_hsv2,
             color=primary_color, alpha=ALPHA, width=WIDTH, label='All other')

# Figure formatting
plt.rcParams["font.family"] = "Times New Roman"  # set legend font
for i in range(2):
    for j in range(2):
        ax[i][j].set_xlabel('Age groups', **{'fontname': 'Times New Roman'}, fontsize=10)  # set x-label font
        handles_hsv1, labels_hsv1 = ax[0][j].get_legend_handles_labels()
        handles_hsv2, labels_hsv2 = ax[1][j].get_legend_handles_labels()
        ax[0][j].legend(handles_hsv1[::-1], labels_hsv1[::-1], loc='upper right')
        ax[1][j].legend(handles_hsv2[::-1], labels_hsv2[::-1], loc='lower right')
        ax[i][j].grid(axis='y', alpha=0.5, linestyle='--')
        ax[i][j].set_ylim([0, 0.07])
        ax[i][j].set_xticklabels(['18-24', '25-29', '30-34', '35-49'],  # set x tick font
                                 fontname='Times New Roman', fontsize=10)
        ax[i][j].set_yticklabels(['0.00', '0.01', '0.02', '0.03', '0.04', '0.05', '0.06', '0.07'],  # set y tick font
                                 fontname='Times New Roman', fontsize=10)
        ax[i][j].set_ylabel('QALYs lost', **{'fontname': 'Times New Roman'}, fontsize=10)  # set y-label font

ax[0][0].set_title('Genital HSV-1 infection (Men)', fontsize=10, fontname='Times New Roman')
ax[0][1].set_title('Genital HSV-1 infection (Women)', fontsize=10, fontname='Times New Roman')
ax[1][0].set_title('HSV-2 infection (Men)', fontsize=10, fontname='Times New Roman')
ax[1][1].set_title('HSV-2 infection (Women)', fontsize=10, fontname='Times New Roman')
plt.suptitle('Breakdowns of QALYs lost per HSV infection by sex and viral type',
             fontsize=10, fontname='Times New Roman', fontweight='bold')

# save file
plt.tight_layout(rect=(0, 0, 1, 0.95))
plt.savefig('{}Fig4-breakdown_utilities.png'.format(output_path), dpi=600,
            bbox_inches='tight')
plt.show()
