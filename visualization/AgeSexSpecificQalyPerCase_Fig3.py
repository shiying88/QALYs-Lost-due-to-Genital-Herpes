from supports.VisualizationSupport import *

# outcomes: age_sex_specific_qaly, non_age_sex_qaly_avg, t_qaly_f_avg, t_qaly_m_avg, t_qaly_hsv_avg
read_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/trees/dics/'
output_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/visualization/graphs/'
# read python dict back from the file
pkl_file = open('{}DicHSV1_dr3_psycho1_noGBD_20221010_dN100.pkl'.format(read_path), 'rb')
hsv1_Dic = pickle.load(pkl_file)
pkl2_file = open('{}DicHSV2_dr3_psycho1_noGBD_20221010_dN100.pkl'.format(read_path), 'rb')
hsv2_Dic = pickle.load(pkl2_file)
pkl_file.close()
pkl2_file.close()

# Age and sex specific QALYs loss per infection
mean_male1, ci_male1, mean_female1, ci_female1 = get_age_sex_qaly_loss(dic=hsv1_Dic)
mean_male2, ci_male2, mean_female2, ci_female2 = get_age_sex_qaly_loss(dic=hsv2_Dic)

# draw plots
width = 0.2
bar_alpha = 1
hsv2_female_c = '#2596be'
hsv2_male_c = lighten_color(hsv2_female_c, 0.6)
hsv1_female_c = '#5D62D5'
hsv1_male_c = lighten_color(hsv1_female_c, 0.6)
# error bars
err_dic_male1 = {'capsize': 3, 'ecolor': lighten_color(hsv1_male_c, 1.3)}
err_dic_female1 = {'capsize': 3, 'ecolor': lighten_color(hsv1_female_c, 1.3)}
err_dic_male2 = {'capsize': 3, 'ecolor': lighten_color(hsv2_male_c, 1.3)}
err_dic_female2 = {'capsize': 3, 'ecolor': lighten_color(hsv2_female_c, 1.3)}
# other format options
x_pos = np.arange(len(age_list))
x_ticks = np.arange(5)
age_text_list = ['18-24', '25-29', '30-34', '35-49']

# figure size
plt.figure(figsize=(5.5, 5.5))
# bar plot
plt.bar(x_pos - 0.2, mean_male1, width=width, color=hsv1_male_c, alpha=bar_alpha, label='HSV-1, men',
        yerr=[ci / 4 for ci in ci_male1], error_kw=err_dic_male1)
plt.bar(x_pos, mean_female1, width=width, color=hsv1_female_c, alpha=bar_alpha, label='HSV-1, women',
        yerr=[ci / 4 for ci in ci_female1], error_kw=err_dic_female1)
plt.bar(x_pos + 0.2, mean_male2, width=width, color=hsv2_male_c, alpha=bar_alpha, label='HSV-2, men',
        yerr=[ci / 4 for ci in ci_male2], error_kw=err_dic_male2)
plt.bar(x_pos + 0.4, mean_female2, width=width, color=hsv2_female_c, alpha=bar_alpha, label='HSV-2, women',
        yerr=[ci / 4 for ci in ci_female2], error_kw=err_dic_female2)
plt.rcParams["font.family"] = "Times New Roman"  # set font for legend
# x ticks text
plt.xticks(x_pos, labels=age_text_list, rotation=0)
plt.yticks(fontname='Times New Roman')
plt.xticks(fontname='Times New Roman')
# legend
plt.legend(loc='upper right', fontsize=10)
# grid
plt.grid(axis='y', alpha=0.5, linestyle='--')
# labels and title
plt.xlabel('Age groups (years)', fontsize=10, **{'fontname': 'Times New Roman'})
plt.ylabel('QALYs lost / infection', fontsize=10, **{'fontname': 'Times New Roman'})

plt.title('Age and sex specific QALYs lost per HSV infection (2018)',
          fontname='Times New Roman',
          fontweight='bold',
          fontsize=10)
# save file
plt.tight_layout()
plt.savefig('{}main_age-sex-specific_per_case_20221010.png'.format(output_path), dpi=350, bbox_inches='tight')
plt.show()
