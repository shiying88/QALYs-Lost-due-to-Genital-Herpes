from supports.VisualizationSupport import *

read_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/trees/dics/'
output_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/visualization/graphs/'

hsv1_main_per_case = 0.01403
hsv1_main_total = 1925.52
hsv1_outcomes = [[0.01359, 0.01452, 0.0083, 0.02593, 0.00633],
                 [1864.9, 1992.89, 1142.59, 3542.85, 879.93]]
hsv2_main_per_case = 0.05214
hsv2_main_total = 31221.58
hsv2_outcomes = [[0.04795, 0.05755, 0.02137, 0.08493, 0.03162],
                 [28710.51, 34466.63, 12827.9, 50830.32, 18674.25]]

# x-axis format
label_text = ['0.06 discount rate', '0 discount rate',
              'GBD symptomatic outbreak disutility',
              '0.02 psychosocial disutility', '0 psychosocial disutility']
x_pos = np.arange(len(label_text))
# bar format
width = 0.8
bar_alpha = 1
hsv1_c = '#5D62D5'
hsv2_c = '#2596be'
base_case_c = 'sienna'
base_case_ls = '--'
# error bar
capsize = 3
hsv1_e_c = lighten_color(hsv1_c, 1.3)
hsv2_e_c = lighten_color(hsv2_c, 1.3)

# figure size
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 6))

# plot 1: HSV-1 per case outcomes
ax[0][0].barh(x_pos, hsv1_outcomes[0], height=width, color=hsv1_c, alpha=bar_alpha)
ax[0][0].axvline(x=hsv1_main_per_case, color=base_case_c, linestyle=base_case_ls)

# plot 2: HSV-2 per case outcomes
ax[0][1].barh(x_pos, hsv2_outcomes[0], height=width, color=hsv2_c, alpha=bar_alpha)
ax[0][1].axvline(x=hsv2_main_per_case, color=base_case_c, linestyle=base_case_ls)

# plot 3: HSV-1 total outcomes
ax[1][0].barh(x_pos, [x / 1000 for x in hsv1_outcomes[1]], height=width, color=hsv1_c, alpha=bar_alpha)
ax[1][0].axvline(x=hsv1_main_total / 1000, color=base_case_c, linestyle=base_case_ls)

# plot 4: HSV-2 total outcomes
ax[1][1].barh(x_pos, [x / 1000 for x in hsv2_outcomes[1]], height=width, color=hsv2_c, alpha=bar_alpha)
ax[1][1].axvline(x=hsv2_main_total / 1000, color=base_case_c, linestyle=base_case_ls)

plt.rcParams["font.family"] = "Times New Roman"  # set legend font

for i in range(2):
    ax[i][0].set_yticks(x_pos)
    ax[i][0].set_yticklabels(label_text, fontname='Times New Roman', fontsize=10)
    ax[i][0].grid(axis='x', alpha=0.5, linestyle='--')

for k in range(2):
    ax[k][1].set_yticks(x_pos)
    ax[k][1].set_yticklabels(label_text, rotation=0, fontname='Times New Roman', fontsize=10)
    ax[k][1].grid(axis='x', alpha=0.5, linestyle='--')

ax[0][0].set_xticklabels(np.arange(0, 0.030, 0.005), fontname='Times New Roman', fontsize=10)
ax[0][1].set_xticklabels(np.arange(0, 0.10, 0.02), fontname='Times New Roman', fontsize=10)
ax[1][0].set_xticks(np.arange(0, 5, 1))
ax[1][0].set_xticklabels(np.arange(0, 5, 1), fontname='Times New Roman', fontsize=10)
ax[1][1].set_xticklabels(np.arange(0, 60, 10), fontname='Times New Roman', fontsize=10)
ax[0][0].set_title('a)', loc='left', fontname='Times New Roman', fontsize=10)
ax[0][1].set_title('b)', loc='left', fontname='Times New Roman', fontsize=10)
ax[1][0].set_title('c)', loc='left', fontname='Times New Roman', fontsize=10)
ax[1][1].set_title('d)', loc='left', fontname='Times New Roman', fontsize=10)
ax[0][0].set_xlabel('QALYs lost per HSV-1 infection', **{'fontname': 'Times New Roman'}, fontsize=10)
ax[0][1].set_xlabel('QALYs lost per HSV-2 infection', **{'fontname': 'Times New Roman'}, fontsize=10)
ax[1][0].set_xlabel('total QALYs lost due to HSV-1 infections (x$10^3$)', **{'fontname': 'Times New Roman'},
                    fontsize=10)
ax[1][1].set_xlabel('total QALYs lost due to HSV-2 infections (x$10^3$)', **{'fontname': 'Times New Roman'},
                    fontsize=10)

plt.suptitle('One-way sensitivity analysis for key parameters',
             fontsize=10, fontname='Times New Roman', fontweight='bold')

# save file
plt.tight_layout(rect=(0, 0, 1, 0.95))
plt.savefig('{}one_way_sensitivity_analysis_20221010.png'.format(output_path), dpi=350, bbox_inches='tight')
plt.show()
