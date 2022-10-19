from supports.VisualizationSupport import *

read_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/trees/dics/'
output_path = '/Users/shiyingyou/PycharmProjects/HSV-QALY/New_HSV_tree/visualization/graphs/'

# read python dict back from the file
pkl_file = open('{}DicHSV1_dr3_psycho1_noGBD_20221010_dN100.pkl'.format(read_path), 'rb')
hsv1_Dic = pickle.load(pkl_file)
pkl2_file = open('{}DicHSV2_dr3_psycho1_noGBD_20221010_dN100.pkl'.format(read_path), 'rb')
hsv2_Dic = pickle.load(pkl2_file)
pkl_file.close()
pkl2_file.close()

mean_list1, ci_list1 = get_total_qaly_loss(dic=hsv1_Dic, scalar=1000)
mean_list2, ci_list2 = get_total_qaly_loss(dic=hsv2_Dic, scalar=1000)

# draw plots
width = 0.4
bar_alpha = 0.8
ci_alpha = 0.8
hsv2_c = '#2596be'
hsv2_female_c = lighten_color(hsv2_c, 0.8)
hsv2_male_c = lighten_color(hsv2_c, 0.64)
hsv1_c = '#5D62D5'
hsv1_female_c = lighten_color(hsv1_c, 0.8)
hsv1_male_c = lighten_color(hsv1_c, 0.64)
x_pos = np.arange(3)
x_text = ['Men', 'Women', 'Total']
# figure size
plt.figure(figsize=(5.5, 5.5))
# bar plot
plt.bar(x_pos - 0.2, mean_list1, width=width, color=hsv1_c, alpha=bar_alpha, label='HSV-1',
        yerr=[ci/4 for ci in ci_list1], ecolor=lighten_color(hsv1_c, 1.2), capsize=3)
plt.bar(x_pos + 0.2, mean_list2, width=width, color=hsv2_c, alpha=bar_alpha, label='HSV-2',
        yerr=[ci/4 for ci in ci_list2], ecolor=lighten_color(hsv2_c, 1.2), capsize=3)

##########
# FORMAT #
##########
plt.rcParams["font.family"] = "Times New Roman"  # font for legend
# x ticks text, y limit
plt.xticks(x_pos, labels=x_text, rotation=0)
plt.xticks(fontname='Times New Roman', fontsize=10)  # font for x ticks
plt.yticks(fontname='Times New Roman', fontsize=10)  # font for y ticks
# grid
plt.grid(axis='y', alpha=0.5, linestyle='--')
# labels and title
plt.legend(loc='upper left')
plt.title('Total QALYs lost by sex and HSV type in the US (2018)',
          fontsize=10, fontname='Times New Roman', fontweight='bold')
plt.xlabel('Groups', **{'fontname': 'Times New Roman'}, fontsize=10)  # font for x label
plt.ylabel('Total QALYs lost (thousand)', **{'fontname': 'Times New Roman'}, fontsize=10)  # font for y label

# save file
plt.tight_layout()
plt.savefig('{}total_loss_20221010.png'.format(output_path), dpi=350, bbox_inches='tight')
plt.show()

