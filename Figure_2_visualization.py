# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:42:36 2017

@author: hhaque
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import PtitPrince
import seaborn as sns
import os
from scipy import stats
from matplotlib.ticker import FixedLocator
import matplotlib.lines as mlines
from decimal import Decimal



def get_CI_limits(condition, sign):
    fpath = 'K:\\palva\\VWM2\\_Documentation\\2019_recordings\\24 - Data sharing\\plot_data\\Figure_2\\Wilcoxon\\jackknife\\'
    condition = condition.replace('0','')
    if condition[-3:] == '4-2':
        condition = '{}4_{}2'.format(condition[:-3], condition[:-3])
    fname = 'cPLV_jackknife_{}_bothTWs_{}.xlsx'.format(condition, sign)
    fpath = fpath + fname

    df = pd.read_excel(fpath)
    
    jackknife_col_names = ['Jackknife_{}'.format(i) for i in range(20)]
    df_jackknife = df.set_index('Frequency')[jackknife_col_names]
    df_jackknife['CI_max'] = df_jackknife.max(axis=1)
    df_jackknife['CI_min'] = df_jackknife.min(axis=1)
    
    return df_jackknife


master_folder = 'path_to_extracted_data' # set to directory of extracted data



'''Plot Fig 2C'''

os.chdir('{}\\plot_data\\Figure_2\\Wilcoxon'.format(master_folder))

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fig = plt.figure(figsize = (7.3, 3.7))

features = ['Shape', 'Color', 'Spatial']
colors = sns.color_palette('muted')

loads = ['02', '04', '04-02']

threshold = 0.00672

count = 1
for i in range(len(features)):
    for load in loads:
        condition = features[i] + load

        pos = 'wilcoxon_%s_pos.xlsx' % condition
        neg = 'wilcoxon_%s_neg.xlsx' % condition

        if load == '02':
            ax = fig.add_subplot(2, 3, count)
            
            ax.tick_params(which='both', top=False, right=False)
            ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
            ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
            ax.set_xscale('log')
            ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
            ax.set_xlim([3, 120])
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.set_ylim([-0.07, 0.09])
            ax.yaxis.set_major_locator(FixedLocator([-0.04, 0.00, 0.04, 0.08]))
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_linewidth(0.6)
            ax.spines['left'].set_linewidth(0.6)
            
            the_title = features[i]
            if features[i] == 'Spatial':
                the_title = 'Location'
    
            plt.title(the_title, size = 8)
            fig.subplots_adjust(hspace= .5)
            
            plt.ylabel('K', size = 7, labelpad = 0.5)
            plt.xlabel('Frequency (Hz)', size = 7)
          
        for excel_file in ['pos', 'neg']:    
            df = pd.read_excel(eval(excel_file))
            df.columns = ['Frequency','Property']
           
            if load == '02':    
                ax.plot(df.Frequency, df.Property, lw = 1, linestyle = '--', color = colors[i], label = '_nolegend_')
            if load == '04':    
                ax.plot(df.Frequency, df.Property, lw = 1, color = colors[i], label = '_nolegend_')
            if load == '04-02':    
                ax.plot(df.Frequency, df.Property, lw = 1, linestyle = (0, (1, 1)), color = colors[i], label = '_nolegend_')
                
            df_jackknife = get_CI_limits(condition, excel_file)
            plt.fill_between(df.Frequency, df_jackknife.CI_min, df_jackknife.CI_max, color = colors[i], alpha = 0.25, linewidth = 0.01)
        
        if load == '04-02':
            count += 1
            
    plt.fill_between(df.Frequency, -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)
    
    x = 0.085
    x_neg = 0.06
    
    plt.plot([11.25, 13.06], [x, x], color = 'black', lw = 0.3) # alpha
    plt.plot([3.10, 4.02], [x, x], color = 'black', lw = 0.3) # delta
    plt.plot([6.00, 8.05], [-x_neg, -x_neg], color = 'black', lw = 0.3) # theta
    plt.plot([17.14, 20.10], [-x_neg, -x_neg], color = 'black', lw = 0.3) # beta
    
    if i == 2:
        dashed_line = mlines.Line2D([], [], color='black', linestyle = '--', lw = 1, label='Load 02')
        solid_line = mlines.Line2D([], [], color='black', lw = 1, label='Load 04')
        dotted_line = mlines.Line2D([], [], color='black', linestyle = (0, (1, 1)), lw = 1, label='Load 04-02')
        plt.legend(handles=[dashed_line, solid_line, dotted_line], frameon = False, prop={'size': 6})
        


'''Plot Fig 2D'''

os.chdir('{}\\plot_data\\Figure_2\\Contrast_betn_features'.format(master_folder))

ax = fig.add_subplot(2, 3, 4)
ax.tick_params(which='both', top=False, right=False)

ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
ax.set_xscale('log')
ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
ax.set_xlim([3, 120])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_ylim([-0.017, 0.07])
ax.yaxis.set_major_locator(FixedLocator([-0.06, -0.04, -0.02, 0.00, 0.02, 0.04, 0.06]))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_linewidth(0.6)
ax.spines['left'].set_linewidth(0.6)

features = ['Shape-Color', 'Shape-Spatial', 'Color-Spatial']
colors = sns.color_palette()[3:10]

for idx, feature in enumerate(features):

    pos = '{}_pos.xlsx'.format(feature)
    neg = '{}_neg.xlsx'.format(feature)

    fig.subplots_adjust(hspace= .5)
    
    plt.ylabel('K', size = 7, labelpad = 0.5)
    plt.xlabel('Frequency (Hz)', size = 7)
  
    for excel_file in ['pos', 'neg']:    
        df = pd.read_excel(eval(excel_file))
        df.columns = ['Frequency','Property']
        
        if excel_file == 'pos':
            ax.plot(df.Frequency, df.Property, lw = 1, color = colors[idx], label = feature)
        else:
            ax.plot(df.Frequency, df.Property, lw = 1, color = colors[idx], label = '_nolegend_')

plt.fill_between(df.Frequency, -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)

plt.legend(frameon = False, prop={'size': 6})



'''Plot Fig 2E'''

fpath = '{}\\plot_data\\Figure_2\\Graph_strength'.format(master_folder)

frequencies = ['Theta', 'Alpha']

count = 5
for idx, frequency in enumerate(frequencies):
    fname = '{}//GS_{}.csv'.format(fpath, frequency)
    
    df_master = pd.read_csv(fname, sep = ';')
            
    ax = fig.add_subplot(2, 3, count)
    PtitPrince.RainCloud(x = 'Condition', y = 'Mean_GS', data = df_master, palette = sns.color_palette('muted'),
                         width_viol = .55, box_showfliers = False, box_linewidth = 1,
                         box_whiskerprops = {'linewidth':1}, point_size = 2)
    ax.tick_params(which='both', top=False, right=False)
    ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.ylabel('Mean GS', size = 7)
    plt.xlabel('Feature', size = 7)
    
    #statistical annotation
    df_shape = df_master.loc[df_master['Condition'] == 'Shape', 'Mean_GS']
    df_color = df_master.loc[df_master['Condition'] == 'Color', 'Mean_GS']
    df_spatial = df_master.loc[df_master['Condition'] == 'Spatial', 'Mean_GS']
    
    #between Shape and Color
    t, prob = stats.ttest_rel(df_shape, df_color)
    if prob < 0.05:
        prob = f"{Decimal(f'{prob:.2g}'):f}"
        x1, x2 = 0, 1
        y, h, col = max(df_shape.max(), df_color.max())+0.001, 0.0005, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.5, c=col)
        plt.text((x1+x2)*.5, y+h, 'p = {}'.format(prob), ha='center', va='bottom', color=col, size = 6)
    
    
    #between Shape and Location
    t, prob = stats.ttest_rel(df_shape, df_spatial)
    if prob < 0.05:
        prob = f"{Decimal(f'{prob:.2g}'):f}"
        x1, x2 = 0, 2
        if frequency == 'Alpha':
            y, h, col = max(df_shape.max(), df_spatial.max())+0.01, 0.0005, 'k'
        if frequency == 'Theta':
            y, h, col = max(df_shape.max(), df_spatial.max())+0.004, 0.0005, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.5, c=col)
        plt.text((x1+x2)*.5, y+h, 'p = {}'.format(prob), ha='center', va='bottom', color=col, size = 6)
    
    
    #between Color and Location
    t, prob = stats.ttest_rel(df_color, df_spatial)
    if prob < 0.05:
        prob = f"{Decimal(f'{prob:.2g}'):f}"
        x1, x2 = 1, 2
        y, h, col = max(df_color.max(), df_spatial.max())+0.001, 0.0005, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.5, c=col)
        plt.text((x1+x2)*.5, y+h, 'p = {}'.format(prob), ha='center', va='bottom', color=col, size = 6)
        
    plt.title(frequency, size = 8)
    
    count += 1

plt.tight_layout()