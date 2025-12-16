#!/usr/bin/env python3

'''
AUTHOR: Vincent Chapdelaine
VERSION: 1.1
YEAR: 2024

Goal:
Generate lollipop plots with annotations (e.g. protein features) and represent mutation data with visual encodings
including log fold change (y), position (x), and statistical significance (marker size, transparency).
'''

# Standard libraries
import random
import sys
import re
import warnings
import textwrap
import math
import argparse

# Scientific/data libraries
import pandas as pd
import numpy as np
from scipy.stats import rankdata, binomtest

# Plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

# ---------------------
# Argument Parser Setup
# ---------------------
parser = argparse.ArgumentParser(description='Create lollipop plots from MaGeCK and BED-like input')
parser.add_argument('-b', '--bed', required=True, dest='bed_file', help='BED file with features (start, end, name, protein)')
parser.add_argument('-i', '--input', required=True, dest='excel_file', help='Excel file from MaGeCK/VEP output')
parser.add_argument('--stat_method', default='quantile',dest='stat_method', choices={'binom_sign', 'quantile', 'sign_test'})
parser.add_argument('--Bio_threshold', dest='Biological_threshold',type=float, required=True, help='Biological threshold (two-sided)')
parser.add_argument('--scheme_location', dest='scheme_loc', default='top', choices={'top', 'bottom', 'middle'})
parser.add_argument('--histogram', dest='histogram', action='store_true', help='Add histogram')
parser.add_argument('--violin', dest='violin', action='store_true', help='Add violin plot')
parser.add_argument('--violin_detail', dest='violin_detail', default='low', choices={'low', 'high'})
parser.add_argument('-p', '--Prob_Threshold', dest='p_thresh', type=float, default=1, help='P-value threshold')
parser.add_argument('--no_stem', dest='no_stem', action='store_true', help='Remove stemlines')
parser.add_argument('-F', dest='fdr', action='store_true', help='Use FDR instead of p-value')
parser.add_argument('--highlight_region', dest='highlight', nargs='+', help='Highlight regions: protein-feature format')

# Set global font size
plt.rcParams.update({'font.size': 38})

# ------------------------
# Color and Marker Mapping
# ------------------------
consequence_mapping = {
    'synonymous': ('g', 'D'),       # green, diamond
    'missense': ('purple', 'o'),   # purple, circle
    'non-sense': ('red', 's'),     # red, square
    'splice': ('gold', '^')        # gold, triangle
}

consequence_mapping_2 = {
    'synonymous': 'g',
    'missense': 'purple',
    'non-sense': 'red',
    'splice': 'gold',
    'Negative Control Gene': 'lightgray',
    'No predicted Mutation': 'darkgray'
}


# ------------------------
# Tick Generation Function
# ------------------------
def generate_ticks_with_labels(min_val, max_val, n_steps):
    """
    Create human-readable tick positions and labels using magnitude rounding.
    Ensures the tick spacing aligns with the data scale.
    """
    base_steps = [1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500]
    raw_step = abs(max_val - min_val) / max(n_steps, 1)

    if raw_step == 0:
        return [min_val], [f"{min_val:.2f}"]

    magnitude = 10 ** math.floor(math.log10(raw_step))
    candidates = [b * magnitude for b in base_steps]
    step = next((s for s in candidates if s >= raw_step), candidates[-1] * 10)

    start_tick = math.floor(min_val / step) * step
    end_tick = math.ceil(max_val / step) * step

    ticks = []
    current = start_tick
    while current <= end_tick + 1e-10:
        ticks.append(round(current, 10))
        current += step

    labels = [f"{t}" for t in ticks]
    return ticks, labels

# ------------------------
# Color Utility Functions
# ------------------------
def is_valid_color(color):
    """Check if a color string is recognized by matplotlib."""
    try:
        mcolors.to_rgba(color)
        return True
    except ValueError:
        return False

def get_color_palette(num_entries):
    """Choose a colormap based on the number of needed unique colors."""
    pastel1 = plt.get_cmap('Pastel1').colors
    tab20 = plt.get_cmap('tab20').colors
    if num_entries <= len(pastel1):
        return pastel1
    elif num_entries <= len(tab20):
        return tab20
    else:
        warnings.warn("Too many entries for color palette; reusing tab20")
        return tab20

def get_color(index, palette, existing_colors):
    """Return a unique color from the palette, avoiding duplicates."""
    color = mcolors.to_hex(palette[index % len(palette)])
    while color in existing_colors:
        index += 1
        color = mcolors.to_hex(palette[index % len(palette)])
    existing_colors.add(color)
    return color

def create_color_dict(df):
    """
    Generate a dictionary of {feature name: color}.
    Valid colors are used as-is; invalid or missing get assigned from a palette.
    Introns are ignored in color assignment.
    """
    df = df[~df['name'].isna()].copy()
    if 'color' in df.columns:
        existing_colors = {c for c in df['color'] if pd.notnull(c) and is_valid_color(c)}
        filtered_df = df[~df['name'].str.contains('intron') & (~df['color'].apply(is_valid_color) | df['color'].isnull())]
    else:
        existing_colors = set()
        filtered_df = df[~df['name'].str.contains('intron')]

    num_entries = len(filtered_df)
    color_palette = get_color_palette(num_entries)
    color_dict = {}

    index = 0
    for _, row in df.iterrows():
        if pd.isna(row['name']) or 'intron' in str(row['name']):
            continue
        if 'color' in df.columns and pd.notnull(row['color']) and is_valid_color(row['color']):
            color_dict[row['name']] = row['color']
        else:
            color_dict[row['name']] = get_color(index, color_palette, existing_colors)
            index += 1
    return color_dict


def add_text_or_legend(ax, start, length,size, text, color, legend_dict):
    """
    Add legend hande to dictionnary for plotting
    
    Parameters:
    ax (matplotlib.axes._subplots.AxesSubplot): Axis to add the legend to.
    start (int) : start position of feature
    length (int) : length of the feature
    size (float) : Width of the rendered feature
    text (str) : Feature name
    color (str) : Color of the feature
    legend_dict (dict) : dictionnary to be updated
    
    Returns:
    None
    """
    renderer = ax.figure.canvas.get_renderer()

    # Create the temporary text object for width measurement
    temp_text = ax.text(0, 0, text,fontsize=8)
    
    # Get the width of the text
    text_width = temp_text.get_window_extent(renderer=renderer).width
    
    # Remove the temporary text object
    temp_text.remove()
    
    # Check if there is enough space to display the text

    if text not in legend_dict:
        legend_dict[text] = {'color': color}



def plot_genomic_regions(df, ax, legend_loc='upper left', title='', legend_title='Legend', color_dict=None, Maximum=None, minumum=None, Custom_Xaxis=False, xlabel='Amino acid position'):
    """
    Plot genomic regions on the given axis.
    
    Parameters:
    bed (pd.DataFrame): BED DataFrame.
    ax (matplotlib.axes._subplots.AxesSubplot): Axis to plot on.
    b (int): Baseline adjustment.
    legend_loc (str): Legend location.
    num_ticks (int): Number of ticks on the X-axis.
    xlabel (str): Label for the X-axis.
    title (str): Title of the plot.
    
    Returns:
    None
    """    
    if  color_dict is None :
        color_dict=create_color_dict(df)
    # Dictionary to store legend entries
    legend_dict = {}
    # Plot each region
    previous_end = None
    #ax.plot([0, 0], [0.75, 0.25], color='grey', linestyle='-')
    #ax.plot([Maximum, Maximum], [0.75, 0.25], color='grey', linestyle='-')
    for index, row in df.iterrows():
        start, end = row['start'], row['end']
        length = end - start
        if previous_end is None and start > 0:
            ax.plot([0, start], [0.5, 0.5], color='grey', linestyle='-')
        # Plot regions not described in the BED file as straight lines
        if previous_end is not None and start > previous_end:
            ax.plot([previous_end, start], [0.5, 0.5], color='grey', linestyle='-')

        # Plot described regions
        if not pd.isna(row['name']):
            color = color_dict[row['name']]
            renderer = ax.figure.canvas.get_renderer()
            #rect = plt.Rectangle((start, 0), end - start, 1, color=color, alpha=0.5)
            rect = patches.FancyBboxPatch((start, 0), end - start, 1, color=color, alpha=0.5,boxstyle="round4")
            box=ax.add_patch(rect)
            rect_extent = rect.get_window_extent(renderer=renderer)
            rect_width = rect_extent.width
            # Add text inside the box or to the legend
            add_text_or_legend(ax, start, end - start,rect_width , row['name'], color, legend_dict)
        else:
            ax.plot([start, end], [0.5, 0.5], color='grey', linestyle='-')
        previous_end = end
    if previous_end < Maximum :
        ax.plot([previous_end, Maximum], [0.5, 0.5], color='grey', linestyle='-')
    # Customize plot
    ax.set_ylim([0, 1])
    ax.grid(False)
    ax.yaxis.set_visible(False)
    if Custom_Xaxis :
        ax.set_xticks(Custom_Xaxis[0])
        ax.set_xticklabels(Custom_Xaxis[1])
        ax.set_xlabel(xlabel)
    else :
        ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_frame_on(False)
    if legend_dict:
        legend_handles = [Patch(color=info['color'], label=f"{name}") for name, info in legend_dict.items()]
        return legend_handles


def create_lollipop_plot(ax, x, y, color='b', marker='o', line_style='-', line_width=2, alpha=1.0, size=6, Custom_Xaxis=False,stemline_remove=False, fdr=False, xlabel='Amino acid position', lines=None, yaxis=None):
    """
    Create a lollipop plot on the specified axes.
    
    Parameters:
    ax (matplotlib.axes.Axes): The axes to plot on.
    x (list or array): The x values of the points.
    y (list or array): The y values of the points.
    color (str or list): Color of the markers. Can be a single color or a list of colors. Default is 'blue'.
    marker (str or list): Marker style. Can be a single marker style or a list of marker styles. Default is 'o'.
    line_style (str): Line style. Default is '-' (solid line).
    line_width (float): Line width. Default is 2.
    alpha (float or list): Transparency of the markers. Can be a single value or a list of values. Default is 1.0.
    size (int or list): Size of the markers. Can be a single value or a list of values. Default is 6.
    Custom_Xaxis ([[int],[str]]) : list of xticks and x labels 
    
    Returns:
    None
    """
    if isinstance(color, str):
        color = [color] * len(x)
    if isinstance(marker, str):
        marker = [marker] * len(x)
    if isinstance(alpha, (int, float)):
        alpha = [alpha] * len(x)
    if isinstance(size, (int, float)):
        size = [size] * len(x)
    (markers, stemlines, baseline) = ax.stem(x.to_numpy(), y.to_numpy(), linefmt='gray', markerfmt=" ", basefmt=" ")
    if stemline_remove :
        stemlines.remove()
        baseline.remove()
    else :
        plt.setp(stemlines, alpha=0.15)
        plt.setp(baseline, alpha=0.15)
    for xi, yi, ci, mi, ai, si in zip(x, y, color, marker, alpha, size):
        ax.plot([xi], [yi], marker=mi, color=ci, linestyle='None', alpha=ai, markersize=si)   
    # Customize the plot appearance
    ax.set_ylabel('Log Fold Change')
    if lines:
        for line in lines :
            if line:
                ax.axhline(y=line, color='red', linestyle='--', linewidth=1)
    # Remove grid
    if yaxis :
        ax.set_ylim(ymin=yaxis[0], ymax=yaxis[1])
    ax.grid(False)
    if Custom_Xaxis :
        ax.set_xticks(Custom_Xaxis[0])
        ax.set_xticklabels(Custom_Xaxis[1])
        ax.set_xlabel(xlabel)
    else :
        ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_frame_on(False)

def add_legend(ax, consequence_mapping, pvalue_mapping,transparency=None, add_legend=False, fdr=False):
    """
    Add a legend to the plot for the given consequence mapping and P-value sizes.
    
    Parameters:
    ax (matplotlib.axes.Axes): The axes to add the legend to.
    consequence_mapping (dict): A dictionary mapping consequence types to their colors and markers.
    pvalue_mapping (dict): A dictionary mapping P-value levels to their sizes.
    
    Returns:
    None
    """
    # Create custom legend handles
    consequence_legend_elements = [
        mlines.Line2D([], [], color=value[0], marker=value[1], linestyle='None', markersize=30, label=key)
        for key, value in consequence_mapping.items()
    ]
    
    size_legend_elements = [
        mlines.Line2D([], [], color='gray', marker='o', linestyle='None', markersize=size, label=f'{pvalue}')
        for pvalue, size in pvalue_mapping.items()
    ]
    if transparency:
        transparency_legend_elements = [
            mlines.Line2D([], [], color='purple', marker='o', linestyle='None', markersize=30, alpha=alpha, label=label)
            for alpha, label in zip([0.4, 1.0], transparency)
        ]
    # Add subtitles
    subtitle_fontsize = 'medium'
    leg_a=ax.legend(handles=consequence_legend_elements,loc='upper left', title='Consequences', handlelength=1,bbox_to_anchor=(1, 1), fontsize=subtitle_fontsize, labelspacing=1.25,frameon=False) #,bbox_transform=fig.transFigure
    if transparency :
        leg_c=ax.legend(handles=transparency_legend_elements,loc='center left', title='Biologically Significance', handlelength=1,bbox_to_anchor=(1, 0.55), fontsize=subtitle_fontsize, labelspacing=1.25, frameon=False) #,bbox_transform=fig.transFigure
    leg_b=ax.legend(handles=size_legend_elements, loc='lower left',title='FDR' if fdr else 'P-value', handlelength=1,bbox_to_anchor=(1, 0.35), fontsize=subtitle_fontsize, labelspacing=1.25,frameon=False) #,bbox_transform=fig.transFigure
    ax.add_artist(leg_a)
    ax.add_artist(leg_b)
    if transparency :
        ax.add_artist(leg_c)
    if add_legend :
        leg_d=ax.legend(handles=add_legend, handlelength=1, loc='upper left',bbox_to_anchor=(1, 0.25), labelspacing=1.25,frameon=False,title='Domains')
        ax.add_artist(leg_b)


def highlight_region(df, ax, name,color_dict=None, negative=False,):
    """
    Add a rectangle to highlight a region in the BED DataFrame by name.
    
    Parameters:
    df (pd.DataFrame): BED DataFrame with at least 'name', 'start', and 'end' columns.
    ax (matplotlib.axes.Axes): Axes object to draw the rectangle on.
    name (str): Name of the region to highlight.
    """
    # Check if the name exists in the DataFrame
    if name not in df['name'].values:
        raise ValueError(f"Name '{name}' not found in the DataFrame")
    if color_dict ==None :
        color_dict=create_color_dict(df)
    # Get the start and end positions
    region = df[df['name'] == name].iloc[0]
    start = region['start']
    end = region['end']
    bottom, top = ax.get_ylim()
    origin= min(bottom, top) if negative else 0
    length= abs(origin) if negative else top
    # Add a rectangle covering the whole y-axis
    rect = plt.Rectangle((start, origin),
                             end - start, 
                             length,
                         color=color_dict[name], alpha=0.3)
    ax.add_patch(rect)

def grouped_violin_plot(ax_input, group1, negative=None, positive=None, ymin=None, ymax=None, hue=None):
    """
    Each group: (list of categories, list of values)
    group1 is required, group2 and group3 are optional.
    """
    figure , gs_outer , sharey = ax_input
    all_dfs = []
    x_tick_labels = []
    x_positions = []
    group_ranges = []
    group_labels = []
    def add_group(group, label, start_pos):
        categories, values = group
        df = pd.DataFrame({'Category': categories, 'Value': values})
        unique_cats = list(pd.unique(categories))
        df['X'] = [start_pos + unique_cats.index(c) for c in categories]
        x_tick_labels.extend(unique_cats)
        x_positions.extend(range(start_pos, start_pos + len(unique_cats)))
        group_ranges.append((start_pos, start_pos + len(unique_cats) - 1))
        group_labels.append(label)
        return df, start_pos + len(unique_cats)

    pos = 0
    df1, pos = add_group(group1, '', pos)
    all_dfs.append(df1)

    if negative and not all(s.empty for s in negative):
        df2, pos = add_group(negative, 'Negative Controls', pos)
        Negative_indexes=[list(hue.keys()).index(i) for i in pd.unique(df2['Category'])]
        all_dfs.append(df2)

    if positive and not all(s.empty for s in positive):
        df3, pos = add_group(positive, 'Positive Controls', pos)
        Positive_indexes=[list(hue.keys()).index(i) for i in pd.unique(df3['Category'])]
        all_dfs.append(df3)

    df_all = pd.concat(all_dfs)

    gs_inner = gs_outer.subgridspec(1, len(hue.keys()), wspace=0.1)

    for ind, item in enumerate(list(hue.keys())):
        current = figure.add_subplot(gs_inner[ind], sharey=sharey)
        df=df_all.loc[df_all['Category']==item, :]
        #sns.kdeplot(ax=current, data=df, y="Value", color=hue[item], fill=True, alpha=.5)
        sns.boxplot(ax=current, data=df, y="Value", color=hue[item], width=0.6, saturation=0.9)
        if ymin or ymax :
            current.set_ylim(ymin=ymin, ymax=ymax)
        current.set_xlabel('')
        current.xaxis.set_visible(False)
        current.set_frame_on(False)
        current.yaxis.set_visible(False)
        for spine in current.spines.values():
            spine.set_visible(False)
        if negative and not all(s.empty for s in negative):
            if ind == min(Negative_indexes) :
                Negative_start=current.get_position()
            if ind == max(Negative_indexes) :
                Negative_end=current.get_position()
        if positive and not all(s.empty for s in positive):
            if ind == min(Positive_indexes) :
                Positive_start=current.get_position()
            if ind == max(Positive_indexes) :
                Positive_end=current.get_position()
    if positive and not all(s.empty for s in positive):
        y_line = Positive_start.y0 + 0.05
        line = mlines.Line2D(
            [Positive_start.x0, Positive_end.x1],   # x0→x1
            [y_line,        y_line],      # constant y
            transform=figure.transFigure, # interpret coords in figure fractions
            color="black",
            linewidth=1.5,
            clip_on=False
        )
        figure.add_artist(line)
        # 5. add text centered under that line
        figure.text(
            (Positive_start.x0 + Positive_end.x1)/2,  # midpoint x
            y_line + 0.05,                  # slightly below the line
            "Positive Controls",
            ha="center", va="top",
            transform=figure.transFigure
        )
    if negative and not all(s.empty for s in negative):
        y_line = Negative_start.y0 - 0.05
        line = mlines.Line2D(
            [Negative_start.x0, Negative_end.x1],   # x0→x1
            [y_line,        y_line],      # constant y
            transform=figure.transFigure, # interpret coords in figure fractions
            color="black",
            linewidth=1.5,
            clip_on=False
        )
        figure.add_artist(line)
        # 5. add text centered under that line
        figure.text(
            (Negative_start.x0 + Negative_end.x1)/2,  # midpoint x
            y_line - 0.02,                  # slightly below the line
            "Negative Controls",
            ha="center", va="top",
            transform=figure.transFigure
        )


# Plot each region
if __name__ == '__main__':
    args = parser.parse_args()
    #### Input Parsing
    bed_full = pd.read_csv(args.bed_file, sep='\t', header=0)
    if not set(['start','end','name','proteins']).issubset(bed_full.columns):
        raise ValueError('Bed-like file does not contain the relevant columns')
    protein_list=set(bed_full.proteins)
    try :
        all_sheets = pd.read_excel(args.excel_file, sheet_name=None, dtype=str)
    except:
        raise ValueError('Main input data (-i) does not have a supported format')

    for sheet_name, data_full in all_sheets.items():
        proteins=set(data_full.proteins)
        data_full['lfc']=pd.to_numeric(data_full['lfc'])
        data_full['Position']=pd.to_numeric(data_full['Position'])
        data_full['p_value']=pd.to_numeric(data_full['p_value'])
        ### trouble shoot user error proteins bed and protein input file
        if not protein_list.issubset(proteins): 
            raise ValueError(f'Not all proteins/regions in the main input file (-i) is present in the Bed-like file:{set(proteins)-set(bed_full.proteins)} \n')
        ### Test highlight list
        if args.highlight:
            if not set(args.highlight).issubset(set(bed_full['proteins'] + '-' + bed_full['name'])):
                raise ValueError(f'Not all proteins/regions to be highlighted are in the bed-like file under the protein-feature format :{set(args.highlight)-set(bed_full["proteins"] + "-" + bed_full["name"])} \n')
        
        for protein in protein_list:
            ### Produce a graph per protein for all protein
            ### Produce 
            data=data_full.loc[data_full['proteins']==protein,:]
            data=data.drop(data[data["Position"].isna()].index)
            data=data.loc[data['Consequence'].isin(consequence_mapping.keys())]
            data_negative_control = data_full.loc[data_full['Controls']=='negative_control',:]
            if args.Biological_threshold :
                line_min=np.quantile(pd.to_numeric(data_negative_control['lfc']),float(args.Biological_threshold) )
                line_max=np.quantile(pd.to_numeric(data_negative_control['lfc']),1-float(args.Biological_threshold) )
                data['Bio_Sig'] = [True if i >= line_max or i<= line_min else False for i in data['lfc']]
                print(f"Biological Significance of {sheet_name}: 5% ({np.quantile(pd.to_numeric(data_negative_control['lfc']),0.05)},{np.quantile(pd.to_numeric(data_negative_control['lfc']),0.95 )}) , 1% : ({np.quantile(pd.to_numeric(data_negative_control['lfc']),0.01)},{np.quantile(pd.to_numeric(data_negative_control['lfc']),0.99 )})")
            else :
                line_min=None
                line_max=None
            data_positive_control = data_full.loc[data_full['Controls']=='positive_control',:]
            if args.violin:
                ymax = np.round(max(pd.to_numeric(pd.concat([data['lfc'], data_negative_control['lfc'],data_positive_control['lfc'] ])))+1)
                ymin = np.round(min(pd.to_numeric(pd.concat([data['lfc'], data_negative_control['lfc'],data_positive_control['lfc'] ])))-1)
            else :
                ymax = np.round(max(pd.to_numeric(data['lfc']))+1)
                ymin = np.round(min(pd.to_numeric(data['lfc']))-1)               
            bed_adjusted=bed_full.loc[bed_full['proteins']==protein,:]
            bed=bed_adjusted.copy()
            maximum=max(pd.concat([bed_adjusted['end'],data['Position']])) # max in graph
            minimum=min(pd.concat([bed_adjusted['start'],data['Position']])) # min in graph
            # Figure Configuration
            if args.scheme_loc == 'top' :
                ratios= [1, 20]
                nfigure = 2
            elif args.scheme_loc == 'middle' :
                ratios= [10, 1, 10]
                nfigure = 3
            else :
                ratios= [20, 1]
                nfigure = 2       
            fig = plt.figure(figsize=(40 if args.violin else 25, 25 if args.histogram else 24,))
            if args.violin:
                gs = gridspec.GridSpec(nfigure + args.histogram, args.violin + 1, hspace=0.05, wspace=0, width_ratios=[1,0.5], height_ratios=[2] + ratios, left=0.05, right=1, top=0.95, bottom=0.05)
            else :
                gs = gridspec.GridSpec(nfigure + args.histogram, args.violin + 1, hspace=0.05, wspace=0, height_ratios=[2] + ratios, left=0.05, right=1, top=0.95, bottom=0.05)
            if args.violin:
                if args.scheme_loc == 'middle' :
                    ax_low = fig.add_subplot(gs[0 + args.histogram,0])
                    ax_scheme = fig.add_subplot(gs[1 + args.histogram, 0])
                    ax_high = fig.add_subplot(gs[2 + args.histogram, 0])
                    ax_violin_input = (fig, gs[:, 1], None)
                else :
                    ax_low = fig.add_subplot(gs[(args.scheme_loc != 'bottom') + args.histogram,0])
                    ax_scheme = fig.add_subplot(gs[(args.scheme_loc == 'bottom') + args.histogram, 0], sharex=ax_low)
                    ax_violin_input = (fig, gs[(args.scheme_loc != 'bottom') + args.histogram, 1],ax_low)
                if args.histogram :
                    ax_histofram = fig.add_subplot(gs[0,0], sharex=ax_low)
            else :
                if args.scheme_loc == 'middle' :
                    ax_low = fig.add_subplot(gs[0 + args.histogram])
                    ax_scheme = fig.add_subplot(gs[1 + args.histogram])
                    ax_high = fig.add_subplot(gs[2 + args.histogram])
                else : 
                    ax_low = fig.add_subplot(gs[(args.scheme_loc != 'bottom') + args.histogram])
                    ax_scheme = fig.add_subplot(gs[(args.scheme_loc == 'bottom') + args.histogram])
                if args.histogram :
                    ax_histofram = fig.add_subplot(gs[0])
            ### Histogram plot and figure grid config (according to histogram yes/no)
            if args.histogram :
                ax_histofram.hist(data['Position'], bins=300, color=plt.cm.Paired(0), edgecolor='none')
                ax_histofram.xaxis.set_visible(False)
                for spine in ax_histofram.spines.values():
                    spine.set_visible(False)

            if args.fdr :
                data['p_value']=data['fdr']

            ### parameter data representation
            colors = data['Consequence'].map(lambda x: consequence_mapping[x][0])
            markers = data['Consequence'].map(lambda x: consequence_mapping[x][1])
            alphas = data['Bio_Sig'].map(lambda x: 0.4 if not x else 1.0)
            if args.stat_method=='quantile':
                Biosig_labels=[fr"$Q_{{{args.Biological_threshold:.2f}}} \leq x \leq Q_{{{1 - args.Biological_threshold:.2f}}}$", fr"$x < Q_{{{args.Biological_threshold:.2f}}}$ or $x > Q_{{{1 - args.Biological_threshold:.2f}}}$"]
            else :
                Biosig_labels = [fr"$p$-value $> {args.Biological_threshold:.2f}$", fr"$p$-value $\leq {args.Biological_threshold:.2f}$"]
            sizes = pd.Series([28 if i < args.p_thresh else 20 for i in data['p_value']], index=data.index)
            pvalue_mapping={fr"p-value $\leq {args.p_thresh:.2f}$" : 28, fr"p-value > ${args.p_thresh:.2f}$" : 20}
            ticks, labels = generate_ticks_with_labels(minimum, maximum, 10)
            if args.scheme_loc == 'middle' :
                create_lollipop_plot(ax_low, data.loc[data['lfc']>=0,'Position'], data.loc[data['lfc']>=0,'lfc'], color=colors[data['lfc']>=0], marker=markers[data['lfc']>=0], size=sizes[data['lfc']>=0], alpha=alphas ,stemline_remove=args.no_stem,lines=line_max,yaxis=(0,ymax))
                if args.highlight:
                    for hl in args.highlight : 
                        protein_hl, feature=hl.split('-')
                        if protein_hl == protein :
                            highlight_region(bed_adjusted,ax_low,feature)
                create_lollipop_plot(ax_high, data.loc[data['lfc']<0,'Position'], data.loc[data['lfc']<0,'lfc'], color=colors[data['lfc']<0], marker=markers[data['lfc']<0], size=sizes[data['lfc']<0],Custom_Xaxis=[ticks, labels], alpha=alphas, stemline_remove=args.no_stem, xlabel = 'Amino acid position', lines=line_min,yaxis=(ymin,0))
                if args.highlight:
                    for hl in args.highlight : 
                        protein_hl, feature=hl.split('-')
                        if protein_hl == protein :
                            highlight_region(bed_adjusted,ax_high,feature ,None,True)
                leg=plot_genomic_regions(bed_adjusted,ax_scheme ,legend_loc='upper left', title='',legend_title='Features',Maximum=maximum)
                add_legend(fig, consequence_mapping, pvalue_mapping,transparency=Biosig_labels ,add_legend=leg,fdr=args.fdr)
            else :
                create_lollipop_plot(ax_low, data['Position'], data['lfc'], color=colors, marker=markers, size=sizes,stemline_remove=args.no_stem, alpha=alphas,Custom_Xaxis=[ticks, labels] if args.scheme_loc == 'top' else None, xlabel = 'Amino acid position',yaxis=(ymin,ymax))
                leg=plot_genomic_regions(bed_adjusted,ax_scheme, legend_loc='upper left', title='',legend_title='Features',Maximum=maximum,Custom_Xaxis=[ticks, labels] if args.scheme_loc == 'bottom' else None ,xlabel = 'Amino acid position')
                add_legend(fig, consequence_mapping, pvalue_mapping,transparency=Biosig_labels,add_legend=leg,fdr=args.fdr)
                if args.highlight:
                    for hl in args.highlight : 
                        protein_hl, feature=hl.split('-')
                        if protein_hl == protein :
                            highlight_region(bed_adjusted,ax_low,feature)
            if args.violin :
                if args.violin_detail == 'high' :
                    grouped_violin_plot(ax_violin_input, (data['Consequence_Detail'],data['lfc']),negative=(data_negative_control['Consequence_Detail'],data_negative_control['lfc']),positive=(data_positive_control['Consequence_Detail'],data_positive_control['lfc']), ymin=ymin,ymax=ymax, hue= consequence_mapping_2)
                else : 
                    grouped_violin_plot(ax_violin_input, (data['Consequence'],data['lfc']),negative=(data_negative_control['Consequence'],data_negative_control['lfc']),positive=(data_positive_control['Consequence'],data_positive_control['lfc']),ymin=ymin,ymax=ymax,hue= consequence_mapping_2)
            fig.savefig(f'BEscreen_{sheet_name}_{protein}.pdf',format="pdf",bbox_inches="tight")
