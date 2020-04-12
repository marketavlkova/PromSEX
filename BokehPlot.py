### Script combining all data for MulitFunPlot and producing the plot(s)

### example to run: python3 BokehPlot.py
### e.g.: python3 BokehPlot.py

import sys
import os.path
import requests
import pandas as pd
from io import StringIO
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot, layout
from bokeh.models import ColumnDataSource, CDSView
from bokeh.models.filters import GroupFilter
from bokeh.models.widgets import Tabs, Panel
from bokeh.palettes import Category10_10

def main():

    if os.path.isfile('output/MultiFunKey.csv'):
        ### load file containing descriptions to BC group IDs
        key = pd.read_csv('output/MultiFunKey.csv', delimiter = ':')
        ### convert the key table into a dictionary
        all_keys = dict([(i, a) for i, a in zip(key.ID, key.function)])
        all_keys['ORFs'] = 'ORFs'
        all_keys['Conserved-Hypothetical-ORFs'] = 'CH-ORFs'
        all_keys['Unclassified-Genes'] = 'Unclassified'

    ### if the file for 'final_table' already exists
    ### load it and skip the part creating it
    if os.path.isfile('output/MultiFunPlotData.csv'):
        final_table = pd.read_csv('output/MultiFunPlotData.csv', delimiter = ',')
        ### define output file the the plot
        output_file("output/BokehPlot.html")
    elif os.path.isfile('output/AllGenes.csv'):
        ### define output file the the plot
        output_file("output/BokehPlot.html")
        ### load data files to check for mismatches in gene names after using EcoCyc
        ecocyc_in = pd.read_csv('output/AllGenes.csv', delimiter = ',')
        ecocyc_out = pd.read_csv('output/AllGenesCompare.txt', delimiter = '\t')

        ### initialize list for possible mismatches in gene names
        mismatch = []
        ### loop through all genes names in the output table from EcoCyc
        for gene_out in ecocyc_out['Gene']:
            ### set the variable that checks whether the same gene name was found in the input table for EcoCyc to 0
            found = 0
            ### loop through all gene names in the input table for EcoCyc
            for gene_in in ecocyc_in['Gene']:
                ### if the same gene name is found, set variable checking this to 1
                if (gene_out == gene_in):
                    found = 1
            ### if the gene name from output file isn't present in the input file
            if (found == 0):
                ### add it to the list of mismatches
                mismatch.append(gene_out)

        ### if any mismatches in gene names were found
        if (len(mismatch) > 0):
            ### print this information to the prompt with the list of the gene names missing in input file
            print('Gene name mismatches between input/output files for/from EcoCyc:\n', mismatch, '\n')
            ### and exit from further executing the script with prompt to fix the mismatches before proceeding
            sys.exit('Fix the mismatches in all output files from EocCyc and run the script again!')

        ### load data table with proportion of segregating sites (Seg) and average pairwise identity (Ide) values
        plot_data = pd.read_csv('output/AllBoth.csv', delimiter = ',')
        ### add column names for further use in indexing
        plot_data = plot_data.rename(columns = {'Unnamed 0':'Promoter', 'V4':'PromIde','V5':'PromSeg', 'V9':'GeneIde', "V10":'GeneSeg'})
        ### add gene names to the table with Ide and Seg
        data = pd.concat([plot_data, ecocyc_in], axis = 1)

        ### load data table with gene assignment to functional groups
        MFinput = pd.read_csv('output/AllGenesMFinput.csv', delimiter = ',')

        ### here I need to create a new table with a separate row for each gene in every BC group
        ### initialize list to store all values for the 'final_table'
        vec = []
        ### loop thorugh all row numbers of table with functional group & genes
        for r in range(0, len(MFinput['ID'])):
            ### split the BC group to extract parent BC groups
            ids = MFinput['ID'][r].split('.')
            ### split genes present in that group, so that each gene is a separable variable
            genes = MFinput['Gene'][r].split(' // ')
            ### loop through all those genes
            for g in genes:
                ### initialize a list to store all values for this gene in
                line = []
                ### remove possibly existing characters: " from the gene name
                if '"' in g:
                    g = g.split('"')[1]
                ### add the BC group ID and gene names to the 'line' list
                line.append(MFinput['ID'][r])
                line.append(g)
                ### initialize a list to store BC subgroups
                bc = []
                ### make a loop that executes 5 times (= number of values to be added to the 'line' list as BC subgroups)
                for id in range(0, 5):
                    ### if going through the loop for the first time
                    if (id == 0):
                        ### add the upper most BC group to the 'line' list
                        line.append(ids[id])
                        ### get the description of the upper most BC group from the key dictionary
                        for ak in all_keys.keys():
                            if ids[id] == ak:
                                line.append(all_keys[ak])
                    ### in all other rounds
                    else:
                        ### if the round number is lower than number of subgroups
                        if (id < len(ids)):
                            ### produce the current subgroup combining subgroup from the previous 'line' value and current extention to it
                            bc = '.'.join([line[2 * id], ids[id]])
                            ### get the description of the current subgroup from the key dictionary
                            for ak in all_keys.keys():
                                if bc == ak:
                                    bckey = all_keys[ak]
                        ### if you run out of subgroup values
                        else:
                            ### set bc variable to NaN
                            bc = 'NaN'
                            bckey = 'None'
                        ### add the BC subgroup value and its description to the 'line' list
                        line.append(bc)
                        line.append(bckey)
                ### set variable to check for multiple Ide & Seg entries for the same gene
                n = 0
                ### loop through all row numbers of the table with Seg & Ide values together with gene names
                for prom in range(0, len(data['Gene'])):
                    ### when you find in the the gene name currently processing and it's the first encounter
                    if all([data['Gene'][prom] == g, n == 0]):
                        ### add its Ide and Seg values to the 'line' list
                        line.append(data['PromIde'][prom])
                        line.append(data['PromSeg'][prom])
                        line.append(data['GeneIde'][prom])
                        line.append(data['GeneSeg'][prom])
                        ### and change the 'n' variable to 1 to avoid processing values for this same gene again if they exist multiple times in the table
                        n += 1
                ### add the whole 'line' variable as a list into the 'vec' list
                vec.append(line)

        ### create dataframe out the the complete 'vec' list with specified column names
        final_table = pd.DataFrame(vec, columns = ['ID', 'Gene', 'BC0', 'BC0-fun', 'BC1', 'BC1-fun', 'BC2', 'BC2-fun', 'BC3', 'BC3-fun', 'BC4', 'BC4-fun', 'PromIde', 'PromSeg', 'GeneIde', 'GeneSeg'])
        ### save this dataframe as csv table
        final_table.to_csv('output/MultiFunPlotData.csv')

    else:
        ### load data file from github
        url_data = requests.get('https://raw.githubusercontent.com/marketavlkova/PromSEX/master/MultiFunPlotData.csv').content
        final_table = pd.read_csv(StringIO(url_data.decode('utf-8')))
        del final_table['Unnamed: 0']

        url_key = requests.get('https://raw.githubusercontent.com/marketavlkova/PromSEX/master/MultiFunKey.csv').content
        key = pd.read_csv(StringIO(url_key.decode('utf-8')), delimiter = ':')
        ### convert the key table into a dictionary
        all_keys = dict([(i, a) for i, a in zip(key.ID, key.function)])
        all_keys['ORFs'] = 'ORFs'
        all_keys['Conserved-Hypothetical-ORFs'] = 'CH-ORFs'
        all_keys['Unclassified-Genes'] = 'Unclassified'
        ### define output file the the plot
        output_file("BokehPlot.html")


    #########################################
    ######## GENERAL PLOT DEFINITION ########
    #########################################
    ### define tools you want to be present in the plots
    tls = "pan, tap, hover, box_zoom, box_select, wheel_zoom, reset, save, crosshair"
    ### define source to enable linked selection of genes between the plots
    src = ColumnDataSource(data = dict(x0 = final_table['PromIde'], y0 = final_table['GeneIde'], x1 = final_table['PromSeg'], y1 = final_table['GeneSeg'], z0 = final_table['Gene'], z1 = final_table['BC0'], z2 = final_table['BC0-fun'], z3 = final_table['BC1-fun'], z4 = final_table['BC2-fun'], z5 = final_table['BC3-fun'], z6 = final_table['BC4-fun']))
    ### add information boxes to an inspection tool (hover)
    tltips = [('Gene', '@z0'),('Main group', '@z2'), ('Subgroups', '@z3'), ('', '@z4'), ('', '@z5'), ('', '@z6')]

    ### extract BC groups' IDs from 'BC0' column
    groups = []
    for zero in final_table['BC0']:
        if zero not in groups:
            groups.append(zero)
    ### sort the groups alphabetically
    groups = sorted(groups)

    ### add function descriptions for the ID extracted from 'BC0' in a dictionary
    plot_key = {}
    for grp in groups:
        for k in all_keys.keys():
            if grp == k:
                plot_key[all_keys[k]] = grp

    ### define views to distinguish BC groups during plotting
    view = {}
    for pk in sorted(plot_key.values()):
        view[pk] = CDSView(source = src, filters = [GroupFilter(column_name = 'z1', group = pk)])

    #########################################
    ########## 1st PLOT DEFINITION ##########
    #########################################
    ### check for maximal values to be used in both x and y axis of the plot
    mprom = round(max(final_table['PromSeg']), 1)
    mgene = round(max(final_table['GeneSeg']), 1)
    if (mprom > mgene):
        mplot = mprom + 0.01
    else:
        mplot = mgene + 0.01
    ### set the figure dimension and other characteristics
    pl1 = figure(tools = tls, tooltips = tltips, plot_width = 750, plot_height = 750, x_range = (-0.012, mplot), y_range = (-0.012, mplot), x_axis_label = 'Promoters', y_axis_label = 'Genes')
    ### set plot main title
    pl1.title.text = 'Proportion of segregating sites'
    ### set the Seg values to be in the plotted
    n = 0
    for v in list(plot_key.keys()):
        pl1.scatter('x1', 'y1', source = src, view = view[plot_key[v]], legend_label = v, color = Category10_10[n])
        n += 1
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl1.legend.location = 'top_left'
    pl1.legend.click_policy = 'hide'

    #########################################
    ########## 2nd PLOT DEFINITION ##########
    #########################################
    ### check for maximal values to be used in both x and y axis of the plot
    mprom = round(max(final_table['PromIde']), 1)
    mgene = round(max(final_table['GeneIde']), 1)
    if (mprom > mgene):
        mplot = mprom + 0.4
    else:
        mplot = mgene + 0.4

    ### set the figure dimension and other characteristics
    pl2 = figure(tools = tls, tooltips = tltips, plot_width = 750, plot_height = 750, x_range = (-0.4, mplot), y_range = (-0.4, mplot), x_axis_label = 'Promoters', y_axis_label = 'Genes')
    ### set plot main title
    pl2.title.text = '100 - average pairwise identity'
    ### set the Ide values to be in the plotted
    n = 0
    for v in list(plot_key.keys()):
        pl2.scatter('x0', 'y0', source = src, view = view[plot_key[v]], legend_label = v, color = Category10_10[n])
        n += 1
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl2.legend.location = 'top_left'
    pl2.legend.click_policy = 'hide'

    ### set layout, tab and panel for plot
    lay0 = layout(gridplot([[pl2, pl1]]))
    tab0 = Panel(child = lay0, title = 'Main function groups')
    see = Tabs(tabs = [tab0])
    ### execute plotting of both defined plots
    show(see)

if __name__ == '__main__':
    main()
