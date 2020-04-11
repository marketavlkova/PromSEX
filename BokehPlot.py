### Script combining all data for MulitFunPlot and producing the plot(s)

### example to run: python3 BokehPlot.py
### e.g.: python3 BokehPlot.py

import sys
import os.path
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, CDSView
from bokeh.models.filters import GroupFilter

def main():

    ### if the file for 'final_table' already exists
    ### load it and skip the part creating it
    if os.path.isfile('output/MultiFunPlotData.csv'):
        final_table = pd.read_csv('output/MultiFunPlotData.csv', delimiter = ',')
    else:
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
                    ### in all other rounds
                    else:
                        ### if the round number is lower than number of subgroups
                        if (id < len(ids)):
                            ### produce the current subgroup combining subgroup from the previous 'line' value and current extention to it
                            bc = '.'.join([line[(id + 1)], ids[id]])
                        ### if you run out of subgroup values
                        else:
                            ### set bc variable to NaN
                            bc = 'NaN'
                        ### add the BC subgroup value to the 'line' list
                        line.append(bc)
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
        final_table = pd.DataFrame(vec, columns = ['ID', 'Gene', 'ParentBC0', 'ParentBC1', 'ParentBC2', 'ParentBC3', 'ParentBC4', 'PromIde', 'PromSeg', 'GeneIde', 'GeneSeg'])
        ### save this dataframe as csv table
        final_table.to_csv('output/MultiFunPlotData.csv')

    ### load file containing descriptions to BC group IDs
    key = pd.read_csv('output/MultiFunKey.csv', delimiter = ':')

    #########################################
    ######## GENERAL PLOT DEFINITION ########
    #########################################
    ### define output file the the plot
    output_file("output/BC0.html")
    ### define tools you want to be present in the plots
    tls = "pan, box_zoom, box_select, wheel_zoom, reset, save, crosshair"
    ### define source to enable linked selection of genes between the plots
    src = ColumnDataSource(data = dict(x0 = final_table['PromIde'], y0 = final_table['GeneIde'], x1 = final_table['PromSeg'], y1 = final_table['GeneSeg'], z0 = final_table['Gene'], z1 = final_table['ParentBC0']))

    ### extract BC groups' IDs from 'ParentBC0' column
    groups = []
    for zero in final_table['ParentBC0']:
        if zero not in groups:
            groups.append(zero)
    ### sort the groups alphabetically
    groups = sorted(groups)

    ### extract descriptions to BC groups IDs
    desc = []
    for grp in groups:
        for k in range(0, len(key['ID'])):
            if (grp == key['ID'][k]):
                desc.append(key['function'][k])

    ### combine BC group IDs and their descriptions in a dictionaty
    plot_key = {}
    for i in range(0, len(groups)):
        if all(['BC-' in groups[i], groups[i] not in plot_key.keys()]):
            plot_key[desc[i]] = groups[i]
    plot_key['ORFs'] = 'ORFs'
    plot_key['CH-ORFs'] = 'Conserved-Hypothetical-ORFs'
    plot_key['Unclassified'] = 'Unclassified-Genes'

    ### define views to distinguish BC groups during plotting
    view = {}
    for pk in sorted(plot_key.values()):
        view[pk] = CDSView(source = src, filters = [GroupFilter(column_name = 'z1', group = pk)])

    ### define colours for plotting
    cols = ['black', 'blue', 'red', 'green', 'orange', 'cyan', 'magenta', 'lime', 'purple', 'sienna']

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
    pl1 = figure(tools = tls, plot_width = 800, plot_height = 800, x_range = (-0.012, mplot), y_range = (-0.012, mplot), x_axis_label = 'Promoters', y_axis_label = 'Genes')
    ### set plot main title
    pl1.title.text = 'Proportion of segregating sites'
    ### set the Seg values to be in the plotted
    n = 0
    for v in list(plot_key.keys()):
        pl1.scatter('x1', 'y1', source = src, view = view[plot_key[v]], legend_label = v, color = cols[n])
        n += 1
    # pl1.text(fileBC0.PromSeg.loc['BC-1'], fileBC0.GeneSeg.loc['BC-1'], text = fileBC0.Gene.loc['BC-1'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Metabolism', color = 'black')
    # pl1.text(fileBC0.PromSeg.loc['BC-2'], fileBC0.GeneSeg.loc['BC-2'], text = fileBC0.Gene.loc['BC-2'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Information transfer', color = 'blue')
    # pl1.text(fileBC0.PromSeg.loc['BC-3'], fileBC0.GeneSeg.loc['BC-3'], text = fileBC0.Gene.loc['BC-3'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Regulation', color = 'red')
    # pl1.text(fileBC0.PromSeg.loc['BC-4'], fileBC0.GeneSeg.loc['BC-4'], text = fileBC0.Gene.loc['BC-4'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Transport', color = 'green')
    # pl1.text(fileBC0.PromSeg.loc['BC-5'], fileBC0.GeneSeg.loc['BC-5'], text = fileBC0.Gene.loc['BC-5'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Cell processes', color = 'orange')
    # pl1.text(fileBC0.PromSeg.loc['BC-6'], fileBC0.GeneSeg.loc['BC-6'], text = fileBC0.Gene.loc['BC-6'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Cell structure', color = 'cyan')
    # pl1.text(fileBC0.PromSeg.loc['BC-8'], fileBC0.GeneSeg.loc['BC-8'], text = fileBC0.Gene.loc['BC-8'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Extrachromosomal', color = 'magenta')
    # pl1.text(fileBC0.PromSeg.loc['ORFs'], fileBC0.GeneSeg.loc['ORFs'], text = fileBC0.Gene.loc['ORFs'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'ORFs', color = 'lime')
    # pl1.text(fileBC0.PromSeg.loc['Conserved-Hypothetical-ORFs'], fileBC0.GeneSeg.loc['Conserved-Hypothetical-ORFs'], text = fileBC0.Gene.loc['Conserved-Hypothetical-ORFs'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'CH-ORFs', color = 'purple')
    # pl1.text(fileBC0.PromSeg.loc['Unclassified-Genes'], fileBC0.GeneSeg.loc['Unclassified-Genes'], text = fileBC0.Gene.loc['Unclassified-Genes'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Unclassified', color = 'sienna')
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
        mplot = mprom + 0.1
    else:
        mplot = mgene + 0.1

    ### set the figure dimension and other characteristics
    pl2 = figure(tools = tls, plot_width = 800, plot_height = 800, x_range = (-0.4, 8.5), y_range = (-0.4, 8.5), x_axis_label = 'Promoters', y_axis_label = 'Genes')
    ### set plot main title
    pl2.title.text = '100 - average pairwise identity'
    ### set the Ide values to be in the plotted
    n = 0
    for v in list(plot_key.keys()):
        pl2.scatter('x0', 'y0', source = src, view = view[plot_key[v]], legend_label = v, color = cols[n])
        n += 1
    # pl2.text(fileBC0.PromIde.loc['BC-1'], fileBC0.GeneIde.loc['BC-1'], text = fileBC0.Gene.loc['BC-1'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Metabolism', color = 'black')
    # pl2.text(fileBC0.PromIde.loc['BC-2'], fileBC0.GeneIde.loc['BC-2'], text = fileBC0.Gene.loc['BC-2'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Information transfer', color = 'blue')
    # pl2.text(fileBC0.PromIde.loc['BC-3'], fileBC0.GeneIde.loc['BC-3'], text = fileBC0.Gene.loc['BC-3'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Regulation', color = 'red')
    # pl2.text(fileBC0.PromIde.loc['BC-4'], fileBC0.GeneIde.loc['BC-4'], text = fileBC0.Gene.loc['BC-4'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Transport', color = 'green')
    # pl2.text(fileBC0.PromIde.loc['BC-5'], fileBC0.GeneIde.loc['BC-5'], text = fileBC0.Gene.loc['BC-5'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Cell processes', color = 'orange')
    # pl2.text(fileBC0.PromIde.loc['BC-6'], fileBC0.GeneIde.loc['BC-6'], text = fileBC0.Gene.loc['BC-6'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Cell structure', color = 'cyan')
    # pl2.text(fileBC0.PromIde.loc['BC-8'], fileBC0.GeneIde.loc['BC-8'], text = fileBC0.Gene.loc['BC-8'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Extrachromosomal', color = 'magenta')
    # pl2.text(fileBC0.PromIde.loc['ORFs'], fileBC0.GeneIde.loc['ORFs'], text = fileBC0.Gene.loc['ORFs'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'ORFs', color = 'lime')
    # pl2.text(fileBC0.PromIde.loc['Conserved-Hypothetical-ORFs'], fileBC0.GeneIde.loc['Conserved-Hypothetical-ORFs'], text = fileBC0.Gene.loc['Conserved-Hypothetical-ORFs'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'CH-ORFs', color = 'purple')
    # pl2.text(fileBC0.PromIde.loc['Unclassified-Genes'], fileBC0.GeneIde.loc['Unclassified-Genes'], text = fileBC0.Gene.loc['Unclassified-Genes'], y_offset = -2.5, x_offset = -7, text_font_size = {'value':'8pt'}, legend_label = 'Unclassified', color = 'sienna')
    ### set legend position and the attribute to hide/show it upon clicking in the legend
    pl2.legend.location = 'top_left'
    pl2.legend.click_policy = 'hide'

    ### execute plotting of both defined plots
    show(gridplot([[pl2, pl1]]))

if __name__ == '__main__':
    main()
