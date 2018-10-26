#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 2.7.0 ###########################################################
# Libraries ####################################################################
import csv
import pprint
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import numpy as np
import itertools

################################################################################
# Description ##################################################################
################################################################################
"""

"""

################################################################################
# Base-level Functions #########################################################
################################################################################
def file_to_datalist(data_filename):
    """
    Starting with a filename, import and convert data to a list.
    """

    with open(data_filename, 'r') as data_file:
        csv_reader = csv.reader(data_file, delimiter = "\t")
        all_data = list(csv_reader)

    return all_data

################################################################################
# Task-specific Functions ######################################################
################################################################################
def data_to_dict():
    """
    Parse array data to dict.
    """

    data = file_to_datalist(expression_data_filename)
    headers = data[128][1:]
    expression_data_dict = {}
    for row in data[129:]:
        array_id = row[0]
        express_values = [float(x) for x in row[1:]]

        if array_id not in expression_data_dict:
            expression_data_dict[array_id] = express_values
        else:
            print array_id, "already in dict"

    return expression_data_dict, headers


def annotations_to_dict():
    """
    Parse array annotations to dict.
    """
    annotations_data = file_to_datalist(annot_filename)
    annotations_data_dict = {}
    for row in annotations_data[21:]:
        array_id = row[0]
        external_ids = row[1:]
        annotations_data_dict[array_id] = external_ids

    return annotations_data_dict


def condition_status():
    """
    Extract condition fusion status for the experiments.
    """

    experiment_descriptions = file_to_datalist(experiment_description_filename)
    condition_exp_dict = {}
    for row in experiment_descriptions[1:]:
        experiment_id = row[0].split(" ")[0]
        condition_status = row[condition_index]

        if condition_status.lower() == "positive" or condition_status.lower() == "negative":
            if condition_status in condition_exp_dict:
                condition_exp_dict[condition_status].append(experiment_id)
            else:
                condition_exp_dict[condition_status] = [experiment_id]

    return condition_exp_dict


def expression_by_condition():
    """
    Compile expression by each gene for pos/neg condition fusion gene.
    """

    target_gene_expression_dict = {}
    for target_gene, ncbi_id in target_genes.iteritems():

        # identify the right array probe id based on known NCBI identifier
        for array_id, external_ids in annotations_data_dict.iteritems():
            if ncbi_id in external_ids[1]:
                target_gene_array_id = array_id
                target_gene_rename = "-".join([target_gene, target_gene_array_id])

        target_gene_expression_dict[target_gene_rename] = {}
        
        # for each condition_status retrieve the expression values for the target gene
        # from the associated experiments
        for condition_status, experiment_ids in condition_exp_dict.iteritems():
            target_gene_expression_dict[target_gene_rename][condition_status] = []
            
            for experiment_id in experiment_ids:
                target_gene_expression_index = headers.index(experiment_id)
                target_gene_expression = expression_data_dict[target_gene_array_id][target_gene_expression_index]
                target_gene_expression_dict[target_gene_rename][condition_status].append(target_gene_expression)

    return target_gene_expression_dict


def stat_test():
    """
    Perform some statistical test on target sets of data.
    """

    # perform independent t-test to determine if expression of genes is different
    # between the e2a-pbx fusion gene pos/neg groups
    for target_gene, condition_status_dict in target_gene_expression_dict.iteritems():
        stat, pval = stats.ttest_ind(condition_status_dict.values()[0], condition_status_dict.values()[1])
        print target_gene
        print stat, pval    
    

def grouped_lmplot_builder():
    """
    Build a boxplot for each individual CAGE peak.
    Error bars solution from here:
    https://gist.github.com/CnrLwlss/22f4abbb11e4ab01ae8bb35e041941bb#file-customerrorbars-py
    """

    ## collate data for figure
    all_data = []
    for target_gene, condition_status_dict in target_gene_expression_dict.iteritems():
        for condition_status, expression_values in condition_status_dict.iteritems():
            for expression_value in expression_values:
                all_data.append([target_gene, condition_status, expression_value])

    target_genes = [x[0] for x in all_data]
    conditions = [x[1] for x in all_data]
    expr_vals = [x[2] for x in all_data]

    all_data_dict = {"target_genes":target_genes, condition:conditions, "expr_vals":expr_vals}
    all_data_df = pd.DataFrame(data = all_data_dict)
    all_data_df.sort_values(["target_genes", condition], ascending=False)

    # plot data
    # figure start
    plt.figure(figsize=(10,10))

    # plot order
    target_genes_ordered = list(set([x[0] for x in all_data]))
    target_genes_ordered.sort()
    condition_ordered = list(set([x[1] for x in all_data]))
    condition_ordered.sort()

    #plot
    g = sns.swarmplot(x="target_genes", y="expr_vals", hue=condition, data=all_data_df, palette="Blues_d", dodge=True, order=target_genes_ordered, hue_order=condition_ordered)

    # error bars
    xcentres=[0,1,2,3]
    delt=0.2
    xneg=[x-delt for x in xcentres]
    xpos=[x+delt for x in xcentres]
    xvals=xneg+xpos
    xvals.sort()
    yvals = all_data_df.groupby(['target_genes', condition]).expr_vals.mean()
    yerr = all_data_df.groupby(['target_genes', condition]).expr_vals.std()
    (_, caps, _)=g.errorbar(x=xvals,y=yvals,yerr=yerr,capsize=4,ecolor="black", barsabove=True)
    for cap in caps:
        cap.set_markeredgewidth(2)

    # plot manip
    g.set(ylim=(-100, 2000))
    sns.despine(offset=10, trim=True)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # save plot
    bp_figure_filename = os.path.join("./", ".".join([condition, "swarmerror", "eps"]))
    plt.savefig(bp_figure_filename, format="eps")


################################################################################
# Initiating Variables #########################################################
################################################################################
expression_data_filename = "./data/GSE11877_series_matrix.txt"
annot_filename = "./data/A-AFFY-44.adf.txt"
experiment_description_filename = "./data/E-GEOD-11877.sdrf.txt"
target_genes = {"ROR1":"BC006374", "WNT5A":"AI968085", "WNT66":"AF169963", "LRP6":"AV725248"}
conditions_index_dict = {"E2A-PBX":5}
################################################################################
# Execution ####################################################################
################################################################################
expression_data_dict, headers = data_to_dict()
annotations_data_dict = annotations_to_dict()

for condition, condition_index in conditions_index_dict.iteritems():
    print condition
    condition_exp_dict = condition_status()
    target_gene_expression_dict = expression_by_condition()
    grouped_lmplot_builder()
    stat_test()
