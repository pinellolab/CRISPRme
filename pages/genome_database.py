#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 17:14:56 2020

@author: francesco
"""
import os
import pandas as pd
from os.path import isfile, isdir, join
from os import listdir
from app import app, current_working_directory
import dash_table
import dash_daq as daq
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import dash
from dash.exceptions import PreventUpdate
import subprocess
from . import genome_database
from . import main_page
# from genome_database import get_genomes


def get_genomes(pathDir):

    #dir_path = os.path.dirname(os.path.realpath(__file__))+"/logGenomes/"
    #logs = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]

    genomes = pd.DataFrame(columns=["Reference Genome", "Enriched Genome",
                                    "Index PAM", "Annotation File", "Samples ID File", "# Bulges"])
    """
    for log in logs:
        with open(join(dir_path, log)) as f:
            refgen = os.path.basename(f.readline().split("\t")[1])
            vcf = os.path.basename(f.readline().split("\t")[1])
            pam = os.path.basename(f.readline().split("\t")[1])
            ann = os.path.basename(f.readline().split("\t")[1])
            samples = os.path.basename(f.readline().split("\t")[1])
            bulges = f.readline().split("\t")[1]
            enrgen = f.readline().split("\t")[1]
    """
    genomes_dirs = [f for f in os.listdir(
        pathDir+"/Genomes/") if isdir(join(pathDir+"/Genomes/", f))]
    genome_library_dirs = [f for f in os.listdir(
        pathDir+"/genome_library/") if isdir(join(pathDir+"/genome_library/", f))]

    for d in genomes_dirs:
        if "+" not in d:  # Reference genome
            genRef = d
            genEnr = "NA"
            #pam = "NA"
            ann = "NA"
            samples = "NA"
            #bulges = "NA"
            pams = []
            bMaxs = []
            for lib in genome_library_dirs:
                libPars = "_".join(lib.split("_")[2:])
                if genRef == libPars:
                    partsLib = lib.split("_")
                    pams.append(partsLib[0])
                    bMaxs.append(partsLib[1])
            for annot in os.listdir(pathDir+"/annotations/"):
                if genRef == annot.split(".")[0]:
                    ann = annot
                    break
            if len(pams) > 0:
                for i in range(len(pams)):
                    genomes = genomes.append({"Reference Genome": genRef, "Enriched Genome": genEnr,
                                              "Index PAM": pams[i], "Annotation File": ann, "Samples ID File": samples, "# Bulges": bMaxs[i]}, ignore_index=True)
            else:
                genomes = genomes.append({"Reference Genome": genRef, "Enriched Genome": genEnr,
                                          "Index PAM": "NA", "Annotation File": ann, "Samples ID File": samples, "# Bulges": "NA"}, ignore_index=True)
        else:  # Enriched Genome
            genEnr = d
            genRef = d.split("+")[0]
            ann = "NA"
            samples = "NA"
            pams = []
            bMaxs = []
            for lib in genome_library_dirs:
                libPars = "_".join(lib.split("_")[2:])
                if genEnr == libPars:
                    partsLib = lib.split("_")
                    pams.append(partsLib[0])
                    bMaxs.append(partsLib[1])
            for annot in os.listdir(pathDir+"/annotations/"):
                if genRef == annot.split(".")[0]:
                    ann = annot
                    break
            for s in os.listdir(pathDir+"/samplesID/"):
                if genEnr == s[8:len(s)-4]:
                    samples = s
                    break
            if len(pams) > 0:
                for i in range(len(pams)):
                    genomes = genomes.append({"Reference Genome": genRef, "Enriched Genome": genEnr,
                                              "Index PAM": pams[i], "Annotation File": ann, "Samples ID File": samples, "# Bulges": bMaxs[i]}, ignore_index=True)
            else:
                genomes = genomes.append({"Reference Genome": genRef, "Enriched Genome": genEnr,
                                          "Index PAM": "NA", "Annotation File": ann, "Samples ID File": samples, "# Bulges": "NA"}, ignore_index=True)

    return genomes

# Change annotation of selected Genome row


@app.callback(
    Output('ann-job', 'children'),
    [Input('change-ann', 'n_clicks')],
    [State('genomes-table', 'derived_virtual_data'),
     State('genomes-table', 'selected_rows'),
     State('tooltip-label-new-annotation-selected', 'children'),
     State('radioitems-new-annotation', 'value')]
)
def change_annotation(nChg, rows, selected_row, selected_new_ann, selected_type):
    if nChg is None:
        raise PreventUpdate
    if len(selected_row) == 0:
        return dbc.Alert("Warning! Please select a Genome!", is_open=True, duration=5000, color='warning')
    selected_new_ann = selected_new_ann.split('Full Path: ')[-1]
    if not isfile(selected_new_ann):
        return dbc.Alert("Error! The selected file does not exist!", is_open=True, duration=7000, color='danger')
    if selected_type is None:
        return dbc.Alert("Warning! Please select an option (Overwrite or Extend)!", is_open=True, duration=7000, color='warning')
    row = pd.DataFrame(rows).iloc[selected_row, :]
    annotationFile = row["Annotation File"].iloc[0]

    if selected_type == 'replace':
        subprocess.run(['cp', selected_new_ann,
                        current_working_directory + 'annotations/' + annotationFile])
    else:
        with open(current_working_directory + 'annotations/' + annotationFile, 'a') as oldAnn:
            with open(selected_new_ann, 'r') as newAnn:
                for line in newAnn:
                    oldAnn.write("\n"+line.strip())

    # from GUI import annotations as ann
    # ann.startChangeAnn(current_working_directory+"/annotations/"+str(annotationFile))

    # Return an alert if the job is completed
    return dbc.Alert(
        "Completed! The annotation file is updated!",
        id="alert-auto",
        is_open=True,
        duration=5000,
    )


@app.callback(
    Output('genomes-table', 'data'),
    [Input('genomes-table', "page_current"),
     Input('genomes-table', "page_size"),
     Input('genomes-table', 'sort_by'),
     Input('genomes-table', 'filter_query')])
def updateGenomePageTable(page_current, page_size, sort_by, filter):
    filtering_expressions = filter.split(' && ')
    dff = get_genomes(current_working_directory)
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = main_page.split_filter_part(
            filter_part)

        if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == 'contains':
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == 'datestartswith':
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    if len(sort_by):
        dff = dff.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    page = page_current
    size = page_size
    return dff.iloc[page * size: (page + 1) * size].to_dict('records')
