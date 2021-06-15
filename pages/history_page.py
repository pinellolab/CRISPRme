#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 18:59:01 2020

@author: francesco
"""
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import base64  # for decoding upload content
import io  # for decoding upload content
import sys
import os
from os import listdir
from os.path import isfile, isdir, join  # for getting directories
import pandas as pd
from app import app
# GUImessage as Gmsg
from . import GUImessage
# from .GUImessage import GUImessage as Gmsg
from .main_page import availableGenomes, availablePAM
import math
from app import app, current_working_directory, URL


def supportFilterHistory(result_df, genome_f, pam_f):
    if genome_f is not None:
        result_df.drop(
            result_df[(result_df['Genome'] != genome_f)].index, inplace=True)
    if pam_f is not None:
        keep_values = []
        for index, row in result_df.iterrows():
            if row.PAM not in pam_f:
                keep_values.append(index)
        result_df.drop(labels=keep_values, inplace=True)

    max_page = len(result_df.index)
    max_page = math.floor(max_page / 1000000) + 1
    return result_df, max_page


# @app.callback(
#     Output('div-history-filter-query', 'children'),
#     [Input('button-filter-history', 'n_clicks')],
#     [State('dropdown-genomes-history', 'value'),
#      State('dropdown-pam-history', 'value')]
# )
# def updateHistoryFilter(n, genome, pam):
#     if n is None:
#         raise PreventUpdate
#     return str(genome) + ',' + str(pam)


@app.callback(
    Output('results-table', 'style_data_conditional'),
    [Input('results-table', 'selected_cells')],
    [State('results-table', 'data')]
)
def highlightRow(sel_cel, all_guides):
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    job_name = all_guides[int(sel_cel[0]['row'])]['Job']
    return [
        {
            'if': {
                'filter_query': '{Job} eq "' + job_name + '"'
            },
            'background-color': 'rgba(0, 0, 255,0.15)'

        }
    ]


def get_results():
    results_dirs = [f for f in listdir(current_working_directory + '/Results/') if isdir(join(
        current_working_directory + '/Results/', f)) and isfile(current_working_directory + '/Results/' + f + '/.Params.txt')]
    col = 'Job\tGenome_Selected\tVariant_Selected\tMismatches\tDNA_bulge\tRNA_bulge\tPAM\tNumber_Guides'
    resultParamDataframe = pd.DataFrame(columns=col.split('\t'))
    for job in results_dirs:
        try:
            if os.path.exists(current_working_directory + '/Results/' + job + '/.Params.txt'):
                with open(current_working_directory + '/Results/' + job + '/.Params.txt') as p:
                    all_params = p.read()
                    mms = (next(s for s in all_params.split('\n')
                                if 'Mismatches' in s)).split('\t')[-1]
                    genome_selected = (next(s for s in all_params.split(
                        '\n') if 'Genome_selected' in s)).split('\t')[-1]
                    with open(current_working_directory + '/Results/' + job + '/log.txt') as lo:
                        all_log = lo.read()
                    job_start = (next(s for s in all_log.split('\n')
                                      if 'Job\tStart' in s)).split('\t')[-1]
                    if '+' in genome_selected:
                        genome_selected = genome_selected.split('+')[0] + '+'
                    dna = (next(s for s in all_params.split(
                        '\n') if 'DNA' in s)).split('\t')[-1]
                    rna = (next(s for s in all_params.split(
                        '\n') if 'RNA' in s)).split('\t')[-1]
                    genome_idx = (next(s for s in all_params.split(
                        '\n') if 'Genome_idx' in s)).split('\t')[-1]
                    if '+' in genome_idx:
                        splitted = genome_idx.split(',')
                        genome_idx = []
                        for ele in splitted:
                            genome_idx.append(ele.split('+')[-1])
                        genome_idx = ','.join(genome_idx)
                    else:
                        genome_idx = 'Reference'
                    pam = (next(s for s in all_params.split(
                        '\n') if 'Pam' in s)).split('\t')[-1]
                    # comparison=(next(s for s in all_params.split(
                    #     '\n') if 'Ref_comp' in s)).split('\t')[-1]
                    if os.path.exists(current_working_directory + '/Results/' + job + '/.guides.txt'):
                        with open(current_working_directory + '/Results/' + job + '/.guides.txt') as g:
                            n_guides = str(len(g.read().strip().split('\n')))
                    else:
                        n_guides = 'NA'
                    resultParamDataframe = resultParamDataframe.append({'Job': job, 'Genome_Selected': genome_selected, 'Variant_Selected': genome_idx, 'Mismatches': mms, 'DNA_bulge': dna,
                                                                        'RNA_bulge': rna, 'PAM': pam, 'Number_Guides': n_guides, 'Start': job_start}, ignore_index=True)
        except:
            continue
    try:
        resultParamDataframe['Start'] = pd.to_datetime(
            resultParamDataframe['Start'])
        resultParamDataframe.sort_values(
            by=['Start'], inplace=True, ascending=False)
    except:
        pass
    # resultParamDataframe = resultParamDataframe.sort_values(
    #     ['Mismatches', 'DNA_bulge', 'RNA_bulge'], ascending=[True, True, True])
    return resultParamDataframe


# # Remove Results and Apply Filtering
# @app.callback(
#     [Output('div-history-table', 'children'),
#      Output('div-current-page-history', 'children')],
#     [Input('button-delete-history-0', 'n_clicks_timestamp'),
#      Input('button-delete-history-1', 'n_clicks_timestamp'),
#      Input('button-delete-history-2', 'n_clicks_timestamp'),
#      Input('button-delete-history-3', 'n_clicks_timestamp'),
#      Input('button-delete-history-4', 'n_clicks_timestamp'),
#      Input('button-delete-history-5', 'n_clicks_timestamp'),
#      Input('button-delete-history-6', 'n_clicks_timestamp'),
#      Input('button-delete-history-7', 'n_clicks_timestamp'),
#      Input('button-delete-history-8', 'n_clicks_timestamp'),
#      Input('button-delete-history-9', 'n_clicks_timestamp'),
#      Input('prev-page-history', 'n_clicks_timestamp'),
#      Input('next-page-history', 'n_clicks_timestamp'),
#      Input('div-history-filter-query', 'children')],
#     [State('button-delete-history-0', 'data-jobid'),
#      State('button-delete-history-1', 'data-jobid'),
#      State('button-delete-history-2', 'data-jobid'),
#      State('button-delete-history-3', 'data-jobid'),
#      State('button-delete-history-4', 'data-jobid'),
#      State('button-delete-history-5', 'data-jobid'),
#      State('button-delete-history-6', 'data-jobid'),
#      State('button-delete-history-7', 'data-jobid'),
#      State('button-delete-history-8', 'data-jobid'),
#      State('button-delete-history-9', 'data-jobid'),
#      State('button-filter-history', 'n_clicks_timestamp'),
#      State('url', 'search'),
#      State('div-current-page-history', 'children')
#      ]
# )
# def removeJobIDandFilter(n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, nPrev, nNext, filter_q, jID0, jID1, jID2, jID3, jID4, jID5, jID6, jID7, jID8, jID9, n, search, current_page):
#     # Get last pressed button
#     if not n0:
#         n0 = 0
#     if not n1:
#         n1 = 0
#     if not n2:
#         n2 = 0
#     if not n3:
#         n3 = 0
#     if not n4:
#         n4 = 0
#     if not n5:
#         n5 = 0
#     if not n6:
#         n6 = 0
#     if not n7:
#         n7 = 0
#     if not n8:
#         n8 = 0
#     if not n9:
#         n9 = 0
#     if not nPrev:
#         nPrev = 0
#     if not nNext:
#         nNext = 0
#     if not n:
#         n = 0
#     btn_group = []
#     btn_group.append(n0)
#     btn_group.append(n1)
#     btn_group.append(n2)
#     btn_group.append(n3)
#     btn_group.append(n4)
#     btn_group.append(n5)
#     btn_group.append(n6)
#     btn_group.append(n7)
#     btn_group.append(n8)
#     btn_group.append(n9)
#     btn_group.append(nPrev)
#     btn_group.append(nNext)
#     btn_group.append(n)

#     genome_filter = filter_q.split(',')[0]
#     pam_filter = filter_q.split(',')[1]
#     if genome_filter == 'None':
#         genome_filter = None
#     if pam_filter == 'None':
#         pam_filter = None
#     current_page = current_page.split('/')[0]
#     current_page = int(current_page)
#     results = get_results()
#     selectedID = ''

#     if max(btn_group) == 0:
#         selectedID = ''
#     elif max(btn_group) == n0:
#         selectedID = jID0
#     elif max(btn_group) == n1:
#         selectedID = jID1
#     elif max(btn_group) == n2:
#         selectedID = jID2
#     elif max(btn_group) == n3:
#         selectedID = jID3
#     elif max(btn_group) == n4:
#         selectedID = jID4
#     elif max(btn_group) == n5:
#         selectedID = jID5
#     elif max(btn_group) == n6:
#         selectedID = jID6
#     elif max(btn_group) == n7:
#         selectedID = jID7
#     elif max(btn_group) == n8:
#         selectedID = jID8
#     elif max(btn_group) == n9:
#         selectedID = jID9
#     elif max(btn_group) == n:  # Filter Button selected, return the first page of the filtered table
#         results, max_page = supportFilterHistory(
#             results, genome_filter, pam_filter)
#         return generate_table_results(results, 1), '1/' + str(max_page)
#     elif max(btn_group) == nNext:  # Next Button of the table
#         current_page = current_page + 1
#         results, max_page = supportFilterHistory(
#             results, genome_filter, pam_filter)
#         if ((current_page - 1) * 10) > len(results):
#             current_page = current_page - 1
#             if current_page < 1:
#                 current_page = 1
#         return generate_table_results(results, current_page), str(current_page) + '/' + str(max_page)
#     elif max(btn_group) == nPrev:  # Go to previous page
#         current_page = current_page - 1
#         if current_page < 1:
#             current_page = 1
#         results, max_page = supportFilterHistory(
#             results, genome_filter, pam_filter)
#         return generate_table_results(results, current_page), str(current_page) + '/' + str(max_page)

#     if selectedID == '':
#         raise PreventUpdate
#     result_removed = GUImessage.deleteResultConfirm(selectedID)
#     if result_removed:  # If the result was removed, update the table
#         results, max_page = supportFilterHistory(
#             get_results(), genome_filter, pam_filter)
#         return generate_table_results(results, 1), '1/' + str(max_page)
#     else:
#         raise PreventUpdate
#     results, max_page = supportFilterHistory(
#         get_results(), genome_filter, pam_filter)
#     return generate_table_results(results, 1), '1/' + str(max_page)


def generate_table_results(dataframe, page, max_rows=1000000):
    '''
    Generate table for History page
    '''
    fl = []
    rows_remaining = len(dataframe) - (page - 1) * max_rows
    header = html.Thead(
        html.Tr([html.Th(col, style={'vertical-align': 'middle', 'text-align': 'center'}) if col != 'Load' and col !=
                 'Delete' else html.Th('', style={'vertical-align': 'middle', 'text-align': 'center'}) for col in dataframe.columns])
    )
    body_history = []
    add_button = 0
    for i in range(min(rows_remaining, max_rows)):
        add_button += 1
        row_hist = []
        for col in dataframe.columns:
            if col == 'Job':
                jobID = str(dataframe.iloc[i + (page - 1)*max_rows][col])
                row_hist.append(html.Td(html.A(jobID, target='_blank', href=URL + '/load?job=' + jobID), style={
                    'vertical-align': 'middle', 'text-align': 'center'}))
            else:
                row_hist.append(
                    html.Td(dataframe.iloc[i + (page - 1)*max_rows][col], style={
                            'vertical-align': 'middle', 'text-align': 'center'})
                )
        body_history.append(html.Tr(row_hist))
    fl.append(
        html.Table([
            header,
            html.Tbody(body_history)
        ], style={'display': 'inline-block'},)
    )

    for i in range(add_button, 10):  # Add hidden buttons for callback removeJobId compatibility
        fl.append(html.Button(
            str(i), id='button-delete-history-'+str(i), **{'data-jobid': 'None'}, style={'display': 'none'}
        ))
    return fl


def historyPage():
    '''
    Create the History page
    '''
    results = get_results()
    final_list = []

    final_list.append(
        html.Div(
            [
                html.H3('Results History'),
                html.P(
                    'List of available results. Click on the link to open the corresponding load page in a new tab.')
            ]
        )
    )
    # final_list.append(
    #     html.Div
    #     (
    #         [
    #             dbc.Row(
    #                 [
    #                     dbc.Col(html.Div(dcc.Dropdown(options=availableGenomes(
    #                     ), id='dropdown-genomes-history', placeholder='Select a Genome'))),
    #                     dbc.Col(html.Div(dcc.Dropdown(options=availablePAM(
    #                     ), id='dropdown-pam-history', placeholder='Select a PAM'))),
    #                     dbc.Col(html.Div(html.Button(
    #                         'Filter', id='button-filter-history')))
    #                 ]
    #             ),
    #         ],
    #     )
    # )
    final_list.append(html.Div(
        'None,None', id='div-history-filter-query', style={'display': 'none'}))

    final_list.append(
        html.Div(
            generate_table_results(results, 1),
            id='div-history-table',
            style={'text-align': 'center'}
        ),
    )
    final_list.append(
        html.Div(id='div-remove-jobid', style={'display': 'none'})
    )
    final_list.append(
        html.Div(
            [
                html.Br(),
                # html.Button('Prev', id='prev-page-history'),
                # html.Button('Next', id='next-page-history')
            ],
            style={'text-align': 'center'}
        )
    )
    max_page = len(results.index)
    max_page = math.floor(max_page / 1000000) + 1
    # final_list.append(html.Div('1/' + str(max_page),
    #                            id='div-current-page-history'))
    page = html.Div(final_list, style={'margin': '1%'})
    return page
