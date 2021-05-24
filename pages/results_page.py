from sqlite3.dbapi2 import Row
import sys
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State
from app import URL, app
# from app import app
import pandas as pd
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_table
from app import current_working_directory, cache, app_main_directory, operators
from PostProcess import CFDGraph
from PostProcess.supportFunctions.loadSample import associateSample
from os.path import isfile, isdir, join  # for getting directories
import os
from os import listdir
import subprocess
import math
import base64  # for decoding upload content
import time
import re
import webbrowser as wb
import sqlite3
from PostProcess import query_manager
import flask
# import send_from_directory


PAGE_SIZE = 10  # number of entries in each page of the table in view report
BARPLOT_LEN = 4  # number of barplots in each row of Populations Distributions
# Columns for dash datatable in REF search
# dff_view_names = ['Bulge_type', 'crRNA', 'Reference_sequence', 'Off_target_motif', 'Chromosome',
#                           'Position', 'Direction', 'Mismatches',
#                           'Bulge_Size', 'PAM_gen', 'Samples', 'SNP',
#                           'CFD', 'CFD_ref', 'Highest_CFD_Risk_Score',
#                           'AF', 'Annotation_Type']
# COL_REF = ['Bulge Type', 'crRNA', 'Off_target_motif', 'Chromosome', 'Position', 'Cluster Position',
#            'Direction', 'Mismatches', 'Bulge Size', 'Total', 'Annotation Type']
# COL_REF_TYPE = ['text', 'text', 'text', 'text', 'numeric',
#                 'numeric', 'text', 'numeric', 'numeric', 'numeric', 'text']
COL_REF = ['Bulge Type', 'crRNA', 'Off target_motif', 'Reference sequence', 'Chromosome',
           'Position', 'Direction', 'Mismatches',
           'Bulge Size', 'PAM gen', 'Samples', 'Variant',
           'CFD', 'CFD ref', 'Highest CFD Risk Score',
           'AF', 'Annotation Type']
COL_REF_TYPE = ['text', 'text', 'text', 'text', 'text', 'numeric',
                'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'text',
                'numeric', 'numeric', 'numeric', 'numeric', 'text']
# COL_REF_RENAME = {0: 'Bulge Type', 1: 'crRNA', 2: 'Off_target_motif', 3: 'Chromosome', 4: 'Position', 5: 'Cluster Position', 6: 'Direction',
#                   7: 'Mismatches', 8: 'Bulge Size', 9: 'Total', 10: 'Correct Guide', 11: 'Annotation Type'}
COL_REF_RENAME = {0: 'Bulge Type', 1: 'crRNA', 2: 'Off target motif', 3: 'Reference sequence', 4: 'Chromosome', 5: 'Position', 6: 'Cluster Position', 7: 'Direction',
                  8: 'Mismatches', 9: 'Bulge Size', 10: 'Total', 11: 'PAM gen', 12: 'Variant Unique', 13: 'Samples', 14: 'Annotation Type', 15: 'Real Guide',
                  16: 'rsID', 17: 'AF', 18: 'Variant', 19: '#Seq in cluster', 20: 'CFD', 21: 'CFD ref', 22: 'Highest CFD Risk Score'}
# Columns for dash datatable in VAR and BOTH search
# COL_BOTH = ['Bulge Type', 'crRNA', 'Off_target_motif', 'Reference_sequence', 'Chromosome', 'Position', 'Cluster Position','Direction', 'Mismatches', 'Bulge Size', 'Total', 'PAM Creation', 'Samples Summary', 'Annotation Type', 'Real_Guide', 'rsID', 'AF', 'SNP', '#Seq_in_cluster', 'CFD', 'CFD_ref']
COL_BOTH = ['Bulge Type', 'crRNA', 'Off target motif', 'Reference sequence', 'Chromosome',
                          'Position', 'Direction', 'Mismatches',
                          'Bulge Size', 'PAM gen', 'Samples', 'Variant',
                          'CFD', 'CFD ref', 'Highest CFD Risk Score',
                          'AF', 'Annotation Type']
COL_BOTH_TYPE = ['text', 'text', 'text', 'text', 'text', 'numeric',
                 'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'text',
                 'numeric', 'numeric', 'numeric', 'numeric', 'text']
COL_BOTH_RENAME = {0: 'Bulge Type', 1: 'crRNA', 2: 'Off target motif', 3: 'Reference sequence', 4: 'Chromosome', 5: 'Position', 6: 'Cluster Position', 7: 'Direction',
                   8: 'Mismatches', 9: 'Bulge Size', 10: 'Total', 11: 'PAM gen', 12: 'Variant Unique', 13: 'Samples', 14: 'Annotation Type', 15: 'Real Guide',
                   16: 'rsID', 17: 'AF', 18: 'Variant', 19: '#Seq in cluster', 20: 'CFD', 21: 'CFD ref', 22: 'Highest CFD Risk Score'}
GENOME_DATABASE = ['Reference', 'Enriched',
                   'Samples', 'Dictionary', 'Annotation']


def resultPage(job_id):
    '''
    La funzione ritorna il layout della pagina risultati (tabella delle guide + eventuali immagini). Nella tabella delle guide
    carico il profile ottenuto dalla ricerca. Carica inoltre l'ACFD, che è il cfd score aggregato per tutti i risultati di una singola guida.
    Crea poi 10 bottoni: un numero pari a mismatches + 2  che sono visibili, il resto con style = {'display':'none'}, così ho sempre il numero
    esatto di bottoni per mismatches in base ai mms dati in input nella ricerca (serve a risolvere problemi con le callback che hanno input
    da elementi non creati. In questo caso io creo tutti i possibili bottoni ma ne rendo visibili/disponibili solo il numero corretto in base
    ai mms).
    '''
    value = job_id
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    warning_message = []
    if (not isdir(job_directory)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))

    count_guides = 0
    with open(current_working_directory + 'Results/' + value + '/guides.txt') as g:
        for line in g:
            count_guides += 1

    # Load mismatches
    with open(current_working_directory + 'Results/' + value + '/Params.txt') as p:
        all_params = p.read()
        real_genome_name = (next(s for s in all_params.split('\n')
                                 if 'Genome_idx' in s)).split('\t')[-1]
        mms = (next(s for s in all_params.split('\n')
                    if 'Mismatches' in s)).split('\t')[-1]
        bulge_dna = (next(s for s in all_params.split(
            '\n') if 'DNA' in s)).split('\t')[-1]
        bulge_rna = (next(s for s in all_params.split(
            '\n') if 'RNA' in s)).split('\t')[-1]
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]
        max_bulges = (next(s for s in all_params.split('\n')
                           if 'Max_bulges' in s)).split('\t')[-1]
        pam_name = (next(s for s in all_params.split('\n')
                         if 'Pam' in s)).split('\t')[-1]

    genome_name = genome_type_f
    if '+' in real_genome_name:
        genome_name = [genome_name]
        splitted_genome_names = real_genome_name.strip().split(',')
        for name in splitted_genome_names:
            name_corrected = name.split("+")[1]
            genome_name.append(name_corrected)
        genome_name = '+'.join(genome_name)
    # genome_name = genome_type_f
    # genome_type = 'ref'
    # if '+' in genome_type_f:
    #     genome_type = 'var'
    #     genome_name = genome_name.split('')[0] + ' Variants'
    # else:
    #     genome_name = genome_name.split('_')[0] + ' Reference'
    if 'True' in ref_comp:
        genome_type = 'both'
    else:
        genome_type = 'ref'
    mms = int(mms[0])

    # load acfd for each guide
    with open(current_working_directory + 'Results/' + value + '/'+value+'.acfd.txt') as a:
        all_scores = a.read().strip().split('\n')

    list_error_guides = []
    if os.path.exists(current_working_directory + 'Results/' + value + '/guides_error.txt'):
        with open(current_working_directory + 'Results/' + value + '/guides_error.txt') as error_g:
            for e_g in error_g:
                list_error_guides.append(e_g.strip())

    col_targetfor = '('
    for i in range(1, mms + int(max_bulges)):
        col_targetfor = col_targetfor + str(i) + '-'
    col_targetfor = col_targetfor + str(mms + int(max_bulges))
    col_targetfor = col_targetfor + ' Mismatches + Bulges)'

    columns_profile_table = [
        {'name': ['', 'gRNA (protospacer+PAM)'], 'id':'Guide', 'type':'text'},
        {'name': ['', 'Nuclease', ''], 'id':'Nuclease', 'type':'text'},
        {'name': ['', 'CFD (0-100)'], 'id':'CFD', 'type':'text'},
        # {'name': ['', 'Doench 2016'], 'id':'Doench 2016',
        #    'type':'text'},  # Doench, only for REF or VAR
        # {'name': ['', 'Doench 2016', ''], 'id':'Doench 2016',
        #     'type':'text'},  # REF Doench, only for Both
        # {'name': ['Doench 2016', 'Enriched'], 'id':'Enriched',
        #    'type':'text'},  # VAR Doench, only for Both
        {'name': [
            'Off-targets for Mismatch (MM) and Bulge (B) Value', 'Total'], 'id':'Total', 'type':'text'}
    ]  # Column of headers . Remove the entries accordingly when checking type of genome

    # for i in range (1, mms + int(max_bulges) + 1):
    #    columns_profile_table.append({'name':['Off-targets for Mismatch (MM) + Bulge (B) Value', str(i) + ' MM + B'], 'id': str(i) + ' MM + B', 'type':'text'})
    columns_profile_table.append(
        {'name': ['Off-targets for Mismatch (MM) and Bulge (B) Value', '# Bulges'], 'id': '# Bulges', 'type': 'text'})
    for i in range(mms + 1):
        columns_profile_table.append(
            {'name': ['Off-targets for Mismatch (MM) and Bulge (B) Value', str(i) + 'MM'], 'id': str(i) + 'MM', 'type': 'text'})

    remove_indices = set()

    if 'NO SCORES' in all_scores:
        # remove_indices.update([1,2,3,4])    #Remove CFD and Doench header
        remove_indices.update('CFD', 'Doench 2016', 'Reference', 'Enriched')

    if genome_type == 'ref':
        # remove_indices.update([3,4,6,7])
        remove_indices.update(
            ['Reference', 'Enriched'])
    # elif genome_type == 'both':
        # remove_indices.update([6])
        # remove_indices.update(['Doench 2016', 'On-Targets Enriched'])
    else:
        # remove_indices.update([3,4,5,7])
        remove_indices.update(
            ['Reference', 'Enriched', 'On-Targets Reference', 'Samples in Class 0 - 0+ - 1 - 1+'])

    # Remove headers not used in selected search result
    columns_profile_table = [i for j, i in enumerate(
        columns_profile_table) if columns_profile_table[j]['id'] not in remove_indices]

    final_list = []
    if list_error_guides:
        final_list.append(
            dbc.Alert(
                [
                    'Warning: Some guides have too many targets! ',
                    html.A("Click here", href=URL + "/data/" + job_id +
                           '/guides_error.txt', className="alert-link"),
                    ' to view them'
                ], color='warning')
        )
    final_list.append(
        html.H3('Result Summary - ' + genome_name + ' - ' + pam_name + ' - Mismatches ' +
                str(mms) + ' - DNA bulges ' + bulge_dna + ' - RNA bulges ' + bulge_rna)
    )

    add_to_description = html.P(
        'General summary for input guides. For each guide, is reported the count of targets in reference and variant genome grouped by mismatches count and bulge size.'
    )
    if genome_type == 'both':
        add_to_description = html.P(
            [
                'General summary for input guides. For each guide, is reported the count of targets in reference and variant genome grouped by mismatches count and bulge size.',
                # html.Span(
                #     "Samples for each Class is provided",
                #     id="tooltip-sample-class",
                #     style={"textDecoration": "underline", "cursor": "pointer"}
                # ),
                # ', along with the number of Off-Targets found for each Mismatch + Bulge value, for both Reference and Enriched Genomes.',
                # dbc.Tooltip(
                #     [
                #         html.Div([html.P([html.B('Class 0:'), ' Samples that does not have any On-Targets']),
                #                   html.P([html.B(
                #                       'Class 0+:'), ' Samples that have a subset of the Reference Genome On-Targets']),
                #                   html.P([html.B(
                #                       'Class 1:'), ' Samples that have the same On-Targets as the Reference Genome']),
                #                   html.P([html.B('Class 1+:'), ' Samples that creates at least a new On-Target, that is not present in the Reference Genome'])],
                #                  style={'display': 'inline-block'})
                #     ],
                #     target="tooltip-sample-class", style={'font-size': '12px'}
                # )
            ]
        )
    final_list.append(add_to_description)
    final_list.append(
        html.Div(
            dbc.Row(
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.P('Generating download link, Please wait...',
                                       id='download-link-general-table'),
                                dcc.Interval(interval=1*1000,
                                             id='interval-general-table'),
                                html.Div(current_working_directory + 'Results/' + job_id + '/' + job_id +
                                         '.general_table.txt', style={'display': 'none'}, id='div-info-general-table')
                            ]
                        ),
                        html.Div(
                            [
                                html.P('Generating download link, Please wait...',
                                       id='download-link-integrated-results'),
                                dcc.Interval(interval=1*1000,
                                             id='interval-integrated-results'),
                                html.Div(current_working_directory + 'Results/' + job_id + '/' + job_id +
                                         '.bestMerge.txt.integrated_results.zip', style={'display': 'none'}, id='div-info-integrated-results')
                            ]
                        )
                    ]
                )
            )
        )
    )
    final_list.append(
        html.Div(
            html.Div(
                dash_table.DataTable(
                    id='general-profile-table',
                    # page_size=PAGE_SIZE,
                    columns=columns_profile_table,
                    merge_duplicate_headers=True,
                    # fixed_rows={ 'headers': True, 'data': 0 },
                    # data = profile.to_dict('records'),
                    selected_cells=[{'row': 0, 'column': 0}],
                    css=[{'selector': '.row', 'rule': 'margin: 0', 'selector': 'td.cell--selected, td.focused', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}, {
                        'selector': 'td.cell--selected *, td.focused *', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}],
                    page_current=0,
                    page_size=10,
                    page_action='custom',
                    # virtualization = True,
                    filter_action='custom',
                    filter_query='',

                    sort_action='custom',
                    sort_mode='multi',
                    sort_by=[],
                    style_table={
                        # 'margin-left': "10%",
                        'max-height': '260px',
                        'overflowY': 'auto',
                        # 'overflowX': 'hidden',
                    },
                    style_data={
                        'whiteSpace': 'pre',
                        'height': 'auto',
                        'font-size': '1.30rem',
                    },
                    # style_cell={
                    #    'width':f'{1/len(columns_profile_table)*100}%'
                    # },
                    style_data_conditional=[
                        {
                            'if': {
                                'column_id': 'Genome'
                            },
                            'font-weight': 'bold',
                            'textAlign': 'center'

                        },
                        # {'if': {'column_id': 'Guide'},
                        #                    'width': '10%',
                        #                    }

                    ],
                    style_cell_conditional=[{'if': {'column_id': 'Guide'},
                                             'width': '20%',
                                             },
                                            {'if': {'column_id': 'Total'},
                                             'width': '15%',
                                             },
                                            {'if': {'column_id': 'Doench 2016'},
                                             'width': '5%',
                                             },
                                            {'if': {'column_id': '# Bulges'},
                                             'width': '5%',
                                             },
                                            {'if': {'column_id': 'Nuclease'},
                                             'width': '5%',
                                             }],
                    #                        {'if': {'column_id': 'Reference'},
                    #                        'width': '10%',
                    #                        }],
                ), id='div-general-profile-table', style={"margin-left": "5%", "margin-right": "5%"})
        )
    )

    final_list.append(html.Br())

    if genome_type == 'ref':
        final_list.append(
            dcc.Tabs(id="tabs-reports", value='tab-query-table', children=[
                dcc.Tab(label='Custom Ranking', value='tab-query-table'),
                dcc.Tab(label='Summary by Mismatches/Bulges',
                        value='tab-summary-by-guide'),
                dcc.Tab(label='Query Genomic Region',
                        value='tab-summary-by-position'),
                dcc.Tab(label='Graphical Reports',
                        value='tab-summary-graphical'),


            ])
        )
    else:
        # Barplot for population distributions
        final_list.append(
            html.Div(
                [
                    dbc.Row(
                        [
                            dbc.Col(html.Button(
                                "Show/Hide Target Distribution in SuperPopulations", id="btn-collapse-populations")),
                            # dbc.Col(html.A('Download full list of targets', target = '_blank', id = 'download-full-list' ))
                        ]
                    ),
                    dbc.Collapse(
                        dbc.Card(dbc.CardBody(
                            html.Div(id='content-collapse-population')
                        )),
                        id="collapse-populations",
                    ),
                ], hidden=True
            )
        )
        final_list.append(html.Br())
        final_list.append(
            dcc.Tabs(id="tabs-reports", value='tab-query-table', children=[
                dcc.Tab(label='Custom Ranking', value='tab-query-table'),
                dcc.Tab(label='Summary by Mismatches/Bulges',
                        value='tab-summary-by-guide'),
                dcc.Tab(label='Summary by Sample',
                        value='tab-summary-by-sample'),
                dcc.Tab(label='Query Genomic Region',
                        value='tab-summary-by-position'),
                dcc.Tab(label='Graphical Reports',
                        value='tab-summary-graphical'),
                dcc.Tab(label='Personal Risk Cards',
                        value='tab-graphical-sample-card'),
            ])
        )
    final_list.append(html.Div(id='div-tab-content'))

    final_list.append(html.Div(genome_type, style={
                      'display': 'none'}, id='div-genome-type'))
    result_page = html.Div(final_list, style={'margin': '1%'})
    return result_page

# Generate download link summary_by_sample


@app.callback(
    [Output('download-link-summary_by_sample', 'children'),
     Output('interval-summary_by_sample', 'disabled')],
    [Input('interval-summary_by_sample', 'n_intervals')],
    [State('div-info-summary_by_sample', 'children'),
     State('url', 'search')]
)
def downloadLinkSample(n, file_to_load, search):  # file to load =
    if n is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    # file_to_load = file_to_load + '.zip'
    file_to_load = file_to_load + '.txt'
    file_to_load = file_to_load.strip().split('/')[-1]
    # print(file_to_load)
    if os.path.exists(current_working_directory + 'Results/' + job_id + '/' + file_to_load):
        return html.A('Download file', href=URL+'/Results/' + job_id + '/' + file_to_load, target='_blank'), True

    return 'Generating download link, Please wait...', False


@app.callback(
    [Output('download-link-general-table', 'children'),
     Output('interval-general-table', 'disabled')],
    [Input('interval-general-table', 'n_intervals')],
    [State('div-info-general-table', 'children'),
     State('url', 'search')]
)
def downloadGeneralTable(n, file_to_load, search):  # file to load =
    if n is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    file_to_load = file_to_load.split('/')[-1]
    # print(file_to_load)
    if os.path.exists(current_working_directory + 'Results/' + job_id + '/' + file_to_load):
        return html.A('Download General Table', href=URL+'/Results/' + job_id + '/' + file_to_load, target='_blank'), True

    return 'Generating download link, Please wait...', False

# downalod integrated results


@app.callback(
    [Output('download-link-integrated-results', 'children'),
     Output('interval-integrated-results', 'disabled')],
    [Input('interval-integrated-results', 'n_intervals')],
    [State('div-info-integrated-results', 'children'),
     State('url', 'search')]
)
def downloadGeneralTable(n, file_to_load, search):  # file to load =
    if n is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    file_to_load = file_to_load.split('/')[-1]
    # print(file_to_load)
    if os.path.exists(current_working_directory + 'Results/' + job_id + '/' + file_to_load):
        return html.A('Download Integrated Results', href=URL+'/Results/' + job_id + '/' + file_to_load, target='_blank'), True

    return 'Generating download link, Please wait...', False

# Generate download link sumbysample


@app.callback(
    [Output('download-link-sumbysample', 'children'),
     Output('interval-sumbysample', 'disabled')],
    [Input('interval-sumbysample', 'n_intervals')],
    [State('div-info-sumbysample-targets', 'children'),
     State('url', 'search')]
)
def downloadLinkSample(n, file_to_load, search):  # file to load = job_id.HG001.guide
    if n is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    file_to_load = file_to_load + '.zip'
    if os.path.exists(current_working_directory + 'Results/' + job_id + '/' + file_to_load):
        return html.A('Download zip', href=URL+'/Results/' + job_id + '/' + file_to_load, target='_blank'), True

    return 'Generating download link, Please wait...', False

# Generate download link sumbyguide


@app.callback(
    [Output('download-link-sumbyguide', 'children'),
     Output('interval-sumbyguide', 'disabled')],
    [Input('interval-sumbyguide', 'n_intervals')],
    [State('div-info-sumbyguide-targets', 'children'),
     State('url', 'search')]
)
def downloadLinkGuide(n, file_to_load, search):  # file to load = job_id.RNA.1.0.guide
    if n is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    file_to_load = file_to_load + '.zip'
    if os.path.exists(current_working_directory + 'Results/' + job_id + '/' + file_to_load):
        return html.A('Download zip', href=URL+'/Results/' + job_id + '/' + file_to_load, target='_blank'), True

    return 'Generating download link, Please wait...', False


@app.server.route('/Results/<path:path>')
def download_file(path):
    # print(current_working_directory)
    # print('test', path)
    return flask.send_from_directory(os.path.join(current_working_directory, 'Results/'), path, as_attachment=True)


# Filter/sort IUPAC decomposition table for cluster page
@app.callback(
    Output('table-scomposition-cluster', 'data'),
    [Input('table-scomposition-cluster', "page_current"),
     Input('table-scomposition-cluster', "page_size"),
     Input('table-scomposition-cluster', 'sort_by'),
     Input('table-scomposition-cluster', 'filter_query')],
    [State('url', 'search'),
     State('url', 'hash')]
)
def update_iupac_scomposition_table_cluster(page_current, page_size, sort_by, filter, search, hash):
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    hash = hash.split('#')[1]
    guide = hash[:hash.find('-Pos-')]
    chr_pos = hash[hash.find('-Pos-') + 5:]
    chromosome = chr_pos.split('-')[0]
    position = chr_pos.split('-')[1]

    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    if genome_type == 'ref':
        raise PreventUpdate

    filtering_expressions = filter.split(' && ')
    dff = global_store_general(current_working_directory + 'Results/' + job_id + '/' +
                               job_id + '.' + chromosome + '_' + position + '.' + guide + '.scomposition.txt')
    if dff is None:
        raise PreventUpdate

    # #Grep annotation
    # file_to_grep = '.Annotation.targets.txt'
    # get_annotation = subprocess.Popen(['LC_ALL=C fgrep ' + guide + ' ' + current_working_directory + 'Results/'+ job_id + '/' + job_id + file_to_grep + ' |  awk \'$6==' + position + ' && $4==\"' + chromosome + '\"\''], shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # out, err = get_annotation.communicate()
    # annotation_type = out.decode('UTF-8').strip().split('\t')[-1]

    # if genome_type == 'var':
    #     dff.rename(columns = {0:'Bulge Type', 1:'crRNA', 2:'DNA', 3:'Chromosome', 4:'Position', 5:'Cluster Position', 6:'Direction',
    #     7:'Mismatches', 8:'Bulge Size', 9:'Total', 10:'Min Mismatches', 11:'Max Mismatches', 12:'Samples', 13:'Correct Guide', 14:'Annotation Type',15:'Top Subcluster'}, inplace = True)
    # else:
    dff.rename(columns=COL_BOTH_RENAME, inplace=True)
    # dff.drop(dff.columns[[-1,]], axis=1, inplace=True)         #NOTE Drop the Correct Guide column
    del dff['Correct Guide']
    # dff['Annotation Type'] = annotation_type
    del dff['Variant Unique']
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

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
            ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
                for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    # Calculate sample count

    data_to_send = dff.iloc[
        page_current*page_size:(page_current + 1)*page_size
    ].to_dict('records')
    if genome_type != 'ref':
        dict_sample_to_pop, dict_pop_to_superpop = associateSample.loadSampleAssociation(
            job_directory + 'sampleID.txt')[:2]
        for row in data_to_send:
            summarized_sample_cell = dict()
            for s in row['Samples'].split(','):
                if s == 'n':
                    break  # If a target have n, it means it's REF, because either all have samples or the single target is REF
                try:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] += 1
                except:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] = 1
            if summarized_sample_cell:
                row['Samples Summary'] = ', '.join(
                    [str(summarized_sample_cell[sp]) + ' ' + sp for sp in summarized_sample_cell])
            else:
                row['Samples Summary'] = 'n'
    return data_to_send


@app.callback(
    Output('table-position-target', 'data'),
    [Input('table-position-target', "page_current"),
     Input('table-position-target', "page_size"),
     Input('table-position-target', 'sort_by'),
     Input('table-position-target', 'filter_query'),
     Input('hide-reference-targets', 'value')],
    [State('url', 'search'),
     State('url', 'hash')]
)
def update_table_cluster(page_current, page_size, sort_by, filter, hide_reference, search, hash):
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    hash = hash.split('#')[1]
    guide = hash[:hash.find('-Pos-')]
    chr_pos = hash[hash.find('-Pos-') + 5:]
    chromosome = chr_pos.split('-')[0]
    position = chr_pos.split('-')[1]

    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    filtering_expressions = filter.split(' && ')
    dff = global_store_general(current_working_directory + 'Results/' + job_id +
                               '/' + job_id + '.' + chromosome + '_' + position + '.' + guide + '.txt')
    if dff is None:
        raise PreventUpdate

    if genome_type == 'ref':
        dff.rename(columns=COL_REF_RENAME, inplace=True)
    else:
        dff.rename(columns=COL_BOTH_RENAME, inplace=True)

    if genome_type != 'ref':
        # add_samples = [dff['Samples'][0]] * dff.shape[0]
        # check_minmms = dff['Min Mismatches']
        # if dff['Variant Unique'][0] != 'y' and dff['PAM Creation'][0] == 'n':
        #     for pos_minmms, minmms in enumerate(check_minmms):
        #         if minmms == '-':
        #             add_samples[pos_minmms] = 'n'
        # dff['Samples'] = add_samples
        del dff['Variant Unique']
        # dff.drop(dff.head(1).index, inplace=True)       #Remove first target, that is the top1 with no iupac (lowest mm of scomposed target) and is
        # needed only for summary by guide, not the show target part
        # NOTE 18/03 removed all the iupac char, the first line is needed to be shown
    # dff['Annotation Type'] = list(dff['Annotation Type'])[0]
    del dff['Correct Guide']

    if 'hide-ref' in hide_reference or genome_type == 'var':
        dff.drop(dff[(dff['Samples'] == 'n')].index, inplace=True)

    if 'hide-cluster' in hide_reference:
        dff = dff.head(1)

    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

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
            ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
                for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    # Calculate sample count

    data_to_send = dff.iloc[
        page_current*page_size:(page_current + 1)*page_size
    ].to_dict('records')
    if genome_type != 'ref':
        dict_sample_to_pop, dict_pop_to_superpop = associateSample.loadSampleAssociation(
            job_directory + 'sampleID.txt')[:2]
        for row in data_to_send:
            summarized_sample_cell = dict()
            for s in row['Samples'].split(','):
                if s == 'n':
                    break
                try:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] += 1
                except:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] = 1
            if summarized_sample_cell:
                row['Samples Summary'] = ', '.join(
                    [str(summarized_sample_cell[sp]) + ' ' + sp for sp in summarized_sample_cell])
            else:
                row['Samples Summary'] = 'n'
    return data_to_send

# Return the targets for the selected cluster


def clusterPage(job_id, hash):
    guide = hash[:hash.find('-Pos-')]
    chr_pos = hash[hash.find('-Pos-') + 5:]
    chromosome = chr_pos.split('-')[0]
    position = chr_pos.split('-')[1]
    if (not isdir(current_working_directory + 'Results/' + job_id)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    style_hide_reference = {'display': 'none'}
    value_hide_reference = []
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'
        style_hide_reference = {}
        value_hide_reference = ['hide-ref', 'hide-cluster']
    final_list = []
    final_list.append(
        html.H3('Selected Position: ' + chromosome + ' - ' + position)
    )

    if genome_type == 'ref':
        cols = [{"name": i, "id": i, 'type': t, 'hideable': True}
                for i, t in zip(COL_REF, COL_REF_TYPE)]
        file_to_grep = '.bestMerge.txt'
    else:
        cols = [{"name": i, "id": i, 'type': t, 'hideable': True}
                for i, t in zip(COL_BOTH, COL_BOTH_TYPE)]
        file_to_grep = '.bestMerge.txt'
    # print('qui cluster before grep')

    cluster_grep_result = current_working_directory + 'Results/' + job_id + \
        '/' + job_id + '.' + chromosome + '_' + position + '.' + guide + '.txt'
    put_header = 'head -1 ' + current_working_directory + 'Results/' + job_id + \
        '/' + job_id + file_to_grep + ' > ' + cluster_grep_result + ' ; '
    # print('esiste cluster?' , str(os.path.exists(cluster_grep_result)) )
    # Example    job_id.chr3_100.guide.txt
    if not os.path.exists(cluster_grep_result):
        # Grep annotation for ref
        if genome_type == 'ref':  # NOTE HEADER NON SALVATO
            get_annotation = subprocess.Popen(['LC_ALL=C fgrep ' + guide + ' ' + current_working_directory + 'Results/' + job_id + '/' + job_id + '.Annotation.targets.txt' +
                                               ' |  awk \'$6==' + position + ' && $4==\"' + chromosome + '\"\''], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = get_annotation.communicate()
            annotation_type = out.decode('UTF-8').strip().split('\t')[-1]
            os.popen(put_header + 'LC_ALL=C fgrep ' + guide + ' ' + current_working_directory + 'Results/' + job_id + '/' + job_id + file_to_grep +
                     ' | awk \'$6==' + position + ' && $4==\"' + chromosome + '\" {print $0\"\\t' + annotation_type + '\"}\' > ' + cluster_grep_result).read()
        else:
            # print('qui cluster in grep')     #NOTE HEADER NON SALVATO
            os.popen(put_header + 'LC_ALL=C fgrep ' + guide + ' ' + current_working_directory + 'Results/' + job_id + '/' + job_id + file_to_grep + ' | awk \'$6==' + position + ' && $4==\"' + chromosome +
                     '\"\' > ' + cluster_grep_result).read()  # NOTE top1 will have sample and annotation, other targets will have '.'-> 18/03 all samples and annotation are already writter for all targets
        os.system('zip '+'-j ' + cluster_grep_result.replace('.txt',
                                                             '.zip') + ' ' + cluster_grep_result+" &")
    final_list.append(
        html.Div(job_id + '.' + chromosome + '_' + position + '.' + guide,
                 style={'display': 'none'}, id='div-info-sumbyposition-targets')
    )

    scomposition_file = current_working_directory + 'Results/' + job_id + '/' + \
        job_id + '.' + chromosome + '_' + position + '.' + guide + '.scomposition.txt'
    file_to_grep = '.bestMerge.txt'

    iupac_scomposition_visibility = {'display': 'none'}
    if genome_type != 'ref':
        iupac_scomposition_visibility = {}
        # Example    job_id.chr_pos.guide.scomposition.txt
        if not os.path.exists(scomposition_file):
            os.popen('LC_ALL=C fgrep ' + guide + ' ' + current_working_directory + 'Results/' + job_id + '/' + job_id + file_to_grep +
                     ' |  awk \'$6==' + position + ' && $4==\"' + chromosome + '\" && $13!=\"n\"\' > ' + scomposition_file).read()

    final_list.append(html.P(
        [
            html.P('List of all the configurations for the target in the selected position.',
                   style=iupac_scomposition_visibility),
            dcc.Checklist(
                options=[{'label': 'Hide Reference Targets', 'value': 'hide-ref'},
                         {'label': 'Show only TOP1 Target', 'value': 'hide-cluster'}],
                id='hide-reference-targets', value=value_hide_reference, style=style_hide_reference
            ),
            html.Div(
                [
                    html.P('Generating download link, Please wait...',
                           id='download-link-sumbyposition'),
                    dcc.Interval(interval=5*1000, id='interval-sumbyposition')
                ]

            )
        ]
    )
    )

    cols_for_scomposition = cols.copy()
    cols_for_scomposition.append(
        {"name": 'Samples', "id": 'Samples', 'type': 'text', 'hideable': True})
    final_list.append(
        html.Div(
            dash_table.DataTable(
                # TABLE that represent scomposition of iupac of selected target, take rows from top_1.samples.txt
                id='table-scomposition-cluster',
                columns=cols_for_scomposition,
                # data = df.to_dict('records'),
                virtualization=True,
                fixed_rows={'headers': True, 'data': 0},
                # fixed_columns = {'headers': True, 'data':1},
                style_cell={'width': '150px'},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                filter_action='custom',
                filter_query='',
                style_table={
                    'max-height': '600px'
                },
                style_data_conditional=[
                    {
                        'if': {
                            'filter_query': '{Variant Unique} eq F',
                        },
                        'background-color': 'rgba(0, 0, 0,0.15)'

                    }
                ],
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}, {'selector': 'td.cell--selected, td.focused', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}, {
                    'selector': 'td.cell--selected *, td.focused *', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}],

            ),
            style=iupac_scomposition_visibility
        )
    )

    final_list.append(html.Hr())

    # Cluster Table
    final_list.append(
        # The rows highlighted in red indicates that the target was found only in the genome with variants.',
        'List of Targets found for the selected position. Other possible configurations of the target are listed in the table above, along with the corresponding samples list.',
    )
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id='table-position-target',
                columns=cols,
                # data = df.to_dict('records'),
                virtualization=True,
                fixed_rows={'headers': True, 'data': 0},
                # fixed_columns = {'headers': True, 'data':1},
                style_cell={'width': '150px'},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                filter_action='custom',
                filter_query='',
                style_table={
                    'max-height': '600px'
                    # 'overflowY': 'scroll',
                },
                style_data_conditional=[
                    {
                        'if': {
                            'filter_query': '{Variant Unique} eq F',
                        },
                        'background-color': 'rgba(0, 0, 0,0.15)'

                    }
                ],
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}, {'selector': 'td.cell--selected, td.focused', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}, {
                    'selector': 'td.cell--selected *, td.focused *', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}],

            ),
            id='div-result-table',
        )
    )
    return html.Div(final_list, style={'margin': '1%'})

# Filter and sorting sample targets


@app.callback(
    Output('table-sample-target', 'data'),
    [Input('table-sample-target', "page_current"),
     Input('table-sample-target', "page_size"),
     Input('table-sample-target', 'sort_by'),
     Input('table-sample-target', 'filter_query')],
    [State('url', 'search'),
     State('url', 'hash')]
)
def update_table_sample(page_current, page_size, sort_by, filter, search, hash):
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    hash = hash.split('#')[1]
    guide = hash[:hash.find('-Sample-')]
    sample = hash[hash.rfind('-') + 1:]
    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    filtering_expressions = filter.split(' && ')
    dff = global_store_general(current_working_directory + 'Results/' +
                               job_id + '/' + job_id + '.' + str(sample) + '.' + guide + '.txt')
    if dff is None:
        raise PreventUpdate

    dff.rename(columns=COL_BOTH_RENAME, inplace=True)
    # del dff['Real_Guide']  # NOTE Drop the Correct Guide column
    # del dff['Variant Unique']
    dff = dff[COL_BOTH]
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == 'contains':
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == 'datestartswith':
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    # if len(sort_by):
    #     dff = dff.sort_values(
    #         ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
    #             for col in sort_by],
    #         ascending=[
    #             col['direction'] == 'asc'
    #             for col in sort_by
    #         ],
    #         inplace=False
    #     )

    # Calculate sample count
    # dict_sample_to_pop, dict_pop_to_superpop = associateSample.loadSampleAssociation(
    #     job_directory + 'sampleID.txt')[:2]
    data_to_send = dff.iloc[page_current *
                            page_size:(page_current + 1)*page_size].to_dict('records')
    # for row in data_to_send:
    #     summarized_sample_cell = dict()
    #     for s in row['Samples'].split(','):
    #         if s == 'n':
    #             break
    #         try:
    #             summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] += 1
    #         except:
    #             summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] = 1
    #     if summarized_sample_cell:
    #         row['Samples Summary'] = ', '.join(
    #             [str(summarized_sample_cell[sp]) + ' ' + sp for sp in summarized_sample_cell])
    #     else:
    #         row['Samples Summary'] = 'n'
    return data_to_send
# Return the targets found for the selected sample


def samplePage(job_id, hash):
    # print("SAMPLE PAGE LOADED FOR", job_id, hash)
    guide = hash[:hash.find('-Sample-')]
    sample = hash[hash.rfind('-') + 1:]
    if (not isdir(current_working_directory + 'Results/' + job_id)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))

    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    final_list = []
    final_list.append(
        # html.P('List of Targets found for the selected Sample - ' + sample + ' - and guide - ' + guide + ' -')
        html.H3('Selected Sample: ' + sample)
    )
    final_list.append(
        html.P(
            [
                # 'The rows highlighted in red indicates that the target was found only in the genome with variants.',
                'List of Targets found for the selected sample.',
                html.Div(
                    [
                        html.P('Generating download link, Please wait...',
                               id='download-link-sumbysample'),
                        dcc.Interval(interval=5*1000,
                                     id='interval-sumbysample')
                    ]
                )
            ]
        )
    )

    file_to_grep = current_working_directory + 'Results/' + \
        job_id + '/' + job_id + '.bestMerge.txt'
    sample_grep_result = current_working_directory + 'Results/' + \
        job_id + '/' + job_id + '.' + str(sample) + '.' + guide + '.txt'
    final_list.append(
        html.Div(job_id + '.' + str(sample) + '.' + guide,
                 style={'display': 'none'}, id='div-info-sumbysample-targets')
    )
    os.system('LC_ALL=C fgrep ' + guide + ' ' + file_to_grep +
              ' | awk \'$14~\"' + str(sample) + '\"\' > ' + sample_grep_result)

    os.system('zip '+'-j ' + sample_grep_result.replace('.txt',
                                                        '.zip') + ' ' + sample_grep_result + " &")

    cols = [{"name": i, "id": i, 'type': t, 'hideable': True}
            for i, t in zip(COL_BOTH, COL_BOTH_TYPE)]

    final_list.append(
        html.Div(
            dash_table.DataTable(
                id='table-sample-target',
                columns=cols,
                # data = df.to_dict('records'),
                virtualization=True,
                fixed_rows={'headers': True, 'data': 0},
                # fixed_columns = {'headers': True, 'data':1},
                style_cell={'width': '150px'},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                filter_action='custom',
                filter_query='',
                style_table={
                    'max-height': '600px'
                    # 'overflowY': 'scroll',
                },
                style_data_conditional=[
                    {
                        'if': {
                            'filter_query': '{Variant Unique} eq F',
                        },
                        'background-color': 'rgba(0, 0, 0,0.15)'

                    }
                ],
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}, {'selector': 'td.cell--selected, td.focused', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}, {
                    'selector': 'td.cell--selected *, td.focused *', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}],
            ),
            id='div-result-table',
        )
    )
    return html.Div(final_list, style={'margin': '1%'})


@cache.memoize()
def global_store_general(path_file_to_load):
    '''
    Caching dei file targets per una miglior performance di visualizzazione
    '''
    if 'scomposition' in path_file_to_load:
        rows_to_skip = 0
    else:
        rows_to_skip = 0  # Skip header
    if path_file_to_load is None:
        return ''
    if os.path.getsize(path_file_to_load) > 0:
        # df = pd.read_csv(path_file_to_load, sep='\t', header=None, skiprows=rows_to_skip, usecols=range(0, 16))
        df = pd.read_csv(path_file_to_load, sep='\t', header=None,
                         skiprows=rows_to_skip, usecols=range(0, 24))
    else:
        df = None
    return df
# Filter etc for second table


@app.callback(
    [Output('second-table-subset-targets', 'data'),  # Table showing iupac scomposition
     Output('second-table-subset-targets', 'style_data_conditional')],
    [Input('second-table-subset-targets', "page_current"),
     Input('second-table-subset-targets', "page_size"),
     Input('second-table-subset-targets', "sort_by"),
     Input('second-table-subset-targets', 'filter_query')],
    [State('url', 'search'),
     State('url', 'hash'),
     State('table-subset-target', 'active_cell'),
     State('table-subset-target', 'data')]
)
def update_table_subsetSecondTable(page_current, page_size, sort_by, filter, search, hash_guide, active_cel, data):
    # NOTE tabella secondaria della scomposizione ora non serve, non cancello il codice ma uso PreventUpdate per non azionare la funzione
    if False:
        raise PreventUpdate
    if active_cel is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    guide = hash_guide[1:hash_guide.find('new')]
    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'
    if search is None:
        raise PreventUpdate

    if genome_type == 'ref':
        raise PreventUpdate

    filtering_expressions = filter.split(' && ')
    bulge_t = data[active_cel['row']]['Bulge Type']
    bulge_s = str(data[active_cel['row']]['Bulge Size'])
    mms = str(data[active_cel['row']]['Mismatches'])
    chrom = str(data[active_cel['row']]['Chromosome'])
    pos = str(data[active_cel['row']]['Cluster Position'])
    # annotation_type = str(data[active_cel['row']]['Annotation Type'])

    scomposition_file = job_directory + job_id + '.' + bulge_t + '.' + bulge_s + \
        '.' + mms + '.' + guide + '.' + chrom + '.' + pos + '.scomposition.txt'
    file_to_grep = job_directory + job_id + '.' + bulge_t + \
        '.' + bulge_s + '.' + mms + '.' + guide + '.txt'
    # file_to_grep_alt = '.altMerge.txt'
    os.system(f'LC_ALL=C fgrep {pos} {file_to_grep} | LC_ALL=C fgrep {chrom} | awk \'$14!=\"n\"' +
              '\' > ' + scomposition_file)  # , shell = True)
    # subprocess.call(['LC_ALL=C fgrep ' + guide + ' ' + current_working_directory + 'Results/'+ job_id + '/' + job_id + file_to_grep_alt + ' |  awk \'$7==' + pos + ' && $5==\"' + chrom + '\" && $10==' + bulge_s + ' && $14!=\"n\"' +'\' >> ' + scomposition_file], shell = True)
    # Check if result grep has at least 1 result
    if os.path.getsize(scomposition_file) > 0:
        df = pd.read_csv(scomposition_file, header=None, sep='\t')
        # df['Annotation Type'] = annotation_type
    else:
        raise PreventUpdate

    df.rename(columns=COL_BOTH_RENAME, inplace=True)
    df.drop(df[(~(df['Cluster Position'] == int(data[active_cel['row']]['Cluster Position']))) | (
        ~(df['Chromosome'] == data[active_cel['row']]['Chromosome']))].index, inplace=True)
    dff = df

    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

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
            ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
                for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    cells_style = [
        {
            'if': {
                'filter_query': '{Variant Unique} eq F',
                                # 'filter_query': '{Direction} eq +',
                                # 'column_id' :'Bulge Type'
            },
            # 'border-left': '5px solid rgba(255, 26, 26, 0.9)',
            # 'rgb(255, 102, 102)'
            'background-color': 'rgba(0, 0, 0,0.15)'
        }
    ]

    # Calculate sample count

    data_to_send = dff.iloc[
        page_current*page_size:(page_current + 1)*page_size
    ].to_dict('records')
    if genome_type != 'ref':
        dict_sample_to_pop, dict_pop_to_superpop = associateSample.loadSampleAssociation(
            job_directory + 'sampleID.txt')[:2]
        for row in data_to_send:
            summarized_sample_cell = dict()
            for s in row['Samples'].split(','):
                if s == 'n':
                    break
                try:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] += 1
                except:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] = 1
            if summarized_sample_cell:
                row['Samples Summary'] = ', '.join(
                    [str(summarized_sample_cell[sp]) + ' ' + sp for sp in summarized_sample_cell])
            else:
                row['Samples Summary'] = 'n'
    return data_to_send, cells_style
# Create second table for subset targets page, and show corresponding samples    -> CHANGED, now show IUPAC scomposition


@app.callback(
    [Output('div-second-table-subset-targets', 'children'),
     Output('table-subset-target', 'style_data_conditional')],
    [Input('table-subset-target', 'active_cell')],
    [State('table-subset-target', 'data'),
     State('table-subset-target', 'columns'),
     State('url', 'search'),
     State('table-subset-target', 'style_data_conditional'),
     State('table-subset-target', 'selected_cells')]
)
def loadFullSubsetTable(active_cel, data, cols, search, style_data, sel_cell):
    # NOTE tabella secondaria della scomposizione ora non serve, non cancello il codice ma uso PreventUpdate per non azionare la funzione
    if False:
        raise PreventUpdate
    if active_cel is None:
        raise PreventUpdate
    fl = []
    job_id = search.split('=')[-1]
    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    if genome_type == 'ref':
        raise PreventUpdate
    fl.append(html.Hr())
    fl.append('List of all the configurations for the selected target.')
    fl.append(html.Br())
    cols.append({"name": 'Samples', "id": 'Samples',
                 'type': 'text', 'hideable': True})

    fl.append(
        html.Div(
            dash_table.DataTable(
                id='second-table-subset-targets',
                columns=cols,
                virtualization=True,
                fixed_rows={'headers': True, 'data': 0},
                style_cell={'width': '150px'},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                filter_action='custom',
                filter_query='',
                style_table={
                    'max-height': '600px'
                },
                style_cell_conditional=[
                    {
                        'if': {'column_id': 'Samples'},
                        'textAlign': 'left'
                    }
                ],
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}, {'selector': 'td.cell--selected, td.focused', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}, {
                    'selector': 'td.cell--selected *, td.focused *', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}],
                style_data_conditional=[
                    {
                        'if': {
                            'filter_query': '{Variant Unique} eq F',
                        },
                        'background-color': 'rgba(0, 0, 0,0.15)'
                    }
                ]
            )
        )
    )

    pos_cluster = data[int(sel_cell[0]['row'])]['Cluster Position']
    chrom = data[int(sel_cell[0]['row'])]['Chromosome']
    cells_style = [
        style_data[0],
        {
            'if': {
                'filter_query': '{Cluster Position} eq "' + str(pos_cluster) + '"',
            },
            'background-color': 'rgba(0, 0, 255,0.15)'  # 'rgb(255, 102, 102)'

        }
    ]
    return fl, cells_style

# Update primary table of 'Show targets' of Summary by Guide


@app.callback(
    Output('table-subset-target', 'data'),
    [Input('table-subset-target', "page_current"),
     Input('table-subset-target', "page_size"),
     Input('table-subset-target', "sort_by"),
     Input('table-subset-target', 'filter_query'),
     Input('hide-reference-targets', 'value')],
    [State('url', 'search'),
     State('url', 'hash')]
)
def update_table_subset(page_current, page_size, sort_by, filter, hide_reference, search, hash_guide):
    '''
    La funzione ritorna uno split dei risultati in base ad un filtering o a un sort da parte dell'utente. Inoltre aggiorna i risultati
    visualizzati quando il bottone next page / prev page è cliccato. (Codice preso dalla pagina dash datatable sul sorting con python)
    Inoltre carica i file targets, o scores se presente, e lo trasforma in un dataframe, cambiando il nome delle colonne per farle corrispondere
    all'id delle colonne della tabella nella pagina.
    Se non ci sono targets ritorna un avviso di errore
    '''
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'
    value = job_id
    if search is None:
        raise PreventUpdate
    filtering_expressions = filter.split(' && ')
    # filtering_expressions.append(['{crRNA} = ' + guide])
    guide = hash_guide[1:hash_guide.find('new')]
    mms = hash_guide[-1:]
    bulge_s = hash_guide[-2:-1]
    if 'DNA' in hash_guide:
        bulge_t = 'DNA'
    elif 'RNA' in hash_guide:
        bulge_t = 'RNA'
    else:
        bulge_t = 'X'
    df = global_store_subset(value, bulge_t, bulge_s, mms, guide)
    dff = df

    # print(dff, 'line 1392')
    # if genome_type == 'ref':
    #    dff.rename(columns = COL_REF_RENAME, inplace = True)
    # else:
    dff.rename(columns=COL_BOTH_RENAME, inplace=True)
    # print(dff, 'line 1397')

    if 'hide-ref' in hide_reference or genome_type == 'var':
        dff.drop(df[(df['Samples'] == 'n')].index, inplace=True)

    try:  # For VAR and BOTH
        del dff['Variant Unique']
    except:  # For REF
        pass

    # print(dff)

    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)
        if col_name == 'Samples Summary':
            col_name = 'Samples'
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
            ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
                for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    cells_style = [
        {
            'if': {
                'filter_query': '{Cluster Position} eq "' + guide + '"',
                                # 'column_id' :'{#Bulge type}',
                                # 'column_id' :'{Total}'
            },
            # 'border-left': '5px solid rgba(255, 26, 26, 0.9)',
            # 'rgb(255, 102, 102)'
            'background-color': 'rgba(0, 0, 255,0.15)'

        },
    ]

    # Calculate sample count

    data_to_send = dff.iloc[
        page_current*page_size:(page_current + 1)*page_size
    ].to_dict('records')
    if genome_type != 'ref':
        dict_sample_to_pop, dict_pop_to_superpop = associateSample.loadSampleAssociation(
            job_directory + 'sampleID.txt')[:2]
        for row in data_to_send:
            summarized_sample_cell = dict()
            for s in row['Samples'].split(','):
                if s == 'n':
                    break
                try:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] += 1
                except:
                    summarized_sample_cell[dict_pop_to_superpop[dict_sample_to_pop[s]]] = 1
            if summarized_sample_cell:
                row['Samples Summary'] = ', '.join(
                    [str(summarized_sample_cell[sp]) + ' ' + sp for sp in summarized_sample_cell])
            else:
                row['Samples Summary'] = 'n'
    return data_to_send  # , cells_style + style_data_table


def guidePagev3(job_id, hash):
    guide = hash[:hash.find('new')]
    mms = hash[-1:]
    bulge_s = hash[-2:-1]
    if 'DNA' in hash:
        bulge_t = 'DNA'
    elif 'RNA' in hash:
        bulge_t = 'RNA'
    else:
        bulge_t = 'X'
    add_header = ' - Mismatches ' + str(mms)
    if bulge_t != 'X':
        add_header += ' - ' + str(bulge_t) + ' ' + str(bulge_s)
    value = job_id
    if (not isdir(current_working_directory + 'Results/' + job_id)):
        return html.Div(dbc.Alert("The selected result does not exist", color="danger"))
    with open(current_working_directory + 'Results/' + value + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    job_directory = current_working_directory + 'Results/' + job_id + '/'
    genome_type = 'ref'
    style_hide_reference = {'display': 'none'}
    value_hide_reference = []
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'
        style_hide_reference = {}
        value_hide_reference = ['hide-ref']
    final_list = []
    final_list.append(html.H3('Selected Guide: ' + guide + add_header))
    final_list.append(
        html.P(
            [
                # 'Select a row to view the target IUPAC character scomposition. The rows highlighted in red indicates that the target was found only in the genome with variants.',
                'List of Targets found for the selected guide. Select a row to view other possible configurations of the target, along with the corresponding samples list.',
                dcc.Checklist(options=[{'label': 'Hide Reference Targets', 'value': 'hide-ref'}],
                              id='hide-reference-targets', value=value_hide_reference, style=style_hide_reference),
                html.Div(
                    [
                        html.P('Generating download link, Please wait...',
                               id='download-link-sumbyguide'),
                        dcc.Interval(interval=5*1000, id='interval-sumbyguide')
                    ]
                )
            ]
        )
    )
    cols = [{"name": i, "id": i, 'type': t, 'hideable': True}
            for i, t in zip(COL_BOTH, COL_BOTH_TYPE)]
    file_to_grep = job_directory + job_id + '.bestMerge.txt'
    # file_to_grep_alt = job_directory + job_id + '.altMerge.txt'

    guide_grep_result = job_directory + job_id + '.' + \
        bulge_t + '.' + bulge_s + '.' + mms + '.' + guide + '.txt'
    # put_header = 'head -1 ' + job_directory + job_id + file_to_grep + ' > ' + guide_grep_result + ' ; '
    final_list.append(
        html.Div(job_id + '.' + bulge_t + '.' + bulge_s + '.' + mms + '.' +
                 guide, style={'display': 'none'}, id='div-info-sumbyguide-targets')
    )

    os.system('LC_ALL=C fgrep ' + guide + ' ' + file_to_grep + ' | LC_ALL=C fgrep ' +
              bulge_t + ' | awk \'$9==' + mms + ' && $10==' + bulge_s + '\'> ' + guide_grep_result)
    os.system('zip '+'-j ' + guide_grep_result.replace('.txt', '.zip') +
              ' ' + guide_grep_result + " &")  # , shell = True)
    global_store_subset(job_id, bulge_t, bulge_s, mms, guide)

    final_list.append(
        html.Div(
            dash_table.DataTable(
                id='table-subset-target',
                columns=cols,
                # data = subset_targets.to_dict('records'),
                virtualization=True,
                fixed_rows={'headers': True, 'data': 0},
                # fixed_columns = {'headers': True, 'data':1},
                style_cell={'width': '150px'},
                page_current=0,
                page_size=PAGE_SIZE,
                page_action='custom',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[],
                filter_action='custom',
                filter_query='',
                style_table={
                    'max-height': '600px'
                    # 'overflowY': 'scroll',
                },
                style_cell_conditional=[
                    {
                        'if': {'column_id': 'Samples'},
                        'textAlign': 'left'
                    }
                ],
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}, {'selector': 'td.cell--selected, td.focused', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}, {
                    'selector': 'td.cell--selected *, td.focused *', 'rule': 'background-color: rgba(0, 0, 255,0.15) !important;'}],
                style_data_conditional=[
                    {
                        'if': {
                            'filter_query': '{Variant Unique} eq F',
                        },
                        'background-color': 'rgba(0, 0, 0,0.15)'

                    },

                ]
            ),
            id='div-result-table',
        )
    )
    final_list.append(html.Br())
    final_list.append(
        html.Div(
            id='div-second-table-subset-targets'
        )
    )

    return html.Div(final_list, style={'margin': '1%'})


@cache.memoize()
def global_store_subset(value, bulge_t, bulge_s, mms, guide):
    '''
    Caching dei file targets per una miglior performance di visualizzazione
    '''
    if value is None:
        return ''
    # Skiprows = 1 to skip header of file
    df = pd.read_csv(current_working_directory + 'Results/' + value + '/' + value + '.' + bulge_t + '.' +
                     bulge_s + '.' + mms + '.' + guide + '.txt', sep='\t', header=None, usecols=range(0, 25))  # , skiprows = 1)
    return df

# Load barplot of population distribution for selected guide


@app.callback(
    Output('content-collapse-population', 'children'),
    [Input('general-profile-table', 'selected_cells')],
    [State('general-profile-table', 'data'),
     State('url', 'search')]
)
def loadDistributionPopulations(sel_cel, all_guides, job_id):
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    guide = all_guides[int(sel_cel[0]['row'])]['Guide']
    job_id = job_id.split('=')[-1]

    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        mms = int((next(s for s in all_params.split('\n')
                        if 'Mismatches' in s)).split('\t')[-1])
        max_bulges = int((next(s for s in all_params.split(
            '\n') if 'Max_bulges' in s)).split('\t')[-1])

    distributions = [dbc.Row(html.P(
        'On- and Off-Targets distributions in the Reference and Variant Genome. For the Variant Genome, the targets are divided into SuperPopulations.', style={'margin-left': '0.75rem'}))]

    for i in range(math.ceil((mms + max_bulges + 1) / BARPLOT_LEN)):
        all_images = []
        for mm in range(i * BARPLOT_LEN, (i + 1) * BARPLOT_LEN):
            if mm < (mms + max_bulges + 1):
                try:
                    all_images.append(
                        dbc.Col(
                            [
                                html.A(
                                    html.Img(
                                        src='data:image/png;base64,{}'.format(base64.b64encode(open(
                                            current_working_directory + 'Results/' + job_id + '/imgs/populations_distribution_' + guide + '_' + str(mm) + 'total.png', 'rb').read()).decode()),
                                        id='distribution-population' + str(mm), width="100%", height="auto"
                                    ),
                                    target="_blank",
                                    href='/Results/' + job_id + '/imgs/' + 'populations_distribution_' +
                                    guide + '_' + str(mm) + 'total.png'
                                ),
                                html.Div(html.P('Distribution ' + str(mm) + ' Mismatches + Bulges ', style={
                                         'display': 'inline-block'}), style={'text-align': 'center'})
                            ]
                        )
                    )
                except:
                    all_images.append(
                        dbc.Col(
                            [
                                html.Div(html.P('No Targets found with ' + str(mm) + ' Mismatches + Bulges', style={
                                         'display': 'inline-block'}), style={'text-align': 'center'}),
                                # html.Div(html.P('Distribution ' + str(mm) + ' Mismatches + Bulges ', style = {'display':'inline-block'} ),style = {'text-align':'center'})
                            ],
                            align='center'
                        )
                    )
            else:
                all_images.append(dbc.Col(html.P('')))

        distributions.append(
            html.Div(
                [
                    dbc.Row(
                        all_images
                    )
                ]
            )
        )
    return distributions


# Open/close barplot for population distribution
@app.callback(
    Output("collapse-populations", "is_open"),
    [Input("btn-collapse-populations", "n_clicks")],
    [State("collapse-populations", "is_open")],
)
def toggleCollapseDistributionPopulations(n, is_open):
    if n:
        return not is_open
    return is_open

# Filtering e sorting per la pagina principale delle guide


@app.callback(
    [Output('general-profile-table', 'data'),
     Output('general-profile-table', 'selected_cells')],
    [Input('general-profile-table', "page_current"),
     Input('general-profile-table', "page_size"),
     Input('general-profile-table', 'sort_by'),
     Input('general-profile-table', 'filter_query')],
    [State('url', 'search')]
)
def update_table_general_profile(page_current, page_size, sort_by, filter, search):
    job_id = search.split('=')[-1]

    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]
        mms = int((next(s for s in all_params.split('\n')
                        if 'Mismatches' in s)).split('\t')[-1])
        max_bulges = int((next(s for s in all_params.split(
            '\n') if 'Max_bulges' in s)).split('\t')[-1])
        nuclease = (next(s for s in all_params.split(
            '\n') if 'Nuclease' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    filtering_expressions = filter.split(' && ')

    # Get error guides
    list_error_guides = []
    if os.path.exists(current_working_directory + 'Results/' + job_id + '/guides_error.txt'):
        with open(current_working_directory + 'Results/' + job_id + '/guides_error.txt') as error_g:
            for e_g in error_g:
                list_error_guides.append(e_g.strip())

    # Get guide from guide.txt
    with open(current_working_directory + 'Results/' + job_id + '/guides.txt') as g:
        guides = g.read().strip().split('\n')
        guides.sort()

    # load acfd for each guide
    with open(current_working_directory + 'Results/' + job_id + '/'+job_id+'.'+'acfd.txt') as a:
        all_scores = a.read().strip().split('\n')

    # Load scores
    if 'NO SCORES' not in all_scores:
        all_scores.sort()
        acfd = [float(a.split('\t')[1]) for a in all_scores if a.split(
            '\t')[0] not in list_error_guides]
        doench = [a.split('\t')[2] for a in all_scores if a.split(
            '\t')[0] not in list_error_guides]
        if genome_type == 'both':
            doench_enr = [a.split('\t')[3] for a in all_scores if a.split('\t')[
                0] not in list_error_guides]
        # acfd = [int(round((100/(100 + x))*100)) for x in acfd]
        acfd = [float("{:.3f}".format(x*100)) for x in acfd]

    # Get target counting from summary by guide
    column_on_target = []
    column_off_target_ref = []
    column_sample_class = []
    column_total = []

    df = []
    table_to_file = list()
    for x, g in enumerate(guides):
        table_to_file.append(g)  # append guide to table
        # append nuclease to table
        table_to_file.append('Nuclease: '+str(nuclease))
        data_general_count = pd.read_csv(current_working_directory + 'Results/' +
                                         job_id + '/' + job_id + '.general_target_count.'+g+".txt", sep='\t')

        data_guides = dict()
        data_guides['Guide'] = g
        data_guides['Nuclease'] = nuclease
        data_general_count_copy = data_general_count.copy()
        count_bulges = list()
        origin_ref = list()
        origin_var = list()
        for the_bulge in range(max_bulges+1):
            origin_ref.append('REF')
            origin_var.append('VAR')
            count_bulges.append(the_bulge)

        count_bulges_concat = count_bulges+count_bulges
        origin_concat = origin_ref+origin_var

        data_general_count_copy.insert(0, 'Genome', origin_concat, True)
        data_general_count_copy.insert(1, 'Bulges', count_bulges_concat, True)

        if 'NO SCORES' not in all_scores:
            data_guides['CFD'] = acfd[x]
            table_to_file.append('CFD: '+str(acfd[x]))  # append CFD to table
            table_to_file.append('\t\t\t\tMismatches')

            # table_to_file.append('IN THE FOLLOWING MATRIX, THE FIRST GROUP OF '+str(max_bulges)+' LINES, ARE REFERED TO REFERENCE TARGET, THE SECOND GROUP OF '+str(max_bulges)+' LINES ARE REFERED TO VARIANT GENOME')

            table_to_file.append(
                data_general_count_copy.to_string(index=False))

            if genome_type == 'both':
                data_guides['Doench 2016'] = doench[x]
                # data_guides['Enriched'] = doench_enr[x]
            else:
                data_guides['Doench 2016'] = doench[x]

        if genome_type == 'both':
            # data_guides['Samples in Class 0 - 0+ - 1 - 1+'] = column_sample_class
            # data_guides['Genome'] = '\nREFERENCE\n-----------\nENRICHED\n'
            tmp = [str(i) for i in range(max_bulges+1)]*2
            tmp.insert(len(tmp)//2, "")
            data_guides["# Bulges"] = "\n".join(tmp)
        else:
            tmp = [str(i) for i in range(max_bulges+1)]
            # tmp.insert(len(tmp)//2, "")
            data_guides["# Bulges"] = "\n".join(tmp)
        # data_guides['Total'] = column_total

        # for i in range (1, mms + int(max_bulges) + 1):  #add count target for each total value
        #    data_guides[str(i) + ' MM + B'] = [str(x[i-1]) for x in column_off_target_ref]      #NOTE i-1 since i starts from 1
        data_guides['Total'] = []
        if genome_type == 'both':
            if max_bulges == 2:
                for i in range(len(data_guides["# Bulges"].split('\n'))-1):
                    if i == 1:
                        data_guides['Total'].append(
                            "REFERENCE\t"+str(sum(data_general_count.iloc[i, :])))
                    elif i == 2:
                        data_guides['Total'].append(
                            "\t"+str(sum(data_general_count.iloc[i, :])))
                        data_guides['Total'].append("\t")
                    elif i == 4:
                        data_guides['Total'].append(
                            "VARIANT\t\t"+str(sum(data_general_count.iloc[i, :])))
                    else:
                        data_guides['Total'].append(
                            "\t"+str(sum(data_general_count.iloc[i, :])))
            elif max_bulges == 1:
                for i in range(len(data_guides["# Bulges"].split('\n'))-1):
                    if i == 1:
                        data_guides['Total'].append(
                            "REFERENCE\t"+str(sum(data_general_count.iloc[i, :])))
                        data_guides['Total'].append("\t")
                    elif i == 3:
                        data_guides['Total'].append(
                            "VARIANT\t\t"+str(sum(data_general_count.iloc[i, :])))
                    else:
                        data_guides['Total'].append(
                            "\t"+str(sum(data_general_count.iloc[i, :])))
            else:
                for i in range(len(data_guides["# Bulges"].split('\n'))-1):
                    if i == 0:
                        data_guides['Total'].append(
                            "REFERENCE\t"+str(sum(data_general_count.iloc[i, :])))
                        data_guides['Total'].append("\t")
                    elif i == 1:
                        data_guides['Total'].append(
                            "VARIANT\t\t"+str(sum(data_general_count.iloc[i, :])))
        else:
            for i in range(len(data_guides["# Bulges"].split('\n'))):
                if i == len(data_guides["# Bulges"].split('\n'))//2:
                    data_guides['Total'].append(
                        "REFERENCE\t"+str(sum(data_general_count.iloc[i, :])))
                else:
                    data_guides['Total'].append(
                        "\t"+str(sum(data_general_count.iloc[i, :])))

        if genome_type == 'both':
            for i in range(mms+1):
                tmp = list(data_general_count.iloc[:, i].values.astype(str))
                tmp.insert(len(tmp)//2, "")
                data_guides[str(i) + 'MM'] = "\n".join(tmp)
        else:
            for i in range(mms+1):
                tmp = list(
                    data_general_count.iloc[:max_bulges+1, i].values.astype(str))
                # tmp.insert(len(tmp)//2, "")
                data_guides[str(i) + 'MM'] = "\n".join(tmp)

        data_guides['Total'] = "\n".join(data_guides['Total'])

        df.append(data_guides)
    dff = pd.DataFrame(df)

    table_to_file_save_dest = current_working_directory + \
        'Results/' + job_id + '/' + job_id + '.general_table.txt'

    outfile = open(table_to_file_save_dest, 'w')
    for elem in table_to_file:
        outfile.write(elem+'\n')
    outfile.close()

    # zip integrated results
    integrated_file = current_working_directory + 'Results/' + \
        job_id + '/' + job_id + '.bestMerge.txt.integrated_results.tsv'
    integrated_to_zip = current_working_directory + 'Results/' + \
        job_id + '/' + job_id + '.bestMerge.txt.integrated_results.zip'
    if not os.path.exists(integrated_to_zip):
        os.system(f"zip -j {integrated_to_zip} {integrated_file} &")

    if 'NO SCORES' not in all_scores:
        try:
            dff = dff.sort_values(['CFD', 'Doench 2016'],
                                  ascending=[False, False])
        except:  # for BOTH
            dff = dff.sort_values(['CFD', 'Enriched'],
                                  ascending=[False, False])
    else:
        try:
            dff = dff.sort_values('On-Targets Reference', ascending=True)
        except:
            dff = dff.sort_values('On-Targets Enriched', ascending=True)

    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

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
            ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
                for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    # Calculate sample count

    data_to_send = dff.iloc[
        page_current*page_size:(page_current + 1)*page_size
    ].to_dict('records')
    return data_to_send, [{'row': 0, 'column': 0}]

# Update color on selected row


@app.callback(
    Output('general-profile-table', 'style_data_conditional'),
    [Input('general-profile-table', 'selected_cells')],
    [State('general-profile-table', 'data')]
)
def colorSelectedRow(sel_cel, all_guides):
    if sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate
    guide = all_guides[int(sel_cel[0]['row'])]['Guide']
    return [
        {
            'if': {
                'filter_query': '{Guide} eq "' + guide + '"',
            },
            'background-color': 'rgba(0, 0, 255,0.15)'  # 'rgb(255, 102, 102)'

        },
        {
            'if': {
                'column_id': 'Genome'
            },
            'font-weight': 'bold',
            'textAlign': 'center'

        }
    ]


@app.callback(
    [Output('div-table-position', 'children'),
     Output('div-current-page-table-position', 'children')],
    [Input('div-position-filter-query', 'children')],
    [State('button-filter-position', 'n_clicks_timestamp'),
     State('url', 'search'),
     State('general-profile-table', 'selected_cells'),
     State('general-profile-table', 'data'),
     State('div-current-page-table-position', 'children'),
     State('div-mms-bulges-position', 'children')]
)
def filterPositionTable(filter_q, n, search, sel_cel, all_guides, current_page, mms_bulge):  # nPrev, nNext,

    if sel_cel is None:
        raise PreventUpdate
    if n is None:
        raise PreventUpdate

    filter_q = filter_q.split(',')
    chrom = filter_q[0]
    if chrom == 'None':
        raise PreventUpdate
        # chrom = None
    # pos = filter_q[1]
    # if pos == 'None':
    #     raise PreventUpdate
    #     # pos = None
    start = filter_q[1]
    if start == 'None':
        raise PreventUpdate

    end = filter_q[2]
    if end == 'None':
        raise PreventUpdate

    current_page = current_page.split('/')[0]
    current_page = int(current_page)
    mms = int(mms_bulge.split('-')[0])
    max_bulges = int(mms_bulge.split('-')[1])

    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    guide = all_guides[int(sel_cel[0]['row'])]['Guide']

    file_to_grep = job_directory + job_id + '.bestMerge.txt'
    # file_to_grep_alt = job_directory + job_id +'.altMerge.txt'
    pos_grep_result = current_working_directory + \
        'Results/' + job_id + '/' + job_id + '.' + start + "." + end + '.txt'

    # os.system(
    #    f'LC_ALL=C fgrep {guide} {file_to_grep} | LC_ALL=C grep -P \"{chrom}\t{pos}\t\" > {pos_grep_result}')
    # os.system(
    #     f'LC_ALL=C fgrep {guide} {file_to_grep} | awk \'$5 == \"{chrom}\" && ($6>={start} && $6<={end})\' | sort -k6,6n > {pos_grep_result}')
    os.system(
        f'awk \'$16 == \"{guide}\" && $5 == \"{chrom}\" && ($6>={start} && $6<={end})\' {file_to_grep} | sort -k6,6n > {pos_grep_result}')
    # pos_grep_result_zip = pos_grep_result.replace('txt', 'zip')
    # os.system(f'zip -j {pos_grep_result_zip} {pos_grep_result}')

    with open(file_to_grep, 'r') as ftg:
        header = ftg.readline().split('\t')[:24]
    try:
        df = pd.read_csv(pos_grep_result, sep='\t',
                         header=None, usecols=range(0, 24))
    except:
        df = pd.DataFrame(columns=header)
    df.rename(columns=COL_BOTH_RENAME, inplace=True)
    # df.columns = header
    # print(df, 'line 2057')
    try:
        df = df[COL_BOTH]
        df_check = True
    except:
        df_check = False #skip df parsing and report no results found
    # df.columns = COL_BOTH
    # df[''] = [''] * df.shape[0]
    # df_cols = df.columns.tolist()
    # df_cols.remove('Samples')
    # df_cols.remove('')
    # df_cols.append('Samples')
    # # df_cols.insert(0, '')
    # df = df[df_cols]
    # print(df, 'position df line 2065')
    if df_check:
        out_1 = [
            dash_table.DataTable(
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}],
                id="table-position",
                export_format="csv",
                columns=[{"name": COL_BOTH[count], "id": i, 'hideable':True}
                         for count, i in enumerate(df.columns)],
                # columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                style_cell_conditional=[{
                    'if': {'column_id': 'Samples'},
                    'textAlign': 'left'
                }],
                style_table={
                    'overflowX': 'scroll',
                },
                page_size=10,

            )
        ]
    else:
        out_1 = [html.P('No results found with this genomic coordinates')]
    os.system(f"rm {pos_grep_result}")
    return out_1, '1/' + str(1)

# Callback to update the hidden div filter position


@app.callback(
    Output('div-position-filter-query', 'children'),
    [Input('button-filter-position', 'n_clicks')],
    [State('dropdown-chr-table-position', 'value'),
     State('input-position-start', 'value'),
     State('input-position-end', 'value')]
)
def updatePositionFilter(n, chrom, pos_start, pos_end):  # , pos_end

    if n is None:
        raise PreventUpdate

    if pos_start == '':
        pos_start = 'None'
    if pos_end == '':
        pos_end = 'None'
    return str(chrom) + ',' + str(pos_start) + ',' + str(pos_end)

# Callback to view next/prev page on sample table


@app.callback(
    [Output('div-table-samples', 'children'),
     Output('div-current-page-table-samples', 'children')],
    [Input('prev-page-sample', 'n_clicks_timestamp'),
     Input('next-page-sample', 'n_clicks_timestamp'),
     Input('div-sample-filter-query', 'children')],
    [State('button-filter-population-sample', 'n_clicks_timestamp'),
     State('url', 'search'),
     State('general-profile-table', 'selected_cells'),
     State('general-profile-table', 'data'),
     State('div-current-page-table-samples', 'children')]
)
def filterSampleTable(nPrev, nNext, filter_q, n, search, sel_cel, all_guides, current_page):
    if sel_cel is None:
        raise PreventUpdate
    if nPrev is None and nNext is None and n is None:
        raise PreventUpdate

    if nPrev is None:
        nPrev = 0
    if nNext is None:
        nNext = 0
    if n is None:
        n = 0

    sup_pop = filter_q.split(',')[0]
    pop = filter_q.split(',')[1]
    samp = filter_q.split(',')[2]
    if sup_pop == 'None':
        sup_pop = None
    if pop == 'None':
        pop = None
    if samp == 'None' or samp == 'NONE':
        samp = None
    current_page = current_page.split('/')[0]
    current_page = int(current_page)
    btn_sample_section = []
    btn_sample_section.append(n)
    btn_sample_section.append(nPrev)
    btn_sample_section.append(nNext)
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    population_1000gp = associateSample.loadSampleAssociation(
        job_directory + 'sampleID.txt')[2]
    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        genome_type_f = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        ref_comp = (next(s for s in all_params.split(
            '\n') if 'Ref_comp' in s)).split('\t')[-1]

    genome_type = 'ref'
    if '+' in genome_type_f:
        genome_type = 'var'
    if 'True' in ref_comp:
        genome_type = 'both'

    guide = all_guides[int(sel_cel[0]['row'])]['Guide']
    if genome_type == 'both':
        col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference',
                            'Targets in Variant', 'Targets in Population', 'Targets in Super Population', 'PAM Creation', 'Class']
    else:
        col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference',
                            'Targets in Variant', 'Targets in Population', 'Targets in Super Population', 'PAM Creation', 'Class']
    # Last button pressed is filtering, return the first page of the filtered table
    if max(btn_sample_section) == n:
        if genome_type == 'both':
            df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                             guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
            df = df.sort_values('Targets in Variant', ascending=False)
            df.drop(['Targets in Reference'], axis=1, inplace=True)
        else:
            df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                             guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
            df = df.sort_values('Targets in Reference', ascending=False)
            df.drop(['Targets in Reference'], axis=1, inplace=True)
            df.drop(['Class'], axis=1, inplace=True)
        more_info_col = []
        for i in range(df.shape[0]):
            more_info_col.append('Show Targets')
        df[''] = more_info_col
        if (sup_pop is None or sup_pop == '') and (pop is None or pop == '') and (samp is None or samp == ''):  # No filter value selected
            max_page = len(df.index)
            max_page = math.floor(max_page / 10) + 1
            return generate_table_samples(df, 'table-samples', 1, guide, job_id), '1/' + str(max_page)
        if samp is None or samp == '':
            if pop is None or pop == '':
                df.drop(
                    df[(~(df['Population'].isin(population_1000gp[sup_pop])))].index, inplace=True)
            else:
                df.drop(df[(df['Population'] != pop)].index, inplace=True)
        else:
            df.drop(df[(df['Sample'] != samp)].index, inplace=True)
        max_page = len(df.index)
        max_page = math.floor(max_page / 10) + 1
        return generate_table_samples(df, 'table-samples', 1, guide, job_id), '1/' + str(max_page)
    else:
        if max(btn_sample_section) == nNext:
            current_page = current_page + 1
            if genome_type == 'both':
                df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                                 guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
                df = df.sort_values('Targets in Variant', ascending=False)
                df.drop(['Targets in Reference'], axis=1, inplace=True)
            else:
                df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                                 guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
                df = df.sort_values('Targets in Reference', ascending=False)
                df.drop(['Targets in Reference'], axis=1, inplace=True)
                df.drop(['Class'], axis=1, inplace=True)
            more_info_col = []
            for i in range(df.shape[0]):
                more_info_col.append('Show Targets')
            df[''] = more_info_col
            # Active filter
            if pop or sup_pop or samp:
                if samp is None or samp == '':
                    if pop is None or pop == '':
                        df.drop(
                            df[(~(df['Population'].isin(population_1000gp[sup_pop])))].index, inplace=True)
                    else:
                        df.drop(df[(df['Population'] != pop)].index,
                                inplace=True)
                else:
                    df.drop(df[(df['Sample'] != samp)].index, inplace=True)

            if ((current_page - 1) * 10) > len(df):
                current_page = current_page - 1
                if current_page < 1:
                    current_page = 1
            max_page = len(df.index)
            max_page = math.floor(max_page / 10) + 1
            return generate_table_samples(df, 'table-samples', current_page, guide, job_id), str(current_page) + '/' + str(max_page)

        else:  # Go to previous page
            current_page = current_page - 1
            if current_page < 1:
                current_page = 1
            if genome_type == 'both':
                df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                                 guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
                df = df.sort_values('Targets in Variant', ascending=False)
                df.drop(['Targets in Reference'], axis=1, inplace=True)
            else:
                df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                                 guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
                df = df.sort_values('Targets in Reference', ascending=False)
                df.drop(['Targets in Reference'], axis=1, inplace=True)
                df.drop(['Class'], axis=1, inplace=True)
            more_info_col = []
            for i in range(df.shape[0]):
                more_info_col.append('Show Targets')
            df[''] = more_info_col
            if pop or sup_pop or samp:
                if samp is None or samp == '':
                    if pop is None or pop == '':
                        df.drop(
                            df[(~(df['Population'].isin(population_1000gp[sup_pop])))].index, inplace=True)
                    else:
                        df.drop(df[(df['Population'] != pop)].index,
                                inplace=True)
                else:
                    df.drop(df[(df['Sample'] != samp)].index, inplace=True)
            max_page = len(df.index)
            max_page = math.floor(max_page / 10) + 1
            return generate_table_samples(df, 'table-samples', current_page, guide, job_id), str(current_page) + '/' + str(max_page)
    raise PreventUpdate

# Callback to update the hidden div filter


@app.callback(
    Output('div-sample-filter-query', 'children'),
    [Input('button-filter-population-sample', 'n_clicks')],
    [State('dropdown-superpopulation-sample', 'value'),
     State('dropdown-population-sample', 'value'),
     State('input-sample', 'value')]
)
def updateSampleFilter(n, superpopulation, population, sample):
    if n is None:
        raise PreventUpdate
    return str(superpopulation) + ',' + str(population) + ',' + str(sample).replace(' ', '').upper()
# Callback to update the sample based on population selected


@app.callback(
    [Output('dropdown-sample', 'options'),
     Output('dropdown-sample', 'value')],
    [Input('dropdown-population-sample', 'value')],
    [State('url', 'search')]
)
def updateSampleDrop(pop, search):
    if pop is None or pop == '':
        return [], None
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    dict_pop = associateSample.loadSampleAssociation(
        job_directory + 'sampleID.txt')[3]
    return [{'label': sam, 'value': sam} for sam in dict_pop[pop]], None

# Callback to update the population tab based on superpopulation selected


@app.callback(
    [Output('dropdown-population-sample', 'options'),
     Output('dropdown-population-sample', 'value')],
    [Input('dropdown-superpopulation-sample', 'value')],
    [State('url', 'search')]
)
def updatePopulationDrop(superpop, search):
    if superpop is None or superpop == '':
        raise PreventUpdate
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    population_1000gp = associateSample.loadSampleAssociation(
        job_directory + 'sampleID.txt')[2]
    return [{'label': i, 'value': i} for i in population_1000gp[superpop]], None


def generate_table_position(dataframe, id_table, page, mms, bulges, guide='', job_id='', max_rows=10):
    rows_remaining = len(dataframe) - (page - 1) * max_rows
    header = [html.Tr([
        html.Th('Chromosome', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Position', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Best Target', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Min Mismatch', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Min Bulge', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Bulge', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Targets in Cluster by Mismatch Value', colSpan=str(
            mms + 1), style={'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('', rowSpan='2'),
    ])
    ]
    mms_header = []
    for mm in range(mms + 1):
        mms_header.append(html.Th(
            str(mm) + ' MM', style={'vertical-align': 'middle', 'text-align': 'center'}))
    header.append(html.Tr(mms_header))

    data = []
    for i in range(min(rows_remaining, max_rows)):
        first_cells = [
            html.Td(dataframe.iloc[i + (page - 1)*max_rows]['Chromosome'], rowSpan=str(
                bulges + 1),  style={'vertical-align': 'middle', 'text-align': 'center'}),
            html.Td(dataframe.iloc[i + (page - 1)*max_rows]['Position'], rowSpan=str(
                bulges + 1),  style={'vertical-align': 'middle', 'text-align': 'center'}),
            html.Td(dataframe.iloc[i + (page - 1)*max_rows]['Best Target'], rowSpan=str(
                bulges+1),  style={'vertical-align': 'middle', 'text-align': 'center'}),
            html.Td(dataframe.iloc[i + (page - 1)*max_rows]['Min Mismatch'], rowSpan=str(
                bulges+1),  style={'vertical-align': 'middle', 'text-align': 'center'}),
            html.Td(dataframe.iloc[i + (page - 1)*max_rows]['Min Bulge'], rowSpan=str(
                bulges+1),  style={'vertical-align': 'middle', 'text-align': 'center'}),
            html.Th('0 Bulge', style={
                    'vertical-align': 'middle', 'text-align': 'center', 'padding-left': '0'})
        ]

        mm_cells = [html.Td(dataframe.iloc[i + (page - 1)*max_rows][col], style={
                            'vertical-align': 'middle', 'text-align': 'center'}) for col in dataframe.columns[5:5+mms+1]]
        data.append(html.Tr(first_cells + mm_cells + [html.Td(
            html.A('Show Targets',  href='result?job=' + job_id + '#' + guide + '-Pos-' + str(dataframe.iloc[i + (
                page - 1)*max_rows]['Chromosome']) + '-' + str(dataframe.iloc[i + (page - 1)*max_rows]['Position']), target='_blank'),
            rowSpan=str(bulges+1), style={'vertical-align': 'middle', 'text-align': 'center'})
        ]))
        for b in range(bulges):
            data.append(
                html.Tr(
                    [html.Th(str(
                        b + 1) + ' Bulge', style={'vertical-align': 'middle', 'text-align': 'center'})]
                    +
                    [html.Td(dataframe.iloc[i + (page - 1)*max_rows][col])
                     for col in dataframe.columns[5 + (b + 1) * (mms + 1): 5 + (b + 1) * (mms+1) + mms + 1]]
                )
            )

    return html.Table(header + data, style={'display': 'inline-block'}, id=id_table)


def generate_table_samples(dataframe, id_table, page, guide='', job_id='', max_rows=10):
    '''
    Per generare una html table. NOTE è diversa da una dash dataTable
    '''
    rows_remaining = len(dataframe) - (page - 1) * max_rows
    return html.Table(
        # Header
        # [html.Tr([html.Th(col) if col != 'Targets in Enriched'  else html.Th(html.Abbr(col, title = 'Counting of targets that are generated by the insertion of a IUPAC nucleotide of a sample',
        # style = {'text-decoration':'underline'})) for col in dataframe.columns]) ] +
        [html.Tr([html.Th(col) for col in dataframe.columns])] +
        # Body
        [html.Tr([
            html.Td(html.A(dataframe.iloc[i + (page - 1)*max_rows][col],  href='result?job=' + job_id + '#' + guide + '-Sample-' + dataframe.iloc[i + (page - 1)*max_rows]['Sample'], target='_blank')) if col == '' else html.Td(dataframe.iloc[i + (page - 1)*max_rows][col]) for col in dataframe.columns
        ]) for i in range(min(rows_remaining, max_rows))],
        style={'display': 'inline-block'},
        id=id_table
    )


def generate_table(dataframe, id_table, genome_type, guide='', job_id='', max_rows=2600):
    '''
    Per generare una html table. NOTE è diversa da una dash dataTable
    '''
    # if genome_type == 'both':
    header = [html.Tr([
        html.Th('Bulge type', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Mismatches', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Bulge Size', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('Targets found in Genome', colSpan=str(3), style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('PAM Creation', rowSpan='2', style={
                'vertical-align': 'middle', 'text-align': 'center'}),
        html.Th('', rowSpan='2'),
    ])
    ]
    # 'Bulge Type' 'Mismatches' 'Bulge Size' 'Targets in Reference' 'Targets in Enriched' 'Combined' 'PAM Creation' ''

    header.append(html.Tr([html.Th(x, style={
                  'vertical-align': 'middle', 'text-align': 'center'}) for x in ['Reference', 'Variant', 'Combined']]))
    # else:
    #    header = [html.Tr([html.Th(col, style = {'vertical-align':'middle', 'text-align':'center'}) for col in dataframe.columns])]
    return html.Table(
        header +
        # Body
        [html.Tr([
            html.Td(html.A(dataframe.iloc[i][col],  href='result?job=' + job_id + '#' + guide + 'new' + dataframe.iloc[i]['Bulge Type'] + str(dataframe.iloc[i]['Bulge Size']) + str(dataframe.iloc[i]['Mismatches']), target='_blank'), style={'vertical-align': 'middle', 'text-align': 'center'}) if col == '' else html.Td(dataframe.iloc[i][col], style={'vertical-align': 'middle', 'text-align': 'center'}) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))],
        style={'display': 'inline-block'},
        id=id_table
    )


def check_existance_sample(job_directory, job_id, sample):
    df = pd.read_csv(job_directory + job_id + 'sampleID.txt', sep='\t')
    samples = df.iloc[:, 0]
    if sample in samples.values:
        return True
    else:
        return False

# Select figures on mms value, sample value


@app.callback(
    # [Output('div-guide-image', 'children'),
    #  Output('div-sample-image', 'children')],
    [Output('div-radar-chart-total', 'children'),
     Output('div-population-barplot', 'children'),
     Output('div-sample-image', 'children'),
     #  Output('div-radar-chart-sample', 'children')],
     Output('row-radar-chart-sample', 'children')],
    [Input('mm-dropdown', 'value'),
     Input('blg-dropdown', 'value'),
     Input('dropdown-superpopulation-sample', 'value'),
     Input('dropdown-population-sample', 'value'),
     Input('dropdown-sample', 'value'),
     Input('general-profile-table', 'selected_cells')],
    [State('url', 'search'),
     State('general-profile-table', 'data')]
)
def updateImagesTabs(mm, bulge, superpopulation, population, sample, sel_cel, search, all_guides):
    # if sel_cel is None:
    #     raise PreventUpdate
    # print('entro update tab')
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    guide = all_guides[int(sel_cel[0]['row'])]['Guide']

    # search for getting job id
    # get guide with sel_cel and all_data
    radar_chart_images = list()
    population_barplots = list()
    guide_images = list()
    sample_images = list()

    try:
        population_barplots.extend(
            [
                html.A(
                    html.Img(
                        src='data:image/png;base64,{}'.format(base64.b64encode(open(
                            current_working_directory + 'Results/' + job_id + '/imgs/populations_distribution_' + guide + '_' + str(int(mm)+int(bulge)) + 'total.png', 'rb').read()).decode()),
                        id='distribution-population' + str(int(mm)+int(bulge)), width="100%", height="auto"
                    ),
                    target="_blank",
                    href='/Results/' + job_id + '/imgs/' + 'populations_distribution_' +
                    guide + '_' +
                    str(int(mm)+int(bulge)) + 'total.png'
                ),
                html.Div(html.P('Distribution ' + str(int(mm)+int(bulge)) + ' Mismatches + Bulges ', style={
                    'display': 'inline-block'}), style={'text-align': 'center'})
            ]
        )
    except:
        population_barplots = [
            html.Div(
                html.H2(
                    "No result found for this combination of mismatches and bulges"
                )
            )
        ]

    # try:
    #     guide_images.extend(  # population barplot
    #         [
    #             html.A(
    #                 html.Img(
    #                     src='data:image/png;base64,{}'.format(base64.b64encode(open(
    #                         current_working_directory + 'Results/' + job_id + '/imgs/populations_distribution_' + guide + '_' + str(int(mm)+int(bulge)) + 'total.png', 'rb').read()).decode()),
    #                     id='distribution-population' + str(int(mm)+int(bulge)), width="100%", height="auto"
    #                 ),
    #                 target="_blank",
    #                 href='/Results/' + job_id + '/imgs/' + 'populations_distribution_' +
    #                 guide + '_' +
    #                 str(int(mm)+int(bulge)) + 'total.png'
    #             ),
    #             html.Div(html.P('Distribution ' + str(int(mm)+int(bulge)) + ' Mismatches + Bulges ', style={
    #                 'display': 'inline-block'}), style={'text-align': 'center'})
    #         ]
    #     )
    # except:
    #     guide_images.append(
    #         html.Div(html.P('No Targets found with ' + str(int(mm)+int(bulge)) + ' Mismatches + Bulges', style={
    #             'display': 'inline-block'}), style={'text-align': 'center'}),
    #         # html.Div(html.P('Distribution ' + str(mm) + ' Mismatches + Bulges ', style = {'display':'inline-block'} ),style = {'text-align':'center'})
    #     )

    # guide = guide.replace("N", "")
    radar_img = '/imgs/summary_single_guide_' + \
        guide + '_' + str(mm) + \
        '.' + str(bulge) + '_TOTAL.png'

    if not os.path.isfile(f"{job_directory}/{radar_img}"):
        # try:
        # print('faccio radar chart')
        os.system(f"python {app_main_directory}/PostProcess/generate_img_radar_chart.py {guide} {job_directory}/guide_dict_{guide}.json {job_directory}/motif_dict_{guide}.json {mm} {bulge} TOTAL {job_directory}/imgs/")
        # except:
        # pass

    img_found = False
    try:
        radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(open(
            current_working_directory + 'Results/' + job_id + '/' + radar_img, 'rb').read()).decode())
        img_found = True
    except:
        pass
        # radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(open(
        #     current_working_directory+'assets/placeholder.png', 'rb').read()).decode())
    try:
        radar_href = '/Results/' + job_id + '/' + radar_img
    except:
        radar_href = ''

    if img_found:
        radar_chart_images.append(
            html.A(
                html.Img(src=radar_src, id='radar-img-guide',
                         width="100%", height="auto"),
                target="_blank",
                href=radar_href
            )
        )
    else:
        radar_chart_images.append(
            html.H2(
                "No result found for this combination of mismatches and bulges"
            )
        )
    class_images = [(sample, 'Samples'), (population, 'Population'),
                    (superpopulation, 'Superpopulation')]

    for c in class_images:
        img_found = False
        if c[0]:
            current_img = job_directory + '/imgs/summary_single_guide_' +\
                guide + '_' + str(mm) + '.'+str(bulge) + '_' + c[0] + '.png'
            if not os.path.isfile(current_img):
                try:
                    os.system(
                        f"python {app_main_directory}/PostProcess/generate_img_radar_chart.py {guide} {job_directory}/guide_dict_{guide}.json {job_directory}/motif_dict_{guide}.json {mm} {bulge} {c[0]} {job_directory}/imgs/")
                except:
                    pass
            try:
                first_img_source = 'data:image/png;base64,{}'.format(
                    base64.b64encode(open(current_img, 'rb').read()).decode())
                img_found = True
            except:
                pass
                # first_img_source = 'data:image/png;base64,{}'.format(base64.b64encode(
                #     open(current_working_directory+'/assets/placeholder.png', 'rb').read()).decode())
            try:
                first_img_href = 'Results/' + job_id + '/imgs/summary_single_guide_' +\
                    guide + '_' + str(mm) + "." + str(bulge) +\
                    '_' + c[0] + '.png'
            except:
                first_img_href = ''
            # sample_images.append(dbc.Row(html.Br()))

            if img_found:
                sample_images.append(
                    dbc.Col(
                        html.A(
                            html.Img(src=first_img_source,
                                     width="100%", height="auto"),
                            target="_blank",
                            href=first_img_href
                        )
                    )
                )
            else:
                sample_images.append(
                    dbc.Col(
                        html.H2(
                            "No result found for this combination of mismatches and bulges"
                        )
                    )
                )
    # reverse list to print plots in correct order since they are append in reverse order into main sample_images list
    reversed_sample_images = sample_images[::-1]
    return radar_chart_images, population_barplots, guide_images, reversed_sample_images


# Open in browser the result directory
# @app.callback(
#     Output('div-open-result-directory', 'children'),
#     [Input('button-open-result-directory', 'n_clicks')],
#     [State('url', 'search'),
#      State('div-open-result-directory', 'children')]
# )
# def openResultDirectory(n, search, guide):
#     if n is None:
#         raise PreventUpdate
#     # TODO decidere se aprire tutta la cartella, solo il file, e nel secondo caso creare copia submit job che non rimuova .targets.GUIDE.txt
#     job_id = search.split('=')[-1]
#     wb.open_new_tab(current_working_directory + 'Results/' +
#                     job_id + '/' + job_id + '.targets.' + guide + '.txt')
#     raise PreventUpdate


@ app.callback(
    [Output('download-link-personal-card', 'children'),
     Output('download-link-personal-card', 'hidden'),
     Output('div-personal-plot', 'children'),
     Output('div-private-plot', 'children'),
     Output('div-table-sample-card', 'children'),
     Output('div-top-target-sample-card', 'children')],
    [Input('button-sample-card', 'n_clicks')],
    [State('dropdown-sample-card', 'value'),
     State('general-profile-table', 'selected_cells'),
     State('general-profile-table', 'data'),
     State('url', 'search')]
)
def generate_sample_card(n, sample, sel_cel, all_guides, search):
    if n is None:
        raise PreventUpdate
    guide = all_guides[int(sel_cel[0]['row'])]['Guide']
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    file_to_grep = job_directory + job_id + '.bestMerge.txt'
    if not os.path.exists(current_working_directory + 'Results/' + job_id + '/' + job_id + '.' + str(sample) + '.' + guide + '.sample_card.txt'):
        df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                         guide+'.txt', sep='\t', skiprows=1, index_col=0, header=None)
        personal = df.loc[sample, 4]
        pam_creation = df.loc[sample, 8]

        file_to_grep = job_directory + job_id + '.bestMerge.txt'
        integrated_to_grep = job_directory+job_id + \
            '.bestMerge.txt.integrated_results.tsv'
        integrated_personal = job_directory + job_id + '.' + \
            str(sample) + '.' + guide + '.integrated.personal.txt'
        integrated_private = job_directory + job_id + '.' + \
            str(sample) + '.' + guide + '.integrated.private.txt'
        # file_to_grep_alt = job_directory + job_id +'.altMerge.txt'
        sample_grep_result = current_working_directory + 'Results/' + \
            job_id + '/' + job_id + '.' + str(sample) + '.' + guide + '.private.txt'

        # copy header from integrated results into sample files
        os.system(f"head -1 {integrated_to_grep} > {integrated_personal}")
        os.system(f"head -1 {integrated_to_grep} > {integrated_private}")
        # grep guide and then sample into personal card data
        os.system(
            f"LC_ALL=C fgrep {guide} {integrated_to_grep} | fgrep {sample} >> {integrated_personal}")
        # grep private targets from personal targets
        os.system(
            f"LC_ALL=C awk \'$32==\"{sample}\"\' {integrated_personal} >> {integrated_private}")
        # grep private targets to generate table and file
        os.system(
            f"LC_ALL=C fgrep {guide} {file_to_grep} | awk \'$14==\"{sample}\"\' > {sample_grep_result}")

        # plot for images in personal card
        os.system(
            f"python {app_main_directory}/PostProcess/CRISPRme_plots_personal.py {integrated_personal} {current_working_directory}/Results/{job_id}/imgs/ {guide}.{sample}.personal > {current_working_directory}/Results/{job_id}/warnings.txt 2>&1")
        os.system(
            f"python {app_main_directory}/PostProcess/CRISPRme_plots_personal.py {integrated_private} {current_working_directory}/Results/{job_id}/imgs/ {guide}.{sample}.private > {current_working_directory}/Results/{job_id}/warnings.txt 2>&1")
        os.system(
            f"rm -f {current_working_directory}/Results/{job_id}/warnings.txt {integrated_private} {integrated_personal}")

        private = 0
        for line in open(sample_grep_result):
            private += 1
        # print(personal, pam_creation, private)
        results_table = pd.DataFrame([[personal, pam_creation, private]], columns=[
            'Personal', 'PAM Creation', 'Private']).astype(str)
        if int(private) > 0:
            tmp_file = current_working_directory + 'Results/' + \
                job_id + '/' + job_id + '.' + str(sample) + ".tmp_card.txt"
            tmp_file_2 = current_working_directory + 'Results/' + \
                job_id + '/' + job_id + '.' + str(sample) + ".tmp_card_2.txt"
            os.system(
                f"LC_ALL=C sort -k21,21rg \"{sample_grep_result}\" > \"{tmp_file}\" ; head -5 {tmp_file} > \"{tmp_file_2}\"")
            ans = pd.read_csv(tmp_file_2, sep='\t',
                              header=None, usecols=range(0, 23))
            with open(file_to_grep) as f_:
                c = f_.readline().strip()
            ans.columns = c.split('\t')[:23]
            ans = ans.astype(str)
            # os.system(f"rm {tmp_file} &") #do not delete temp file until zip is created
            os.system(f"rm {tmp_file_2} &")
            # create zip file to download result card /blocking operation on the system to avoid updating the page before the zip is created
            os.system('zip '+'-j ' + tmp_file.replace('.txt',
                                                      '.zip') + ' ' + tmp_file)
            # do not delete temp file until zip is created
            os.system(f"rm {tmp_file} &")

        with open(current_working_directory + 'Results/' + job_id + '/' + job_id + '.' + str(sample) + '.' + guide + '.sample_card.txt', "w") as file_out:
            file_out.write(
                '\t'.join(results_table.iloc[0, :].values.tolist()) + '\n')
            if int(private) > 0:
                up_to = min([ans.shape[0], 5])
                for row in range(up_to):
                    file_out.write(
                        '\t'.join(ans.iloc[row, :].values.tolist()) + '\n')

        # os.system(f"{app_main_directory}/PostProcess/personal_cards.py {current_working_directory}/Results/{job_id}/{job_id}.{sample}.{guide}.sample_card.txt {current_working_directory}/Results/{job_id}")
        os.system(f"rm {sample_grep_result}")
    else:
        with open(current_working_directory + 'Results/' + job_id + '/' + job_id + '.' + str(sample) + '.' + guide + '.sample_card.txt', "r") as file_in:
            infos = file_in.readline().strip().split('\t')
            results_table = pd.DataFrame([[infos[0], infos[1], infos[2]]], columns=[
                'Personal', 'PAM Creation', 'Private'])
            if int(infos[2]) > 0:
                targets = []
                for line in file_in:
                    targets.append(line.strip().split('\t'))

                ans = pd.DataFrame(targets)
                # ans = ans.iloc[:,:23]
                with open(file_to_grep) as f_:
                    c = f_.readline().strip()
                ans.columns = c.split('\t')[:23]

    # image for personal and private
    try:
        image_personal_top = 'data:image/png;base64,{}'.format(base64.b64encode(open(
            current_working_directory + 'Results/' + job_id + f'/imgs/CRISPRme_top_1000_log_for_main_text_{guide}.{sample}.personal.png', 'rb').read()).decode())
        image_private_top = 'data:image/png;base64,{}'.format(base64.b64encode(open(
            current_working_directory + 'Results/' + job_id + f'/imgs/CRISPRme_top_1000_log_for_main_text_{guide}.{sample}.private.png', 'rb').read()).decode())
    except:
        sys.stderr.write('PERSONAL AND PRIVATE LOLLIPOP PLOTS NOT GENERATED')

    # print(ans, 'table sample line 2784')
    # ans.columns(COL_BOTH_RENAME, inplace=True)
    # print(ans, 'table sample post rename line 2786')
    
    # print(ans, 'table sample line 2792')

    try:
        file_to_load = job_id + '.' + sample + '.tmp_card.zip'
        ans = ans[['#Bulge_type', 'crRNA', 'DNA', 'Reference', 'Chromosome',
                   'Position', 'Direction', 'Mismatches',
                   'Bulge_Size', 'PAM_gen', 'Samples', 'SNP',
                   'CFD', 'CFD_ref', 'Highest_CFD_Risk_Score',
                   'AF', 'Annotation_Type']]
        # ans[''] = [''] * ans.shape[0]  # taaaaaaaaaac
        # ans_cols = ans.columns.tolist()
        # ans_cols.remove('Samples')
        # ans_cols.append('Samples')
        # ans_cols.remove('')
        # # ans_cols.insert(0, '')
        # ans = ans[ans_cols]
        out_1 = [
            html.A('Download private targets', href=URL+'/Results/' +
                   job_id + '/' + file_to_load, target='_blank'),
            False,
            [
                html.P('Top 100 Personal Targets per CFD score'),
                html.A(
                    html.Img(src=image_personal_top, id='sample-personal-top',
                             width="100%", height="auto"),
                    target="_blank"
                )
            ],
            [
                html.P('Top 100 Private Targets per CFD score'),
                html.A(
                    html.Img(src=image_private_top, id='sample-private-top',
                             width="100%", height="auto"),
                    target="_blank"
                )
            ],
            dash_table.DataTable(
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}],
                id="results-table",
                columns=[{"name": i, "id": i} for i in results_table.columns],
                data=results_table.to_dict('records'),
                # style_table={
                #     'overflowX': 'scroll'
                # }
            ),
            dash_table.DataTable(
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}],
                id="results-table-risk",
                columns=[{"name": COL_BOTH[count], "id": i, 'hideable':True}
                         for count, i in enumerate(ans.columns)],
                data=ans.to_dict('records'),
                # style_cell_conditional=[
                #     {'if': {'column_id': 'Bulge_type'},
                #      'width': '250px'},
                # ],
                style_table={
                    'overflowX': 'scroll'
                }
            )
        ]
    except:
        out_1 = [
            # dbc.Col(
            #    html.A(
            #        html.Img(src=image_sample_card, id='sample-card-img',
            #                 width="100%", height="auto"),
            #        target="_blank",
            #    ),
            #    width=10
            # ),
            html.A('Download private targets', href=URL+'/Results/' +
                   job_id + '/' + file_to_load, target='_blank'),
            True,
            [
                html.P('Top 100 Personal Targets per CFD score'),
                html.A(
                    html.Img(src=image_personal_top, id='sample-personal-top',
                             width="100%", height="auto"),
                    target="_blank"
                )
            ],
            [
                html.P('Top 100 Private Targets per CFD score'),
                html.A(
                    html.Img(src=image_private_top, id='sample-private-top',
                             width="100%", height="auto"),
                    target="_blank"
                )
            ],
            dash_table.DataTable(
                css=[{'selector': '.row',
                      'rule': 'margin: 0'}],
                id="results-table",
                columns=[{"name": i, "id": i} for i in results_table.columns],
                data=results_table.to_dict('records')
            ),
            []
        ]

    return list(out_1)

# Load the table/children under the tab value


@ app.callback(
    Output('div-tab-content', 'children'),
    [Input('tabs-reports', 'value'),
     Input('general-profile-table', 'selected_cells')],
    [State('general-profile-table', 'data'),
     State('url', 'search'),
     State('div-genome-type', 'children')]


)
def updateContentTab(value, sel_cel, all_guides, search, genome_type):
    if value is None or sel_cel is None or not sel_cel or not all_guides:
        raise PreventUpdate

    guide = all_guides[int(sel_cel[0]['row'])]['Guide']
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'

    with open(current_working_directory + 'Results/' + job_id + '/Params.txt') as p:
        all_params = p.read()
        mms = (next(s for s in all_params.split('\n')
                    if 'Mismatches' in s)).split('\t')[-1]
        genome_selected = (next(s for s in all_params.split(
            '\n') if 'Genome_selected' in s)).split('\t')[-1]
        max_bulges = (next(s for s in all_params.split('\n')
                           if 'Max_bulges' in s)).split('\t')[-1]

    fl = []
    fl.append(html.Br())
    fl.append(html.H5('Focus on: ' + guide))

    if value == 'tab-summary-by-guide':  # BUG se cambio guida selezionata due volte mi cambia il mms mettendo a 0, provare con un div nascosto
        # Show Summary by Guide table
        fl.append(
            html.P(
                ['Summary table counting the number of targets found in the Reference and Variant Genome for each combination of Bulge Type, Bulge Size and Mismatch. Select \'Show Targets\' to view the corresponding list of targets. ',
                 # html.A('Click here', href = URL + '/data/' + job_id + '/' + job_id + '.targets.' + guide + '.zip' ,target = '_blank', id = 'download-full-list' ), ' to download the full list of targets.'
                 ]
            )
        )
        # fl.append(html.Button(html.A('Download Full list of targets', href = URL + '/data/' + job_id + '/' + job_id + '.targets.' + guide + '.zip' ,target = '_blank', style = {'text-decoration':'none', 'color':'black'} )))

        # fl.append(html.Button('Open Result Directory',
        #                       id='button-open-result-directory'))
        # fl.append(html.Div(guide, id='div-open-result-directory',
        #                    style={'display': 'none'}))
        fl.append(html.Br())
        df = pd.read_csv(job_directory + job_id +
                         '.summary_by_guide.' + guide + '.txt', sep='\t')
        more_info_col = []
        total_col = []
        for i in range(df.shape[0]):
            more_info_col.append('Show Targets')
            total_col.append(df['Bulge Size'])

        # #Swap pam creaation and Combined column
        # if genome_type == 'both':   #TODO fixare per var e ref
            # df['Combined'] = df['Targets in Reference'] + df['Targets in Enriched']
            # #['Guide' 'Bulge Type' 'Bulge Size' 'Mismatches' 'Targets in Reference', 'Targets in Enriched' 'PAM Creation' 'Combined']
            # df = df[['Guide' , 'Bulge Type', 'Mismatches', 'Bulge Size' , 'Targets in Reference', 'Targets in Enriched', 'Combined', 'PAM Creation']]
        df[''] = more_info_col
        # df['Total'] = df['Bulge Size'] + df['Mismatches']
        # if genome_type == 'both' and genome_type == 'var':
        # df = df.sort_values(['Total', 'Targets in Enriched'], ascending = [True, False])
        # else:
        # df = df.sort_values('Total', ascending = True)
        # del df['Total']
        # del df['Guide']
        fl.append(html.Div(
            generate_table(df, 'table-summary-by-guide', genome_type, guide, job_id), style={'text-align': 'center'}
        )
        )
        return fl
    elif value == 'tab-summary-by-sample':
        # Show Summary by Sample table
        fl.append(
            html.P('Summary table counting the number of targets found in the Variant Genome for each sample. Filter the table by selecting the Population or Superpopulation desired from the dropdowns.')
        )
        if genome_type == 'both':
            # col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference', 'Targets in Enriched', 'Targets in Population', 'Targets in Super Population', 'PAM Creation', 'Class']
            col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference',
                                'Targets in Variant', 'Targets in Population', 'Targets in Super Population', 'PAM Creation']
            df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                             guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
            df = df.sort_values('Targets in Variant', ascending=False)
            df.drop(['Targets in Reference'], axis=1, inplace=True)
        else:
            # col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference', 'Targets in Enriched', 'Targets in Population', 'Targets in Super Population', 'PAM Creation', 'Class']
            col_names_sample = ['Sample', 'Sex', 'Population', 'Super Population',  'Targets in Reference',
                                'Targets in Variant', 'Targets in Population', 'Targets in Super Population', 'PAM Creation']
            df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                             guide + '.txt', sep='\t', names=col_names_sample, skiprows=1)
            df = df.sort_values('Targets in Variant', ascending=False)
            df.drop(['Targets in Reference'], axis=1, inplace=True)
            df.drop(['Class'], axis=1, inplace=True)
        more_info_col = []
        for i in range(df.shape[0]):
            more_info_col.append('Show Targets')
        df[''] = more_info_col

        population_1000gp = associateSample.loadSampleAssociation(
            job_directory + 'sampleID.txt')[2]
        super_populations = [{'label': i, 'value': i}
                             for i in population_1000gp.keys()]
        populations = []
        for k in population_1000gp.keys():
            for i in population_1000gp[k]:
                populations.append({'label': i, 'value': i})
        fl.append(
            html.Div
            (
                [
                    html.Div(job_directory + job_id + '.summary_by_samples.' + guide,
                             style={'display': 'none'}, id='div-info-summary_by_sample'),
                    dbc.Row(
                        [
                            dbc.Col(html.Div(dcc.Dropdown(
                                options=super_populations, id='dropdown-superpopulation-sample', placeholder='Select a Super Population'))),
                            dbc.Col(html.Div(dcc.Dropdown(
                                options=populations, id='dropdown-population-sample', placeholder='Select a Population'))),
                            # dbc.Col(html.Div(dcc.Dropdown( id = 'dropdown-sample', placeholder = 'Select a Sample'))),
                            dbc.Col(
                                html.Div(dcc.Input(id='input-sample', placeholder='Select a Sample'))),
                            dbc.Col(html.Div(html.Button(
                                    'Filter', id='button-filter-population-sample')))
                        ]
                    ),
                    dbc.Row(
                        dbc.Col(
                            html.Div(
                                [
                                    html.P('Generating download link, Please wait...',
                                           id='download-link-summary_by_sample'),
                                    dcc.Interval(
                                        interval=1*1000, id='interval-summary_by_sample'),
                                ]
                            )
                        )
                    )
                ],
                style={'width': '50%'}
            )
        )
        fl.append(html.Div('None,None,None', id='div-sample-filter-query',
                           style={'display': 'none'}))  # Folr keep current filter:  Superpop,Pop
        fl.append(html.Div(
            generate_table_samples(df, 'table-samples', 1, guide, job_id), style={'text-align': 'center'}, id='div-table-samples'
        )
        )
        fl.append(
            html.Div(
                [
                    html.Button('Prev', id='prev-page-sample'),
                    html.Button('Next', id='next-page-sample')
                ],
                style={'text-align': 'center'}
            )
        )
        max_page = len(df.index)
        max_page = math.floor(max_page / 10) + 1
        fl.append(html.Div('1/' + str(max_page),
                           id='div-current-page-table-samples'))
        return fl
    elif value == 'tab-summary-by-position':
        # Show Summary by position table
        fl.append(
            html.P(
                'Summary table containing all the targets found in a specific range of positions (chr, start, end) of the genome.')
        )

        fl.append(
            html.P('Filter the table by selecting the chromosome of interest and writing the start and end position of the region to view.')
        )
        # Dropdown chromosomes
        try:
            onlyfile = [f for f in listdir(current_working_directory + 'Genomes/' + genome_selected) if (isfile(join(
                current_working_directory + 'Genomes/' + genome_selected, f)) and (f.endswith('.fa') or f.endswith('.fasta')))]
        except:
            onlyfile = ['chr' + str(i) + '.fa' for i in range(1, 23)]
            onlyfile.append('chrX.fa')
            # NOTE in case no chr in GENOMES/ i put 22 chr + X Y M
            onlyfile.append('chrY.fa')
            onlyfile.append('chrM.fa')
        # removed .fa for better visualization
        onlyfile = [x[:x.rfind('.')] for x in onlyfile]
        chr_file = []
        chr_file_unset = []
        for chr_name in onlyfile:
            chr_name = chr_name.replace('.enriched', '')
            if '_' in chr_name:
                chr_file_unset.append(chr_name)
            else:
                chr_file.append(chr_name)
        chr_file.sort(key=lambda s: [int(t) if t.isdigit(
        ) else t.lower() for t in re.split(r'(\d+)', s)])
        chr_file_unset.sort(key=lambda s: [int(t) if t.isdigit(
        ) else t.lower() for t in re.split(r'(\d+)', s)])
        chr_file += chr_file_unset
        chr_file = [{'label': chr_name, 'value': chr_name}
                    for chr_name in chr_file]

        # Colonne tabella: chr, pos, target migliore, min mm, min bulges, num target per ogni categoria di mm e bulge, show targets; ordine per total, poi mm e poi bulge
        # start_time = time.time()
        # df = pd.read_csv( job_directory + job_id + '.summary_by_position.' + guide +'.txt', sep = '\t')
        # df.rename(columns = {'#Chromosome':'Chromosome'}, inplace = True)
        # more_info_col = []
        # for i in range(df.shape[0]):
        #    more_info_col.append('Show Targets')
        # df[''] = more_info_col
        # TODO inserire failsafe se non ci sono chr, esempio elenco chr da 1 a 22
        fl.append(
            html.Div
            (
                [
                    dbc.Row(
                        [
                            dbc.Col(html.Div(dcc.Dropdown(
                                options=chr_file, id='dropdown-chr-table-position', placeholder='Select a chromosome'))),
                            # dbc.Col(
                            # html.Div(dcc.Input(placeholder='Position', id='input-position'))),
                            dbc.Col(
                                html.Div(dcc.Input(placeholder='Start Position', id='input-position-start'))),
                            dbc.Col(
                                html.Div(dcc.Input(placeholder='End Position', id='input-position-end'))),
                            dbc.Col(html.Div(html.Button(
                                    'Filter', id='button-filter-position'))),
                            html.Br()
                            # )
                        ]
                    ),
                ],
                style={'width': '50%'}
            )
        )
        # print('Position dataframe ready', time.time() - start_time)
        # Folr keep current filter:  chr,pos_start,pos_end
        fl.append(html.Div('None,None,None',
                           id='div-position-filter-query', style={'display': 'none'}))
        # start_time = time.time()
        fl.append(html.Br())
        fl.append(html.Div(
            style={'text-align': 'center'}, id='div-table-position'
        )
        )
        max_page = 1
        fl.append(html.Div('1/' + str(max_page),
                           id='div-current-page-table-position'))
        fl.append(html.Div(mms + '-' + max_bulges,
                           id='div-mms-bulges-position', style={'display': 'none'}))
        return fl
    elif value == 'tab-graphical-sample-card':
        df = pd.read_csv(job_directory + job_id + '.summary_by_samples.' +
                         guide+'.txt', skiprows=1, sep='\t', header=None)
        samples = df.iloc[:, 0]
        fl.append(
            html.P(
                'Summary page containing the single Personal Risk card to be inspected and downloaded')
        )
        fl.append(
            html.Div
            (
                [
                    dbc.Row(
                        [
                            dbc.Col(html.Div(dcc.Dropdown(id='dropdown-sample-card', options=[
                                    {'label': sam, 'value': sam} for sam in samples], placeholder='Select a Sample'))),
                            dbc.Col(html.Div(html.Button(
                                    'Generate', id='button-sample-card'))),
                            dbc.Col(
                                html.Div(id='download-link-personal-card', hidden=True))
                        ]
                    ),
                ],
                style={'width': '50%'}
            )
        )
        fl.append(html.Div(
            [
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(html.Div('', id='div-personal-plot')),
                        dbc.Col(html.Div('', id='div-private-plot'))
                    ]
                )
            ]
        ))
        fl.append(html.Div(
            '', id='div-table-sample-card', style={'text-align': 'center', 'margin-left': '1%', 'margin-right': '1%'}
        )
        )
        fl.append(html.Div(
            '', id='div-top-target-sample-card', style={'text-align': 'center', 'margin-left': '1%', 'margin-right': '1%'}
        )
        )
        fl.append(html.Div('', id='div-sample-card'))
        return fl
    elif value == 'tab-query-table':
        fl.append(html.P(
            'Summary page to query the final result file selecting one/two column to group by the table and extract requested targets'))
        path = current_working_directory+"/Results/"+job_id+"/"+job_id+".db"
        conn = sqlite3.connect(path)
        c = conn.cursor()
        dff_view_names = COL_BOTH
        # dff_view_names = ['Bulge type', 'crRNA', 'Off target motif', 'Reference sequence', 'Chromosome',
        #                   'Position', 'Direction', 'Mismatches',
        #                   'Bulge Size', 'PAM gen', 'Samples', 'SNP',
        #                   'CFD', 'CFD ref', 'Highest CFD Risk Score',
        #                   'AF', 'Annotation Type']
        dff = pd.DataFrame(columns=['Bulge_type', 'crRNA', 'DNA', 'Reference', 'Chromosome',
                                    'Position', 'Direction', 'Mismatches',
                                    'Bulge_Size', 'PAM_gen', 'Samples', 'SNP',
                                    'CFD', 'CFD_ref', 'Highest_CFD_Risk_Score',
                                    'AF', 'Annotation_Type'])  # pd.read_sql_query("SELECT * FROM final_table LIMIT 0",conn) #to define column names in the first empty table
        all_value = {'Target1 :with highest CFD': ['Mismatches', 'Bulge_Size', 'Total', 'CFD', 'Highest_CFD_Risk_Score'],  # , 'Highest_CFD_Absolute_Risk_Score'
                     'Target2 :with lowest Mismatches + Bulge Count': ['Mismatches', 'Bulge_Size', 'Total', 'CFD', 'CFD_Risk_Score']}  # , 'CFD_Absolute_Risk_Score'
    # target_options = {'Mismatches': ['Bulge_Size', 'Total', 'CFD'], 'Bulge_Size': ['Mismatches', 'Total', 'CFD'], 'Total': ['Mismatches', 'Bulge_Size', 'CFD'], 'CFD': [
    #     'Mismatches', 'Bulge_Size', 'Total'], 'Highest_CFD_Risk_Score': [], 'Highest_CFD_Absolute_Risk_Score': [], 'CFD_Risk_Score': [], 'CFD_Absolute_Risk_Score': []}
        all_options = {'Target1 :with highest CFD': [' Mismatches', ' Bulges', ' Mismatch+Bulges', ' CFD', ' Risk Score'],  # , ' Absolute Risk Score'
                       'Target2 :with lowest Mismatches + Bulges Count': [' Mismatches', ' Bulges', ' Mismatch+Bulges', ' CFD', ' Risk Score']}  # , ' Absolute Risk Score'
    # target_o

        # all_options = {'Target1 :with highest CFD': [' Mismatches', ' Bulges', ' Mismatch+Bulges', ' CFD'],
        #                'Target2 :with lowest Mismatches + Bulge Count': [' Mismatches', ' Bulges', ' Mismatch+Bulges', ' CFD']}
        # target_options = {'Mismatches': [' Bulges', ' Mismatch+Bulges', ' CFD'], ' Bulges': [' Mismatches', ' Mismatch+Bulges', ' CFD'], ' Mismatch+Bulges': [
        #     ' Mismatches', ' Mismatch+Bulges', ' CFD'], ' CFD': [' Mismatches', ' Bulges', ' Mismatch+Bulges']}  # parameters to be set as default
        label = [{'label': lab} for lab in all_options.keys()]
        value = [{'value': val} for val in all_value.keys()]
        target_opt = [label, value]
        # try:
        #     img_panel = dbc.Row(  # row with plot
        #         [
        #             dbc.Col(
        #                 [
        #                     html.A(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
        #                         current_working_directory + 'Results/' + job_id + f'/imgs/CRISPRme_top_1000_log_for_main_text_{guide}.png', 'rb').read()).decode()),
        #                         id='top-1000-score', width="80%", height="auto"),
        #                         target="_blank")
        #                 ], width={"size": 10, "offset": 2}
        #             )
        #         ]
        #     )
        # except:
        #     img_panel = dbc.Row(  # row with plot
        #         [
        #             dbc.Col(
        #                 [html.Div()]
        #             )
        #         ]
        #     )
        query_tab_content = html.Div(
            [
                # img_panel,
                dbc.Row(  # row with main group by, secondo group by and thresholds
                    [
                        dbc.Col(  # col0 phantom target select
                            [
                                html.Div(
                                    [
                                        html.H4('Order by'),
                                        dcc.RadioItems(
                                            id='target',
                                            options=target_opt,
                                            # options=[{'label': k, 'value': k} for k in all_options.keys()],
                                            value='Target1 :with highest CFD'
                                        )
                                    ]
                                )
                            ],
                            style={'display': 'none'}
                        ),
                        dbc.Col(  # col1 main group by
                            html.Div(
                                [
                                    html.H4('Group by'),
                                    dcc.RadioItems(
                                        id='order')
                                ]
                            ), width=3
                        ),
                        dbc.Col(  # col2 second group by
                            html.Div(
                                [
                                    html.H4('And group by'),
                                    html.P('First select the left group by value',
                                           id="secondtext"),
                                    dcc.RadioItems(id='multiorder'),
                                ]
                            ), width=3
                        ),
                        dbc.Col(  # select threshold
                            html.Div(
                                [
                                    html.H4(
                                        'Select thresholds'),
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.Div(
                                                        [
                                                            html.H6(
                                                                'Min'),
                                                            dcc.Dropdown(
                                                                id='sholddrop')
                                                        ]
                                                    ),
                                                ]
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Div(
                                                        [
                                                            html.H6(
                                                                'Max'),
                                                            dcc.Dropdown(
                                                                id='maxdrop'
                                                            )
                                                        ]
                                                    )
                                                ]
                                            )
                                        ]
                                    )
                                ]
                            ), width=3
                        ),
                        dbc.Col(
                            html.Div(
                                [
                                    html.H4(
                                        'Select ordering'),
                                    dcc.RadioItems(id='Radio-asc-1',
                                                   options=[
                                                       {'label': ' Ascending',
                                                        'value': 'ASC'},
                                                       {'label': ' Descending',
                                                        'value': 'DESC'}
                                                   ], value='ASC',
                                                   labelStyle={
                                                       'display': 'inline-block', 'margin': '10px'},
                                                   )
                                ]
                            ), width=3
                        )
                    ], justify='center'
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.Div(
                                    html.Button('Submit', id='submit-val', n_clicks=0,
                                                # style={
                                                #     'position': 'absolute', 'center': '50%'}
                                                ),
                                )
                            ],  width={"size": 1}
                        ),
                        dbc.Col(
                            [
                                html.Div(
                                    html.Button('Reset', id='reset-val', n_clicks=0,
                                                # style={'position': 'absolute',
                                                #        'center': '50%'}
                                                )
                                )
                            ],  width={"size": 1, "offset": 1},
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Br(),
                                        html.Hr(),
                                    ]
                                )
                            ]
                        )
                    ], justify='center'
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                # html.Br(),
                                html.P(
                                    'Export will download 1000 lines contained in the current view of the table'),
                                html.Div(dash_table.DataTable(
                                    css=[{'selector': '.row',
                                          'rule': 'margin: 0'}],
                                    export_format="csv",
                                    id='live_table',
                                    columns=[{"name": dff_view_names[count], "id": i, 'hideable':True}
                                             for count, i in enumerate(dff.columns)],
                                    # tooltip_data=[
                                    #     {
                                    #         column: {'value': str(value), 'type': 'markdown'}
                                    #         for column, value in row.items()
                                    #     } for row in dff.to_dict('records')
                                    # ],

                                    # style_cell=dict(textAlign='left'),
                                    # style_header=dict(backgroundColor="white"),
                                    style_data={
                                        'backgroundColor': "white",
                                        # 'whiteSpace': 'normal',
                                        # 'height': 'auto'
                                    },
                                    # , 'overflowX': 'auto'
                                    style_table={
                                        'overflowX': 'scroll', 'overflowY': 'scroll', 'max-height': '300px'},
                                    # style_cell=[{
                                    #     # 'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
                                    #     # 'minWidth': f'{1./len(dff.columns)*100}%', 'width': f'{1./len(dff.columns)*100}%', 'maxWidth': f'{1./len(dff.columns)*100}%'
                                    #     'width': '{}%'.format(len(dff.columns)),
                                    #     # 'whiteSpace': 'normal'
                                    # }],
                                    style_cell_conditional=[
                                        {'if': {'column_id': 'SNP'},
                                         # 'overflow': 'hidden',
                                         'textOverflow': 'ellipsis',
                                         'maxWidth': 200,
                                         },
                                        {'if': {'column_id': 'Samples'},
                                         # 'overflow': 'hidden',
                                         'textOverflow': 'ellipsis',
                                         'maxWidth': 300,
                                         }
                                    ],
                                    # style_cell_conditional=[{'if': {'column_id': 'SNP'},
                                    #                          'overflow': 'hidden',
                                    #                           'textOverflow': 'ellipsis'
                                    #                        }],
                                    #                        {'if': {'column_id':'MMBLG_Samples_2'},
                                    #                          'textAlign': 'left',
                                    #                          'overflow': 'hidden',
                                    #                          'textOverflow': 'ellipsis',
                                    #                          'maxWidth': 100,
                                    #                        },
                                    #                        {'if': {'column_id':'Annotation_Type_1'},
                                    #                          'textAlign': 'left',
                                    #                          'overflow': 'hidden',
                                    #                          'textOverflow': 'ellipsis',
                                    #                          'maxWidth': 200,
                                    #                        },
                                    #                        {'if': {'column_id':'MMBLG_Annotation_Type_2'},
                                    #                          'textAlign': 'left',
                                    #                          'overflow': 'hidden',
                                    #                          'textOverflow': 'ellipsis',
                                    #                          'maxWidth': 200,
                                    #                        },
                                    #                        ],
                                    page_current=0,
                                    page_size=1000,
                                    page_action='custom',
                                    tooltip_delay=0,
                                    tooltip_duration=None
                                ), id='div-query-table'),
                            ],
                        ),
                    ],
                ),
                html.Div(
                    [
                        dbc.Row(
                            dbc.Col(
                                [
                                    dbc.Alert(
                                        "Select a main order before submitting the query",
                                        id="message-alert",
                                        color="danger",
                                        dismissable=True,
                                        fade=True,
                                        is_open=False,
                                        duration=4000,
                                    ),

                                ]
                            ),
                        )
                    ],        style={'display': 'inline-block'}
                )
            ]
        )
        fl.append(query_tab_content)

        # fl.append(

        return fl
    else:  # tab-graphical
        # Show Report images
        samp_style = {}
        if genome_type == 'ref':
            samp_style = {'display': 'none'}

        # fl.append(html.Br())
        fl.append(html.P(
            'Summary Graphical report collecting all the plots and images produced during the search'))

        opt_mm = []
        for i in range(int(mms)+1):
            opt_mm.append({'label': str(i), 'value': str(i)})
        opt_blg = []
        for i in range(int(max_bulges)+1):
            opt_blg.append({'label': str(i), 'value': str(i)})

        # fl.append(html.Br())

        # guide = guide.strip()
        # radar_img = '/imgs/summary_single_guide_' + guide + '_' + str(
        #     0) + '_total_0.0_TOTAL.png'  # TODO choose between 0 mm and max used mms
        # if not os.path.isfile(f"{job_directory}/{radar_img}"):
        #     try:
        #         os.system(f"python {app_main_directory}/PostProcess/generate_img_radar_chart.py {guide} {job_directory}/guide_dict_{guide}.json {job_directory}/motif_dict_{guide}.json 0 0 TOTAL {job_directory}/imgs/")
        #     except:
        #         pass
        # barplot_img = 'summary_histogram_' + \
        #     guide + '_' + str(0) + '_total.png'
        # try:  # NOTE serve per non generare errori se il barplot non è stato fatto
        #     barplot_src = 'data:image/png;base64,{}'.format(base64.b64encode(open(
        #         current_working_directory + 'Results/' + job_id + '/' + barplot_img, 'rb').read()).decode())
        # except:
        #     barplot_src = ''
        # try:
        #     barplot_href = 'assets/Img/' + job_id + '/' + barplot_img
        # except:
        #     barplot_href = ''

        # img_found = False
        # try:
        #     radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(open(
        #         current_working_directory + 'Results/' + job_id + '/' + radar_img, 'rb').read()).decode())
        #     img_found = True
        # except:
        #     radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(
        #         open(current_working_directory+'/assets/placeholder.png', 'rb').read()).decode())
        # try:
        #     radar_href = 'assets/Img/' + job_id + '/' + radar_img
        # except:
        #     radar_href = ''

        if genome_type != 'ref':
            population_1000gp = associateSample.loadSampleAssociation(
                job_directory + 'sampleID.txt')[2]

            super_populations = [{'label': i, 'value': i}
                                 for i in population_1000gp.keys()]
            populations = []
            for k in population_1000gp.keys():
                for i in population_1000gp[k]:
                    populations.append({'label': i, 'value': i})
        else:
            super_populations = []
            populations = []
        # fl.append(html.P('Select Mismatch Value'))

        top1000_image = html.Div(
            html.A(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                            current_working_directory + 'Results/' + job_id + f'/imgs/CRISPRme_top_1000_log_for_main_text_{guide}.png', 'rb').read()).decode()),
                            id='top-1000-score', width="80%", height="auto"),
                   target="_blank")
        )

        total_buttons = [
            dbc.Col(
                html.Div(
                    [
                        html.P("Select Mismatches"),
                        dcc.Dropdown(id='mm-dropdown',
                                     options=opt_mm,
                                     value='0',
                                     clearable=False,
                                     )
                    ]
                ), width=4
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select Bulges"),
                        dcc.Dropdown(id='blg-dropdown',
                                     options=opt_blg,
                                     value='0',
                                     clearable=False,
                                     )
                    ]
                ), width=4
            )
        ]
        sample_buttons = [
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Superpopulation", style=samp_style),
                        html.Div(
                            dcc.Dropdown(
                                options=super_populations,
                                id='dropdown-superpopulation-sample',
                                placeholder='SuperPopulation',
                                style=samp_style),
                        )
                    ]
                ), md=4
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Population", style=samp_style),
                        html.Div(dcc.Dropdown(
                            options=populations,
                            id='dropdown-population-sample',
                            placeholder='Population',
                            style=samp_style),
                        )
                    ]
                ), md=4
            ),
            dbc.Col(
                html.Div(
                    [
                        html.P("Select a Sample", style=samp_style),
                        html.Div(dcc.Dropdown(
                            id='dropdown-sample',
                            placeholder='Sample',
                            style=samp_style),
                        )
                    ]
                ), md=4
            )
        ]
        fl.append(
            html.Div(
                [
                    dbc.Row(
                        dbc.Col(top1000_image, width={"size": 10, "offset": 2})
                    ),
                    dbc.Row(
                        total_buttons, justify='center'
                    ),
                    html.Br()
                ]
            )
        )

        radar_chart_total_content = html.Div(id='div-radar-chart-total')
        populations_barplots = html.Div(id='div-population-barplot')
        # radar_chart_sample_content = html.Div(id='div-radar-chart-sample')
        radar_chart_sample_content = dbc.Row(id='row-radar-chart-sample')
        sample_image_content = html.Div(id='div-sample-image')

        if genome_type != 'ref':
            graph_summary_both = [
                dbc.Col(populations_barplots), dbc.Col(radar_chart_total_content)]
        else:
            graph_summary_both = dbc.Col(radar_chart_total_content, width={
                                         "size": 8, "offset": 2})

        fl.append(
            html.Div(
                [
                    dbc.Row(graph_summary_both)
                ]
            )
        )
        fl.append(
            html.Div(
                [
                    html.Br(),
                    dbc.Row(sample_buttons),
                    radar_chart_sample_content
                ]
            )
        )

        # TODO codice per l'integrazione del CFD graph. When .CFDGraph.txt will be integrated, remove the try/except

        cfd_path = job_directory + job_id + '.CFDGraph.txt'
        if not isfile(cfd_path):  # No file found
            return fl

        fl.extend(
            CFDGraph.CFDGraph(cfd_path)
        )

        return fl
    # guide = all_guides[int(sel_cel[0]['row'])]['State']
        # return guide + value
    raise PreventUpdate

# Read the uploaded file and converts into bit


def parse_contents(contents):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    return decoded
# For filtering


def split_filter_part(filter_part):
    '''
    Preso dal sito di dash sul filtering datatables con python
    '''
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find('{') + 1: name_part.rfind('}')]

                value_part = value_part.strip()
                v0 = value_part[0]
                if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
                    value = value_part[1: -1].replace('\\' + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3

# Perform expensive loading of a dataframe and save result into 'global store'
# Cache are in the Cache directory


@ cache.memoize()
def global_store(value):
    '''
    Caching dei file targets per una miglior performance di visualizzazione
    '''
    if value is None:
        return ''
    target = [f for f in listdir(current_working_directory + 'Results/' + value) if isfile(
        join(current_working_directory + 'Results/'+value, f)) and f.endswith('scores.txt')]
    if not target:
        target = [f for f in listdir(current_working_directory + 'Results/' + value) if isfile(
            join(current_working_directory + 'Results/'+value, f)) and f.endswith('targets.txt')]

    df = pd.read_csv(current_working_directory + 'Results/' +
                     value + '/' + target[0], sep='\t', usecols=range(0, 16))
    df.rename(columns={"#Bulge type": 'BulgeType', '#Bulge_type': 'BulgeType', 'Bulge Size': 'BulgeSize',
                       'Bulge_Size': 'BulgeSize', 'Doench 2016': 'Doench2016', 'Doench_2016': 'Doench2016'}, inplace=True)
    return df


@ app.callback(
    Output('result-table', 'data'),
    [Input('result-table', "page_current"),
     Input('result-table', "page_size"),
     Input('result-table', "sort_by"),
     Input('result-table', 'filter_query')],
    [State('url', 'search'),
     State('url', 'hash')]
)
def update_table(page_current, page_size, sort_by, filter, search, hash_guide):
    '''
    La funzione ritorna uno split dei risultati in base ad un filtering o a un sort da parte dell'utente. Inoltre aggiorna i risultati
    visualizzati quando il bottone next page / prev page è cliccato. (Codice preso dalla pagina dash datatable sul sorting con python)
    Inoltre carica i file targets, o scores se presente, e lo trasforma in un dataframe, cambiando il nome delle colonne per farle corrispondere
    all'id delle colonne della tabella nella pagina.
    Se non ci sono targets ritorna un avviso di errore
    '''
    job_id = search.split('=')[-1]
    job_directory = current_working_directory + 'Results/' + job_id + '/'
    guide = hash_guide.split('#')[1]
    value = job_id
    if search is None:
        raise PreventUpdate

    filtering_expressions = filter.split(' && ')
    # filtering_expressions.append(['{crRNA} = ' + guide])
    df = global_store(value)
    dff = df[df['crRNA'] == guide]

    sort_by.insert(0, {'column_id': 'Mismatches', 'direction': 'asc'})
    sort_by.insert(1, {'column_id': 'BulgeSize', 'direction': 'asc'})
    # sort_by.insert(2, {'column_id': 'CFD', 'direction':'desc'})
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)].sort_values([col['column_id'] for col in sort_by],
                                                                                      ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
                inplace=False)
        elif operator == 'contains':
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == 'datestartswith':
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    # NOTE sort_by: [{'column_id': 'BulgeType', 'direction': 'asc'}, {'column_id': 'crRNA', 'direction': 'asc'}]
    # sort_by.insert(0, {'column_id' : 'Mismatches', 'direction': 'asc'})
    # sort_by.insert(0, {'column_id' : 'BulgeSize', 'direction': 'asc'})
    if len(sort_by):
        dff = dff.sort_values(
            ['Samples' if col['column_id'] == 'Samples Summary' else col['column_id']
                for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    # Check if results are not 0
    warning_no_res = ''
    with open(job_directory + job_id + '.targets.txt') as t:
        no_result = False
        t.readline()
        last_line = t.readline()
        if (last_line == '' or last_line == '\n'):
            no_result = True

    if (no_result):
        warning_no_res = dbc.Alert(
            "No results were found with the given parameters", color="warning")

    return dff.iloc[
        page_current*page_size:(page_current + 1)*page_size
    ].to_dict('records')


# Callbacks for querying part--------------------------------------------------------------


# Return the table with the result of the query
@ app.callback(
    # [Output('live_table', 'data'),
    [Output('live_table', 'data'),
     Output('live_table', 'tooltip_data'),
     Output("message-alert", "is_open"), ],
    [Input('submit-val', 'n_clicks'),
     Input('live_table', "page_current")],
    [State('live_table', "page_size"),
     State('general-profile-table', 'selected_cells'),
     State('target', 'value'),
     State('order', 'value'),
     State('general-profile-table', 'data'),
     State('multiorder', 'value'),
     State('sholddrop', 'value'),
     State('Radio-asc-1', 'value'),
     State('maxdrop', 'value'),
     State('url', 'search'),
     State("message-alert", "is_open"),
     ]
)
def update_output(n_clicks, page_current, page_size, sel_cel, target, radio_order, all_guides, orderdrop, sholddrop, asc1, maxdrop, url, alert):
    guide = all_guides[int(sel_cel[0]['row'])]['Guide']

    target = target[0:7]
    """
    print(target)
    # temporal guide for test file
    print(n_clicks)
    print(type(n_clicks))
    """
    if n_clicks > 0:
        if radio_order == None:
            data = []
            tooltip_data = []
            return data, tooltip_data, not alert
        else:
            if sholddrop != None:
                alert = False
                data = query_manager.shold(target, n_clicks, page_current, page_size, radio_order,
                                           orderdrop, sholddrop, maxdrop, asc1, url, guide, current_working_directory)
            else:
                data = query_manager.noshold(target, n_clicks, page_current, page_size,
                                             radio_order, orderdrop, asc1, url, guide, current_working_directory)
            if target[-1] == '1':

                sub_cols = ['Bulge_type_1', 'crRNA_1', 'DNA_1', 'Reference_1', 'Chromosome_1',
                            'Position_1', 'Direction_1', 'Mismatches_1',
                            'Bulge_Size_1', 'PAM_gen_1', 'Samples_1', 'SNP_1',
                            'CFD_1', 'CFD_ref_1', 'Highest_CFD_Risk_Score_1',
                            'AF_1', 'Annotation_Type_1']

                data = data[sub_cols]
                data.columns = [x[:-2] for x in sub_cols]
                data[' '] = [' '] * data.shape[0]
                data_cols = data.columns.tolist()
                data_cols.remove(' ')
                data_cols.insert(0, ' ')

            snps = pd.DataFrame(data['SNP']).to_dict('records')
            data = data.to_dict('records')
            tooltip_data = [{
                            column: {'value': str(value), 'type': 'markdown'}
                            for column, value in row.items()
                            } for row in snps]
    else:
        raise PreventUpdate

    return data, tooltip_data, alert


# to get correct number of page
@ app.callback(
    Output('live_table', 'page_current'),
    [Input('submit-val', 'n_clicks')]
)
def reset_pagenumber(n):
    if n > 0:
        a = 0
        return a
    else:
        raise PreventUpdate


@ app.callback(
    Output('order', 'options'),
    [Input('target', 'value')])
def set_columns_options(selected_target):
    all_value = {'Target1 :with highest CFD': ['Mismatches', 'Bulge_Size', 'Total', 'CFD', 'Highest_CFD_Risk_Score'],  # , 'Highest_CFD_Absolute_Risk_Score'
                 'Target2 :with lowest Mismatches + Bulge Count': ['Mismatches', 'Bulge_Size', 'Total', 'CFD', 'CFD_Risk_Score']}  # , 'CFD_Absolute_Risk_Score'
    # target_options = {'Mismatches': ['Bulge_Size', 'Total', 'CFD'], 'Bulge_Size': ['Mismatches', 'Total', 'CFD'], 'Total': ['Mismatches', 'Bulge_Size', 'CFD'], 'CFD': [
    #     'Mismatches', 'Bulge_Size', 'Total'], 'Highest_CFD_Risk_Score': [], 'Highest_CFD_Absolute_Risk_Score': [], 'CFD_Risk_Score': [], 'CFD_Absolute_Risk_Score': []}
    all_options = {'Target1 :with highest CFD': [' Mismatches', ' Bulges', ' Mismatch+Bulges', ' CFD', ' Risk Score'],  # , ' Absolute Risk Score'
                   'Target2 :with lowest Mismatches + Bulges Count': [' Mismatches', ' Bulges', ' Mismatch+Bulges', ' CFD', ' Risk Score']}  # , ' Absolute Risk Score'
    # target_options = {' Mismatches': [' Bulges', ' Mismatch+Bulges', ' CFD'], ' Bulges': [' Mismatches', ' Mismatch+Bulges', ' CFD'], ' Mismatch+Bulges': [' Mismatches', ' Bulges', ' CFD'], ' CFD': [
    #     ' Mismatches', ' Bulges', ' Mismatch+Bulges'], ' Risk_Score': [], ' Absolute_Risk_Score': [], ' Risk_Score': [], ' Risk_Score': []}
    # main_order_dict = dict()
    # main_order_dict['label'] = [lab for lab in all_options[selected_target]]
    # main_order_dict['value'] = [val for val in all_value[selected_target]]
    # label = [{'label': lab} for lab in all_options[selected_target]]
    # value = [{'value': val} for val in all_value[selected_target]]
    gi = []
    for count in range(0, len(all_value[selected_target])):
        gi.append({'label': all_options[selected_target][count],
                   'value': all_value[selected_target][count]})
    # gi = [{'label': i, 'value': i} for i in all_options[selected_target]]
    # print(gi)
    # return gi
    # print(main_order_dict)
    return gi


'''
@app.callback(
    [Output('sholddrop', 'value'),
    Output('order', 'value'),
    Output('multiorder', 'value'),
    Output('maxdrop', 'value'),
    Output('Radio-asc-1', 'value')],
    [Input('order', 'options')])
def set_columns_value(available_options):
    return None,None,None,None,None


@app.callback(
    [Output('sholddrop', 'value'),
    Output('multiorder', 'value'),
    Output('maxdrop', 'value'),
    Output('Radio-asc-1', 'value')],
    [Input('multiorder', 'options')])
def set_columns_value(available_options):
    return None,None,None,None
'''

# callback to return the parameters in the various cases


@ app.callback(
    [Output('multiorder', 'options'),
     Output('sholddrop', 'options'),
     Output(component_id='secondtext', component_property='style')],
    [Input('order', 'value')]
)
def set_display_children(selected_order):
    # all_options = {'Target1': ['Mismatches', 'Bulge_Size', 'Total','CFD'],'Target2': ['Mismatches', 'Bulge_Size', 'Total','CFD']}
    target_value = {'Mismatches': ['Bulge_Size', 'Total', 'CFD'], 'Bulge_Size': ['Mismatches', 'Total', 'CFD'], 'Total': ['Mismatches', 'Bulge_Size', 'CFD'], 'CFD': [
        'Mismatches', 'Bulge_Size', 'Total'], 'Highest_CFD_Risk_Score': [], 'Highest_CFD_Absolute_Risk_Score': [], 'CFD_Risk_Score': [], 'CFD_Absolute_Risk_Score': []}
    target_label = {'Mismatches': [' Bulges', ' Mismatch+Bulges', ' CFD'], 'Bulge_Size': [' Mismatches', ' Mismatch+Bulges', ' CFD'], 'Total': [' Mismatches', ' Bulges', ' CFD'],
                    'CFD': [' Mismatches', ' Bulges', ' Mismatch+Bulges'], 'Highest CFD Risk Score': [], 'Highest CFD Absolute Risk Score': [], 'CFD Risk Score': [], 'CFD Absolute Risk Score': []}

    gi = []
    if selected_order is not None:
        for count in range(0, len(target_value[selected_order])):
            gi.append({'label': target_label[selected_order][count],
                       'value': target_value[selected_order][count]})

    if selected_order == None:
        return [], [], {'display': 'block'}
    elif selected_order == "Mismatches":
        data = [{'label': '0', 'value': '0'}, {'label': '1', 'value': '1'}, {'label': '2', 'value': '2'}, {
            'label': '3', 'value': '3'}, {'label': '4', 'value': '4'}, {'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {'display': 'none'}
    elif selected_order == "CFD":
        data = [{'label': '0.001', 'value': '0.001'}, {'label': '0.01', 'value': '0.01'}, {'label': '0.1', 'value': '0.1'}, {'label': '0.2', 'value': '0.2'}, {'label': '0.3', 'value': '0.3'}, {
            'label': '0.4', 'value': '0.4'}, {'label': '0.5', 'value': '0.5'}, {'label': '0.6', 'value': '0.6'}, {'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {'display': 'none'}
    elif selected_order == "Total":
        data = [{'label': '0', 'value': '0'}, {'label': '1', 'value': '1'}, {'label': '2', 'value': '2'}, {'label': '3', 'value': '3'}, {
            'label': '4', 'value': '4'}, {'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}, {'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {'display': 'none'}
    elif selected_order == "Bulge_Size":
        data = [{'label': '1', 'value': '1'}, {'label': '2', 'value': '2'}]
        # return [{'label': i, 'value': i} for i in target_options[selected_order]], data, {'display': 'none'}
        return gi, data, {'display': 'none'}
    else:
        return [], [], {'display': 'none'}


@ app.callback(
    Output('maxdrop', 'options'),
    [Input('sholddrop', 'value'),
        Input('order', 'value')]
)
def maxdrop(sholddrop, order):
    if order == 'Mismatches':
        if sholddrop == '0':
            data = [{'label': '1', 'value': '1'}, {'label': '2', 'value': '2'}, {'label': '3', 'value': '3'}, {
                'label': '4', 'value': '4'}, {'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}]
        elif sholddrop == '1':
            data = [{'label': '2', 'value': '2'}, {'label': '3', 'value': '3'}, {
                'label': '4', 'value': '4'}, {'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}]
        elif sholddrop == '2':
            data = [{'label': '3', 'value': '3'}, {'label': '4', 'value': '4'}, {
                'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}]
        elif sholddrop == '3':
            data = [{'label': '4', 'value': '4'}, {
                'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}]
        elif sholddrop == '4':
            data = [{'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}]
        elif sholddrop == '5':
            data = [{'label': '6', 'value': '6'}]
        else:
            data = []

    elif order == 'CFD':
        if sholddrop == '0.001':
            data = [{'label': '0.01', 'value': '0.01'}, {'label': '0.1', 'value': '0.1'}, {'label': '0.2', 'value': '0.2'}, {'label': '0.3', 'value': '0.3'}, {'label': '0.4', 'value': '0.4'}, {
                'label': '0.5', 'value': '0.5'}, {'label': '0.6', 'value': '0.6'}, {'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.01':
            data = [{'label': '0.1', 'value': '0.1'}, {'label': '0.2', 'value': '0.2'}, {'label': '0.3', 'value': '0.3'}, {'label': '0.4', 'value': '0.4'}, {
                'label': '0.5', 'value': '0.5'}, {'label': '0.6', 'value': '0.6'}, {'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.1':
            data = [{'label': '0.2', 'value': '0.2'}, {'label': '0.3', 'value': '0.3'}, {'label': '0.4', 'value': '0.4'}, {'label': '0.5', 'value': '0.5'}, {
                'label': '0.6', 'value': '0.6'}, {'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.2':
            data = [{'label': '0.3', 'value': '0.3'}, {'label': '0.4', 'value': '0.4'}, {'label': '0.5', 'value': '0.5'}, {
                'label': '0.6', 'value': '0.6'}, {'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.3':
            data = [{'label': '0.4', 'value': '0.4'}, {'label': '0.5', 'value': '0.5'}, {'label': '0.6', 'value': '0.6'}, {
                'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.4':
            data = [{'label': '0.5', 'value': '0.5'}, {'label': '0.6', 'value': '0.6'}, {
                'label': '0.7', 'value': '0.7'}, {'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.5':
            data = [{'label': '0.6', 'value': '0.6'}, {'label': '0.7', 'value': '0.7'}, {
                'label': '0.8', 'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.6':
            data = [{'label': '0.7', 'value': '0.7'}, {'label': '0.8',
                                                       'value': '0.8'}, {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.7':
            data = [{'label': '0.8', 'value': '0.8'},
                    {'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.8':
            data = [{'label': '0.9', 'value': '0.9'}]
        elif sholddrop == '0.9':
            data = [{'label': '1', 'value': '1'}]
        else:
            data = []

    elif order == 'Bulge_Size':
        if sholddrop == '1':
            data = [{'label': '2', 'value': '2'}]
        else:
            data = []

    elif order == 'Total':
        if sholddrop == '0':
            data = [{'label': '1', 'value': '1'}, {'label': '2', 'value': '2'}, {'label': '3', 'value': '3'}, {'label': '4', 'value': '4'}, {
                'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}, {'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '1':
            data = [{'label': '2', 'value': '2'}, {'label': '3', 'value': '3'}, {'label': '4', 'value': '4'}, {
                'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}, {'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '2':
            data = [{'label': '3', 'value': '3'}, {'label': '4', 'value': '4'}, {'label': '5', 'value': '5'}, {
                'label': '6', 'value': '6'}, {'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '3':
            data = [{'label': '4', 'value': '4'}, {'label': '5', 'value': '5'}, {
                'label': '6', 'value': '6'}, {'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '4':
            data = [{'label': '5', 'value': '5'}, {'label': '6', 'value': '6'}, {
                'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '5':
            data = [{'label': '6', 'value': '6'}, {
                'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '6':
            data = [{'label': '7', 'value': '7'}, {'label': '8', 'value': '8'}]
        elif sholddrop == '7':
            data = [{'label': '8', 'value': '8'}]
        else:
            data = []

    else:
        data = []
    return data


'''
@app.callback(
    [Output('sholddrop', 'value'),
    Output('multiorder', 'value'),
    Output('maxdrop', 'value'),
    Output('Radio-asc-1', 'value')],
    [Input('multiorder', 'options')])
def set_columns_value(available_options):
    return None,None,None,None
'''


@ app.callback(
    [Output('order', 'value'),
     Output('maxdrop', 'value'),
     Output('sholddrop', 'value '),
     Output('Radio-asc-1', 'value')],
    [Input('reset-val', 'n_clicks')])
def resetbutton(n_clicks):
    # print(n_clicks)
    # print(type(n_clicks))
    if n_clicks > 0:
        return None, None, None, None
    else:
        return None, None, None, None
