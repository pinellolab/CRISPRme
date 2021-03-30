# IMPORT DASH
import dash_table
import dash_daq as daq
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import dash
import os
from os import listdir
import shutil
from app import current_working_directory, app_main_directory, app
from dash.exceptions import PreventUpdate
import subprocess
from urllib.parse import quote as urlquote
import base64
from os.path import isfile, isdir, join  # for getting directories
import pandas as pd


def get_more_VCF():
    onlydir = [f for f in listdir(current_working_directory + 'VCF')
               if isdir(join(current_working_directory + 'VCF', f))]
    # onlydir = [x.replace('_', ' ') for x in onlydir]
    vcf_dir = []
    for dir in onlydir:
        if 'VCFs_1000_genome_project' not in dir and 'None' not in dir:
            vcf_dir.append(dir)
    return pd.DataFrame(vcf_dir, columns=["Available VCFs"])


def get_personal_genomes():
    onlydir = [f for f in listdir(current_working_directory + 'Genomes')
               if isdir(join(current_working_directory + 'Genomes', f))]
    # onlydir = [x.replace('_', ' ') for x in onlydir]
    gen_dir = []
    for dir in onlydir:
        if '+' not in dir:
            gen_dir.append(dir)

    return pd.DataFrame(gen_dir, columns=["Available genomes"])


def save_file(name, content, directory):
    """Decode and store a file uploaded with Plotly Dash."""
    # check if directory exist, if not, create it
    if not os.path.exists(directory):
        os.makedirs(directory)

    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(directory, name), "wb") as fp:
        fp.write(base64.decodebytes(data))


def uploaded_files(directory):
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(directory):
        path = os.path.join(directory, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files


@app.callback(
    Output("file-annotation-list", "children"),
    [Input("upload-annotation", "filename"),
     Input("upload-annotation", "contents")],
)
def upload_annotation(uploaded_filenames, uploaded_file_contents):
    """Save uploaded files and regenerate the file list."""

    UPLOAD_DIRECTORY = current_working_directory + "annotations/"

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data, UPLOAD_DIRECTORY)

    files = uploaded_files(UPLOAD_DIRECTORY)
    if len(files) == 0:
        return [html.Li("No annotations file found in the directory")]
    else:
        return [html.Li(filename) for filename in files]


@app.callback(
    Output("file-PAM-list", "children"),
    [Input("upload-PAM", "filename"),
     Input("upload-PAM", "contents")],
)
def upload_PAM(uploaded_filenames, uploaded_file_contents):
    """Save uploaded files and regenerate the file list."""

    UPLOAD_DIRECTORY = current_working_directory + "pam/"

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data, UPLOAD_DIRECTORY)

    files = uploaded_files(UPLOAD_DIRECTORY)
    if len(files) == 0:
        return [html.Li("No PAM file found in the directory")]
    else:
        return [html.Li(filename) for filename in files]


@app.callback(
    Output("file-sample-list", "children"),
    [Input("upload-SampleID", "filename"),
     Input("upload-SampleID", "contents"),
     Input('genomes-table-vcf-content', 'active_cell')],
    [State('genomes-table-vcf-content', 'data'),
     State('VCF-name', 'value')],
)
def upload_sampleid(uploaded_filenames, uploaded_file_contents, genome_active_cell, genome_data, name_to_save):
    """Save uploaded files and regenerate the file list."""

    col = genome_active_cell['column_id']
    row = genome_active_cell['row']
    genome_selected = genome_data[row][col]

    UPLOAD_DIRECTORY = current_working_directory + "samplesID/" + \
        str(genome_selected)+'_'+str(name_to_save).strip()

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data, UPLOAD_DIRECTORY)

    files = uploaded_files(UPLOAD_DIRECTORY)
    if len(files) == 0:
        return [html.Li("No PAM file found in the directory")]
    else:
        return [html.Li(filename) for filename in files]


@ app.callback(
    [Output("file-VCF-list", "children"),
     Output("more-VCF-table", "columns"),
     Output("more-VCF-table", "data")],
    [Input("upload-VCF", "filename"),
     Input("upload-VCF", "contents"),
     Input('genomes-table-vcf-content', 'active_cell')],
    [State('genomes-table-vcf-content', 'data'),
     State('VCF-name', 'value')],
)
def upload_VCF(uploaded_filenames, uploaded_file_contents, genome_active_cell, genome_data, name_to_save):
    """Save uploaded files and regenerate the file list."""

    col = genome_active_cell['column_id']
    row = genome_active_cell['row']
    genome_selected = genome_data[row][col]

    UPLOAD_DIRECTORY = current_working_directory + "VCF/" + \
        str(genome_selected)+'_'+str(name_to_save).strip()

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data, UPLOAD_DIRECTORY)

    files = uploaded_files(UPLOAD_DIRECTORY)
    personal_VCF = get_more_VCF()
    if len(files) == 0:
        return [html.Li("No VFCs found in the directory")], [{"name": i, "id": i} for i in personal_VCF.columns], personal_VCF.to_dict('records')
    else:
        return [html.Li(filename) for filename in files], [{"name": i, "id": i} for i in personal_VCF.columns], personal_VCF.to_dict('records')


@app.callback(
    Output('genome-name', 'value'),
    [Input('genomes-table-genome-content', 'active_cell')],
    [State('genomes-table-genome-content', 'data')],
)
def update_genome_name(genome_active_cell, genome_data):
    col = genome_active_cell['column_id']
    row = genome_active_cell['row']
    genome_selected = genome_data[row][col]

    UPLOAD_DIRECTORY = current_working_directory +\
        "Genomes/" + str(genome_selected).strip()

    return genome_selected

# @app.callback(
#     Output('genome-name', 'value'),
#     Output("file-annotation-list", "children"),
#     [Input('genomes-table-genome-content', 'active_cell')],
#     [State('genomes-table-genome-content', 'data')],
# )
# def update_genome_name(genome_active_cell, genome_data):
#     col = genome_active_cell['column_id']
#     row = genome_active_cell['row']
#     genome_selected = genome_data[row][col]

#     UPLOAD_DIRECTORY = current_working_directory +\
#         "Genomes/"+str(genome_selected).strip()

#     files = uploaded_files(UPLOAD_DIRECTORY)
#     if len(files) == 0:
#         return genome_selected, [html.Li("No annotations file found in the directory")]
#     else:
#         return genome_selected, [html.Li(filename) for filename in files]


@ app.callback(
    [Output("file-genome-list", "children"),
     Output("genomes-table-genome-content", "columns"),
     Output("genomes-table-genome-content", "data"),
     Output("genomes-table-vcf-content", "columns"),
     Output("genomes-table-vcf-content", "data")],
    [Input("upload-genome", "filename"),
     Input("upload-genome", "contents")],
    [State('genome-name', 'value')],
)
def upload_genome(uploaded_filenames, uploaded_file_contents, name_to_save):
    """Save uploaded files and regenerate the file list."""

    UPLOAD_DIRECTORY = current_working_directory +\
        "Genomes/"+str(name_to_save).strip()

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, data, UPLOAD_DIRECTORY)

    files = uploaded_files(UPLOAD_DIRECTORY)
    genomes = get_personal_genomes()
    if len(files) == 0:
        return [html.Li("No fasta files found in the directory")], [{"name": i, "id": i} for i in genomes.columns], genomes.to_dict('records'), [{"name": i, "id": i} for i in genomes.columns], genomes.to_dict('records')
    else:
        return [html.Li(filename) for filename in files], [{"name": i, "id": i} for i in genomes.columns], genomes.to_dict('records'), [{"name": i, "id": i} for i in genomes.columns], genomes.to_dict('records')


def genomeAndDictionaryManagement():
    final_list = []

    personal_genomes = get_personal_genomes()
    personal_VCF = get_more_VCF()

    new_genome_content = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H3('Upload Genome'),
                            html.P(
                                'Insert a genome identifier'),
                            dcc.Textarea(id='genome-name', placeholder='hg19', required='required', style={
                                'width': '350px', 'height': '30px', 'font-family': 'monospace', 'font-size': 'large'}),
                            dcc.Upload(
                                id="upload-genome",
                                children=html.Button(
                                    ["Select fasta files (one per chromosome)"], style={'width': '350px'}
                                ),
                                multiple=True,
                            ),
                            html.Br(),
                            html.Div(
                                [
                                    dash_table.DataTable(id="genomes-table-genome-content", style_cell={'textAlign': 'left'},
                                                         columns=[
                                                         {"name": i, "id": i} for i in personal_genomes.columns],
                                                         data=personal_genomes.to_dict(
                                                         'records')
                                                         )
                                ],
                                id='div-genome-table', style={'margin-left': '2%', 'width': '325px'}
                            )
                        ]
                    ),
                    dbc.Col(
                        [
                            html.H4(
                                "List of uploaded files for genome:"),
                            html.Ul(id="file-genome-list")
                        ]
                    )
                ]
            )
        ]
    )
    new_VCF_content = html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H3(
                                'Select a genome to bind the VCFs to'),
                            html.Div(
                                [
                                    dash_table.DataTable(id="genomes-table-vcf-content", style_cell={'textAlign': 'left'},
                                                         columns=[
                                                         {"name": i, "id": i} for i in personal_genomes.columns],
                                                         data=personal_genomes.to_dict(
                                                         'records')
                                                         )
                                ],
                                id='div-browse-annotation', style={'margin-left': '2%', 'width': '325px'}
                            ),
                            html.Br(),
                            html.Br(),
                            html.Div(
                                [
                                    dash_table.DataTable(id="more-VCF-table", style_cell={'textAlign': 'left'},
                                                         columns=[
                                        {"name": i, "id": i} for i in personal_VCF.columns],
                                        data=personal_VCF.to_dict(
                                        'records')
                                    )
                                ],
                                id='div-table-vcf', style={'margin-left': '2%', 'width': '325px'}
                            )
                        ]
                    ),
                    dbc.Col(
                        [
                            html.H3('Upload VCFs'),
                            html.P(
                                'Insert a VCF identifier'),
                            dcc.Textarea(id='VCF-name', placeholder='1000G', style={
                                'width': '350px', 'height': '30px', 'font-family': 'monospace', 'font-size': 'large'}),
                            dcc.Upload(
                                id="upload-VCF",
                                children=html.Button(
                                    ["Select VCFs (one per chromosome)"], style={'width': '350px'}
                                ),
                                multiple=True,
                            ),
                            html.H4("List of uploaded files for VCF:"),
                            html.Ul(id="file-VCF-list"),
                        ]
                    ),
                    dbc.Col(
                        [
                            html.H3('Upload sampleID'),
                            html.P(
                                'Insert a sampleID file'),
                            dcc.Upload(
                                id="upload-SampleID",
                                children=html.Button(
                                    ["Select one Sample ID file"], style={'width': '350px'}
                                ),
                                multiple=True,
                            ),
                            html.H4(
                                "Uploaded file for SampleID:"),
                            html.Ul(id="file-sample-list")
                        ]
                    )
                ]
            )
        ]
    )
    new_PAM_content = html.Div(
        [
            html.H3('Upload PAM'),
            html.P(
                'Insert PAM file (the filename must be in this format e.g. 20bp-NGG-SpCas9.txt)'),
            dcc.Upload(
                id="upload-PAM",
                children=html.Button(
                    ["Select PAM file"]
                ),
                multiple=True,
            ),
            html.H4("List of uploaded PAMs:"),
            html.Ul(id="file-PAM-list"),
        ]
    )
    new_annotation_content = html.Div(
        [
            html.H3('Upload Annotation file'),
            html.P(
                'Insert annotation file (the file must be in bed format, e.g. chr1\t100\t1000\tintron)'),
            dcc.Upload(
                id="upload-annotation",
                children=html.Button(
                    ["Select annotation file"]
                ),
                multiple=True,
            ),
            html.H4("List of uploaded annotations:"),
            html.Ul(id="file-annotation-list"),
        ]
    )

    final_list.append(
        dbc.Tabs(
            [
                dbc.Tab(new_genome_content, label='Genomes',
                        tab_id='genome-tab'),
                dbc.Tab(new_VCF_content,
                        label='Personal Variants', tab_id='vcf-tab'),
                dbc.Tab(new_PAM_content,
                        label='PAM Sequences', tab_id='PAM-tab'),
                dbc.Tab(new_annotation_content,
                        label='Annotation Files', tab_id='annotation-tab')
            ],
            active_tab='genome-tab',
            id='management-tabs'
        )
    )
    final_page = html.Div(
        final_list, style={'margin-left': '1%', 'margin-right': '1%'})
    return final_page
