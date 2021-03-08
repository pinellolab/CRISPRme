#NEW: 
#-sequence input
#-extract guides from sequence
#-General guide result table
#-Spcific table with ordered offtargets for each guide
#-Download results

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import dash_table
from dash.exceptions import PreventUpdate
from os import listdir                      #for getting directories
from os.path import isfile, isdir,join      #for getting directories
import subprocess
import base64                               #for decoding upload content
import io                                   #for decoding upload content
import pandas as pd                         #for dash table
import json                                 #for getting and saving report images list
from os import getcwd
import time                                 #measure time for loading df table
from flask_caching import Cache             #for cache of .targets or .scores
import os
import string                               #for job id
import random                               #for job id
import sys                                  #for sys.exit()
import filecmp                              #check if Params files are equals
import dash_bootstrap_components as dbc
import collections                          #For check if guides are the same in two results
from datetime import datetime               #For time when job submitted
from seq_script import extract_seq, convert_pam

PAGE_SIZE = 10                     #number of entries in each page of the table in view report

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css', dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.title = 'CRISPRitz'
app.config['suppress_callback_exceptions'] = True       #necessary if update element in a callback generated in another callback
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

CACHE_CONFIG = {
    # try 'filesystem' if you don't want to setup redis
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': ('Cache')#os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)
app_location = os.path.dirname(os.path.abspath(__file__)) + '/'
operators = [['ge ', '>='],
             ['le ', '<='],
             ['lt ', '<'],
             ['gt ', '>'],
             ['ne ', '!='],
             ['eq ', '='],
             ['contains ']]     #for filtering

#Dropdown available genomes
onlydir = [f for f in listdir('Genomes') if isdir(join('Genomes', f))]
onlydir = [x.replace('_', ' ') for x in onlydir]
gen_dir = []
for dir in onlydir:
    gen_dir.append({'label': dir, 'value' : dir})

#Dropdown available PAM
onlyfile = [f for f in listdir('pam') if isfile(join('pam', f))]
onlyfile = [x.replace('.txt', '') for x in onlyfile]            #removed .txt for better visualization
pam_file = []
for pam_name in onlyfile:
    if 'NGG' in pam_name:
        pam_file.append({'label':pam_name, 'value':pam_name})
    else:
        pam_file.append({'label': pam_name, 'value' : pam_name, 'disabled':True})

#Dropdown available Variants
onlydir = [f for f in listdir('Variants') if isdir(join('Variants', f))]
var_dir = []
for dir in onlydir:
    var_dir.append({'label': dir, 'value' : dir})

#Available mismatches and bulges
av_mismatches = [{'label': i, 'value': i} for i in range(0, 8)]
av_bulges = [{'label': i, 'value': i} for i in range(0, 6)]
av_guide_sequence = [{'label': i, 'value': i} for i in range(15, 26)]
search_bar = dbc.Row(
    [
        #dbc.Col(dbc.Input(type="search", placeholder="Search")),
        dbc.Col(dbc.NavLink('HOME', active = True, href = 'http://127.0.0.1:8050', className= 'testHover', style = {'text-decoration':'none', 'color':'white', 'font-size':'1.5rem'})),
        dbc.Col(dbc.NavLink('ABOUT', active = True, href = 'http://127.0.0.1:8050', className= 'testHover', style = {'text-decoration':'none', 'color':'white', 'font-size':'1.5rem'})),
        dbc.Col(
            dbc.DropdownMenu(
                children=[
                    dbc.DropdownMenuItem("Github", header=True),
                    dbc.DropdownMenuItem("InfOmics/CRISPRitz", href='https://github.com/InfOmics/CRISPRitz'),
                    dbc.DropdownMenuItem("Pinellolab/CRISPRitz", href='https://github.com/pinellolab/CRISPRitz'),
                ],
                #nav=True,
                in_navbar=True,
                label="Downloads",
                style = {'width': '300px !important' } #'height': '400px !important' 
            ),
        ),
        dbc.Col(dbc.NavLink('CONTACTS', active = True, href = 'http://127.0.0.1:8050', className= 'testHover', style = {'text-decoration':'none', 'color':'white', 'font-size':'1.5rem'}))
    ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)
PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"


navbar = dbc.Navbar(
    [
        html.A(
            # Use row and col to control vertical alignment of logo / brand
            dbc.Row(
                [
                    dbc.Col(html.Img(src=PLOTLY_LOGO, height="30px")),
                    dbc.Col(dbc.NavbarBrand("CRISPRitz Web App", className="ml-2", style = {'font-size': '30px'}))
                ],
                align="center",
                no_gutters=True,
            ),
            href='http://127.0.0.1:8050',
        ),
        dbc.NavbarToggler(id="navbar-toggler"),
        dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
    ],
    color="dark",
    dark=True,
)

#For multipage
app.layout = html.Div([
    navbar,
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content'),
    html.P(id = 'signal', style = {'visibility':'hidden'})
])



#new final_list
final_list = []
final_list.extend([#html.H1('CRISPRitz Web Application'),
    html.Div(children='''
        CRISPRitz is a software package containing 5 different tools dedicated to perform predictive analysis and result assessement on CRISPR/Cas experiments. 
    '''),
    html.P()])

final_list.append(
    html.Div(
        [
            html.P(['Download the offline version here: ', html.A('InfOmics/CRISPRitz', href = 'https://github.com/InfOmics/CRISPRitz', target="_blank"), ' or ', html.A('Pinellolab/CRISPRitz', href = 'https://github.com/pinellolab/CRISPRitz', target="_blank") ])
        ]
    )
)
checklist_div = html.Div(
    [
        dbc.FormGroup(
            [
                dbc.Checkbox(
                    id="checkbox-gecko", className="form-check-input"
                ),
                dbc.Label(
                    #html.P(['Activate Gecko ', html.Abbr('comparison', title ='The results of your test guides will be compared with results obtained from a previous computed analysis on gecko library')]) ,
                    html.P('Compare your results with the Gecko library'),
                    html_for="checkbox-gecko",
                    className="form-check-label",
                ),
                dbc.Checkbox(
                    id="checkbox-ref-comp", className="form-check-input"
                ),
                dbc.Label(
                    #html.P(['Activate Reference genome ', html.Abbr('comparison', title ='The results of your test guides will be compared with the results obtained from a computed analysis on the corresponding reference genome. Note: this may increase computational time')]) ,
                    html.P('Compare your results with the corresponding reference genome'),
                    html_for="checkbox-ref-comp",
                    className="form-check-label",
                ),
                dbc.Checkbox(
                    id="checkbox-example-input", className="form-check-input"
                ),
                dbc.Label(
                    #html.P(['Activate Reference genome ', html.Abbr('comparison', title ='The results of your test guides will be compared with the results obtained from a computed analysis on the corresponding reference genome. Note: this may increase computational time')]) ,
                    html.P('Insert example parameters'),
                    html_for="checkbox-example-input",
                    className="form-check-label",
                ),
                # dbc.Checkbox(
                #     id="checkbox-email", className="form-check-input"
                # ),
                # dbc.Label(
                #     'Notify me by email',
                #     html_for="checkbox-email",
                #     className="form-check-label",
                # )
            ],
            check = True
        )
    ],
    id = 'checklist-test-div'
)

modal = html.Div(
    [
        dbc.Modal(
            [
                dbc.ModalHeader("WARNING! Missing inputs"),
                dbc.ModalBody('The following inputs are missing, please select values before submitting the job', id = 'warning-list'),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close" , className="modal-button")
                ),
            ],
            id="modal",
            centered=True
        ),
    ]
)

tab_guides_content = html.Div(
    [
        html.P([
            'Insert crRNA sequence(s), one per line.', 
            html.P('Sequences must have the same length and be provided without the PAM sequence', id = 'testP') ,
        ],
        style = {'word-wrap': 'break-word'}), 

        dcc.Textarea(id = 'text-guides', placeholder = 'GAGTCCGAGCAGAAGAAGAA\nCCATCGGTGGCCGTTTGCCC', style = {'width':'450px', 'height':'160px', 'font-family':'monospace', 'font-size':'large'}),
        #html.P('Note: a maximum number of 1000 sequences can be provided'),
        dbc.FormText('Note: a maximum number of 1000 sequences can be provided', color = 'secondary')
    ],
    style = {'width':'450px'} #same as text-area
)
tab_sequence_content = html.Div(
    [
        html.P(['Search crRNAs by inserting one or more genomic sequences.', html.P('Chromosome ranges can also be supplied')],
        style = {'word-wrap': 'break-word'}), 

        dcc.Textarea(id = 'text-sequence', placeholder = '>sequence 1\nAAGTCCCAGGACTTCAGAAGagctgtgagaccttggc\n>sequence2\nchr1:11,130,540-11,130,751', style = {'width':'450px', 'height':'160px', 'font-family':'monospace', 'font-size':'large'}),
        #html.P('Note: a maximum number of 1000 sequences can be provided'),
        dbc.FormText('Note: a maximum number of 1000 characters can be provided', color = 'secondary')
    ],
    style = {'width':'450px'} #same as text-area
)
final_list.append(
    html.Div(
        html.Div(
            [
                modal,
                html.Div(
                    [
                        html.H3('STEP 1', style = {'margin-top':'0'}), 
                        html.P('Select a genome'),
                        html.Div(
                            dcc.Dropdown(options = gen_dir, clearable = False, id = "available-genome",) #style = {'width':'75%'})
                        ),
                        dbc.FormText('Note: Genomes enriched with variants are indicated with a \'+\' symbol', color='secondary'),
                        # html.P('Add a genome variant', style = {'visibility':'hidden'}),
                        # html.Div(
                        #     dcc.Dropdown(options = var_dir,clearable = False, id = 'available-variant', style = {'width':'75%', 'visibility':'hidden'})
                        # ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.P('Select PAM'),
                                        html.Div(
                                            dcc.Dropdown(options = pam_file, clearable = False, id = 'available-pam')
                                        )
                                    ],
                                    style = {'flex':'0 0 50%', 'margin-top': '10%'}
                                ),
                                #html.P('or'),
                                # html.Div(
                                #     [
                                #         html.P('Insert custom PAM'),
                                #         dcc.Input(type = 'text', id = 'custom-pam', placeholder = 'NGG', disabled = True)
                                #     ]
                                # )
                            ],
                            id = 'div-pam',
                            className = 'flex-div-pam'
                        ),
                        html.Div(
                            [
                                # html.P(#'Send us a request to add a specific genome sequence or a variant, or download the offline version'
                                # [html.A('Contact us', href = 'http://127.0.0.1:8050', target="_blank"),' to request new genomes availability in the dropdown list', html.P('or'), html.P('Download the offline version'),], style = {'margin-top':'10px', 'text-align':'-webkit-center', 'position': 'relative', 'top': '25%'}),
                                html.Ul(
                                    [html.Li(
                                        [html.A('Contact us', href = 'http://127.0.0.1:8050', target="_blank"),' to request new genomes availability in the dropdown list'],
                                        style = {'margin-top':'5%'}
                                    ),
                                    html.Li(
                                        [html.A('Download', href = 'https://github.com/InfOmics/CRISPRitz'), ' the offline version for more custom parameters']
                                    )
                                    ],
                                    style = {'list-style':'inside'}
                                ),
                                html.Div(
                                    html.Button('Insert example parameters', id = 'example-parameters', style={'display':'inline-block'}),
                                    style = {'text-align':'center'}
                                )
                            ],
                            style = {'height':'50%'}
                        ),
                        
                    ],
                    id = 'step1',
                    style = {'flex':'0 0 30%', 'tex-align':'center'}
                ),
                html.Div(style = {'border-right':'solid 1px white'}),
                html.Div(
                    [
                        html.H3('STEP 2', style = {'margin-top':'0'}),
                        html.Div(
                            [
                                html.Div(
                                    [   html.P('Select the input type'),
                                        # dbc.RadioItems(
                                        #     options=[
                                        #         {'label': 'Guides', 'value': 'Guides'},
                                        #         {'label': 'Sequence', 'value': 'Sequence'}
                                        #     ],
                                        #     value='Guides',
                                        #     inline = True,
                                        #     id = 'guide-or-sequence'
                                        # ),
                                        dbc.Tabs(
                                            [
                                                dbc.Tab(tab_guides_content, label='Guides', tab_id= 'guide-tab'),
                                                dbc.Tab(tab_sequence_content, label='Sequence', tab_id = 'sequence-tab')
                                            ],
                                            active_tab='guide-tab',
                                            id = 'tabs'
                                        ),
                                        
                                        # html.Div(
                                        #     [
                                        #         html.P([
                                        #             'Insert crRNA sequence(s), one per line.', 
                                        #             html.P('Sequences must have the same length and be provided without the PAM sequence') ,
                                        #             #html.Abbr('\uD83D\uDEC8', style = {'text-decoration':'none'} ,title = 'One sequence per line. All sequences must have the same lenght and PAM characters are not required')
                                        #         ],
                                        #         style = {'word-wrap': 'break-word'}), 
                            
                                        #         dcc.Textarea(id = 'text-guides', placeholder = 'GAGTCCGAGCAGAAGAAGAA\nCCATCGGTGGCCGTTTGCCC', style = {'width':'450px', 'height':'160px'}),
                                        #         #html.P('Note: a maximum number of 1000 sequences can be provided'),
                                        #         dbc.FormText('Note: a maximum number of 1000 sequences can be provided', color = 'secondary')
                                        #     ],
                                        #     style = {'width':'450px'} #same as text-area
                                        # )
                                    ],
                                    id = 'div-guides'
                                ),
                                html.Div(
                                    [
                                        html.P('Allowed mismatches'),
                                        dcc.Dropdown(options = av_mismatches, clearable = False, id = 'mms', style = {'width':'60px'}),
                                        html.P('Bulge DNA size'),
                                        dcc.Dropdown(options = av_bulges, clearable = False, id = 'dna', style = {'width':'60px'}),
                                        html.P('Bulge RNA size'),
                                        dcc.Dropdown(options = av_bulges, clearable = False, id = 'rna', style = {'width':'60px'}),
                                        dbc.Fade(
                                            [
                                                html.P('crRNA length (without PAM)'),
                                                dcc.Dropdown(options = av_guide_sequence, clearable = False, id = 'len-guide-sequence-ver', style = {'width':'60px'})
                                            ],
                                            id = 'fade-len-guide', is_in= False, appear= False
                                        )
                                    ]
                                )
                            ],
                            className = 'flex-step2'
                        )

                    ],
                    id = 'step2',
                    style = {'flex':'0 0 40%'}
                    
                ),
                html.Div(style = {'border-right':'solid 1px white'}),
                html.Div(
                    [
                        html.H3('Advanced Options'),
                        checklist_div,
                        dcc.Checklist(
                            options = [
                            #{'label':'Gecko comparison', 'value':'GC', 'disabled':False},
                            #{'label':'Reference genome comparison', 'value':'RGC', 'disabled':False},
                            {'label':'Notify me by email','value':'email', 'disabled':False}], 
                            id = 'checklist-advanced',
                        ),
                        dbc.Fade(
                            [
                                dbc.FormGroup(
                                    [
                                        dbc.Label("Email", html_for="example-email"),
                                        dbc.Input(type="email", id="example-email", placeholder="Enter email", className='exampleEmail'),
                                        # dbc.FormText(
                                        #     "Are you on email? You simply have to be these days",
                                        #     color="secondary",
                                        # ),
                                    ]
                                )
                            ],
                            id = 'fade', is_in= False, appear= False
                        ),
                        #html.H3('Submit', style = {'margin-top':'0'}),
                        html.Div(
                            [
                                html.Button('Submit', id = 'check-job'),
                                html.Button('', id = 'submit-job', style = {'visibility':'hidden'})
                            ],
                            style = {'display':'inline-block', 'margin':'0 auto'}   #style="height:55px; width:150px"
                        )
                    ],
                    id = 'step3',
                    style = {'tex-align':'center'},
                    className = 'flex-step3'
                )
            ],
            id = 'div-steps',
            style = {'margin':'1%'},
            className = 'flex-div-steps'
        ),
        style = {'background-color':'rgba(154, 208, 150, 0.39)', 'border-radius': '10px', 'border':'1px solid black'},
        id = 'steps-background'
    )
)
index_page = html.Div(final_list, style = {'margin':'1%'})

#Load Page
final_list = []
#final_list.append(html.H1('CRISPRitz Web Application'))
final_list.append(
    html.Div(
        html.Div(
            html.Div(
                [
                    html.P('Job submitted. Copy this link to view the status and the result page '),
                    html.Div(
                        html.P('link', id = 'job-link'),
                        style = {'border':'2px solid', 'border-color':'blue' ,'width':'70%','display':'inline-block', 'margin':'5px'}
                    )
                ],
                style = {'display':'inline-block'}
            ),
            style = {'display':'inline-block','background-color':'rgba(154, 208, 150, 0.39)', 'border-radius': '10px', 'border':'1px solid black', 'width':'70%'}
        ),
        style = {'text-align':'center'}
    )
)

final_list.append(
    html.Div(
        [
            html.H4('Status report'),
            html.Div(
                [
                    html.Div(
                        html.Ul(
                            [
                                html.Li('Searching crRNA'),
                                html.Li('Annotating result'),
                                html.Li('Generating report')
                            ]
                        ),
                        style = {'flex':'0 0 20%'}
                    ),
                    html.Div(
                        html.Ul(
                            [
                                html.Li('To do', style = {'color':'red'}, id = 'search-status'),
                                html.Li('To do', style = {'color':'red'}, id = 'annotate-result-status'),
                                html.Li('To do', style = {'color':'red'}, id = 'generate-report-status')
                            ],
                            style = {'list-style-type':'none'}
                        )
                    )
                ],
                className = 'flex-status'
            ),
            html.Div(
                [
                    dcc.Link('View Results', style = {'visibility':'hidden'}, id = 'view-results'),
                    html.Div(id = 'no-directory-error')
                ]
            )
        ],
        id = 'div-status-report'
    )
)

final_list.append(html.P('', id = 'done'))

final_list.append(dcc.Interval(id = 'load-page-check', interval=3*1000))
load_page = html.Div(final_list, style = {'margin':'1%'})


#Test bootstrap page
final_list = []
final_list.append(
    html.Div(
    [
        html.P('Test P', id= 'test-P'),
        html.Button('AA', id = 'button-test1'),
        html.Button('BB', id = 'button-test2'),
        html.Button('CC', id = 'button-test3'),
        html.Button('CC10', id = 'button-test10'),
        html.Button ('Gen callback', id = 'gen-callback')
    ]
)
)
final_list.append(html.Div(id='test-div-for-button'))

test_page = html.Div(final_list, style = {'margin':'1%'})
##################################################CALLBACKS##################################################
@app.callback(
    Output('test-P', 'children'), [Input('button-test1', 'n_clicks')]
)
def test1(n):
    if n is None:
        raise PreventUpdate
    print('btn1')
    return ''

@app.callback(
    Output('test-div-for-button', 'children'), [Input('button-test10', 'n_clicks')]
)
def test1(n):
    if n is None:
        raise PreventUpdate
    print('btn10')
    return ''

#################################################
#Fade in/out email
@app.callback(
    Output("fade", "is_in"),
    [Input("checklist-advanced", "value")],
    [State("fade", "is_in")],
)
def toggle_fade(selected_options, is_in):
    if  selected_options is None:
        return False
    if 'email' in selected_options:
        return True
    return False

#Insert/Delete example input
@app.callback(
    [Output('available-genome', 'value'),
    Output('available-pam', 'value'),
    Output('text-guides', 'value'),
    Output('mms', 'value'),
    Output('dna', 'value'),
    Output('rna', 'value'),
    Output('len-guide-sequence-ver', 'value'),
    Output('text-sequence','value')],
    [Input('example-parameters', 'n_clicks')]
)
def inDelExample(n):
    if n is None:
        raise PreventUpdate
    
    #TODO mettere esempio già svolto per non caricare troppo il server
    return gen_dir[0]['value'], '5\'-NGG-3\'', 'GAGTCCGAGCAGAAGAAGAA', '4', '0', '0', '20','>sequence\nTACCCCAAACGCGGAGGCGCCTCGGGAAGGCGAGGTGGGCAAGTTCAATGCCAAGCGTGACGGGGGA'

#Email validity
@app.callback(
    Output('example-email', 'style'),
    [Input('example-email', 'value')]
)
def checkEmailValidity(val):
    if val is None:
        raise PreventUpdate

    if '@' in val:
        return {'border':'1px solid #94f033', 'outline':'0'}
    return {'border':'1px solid red'}

#Fade in guide len dropdown for sequence tabs version
@app.callback(
    Output('fade-len-guide', 'is_in'),
    [Input('tabs', 'active_tab')],
    [State('fade-len-guide', 'is_in')]
)
def resetTab(current_tab, is_in):
    if current_tab is None:
        raise PreventUpdate

    if current_tab == 'guide-tab':
        return False
    return True


#Check input presence
@app.callback(
    [Output('submit-job', 'n_clicks'),
    Output('modal', 'is_open'),
    Output('available-genome', 'className'),
    Output('available-pam', 'className'),
    Output('text-guides', 'style'),
    Output('mms', 'className'),
    Output('dna', 'className'),
    Output('rna', 'className'),
    Output('len-guide-sequence-ver', 'className'),
    Output('warning-list', 'children')],
    [Input('check-job','n_clicks'),
    Input('close','n_clicks')],
    [State('available-genome', 'value'),
    State('available-pam','value'),
    State('text-guides', 'value'),
    State('mms','value'),
    State('dna','value'),
    State('rna','value'),
    State('len-guide-sequence-ver','value'),
    State('tabs','active_tab'),
    State("modal", "is_open")]
)
def checkInput(n, n_close, genome_selected, pam, text_guides, mms, dna, rna, len_guide_seq, active_tab ,is_open):
    if n is None:
        raise PreventUpdate
    if is_open is None:
        is_open = False
    
    classname_red = 'missing-input'
    genome_update = None
    pam_update = None
    text_update = {'width':'450px', 'height':'160px'}
    mms_update = None
    dna_update = None
    rna_update = None
    len_guide_update = None
    update_style = False
    miss_input_list = []
    
    if genome_selected is None or genome_selected is '':
        genome_update = classname_red
        update_style = True
        miss_input_list.append('Genome')
    if pam is None or pam is '':
        pam_update = classname_red
        update_style = True
        miss_input_list.append('PAM')
    # if text_guides is None or text_guides is '':
        # text_update = {'width':'450px', 'height':'160px','border': '1px solid red'}
        # update_style = True
        # miss_input_list.append('crRNA sequence(s)')
    if mms is None or str(mms) is '':
        mms_update = classname_red
        update_style = True
        miss_input_list.append('Allowed Mismatches')
    if dna is None or str(dna) is '':
        dna_update = classname_red
        update_style = True
        miss_input_list.append('Bulge DNA size')
    if rna is None or str(rna) is '':
        rna_update = classname_red
        update_style = True
        miss_input_list.append('Bulge RNA size')
    if (len_guide_seq is None or str(len_guide_seq) is '') and ('sequence-tab' in active_tab):
        len_guide_update = classname_red
        update_style = True
        miss_input_list.append('crRNA length')
    miss_input = html.Div(
        [
            html.P('The following inputs are missing:'),
            html.Ul([html.Li(x) for x in miss_input_list]),
            html.P('Please fill in the values before submitting the job')
        ]
    )
    
    if not update_style:
        return 1, False, genome_update, pam_update, text_update, mms_update, dna_update, rna_update, len_guide_update, miss_input
    return None, not is_open, genome_update, pam_update, text_update, mms_update, dna_update, rna_update, len_guide_update, miss_input

#Submit Job, change url
@app.callback(
    [Output('url', 'pathname'),
    Output('url','search')],
    [Input('submit-job','n_clicks')],
    [State('url', 'href'),
    State('available-genome', 'value'),
    State('available-pam','value'),
    State('text-guides', 'value'),
    State('mms','value'),
    State('dna','value'),
    State('rna','value'),
    State('checkbox-gecko','checked'),
    State('checkbox-ref-comp', 'checked'),
    State('checklist-advanced', 'value'),
    State('example-email','value'),
    State('tabs','active_tab'),
    State('text-sequence','value'),
    State('len-guide-sequence-ver', 'value')]
)
def changeUrl(n, href, genome_selected, pam, text_guides, mms, dna, rna, gecko_opt, genome_ref_opt, adv_opts,dest_email, active_tab, text_sequence, len_guide_sequence):      #NOTE startJob
    '''
    genome_selected can be Human genome (hg19), or Human Genome (hg19) + 1000 Genome Project, the '+' character defines the ref or enr version.
    Note that pam parameter can be 5'-NGG-3', but the corresponding filename is 5'-NGG-3'.txt
    Pam file (5'-NGG-3'.txt) is structured as NGG 3, or TTTN -4. The created pam.txt inside the result directory add the corresponding N's
    Annotations path file is named genome_name_annotationpath.txt, where genome_name is the reference genome name
    '''
    if n is None:
        raise PreventUpdate
    
    #Check input, else give simple input
    if genome_selected is None or genome_selected is '':
        genome_selected = 'hg19_ref'
    if pam is None or pam is '':
        pam = '5\'-NGG-3\''
    if text_guides is None or text_guides is '':
        text_guides = 'GAGTCCGAGCAGAAGAAGAA'
    else:
        text_guides = text_guides.strip()
        if len(text_guides.split('\n')) > 1000:
            text_guides = '\n'.join(text_guides.split('\n')[:1000]).strip()
        if ( not all(len(elem) == len(text_guides.split('\n')[0]) for elem in text_guides.split('\n'))):
            text_guides = selectSameLenGuides(text_guides)
    if (len_guide_sequence is None or str(len_guide_sequence) is '') and ('sequence-tab' in active_tab):
        len_guide_sequence = 20
    if (text_sequence is None or text_sequence is '') and ('sequence-tab' in active_tab):
        text_sequence = '>sequence\nTACCCCAAACGCGGAGGCGCCTCGGGAAGGCGAGGTGGGCAAGTTCAATGCCAAGCGTGACGGGGGA'

    job_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k = 10))
    result_dir = 'Results/' + job_id
    subprocess.run(['mkdir ' + result_dir], shell = True)
    
    search_index = True
    search = True
    annotation = True
    report =  True
    gecko_comp = False
    ref_comparison = False
    send_email = False
    if adv_opts is None:
        adv_opts = []
    if gecko_opt:
        gecko_comp = True
    if genome_ref_opt:
        ref_comparison = True
    if 'email' in adv_opts and dest_email is not None and len(dest_email.split('@')) > 1 and dest_email.split('@')[-1] is not '':
        send_email = True
        with open(result_dir + '/email.txt', 'w') as e:
            e.write(dest_email + '\n')
            e.write('http://127.0.0.1:8050/load?job=' + job_id + '\n')
            e.write(datetime.utcnow().strftime("%m/%d/%Y, %H:%M:%S") + '\n')
            e.write('Job done. Parameters: etc etc')
            e.close()
    
    
    #Set parameters
    genome_selected = genome_selected.replace(' ', '_')
    genome_ref = genome_selected.split('+')[0]              #+ char to separate ref and vcf, eg Human_genome+1000_genome_project
    if genome_ref == genome_selected:
        ref_comparison = False
    #NOTE Indexed genomes names are PAM + _ + bMax + _ + genome_selected
    
    pam_len = 0
    custom_pam = None

    with open('pam/' + pam + '.txt') as pam_file:
        pam_char = pam_file.readline()
        index_pam_value = pam_char.split(' ')[-1]
        if int(pam_char.split(' ')[-1]) < 0:
            end_idx = int(pam_char.split(' ')[-1]) * (-1)
            pam_char = pam_char.split(' ')[0][0 : end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(' ')[-1])
            pam_char = pam_char.split(' ')[0][end_idx * (-1):]
            pam_len = end_idx
            pam_begin = False
    
    if 'sequence-tab' in active_tab:
        #Extract sequence and create the guides
        guides = []
        for name_and_seq in text_sequence.split('>'):
            if '' == name_and_seq:
                continue
            name, seq = name_and_seq.strip().split('\n')
            if 'chr' in seq:
                extracted_seq = extract_seq.extractSequence(name, seq, genome_selected.replace(' ', '_'))
            else:
                extracted_seq = seq.strip()
            guides.extend(convert_pam.getGuides(extracted_seq, pam_char, len_guide_sequence))
            text_guides = '\n'.join(guides).strip()
    

    len_guides = len(text_guides.split('\n')[0])
    if (pam_begin):
        pam_to_file = pam_char + ('N' * len_guides) + ' ' + index_pam_value
    else:
        pam_to_file = ('N' * len_guides) + pam_char + ' ' + index_pam_value

    save_pam_file = open(result_dir + '/pam.txt', 'w')
    save_pam_file.write(pam_to_file)
    save_pam_file.close()
    pam = result_dir + '/pam.txt'
        
    guides_file = result_dir + '/guides.txt'
    if text_guides is not None and text_guides is not '':
        save_guides_file = open(result_dir + '/guides.txt', 'w')
        if (pam_begin):
            text_guides = 'N' * pam_len + text_guides.replace('\n', '\n' + 'N' * pam_len)
        else:
            text_guides = text_guides.replace('\n', 'N' * pam_len + '\n') + 'N' * pam_len
        save_guides_file.write(text_guides)
        save_guides_file.close()     

    if (int(dna) == 0 and int(rna) == 0):
        search_index = False
    max_bulges = rna
    if (int(dna) > int(rna)):
        max_bulges = dna

    if (search_index):
        search = False

    if int(max_bulges) <= 2:
        genome_idx = pam_char + '_' + '2' + '_' + genome_selected
    else:
        genome_idx = pam_char + '_' + '5' + '_' + genome_selected
    genome_idx_ref = genome_idx.split('+')[0]

    #Create Params.txt file
    with open(result_dir + '/Params.txt', 'w') as p:            #NOTE if modified, chenge also mms value in update_table function
        p.write('Genome_selected\t' + genome_selected + '\n')         
        p.write('Genome_ref\t' + genome_ref + '\n')
        if search_index:
            p.write('Genome_idx\t' + genome_idx + '\n')
        else:
            p.write('Genome_idx\t' + 'None\n')
        p.write('Pam\t' + pam_char + '\n')
        p.write('Max_bulges\t' + str(max_bulges) + '\n')
        p.write('Mismatches\t' + str(mms) + '\n')
        p.write('DNA\t' + str(dna) + '\n')
        p.write('RNA\t' + str(rna) + '\n')
        p.write('Gecko\t' + str(gecko_comp) + '\n')
        p.write('Ref_comp\t' + str(ref_comparison) + '\n')
        p.close()

    #Check if input parameters (mms, bulges, pam, guides, genome) are the same as a previous search
    all_result_dirs = [f for f in listdir('Results') if isdir(join('Results', f))]
    all_result_dirs.remove(job_id)
    #all_result_dirs.remove('test')
    for check_param_dir in all_result_dirs:
        if os.path.exists('Results/' + check_param_dir + '/Params.txt'):
            if os.path.exists('Results/' + check_param_dir + '/log.txt'):
                with open('Results/' + check_param_dir + '/log.txt') as log:
                    if ('Job\tDone' in log.read()):
                        if (filecmp.cmp('Results/' + check_param_dir + '/Params.txt', result_dir + '/Params.txt' )):
                                guides1 = open('Results/' + check_param_dir + '/guides.txt').read().split('\n')
                                guides2 = open('Results/' + job_id + '/guides.txt').read().split('\n')
                                if (collections.Counter(guides1) == collections.Counter(guides2)):
                                    search = False
                                    search_index = False
                                    subprocess.run(['cp $PWD/Results/' + check_param_dir + '/' + check_param_dir + '* ' + result_dir + '/'], shell = True)
                                    subprocess.run(['cp $PWD/Results/' + check_param_dir + '/*.png ' + result_dir + '/'], shell = True)
                                    subprocess.run(['rename \'s/' + check_param_dir + '/' + job_id + '/g\' ' + result_dir + '/*'], shell = True)
                                    break           
    
    #Annotation
    if (not search and not search_index):
        annotation = False      

    #Generate report
    if (not search and not search_index):
        report = False         
    
    annotation_filepath = [f for f in listdir('./') if isfile(join('./', f)) and f.startswith(genome_ref)]

    
    subprocess.Popen(['assets/./submit_job.sh ' + 'Results/' + job_id + ' ' + 'Genomes/' + genome_selected + ' ' + 'Genomes/' + genome_ref + ' ' + 'genome_library/' + genome_idx + (
        ' ' + pam + ' ' + guides_file + ' ' + str(mms) + ' ' + str(dna) + ' ' + str(rna) + ' ' + str(search_index) + ' ' + str(search) + ' ' + str(annotation) + (
            ' ' + str(report) + ' ' + str(gecko_comp) + ' ' + str(ref_comparison) + ' ' + 'genome_library/' + genome_idx_ref + ' ' + str(send_email) + ' ' + annotation_filepath[0]
        )
    )], shell = True)
    return '/load','?job=' + job_id

#When url changed, load new page
@app.callback(
    [Output('page-content', 'children'),
    Output('job-link', 'children')],
    [Input('url', 'href'), Input('url','pathname'), Input('url','search')],[State('url','hash')]
    # [State('url','pathname'), 
    # State('url','search')]
)
def changePage( href, path, search, hash_guide):
    # print('href', href)
    # print('hash', hash_guide)
    # print('pathname', path)
    # print('search', search)
    print('hash', hash_guide)
    if path == '/load':
        return load_page, 'http://127.0.0.1:8050/load' + search #NOTE change the url part when DNS are changed
    if path == '/result':
        job_id = search.split('=')[-1]
        if hash_guide is None or hash_guide is '':
            return resultPage(job_id), 'http://127.0.0.1:8050/load' + search
        return guidePage(job_id, hash_guide.split('#')[1]), 'http://127.0.0.1:8050/load' + search
    if path == '/test-page':
        return test_page, 'http://127.0.0.1:8050/load' + search
   
    return index_page, ''

#Check end job
@app.callback(
    [Output('view-results', 'style'),
    Output('annotate-result-status', 'children'),
    Output('search-status', 'children'),
    Output('generate-report-status', 'children'),
    Output('view-results','href'),
    Output('no-directory-error', 'children')],
    [Input('load-page-check', 'n_intervals')],
    [State('url', 'search')]
)
def refreshSearch(n, dir_name):
    if n is None:
        raise PreventUpdate     #TODO fa un controllo subito, così l'utente non deve aspettare 3 secondi per l'update
    
    onlydir = [f for f in listdir('Results') if isdir(join('Results', f))]
    current_job_dir = 'Results/' +  dir_name.split('=')[-1] + '/'
    if dir_name.split('=')[-1] in onlydir:
        onlyfile = [f for f in listdir(current_job_dir) if isfile(join(current_job_dir, f))]
        if 'log.txt' in onlyfile:
            with open(current_job_dir + 'log.txt') as log:
                all_done = 0
                annotate_res_status = html.P('To do', style = {'color':'red'})
                search_status = html.P('To do', style = {'color':'red'})
                report_status = html.P('To do', style = {'color':'red'})
                current_log = log.read()
                if ('Annotation\tDone' in current_log):
                    annotate_res_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                if ('Search-index\tDone' in current_log or 'Search\tDone' in current_log):
                    search_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                if ('Report\tDone' in current_log):
                    report_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                if all_done == 3:
                    return {'visibility':'visible'}, annotate_res_status, search_status, report_status, '/result?job=' + dir_name.split('=')[-1], ''
                else:
                    return {'visibility':'hidden'}, annotate_res_status, search_status, report_status,'', ''
    return {'visibility':'hidden'}, html.P('Not available', style = {'color':'red'}), html.P('Not available', style = {'color':'red'}), html.P('Not available', style = {'color':'red'}), '', dbc.Alert("The selected result does not exist", color = "danger")

#Perform expensive loading of a dataframe and save result into 'global store'
#Cache are in the Cache directory
@cache.memoize()
def global_store(value):
    
    if value is None:
        return ''
    target = [f for f in listdir('Results/' + value) if isfile(join('Results/'+value, f)) and f.endswith('scores.txt') ]
    if not target:
        target = [f for f in listdir('Results/' + value) if isfile(join('Results/'+value, f)) and f.endswith('targets.txt') ]
    
    df = pd.read_csv('Results/' +value + '/' + target[0], sep = '\t')
    df.rename(columns = {"#Bulge type":'BulgeType', '#Bulge_type':'BulgeType','Bulge Size': 'BulgeSize', 'Bulge_Size': 'BulgeSize', 'Doench 2016':'Doench2016','Doench_2016':'Doench2016'}, inplace = True)
    return df

# #Callback to populate the tab, note that it's called when the result_page is loaded (dash implementation), so we do not use raise update to block this first callback
# @app.callback(
#     [Output('signal','children'),
#     Output('result-table','page_current'),
#     Output('result-table', "sort_by"), 
#     Output('result-table','filter_query')],
#     [Input('url', 'pathname')],
#     [State('url', 'search')]
# )
# def populateTable(pathname, search):
#     print('pathname', pathname)
#     if pathname != '/result':
#         raise PreventUpdate

#     job_id = search.split('=')[-1]
#     job_directory = 'Results/' + job_id + '/'
#     print('job dir', job_directory)
#     if(not isdir(job_directory)):
#         return 'not_exists', 0, [], ''
#     #global_store(job_id)
#     print('ok')
#     return job_id, 0, [], ''

#Send the data when next or prev button is clicked on the result table
@app.callback(
    Output('result-table', 'data'),
    [Input('result-table', "page_current"),
     Input('result-table', "page_size"),
     Input('result-table', "sort_by"),
     Input('result-table', 'filter_query')],
     [State('url', 'search'),
     State('url', 'hash')]
)
def update_table(page_current, page_size, sort_by, filter, search, hash_guide):
    job_id = search.split('=')[-1]
    job_directory = 'Results/' + job_id + '/'
    guide = hash_guide.split('#')[1]
    value = job_id
    if search is None:
        raise PreventUpdate

        
    filtering_expressions = filter.split(' && ')
    #filtering_expressions.append(['{crRNA} = ' + guide])     
    df = global_store(value)
    dff = df[df['crRNA'] == guide]
    print('len index before', len(dff.index))
    sort_by.insert(0, {'column_id' : 'Mismatches', 'direction': 'asc'})
    sort_by.insert(1, {'column_id' : 'BulgeSize', 'direction': 'asc'})
    sort_by.insert(2, {'column_id': 'CFD', 'direction':'desc'})
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

    print('len index after', len(dff.index))
    #NOTE sort_by: [{'column_id': 'BulgeType', 'direction': 'asc'}, {'column_id': 'crRNA', 'direction': 'asc'}]
    #sort_by.insert(0, {'column_id' : 'Mismatches', 'direction': 'asc'})
    #sort_by.insert(0, {'column_id' : 'BulgeSize', 'direction': 'asc'})
    if len(sort_by):
        dff = dff.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    #Check if results are not 0
    warning_no_res = ''
    with open(job_directory + job_id + '.targets.txt') as t:
        no_result = False
        t.readline()
        last_line = t.readline()
        if (last_line is '' or last_line is '\n'):
            no_result = True

    if (no_result):
        warning_no_res = dbc.Alert("No results were found with the given parameters", color = "warning")

     
    return dff.iloc[
        page_current*page_size:(page_current+ 1)*page_size
    ].to_dict('records')


#For filtering
def split_filter_part(filter_part):
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


#Read the uploaded file and converts into bit
def parse_contents(contents):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    return decoded

#Show image: Barplot
@app.callback(
    [Output('barplot-img', 'src'),
    Output('link-barplot', 'href')],
    [Input('mms-dropdown','value')],
    [State('url', 'search')]
)
def showImages(mms, search):
    if mms is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    job_directory = 'Results/' + job_id + '/'
    barplot_img = 'summary_histogram_' + str(mms) + 'mm.png'
    try:            #NOTE serve per non generare errori se il barplot non è stato fatto
        barplot_src = 'data:image/png;base64,{}'.format(base64.b64encode(open('Results/' + job_id + '/' + barplot_img, 'rb').read()).decode())
        barplot_href = 'assets/Img/' + job_id + '/' + barplot_img
    except:
        barplot_src = ''
        barplot_href = ''
    # guide = all_guides[int(sel_cel[0]['row'])]['Guide'] 
    # radar_img = 'summary_single_guide_' + guide + '_' + str(mms) + 'mm.png'
    # try:
    #     radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(open('Results/' + job_id + '/' + radar_img, 'rb').read()).decode())
    #     radar_href = 'assets/Img/' + job_id + '/' + radar_img
    # except:
    #     radar_src = ''
    #     radar_href = ''
    return barplot_src, barplot_href

#Show image: Radar chart
@app.callback(
    [Output('radar-img', 'src'),
    Output('link-radar', 'href')],
    [Input('mms-dropdown-guide-specific','value')],
    [State('url', 'search'), State('url','hash')]
)
def showImages(mms, search, hash_guide):
    if mms is None:
        raise PreventUpdate
    job_id = search.split('=')[-1]
    job_directory = 'Results/' + job_id + '/'
    # barplot_img = 'summary_histogram_' + str(mms) + 'mm.png'
    # try:            #NOTE serve per non generare errori se il barplot non è stato fatto
    #     barplot_src = 'data:image/png;base64,{}'.format(base64.b64encode(open('Results/' + job_id + '/' + barplot_img, 'rb').read()).decode())
    #     barplot_href = 'assets/Img/' + job_id + '/' + barplot_img
    # except:
    #     barplot_src = ''
    #     barplot_href = ''
    # guide = all_guides[int(sel_cel[0]['row'])]['Guide'] 
    guide = hash_guide.split('#')[1]
    radar_img = 'summary_single_guide_' + guide + '_' + str(mms) + 'mm.png'
    try:
        radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(open('Results/' + job_id + '/' + radar_img, 'rb').read()).decode())
        radar_href = 'assets/Img/' + job_id + '/' + radar_img
    except:
        radar_src = ''
        radar_href = ''
    return radar_src, radar_href

def generate_table(dataframe, id_table, max_rows=26):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns]) ] +
        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))],
        style = {'display':'inline-block'},
        id = id_table
    )

#FOR BUTTON IN TABLE
# element.style {
#     background: none;
#     border: none;
#     margin: 0;
#     padding: 0;
#     cursor: pointer;
#     font-family: monospace;
#     font-size: large;
#     font-weight: normal;
# }




#If the input guides are different len, select the ones with same length as the first
def selectSameLenGuides(list_guides):
    selected_length = len(list_guides.split('\n')[0])
    same_len_guides_list = []
    for guide in list_guides.split('\n'):
        if len(guide) == selected_length:
            same_len_guides_list.append(guide)
    same_len_guides = '\n'.join(same_len_guides_list).strip()
    return same_len_guides

def resultPage(job_id):
    value = job_id
    job_directory = 'Results/' + job_id + '/'
    warning_message = []
    if (not isdir(job_directory)):
        return html.Div(dbc.Alert("The selected result does not exist", color = "danger"))

    #Load mismatches
    with open('Results/' + value + '/Params.txt') as p:
       mms = (next(s for s in p.read().split('\n') if 'Mismatches' in s)).split('\t')[-1]

    mms = int(mms[0])
    mms_values = [{'label':i, 'value':i} for i in range(mms + 1) ]      
    
    col_profile_general = ['Total On-Targets', 'Total Off-Targets']
    for i in range(mms):
        col_profile_general.append(str(i+1) + ' Mismatches')
    col_type = ['numeric' for i in col_profile_general]
    
    
    #Load profile
    try:
        profile = pd.read_csv('Results/' + value + '/' + value + '.profile_complete.xls')   #NOTE profile_complete has ',' as separator
        if len(profile.columns) == 1:
            profile = pd.read_csv('Results/' + value + '/' + value + '.profile.xls', sep='\t')
    except:
        profile = pd.read_csv('Results/' + value + '/' + value + '.profile.xls', sep = '\t')    #NOTE profile has \t as separator or ','
        if len(profile.columns) == 1:
            profile = pd.read_csv('Results/' + value + '/' + value + '.profile.xls')
    
    columns_profile_table = [{'name':'Guide', 'id' : 'Guide', 'type':'text'}, {'name':'CFD', 'id': 'CFD', 'type':'numeric'}, {'name':'Total On-Targets', 'id' : 'Total On-Targets', 'type':'numeric'}, {'name':'Total Off-targets', 'id' : 'Total Off-Targets', 'type':'numeric'}]
    keep_column = ['GUIDE', 'ONT', 'OFFT']
    for i in range (mms):
        columns_profile_table.append({'name': str(i+1) + ' Mismatches', 'id':str(i+1) + ' Mismatches', 'type':'numeric'})
        keep_column.append(str(i+1) + 'MM')
    columns_profile_table.append({'name':'More Info', 'id':'More Info', 'type':'text'})
    print(profile.columns)
    profile = profile[keep_column]
    rename_columns = {'GUIDE':'Guide',"ONT":'Total On-Targets', 'OFFT':'Total Off-Targets'}
    for i in range(mms):
        rename_columns[str(i+1) + 'MM'] = str(i+1) + ' Mismatches'

    profile.rename(columns = rename_columns, inplace = True)    #Now profile is Guide, Total On-targets, ...
    #link_for_guides = [html.A('Show all...', href = 'http://127.0.0.1:8050/result?job=' + job_id + '#' + i, target = '_blank') for i in profile['Guide']]
    #profile['More Info'] = link_for_guides
    final_list = []

    
    #final_list.append(html.H1('CRISPRitz Web Application'))
    final_list.append(
        html.H3('Result Summary')
    )
    #final_list.append(html.P('Select a Guide to view more details'))
    # final_list.append(html.Div(
    #         generate_table(profile, 'result-page-guide-table'),
    #         style = {'text-align':'center'}
    #     )
    # )

    #load acfd for each guide   #TODO sistemare e controllare
    with open('Results/' + value + '/acfd.txt') as a:
        acfd = a.read().strip().split('\n')
    acfd.remove('crRNA 0')
    acfd.sort()
    acfd = [float(a.split(' ')[-1]) for a in acfd]
    acfd  = [int(round((100/(100 + x))*100)) for x in acfd]
    profile = profile.sort_values('Guide')
    profile['CFD'] = acfd
    profile = profile.sort_values('CFD', ascending = False)
    final_list.append(
        html.Div(
            dash_table.DataTable(
                id = 'general-profile-table',
                page_size=PAGE_SIZE,
                columns = columns_profile_table,
                data = profile.to_dict('records')
            )
            ,id = 'div-general-profile-table')
    )

    final_list.append(html.Br())
    final_list.append(
        html.Div(
            [
                html.Div(
                    [
                        html.H5('Comparison with Reference Genome'),
                        html.P('Select the mismatch value'),
                        dcc.Dropdown(options = mms_values, id = 'mms-dropdown', style = {'flex':'0 0 5%', 'width':'100px'}, clearable = False)
                    ]
                ),
                html.Div(
                    html.A(
                        html.Img(id = 'barplot-img', width="100%", #height="30%", 
                        
                        ),
                        
                        target="_blank",
                        id = 'link-barplot'
                        
                    ),
                    style = {'flex':'0 0 30%'}
                ),
                html.Div(
                html.A(
                    html.Img( width="100%", #height="30%", 
                    
                    ),
                    
                    target="_blank",
                    
                ),
                style = {'flex':'0 0 30%'}
            ),
                
            ],
            className = 'flex-view-images'
        )
    )

    result_page = html.Div(final_list, style = {'margin':'1%'})
    return result_page


def guidePage(job_id, guide):
    value = job_id
    final_list = []
    final_list.append(html.P('List of Targets found for the selected guide'))
    col_list = ['BulgeType', 'crRNA', 'DNA', 'Chromosome', 'Position', 'Direction', 'Mismatches', 'BulgeSize', 'CFD', 'Doench2016']
    col_type = ['text','text','text','text','numeric','text','numeric', 'numeric', 'numeric', 'numeric', 'numeric']
    cols = [{"name": i, "id": i, 'type':t} for i,t in zip(col_list, col_type)]
    job_directory = 'Results/' + job_id + '/'
    #Load mismatches
    with open('Results/' + value + '/Params.txt') as p:
       mms = (next(s for s in p.read().split('\n') if 'Mismatches' in s)).split('\t')[-1]

    mms = int(mms[0])
    mms_values = [{'label':i, 'value':i} for i in range(mms + 1) ]
    global_store(job_id)    #TODO controllare se carica ogni volta o solo la prima
    #NOTE the filtering is done automatically when the page is loaded due to the function update_table since it's triggered when the table is created, putting page_current, sort_by etc
    #at initial values

    # df = global_store(job_id)
    # dff = df
    # filtering_expressions = ['{crRNA} = ' + guide]
    # for filter_part in filtering_expressions:
    #     col_name, operator, filter_value = split_filter_part(filter_part)

    #     if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
    #         # these operators match pandas series operator method names
    #         dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
    #     elif operator == 'contains':
    #         dff = dff.loc[dff[col_name].str.contains(filter_value)]
    #     elif operator == 'datestartswith':
    #         # this is a simplification of the front-end filtering logic,
    #         # only works with complete fields in standard format
    #         dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    final_list.append(
        html.Div(
            dash_table.DataTable(
                id='result-table', 
                columns=cols, 
                #data = dff.to_dict('records'),
                virtualization = True,
                fixed_rows={ 'headers': True, 'data': 0 },
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
                    'height': '300px',
                    #'overflowY': 'scroll',
                },
                # style_data_conditional=[{
                #     "if": {'column_id':'BulgeType', 'filter_query' : 'BulgeType eq "RNA"'}, #{'filter_query' : 'BulgeType eq "RNA"'},
                #     "backgroundColor": "lightblue",
                #     'color': 'white'
                # }],
            ),
            id = 'div-result-table',
        )
    )
    final_list.append(html.Br())
    final_list.append(
        html.Div(
            [
                html.Div(
                    [
                        html.P('Select the mismatch value'),
                        dcc.Dropdown(options = mms_values, id = 'mms-dropdown-guide-specific', style = {'flex':'0 0 5%'}, clearable = False)
                    ]
                ),
                
                html.Div(
                    html.A(
                        html.Img(id = 'radar-img', width="100%", #height="30%", 
                        
                        ),
                        
                        target="_blank",
                        id = 'link-radar'
                        
                    ),
                    style = {'flex':'0 0 30%'}
                ),
                html.Div(
                    html.A(
                        html.Img( width="100%", #height="30%", 
                        
                        ),
                        
                        target="_blank",
                        
                    ),
                    style = {'flex':'0 0 30%'}
                )
                
            ],
            className = 'flex-view-images'
        )
    )
    return html.Div(final_list, style = {'margin':'1%'})

if __name__ == '__main__':
    app.run_server(debug=True)
    cache.clear()       #delete cache when server is closed

    #BUG quando faccio scores, se ho dei char IUPAC nei targets, nel terminale posso vedere 150% 200% etc perche' il limite massimo e' basato su wc -l dei targets, ma possono aumentare se ho molti
    #Iupac
