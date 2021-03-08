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

PAGE_SIZE = 10                     #number of entries in each page of the table in view report

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
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
gen_dir = []
for dir in onlydir:
    gen_dir.append({'label': dir, 'value' : dir})

#Dropdown available PAM
onlyfile = [f for f in listdir('pam') if isfile(join('pam', f))]
pam_file = []
for pam_name in onlyfile:
    pam_file.append({'label': pam_name, 'value' : pam_name})

#Dropdown available Variants
onlydir = [f for f in listdir('Variants') if isdir(join('Variants', f))]
var_dir = []
for dir in onlydir:
    var_dir.append({'label': dir, 'value' : dir})

#For multipage
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content'),
    html.P(id = 'signal', style = {'visibility':'hidden'})
])

# final_list.append(
#     html.Div(
#         [
#             html.P('Insert Job title: ', style = {'padding-top':'5px'}),
#             dcc.Input(id = 'job-name', size = '30')
#         ],
#         className = 'flex-job-title',
#         style = {'margin':'1%', 'width':'23%'}
#     )
# )
# final_list.append(
#     html.Div(
#         [
#             html.P('Insert Email: ', style = {'padding-top':'5px'}),
#             dcc.Input(id = 'email', size = '30')
#         ],
#         className = 'flex-job-title',
#         style = {'margin':'1%', 'width':'23%'}
#     )
# )



#new final_list
final_list = []
final_list.extend([html.H1('CRISPRitz Web Application'),
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

final_list.append(
    html.Div(
        html.Div(
            [
                html.Div(
                    [
                        html.H3('STEP 1', style = {'margin-top':'0'}), 
                        html.P('Select a genome'),
                        html.Div(
                            dcc.Dropdown(options = gen_dir, clearable = False, id = "available-genome", style = {'width':'75%'}),
                            #style = {'width':'50%'}
                        ),
                        html.P('Add a genome variant'),
                        html.Div(
                            dcc.Dropdown(options = var_dir, clearable = False, id = 'available-variant', style = {'width':'75%'})
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.P('Select PAM'),
                                        html.Div(
                                            dcc.Dropdown(options = pam_file, clearable = False, id = 'available-pam', style = {'width':'75%'})
                                        )
                                    ],
                                    style = {'flex':'0 0 50%'}
                                ),
                                #html.P('or'),
                                html.Div(
                                    [
                                        html.P('Insert custom PAM'),
                                        dcc.Input(type = 'text', id = 'custom-pam', placeholder = 'NGG')
                                    ]
                                )
                            ],
                            id = 'div-pam',
                            className = 'flex-div-pam'
                        )
                    ],
                    id = 'step1',
                    style = {'flex':'0 0 40%'}
                ),
                html.Div(style = {'border-right':'solid 1px white'}),
                html.Div(
                    [
                        html.H3('STEP 2', style = {'margin-top':'0'}),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        
                                        html.Div(
                                            [
                                                html.P(['Insert crRNA sequence(s)', html.Abbr('\uD83D\uDEC8', style = {'text-decoration':'none'} ,title = 'One sequence per line. All sequences must have the same lenght and PAM characters are not required')], style = {'word-wrap': 'break-word'}), 
                            
                                                dcc.Textarea(id = 'text-guides', placeholder = 'GAGTCCGAGCAGAAGAAGAA\nCCATCGGTGGCCGTTTGCCC', style = {'width':'275px', 'height':'160px'}),
                                                html.P('or', style = {'position': 'relative', 'left':'50%'}),
                                                html.Div(
                                                    [
                                                        html.Div(
                                                            [
                                                                dcc.Upload('Upload file with crRNA sequences', id = 'upload-guides')
                                                            ],
                                                            style={
                                                                'width': '100%',
                                                                'height': '60px',
                                                                'lineHeight': '60px',
                                                                'borderWidth': '1px',
                                                                'borderStyle': 'dashed',
                                                                'borderRadius': '5px',
                                                                'textAlign': 'center',
                                                                #'margin': '10px'
                                                            }
                                                        ),
                                                        html.P('', id = 'uploaded-filename',style = {'width':'275px', 'visibility':'hidden', 'display':'inline-block'})
                                                    ],
                                                    style = {'text-align':'center'}
                                                )
                                            ],
                                            style = {'width':'275px'} #same as text-area
                                        )
                                    ],
                                    id = 'div-guides'
                                ),
                                html.Div(
                                    [
                                        html.P('Allowed mismatches'),
                                        dcc.Input(value = '0', id = 'mms', type = 'number', min = '0', style = {'width':'60px'}),
                                        html.P('Bulge DNA size'),
                                        dcc.Input(value = '0', id = 'dna', type = 'number', min = '0', style = {'width':'60px'}),
                                        html.P('Bulge RNA size'),
                                        dcc.Input(value = '0', id = 'rna', type = 'number', min = '0', style = {'width':'60px'})
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
                        html.H3('Submit', style = {'margin-top':'0'}),
                        html.Div(
                            [
                                html.Button('Submit', id = 'submit-job')
                            ],
                            style = {'display':'inline-block', 'margin':'0 auto'}
                        )
                    ],
                    id = 'step3',
                    style = {'tex-align':'center'}
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
final_list.append(html.H1('CRISPRitz Web Application'))
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
                                html.Li('Adding variants'),
                                html.Li('Searching crRNA'),
                                html.Li('Generating report')
                            ]
                        ),
                        style = {'flex':'0 0 20%'}
                    ),
                    html.Div(
                        html.Ul(
                            [
                                html.Li('To do', style = {'color':'red'}, id = 'add-variants-status'),
                                html.Li('To do', style = {'color':'red'}, id = 'search-status'),
                                html.Li('To do', style = {'color':'red'}, id = 'generate-report-status')
                            ],
                            style = {'list-style-type':'none'}
                        )
                    )
                ],
                className = 'flex-status'
            ),
            html.Div(
                dcc.Link('View Results', style = {'visibility':'hidden'}, id = 'view-results')
            )
        ],
        id = 'div-status-report'
    )
)

final_list.append(html.P('', id = 'done'))

final_list.append(dcc.Interval(id = 'load-page-check', interval=3*1000))
load_page = html.Div(final_list, style = {'margin':'1%'})

#Result page
final_list = []
final_list.append(html.H1('CRISPRitz Web Application'))

col_list = ['BulgeType', 'crRNA', 'DNA', 'Chromosome', 'Position', 'Direction', 'Mismatches', 'BulgeSize', 'CFD', 'Doench2016']
col_type = ['text','text','text','text','numeric','text','numeric', 'numeric', 'numeric', 'numeric', 'numeric']
cols = [{"name": i, "id": i, 'type':t} for i,t in zip(col_list, col_type)]
final_list.append(
    html.Div(
        dash_table.DataTable(
            id='result-table', 
            columns=cols, 
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
            filter_query=''
        ),
        id = 'div-result-table'
    )
)

final_list.append(html.Br())
final_list.append(
    html.Div(
        [
            dash_table.DataTable(
                id = 'guide-table',
                columns = [{'name':'Available Guides', 'id':'Guides', 'type':'text'}],
                page_size=PAGE_SIZE
                
            ),
            html.Div(
                [
                    html.P('Select the mismatch value'),
                    dcc.Dropdown(id = 'mms-dropdown', style = {'flex':'0 0 5%'}, clearable = False)
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
                    html.Img(id = 'barplot-img', width="100%", #height="30%", 
                    
                    ),
                    
                    target="_blank",
                    id = 'link-barplot'
                    
                ),
                style = {'flex':'0 0 30%'}
            )
            
        ],
        className = 'flex-view-images'
    )
)

result_page = html.Div(final_list, style = {'margin':'1%'})
##################################################CALLBACKS##################################################

#Submit Job, change url
@app.callback(
    [Output('url', 'pathname'),
    Output('url','search')],
    [Input('submit-job', 'n_clicks')],
    [State('url', 'href'),
    State('available-genome', 'value'),
    State('available-variant', 'value'),
    State('available-pam','value'),
    State('custom-pam','value'),
    State('text-guides', 'value'),
    State('upload-guides','contents'),
    State('mms','value'),
    State('dna','value'),
    State('rna','value')]
)
def changeUrl(n, href, genome_ref, variant, pam, custom_pam, text_guides, file_guides, mms, dna, rna):      #NOTE startJob
    if n is None:
        raise PreventUpdate
    
    #TODO check se input è corretto
   
    
    job_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k = 10))
    result_dir = 'Results/' + job_id
    subprocess.run(['mkdir ' + result_dir], shell = True)
    
    enrich = True
    index = True
    search_index = True
    search = True
    annotation = True
    report =  True

    #Check which tool to use
    #TODO controllare se add-variants mi crea una cartella all'interno di SNP genomes, nel caso fare in modo che lo faccia

    #Enrichment
    if variant is None:
        variant = 'None'
        enrich = False
        genome_enr = genome_ref
    else:
        all_genome_enr = [f for f in listdir('variants_genome/SNPs_genome') if isdir(join('variants_genome/SNPs_genome', f))]
        if variant is 'None':
            genome_enr = genome_ref
        else:
            genome_enr = genome_ref + '_' + variant
        if (genome_ref + '_' + variant in all_genome_enr):          #NOTE Enriched genomes dir names are genome_ref name + symbol + variant_dir name
            enrich = False
        
    #subprocess.run(['crispritz.py add-variants ' + 'Variants/' + variant + ' ' + 'Genomes/' + genome_ref], shell = True)
    #Indexing and search
    #NOTE Indexed genomes names are PAM + _ + bMax + _ + genome_ref/genome_enr
    all_genomes_idx = [f for f in listdir('genome_library') if isdir(join('genome_library', f))]
    
    #TODO if pam input area, find way to select left or right of guide
    pam_len = 0
    if custom_pam is not None and custom_pam is not '':
        pam_char = custom_pam
        pam_len = len(pam_char)
        #Save to file as NNNN...PAM, but first check what guides have been inserted
        if text_guides is not None and text_guides is not '':
            n_number = len(text_guides.split('\n')[0])
        else:
            decoded_guides = parse_contents(file_guides).decode('utf-8')
            n_number = len(decoded_guides.split('\n')[0]) - decoded_guides.split('\n')[0].count('N')
        with open(result_dir + '/pam.txt', 'w') as pam_file:
            pam_file.write(('N' * n_number) + pam_char)
            pam = result_dir + '/pam.txt'
    else:
        with open('pam/' + pam) as pam_file:
            pam_char = pam_file.readline()
            
            if int(pam_char.split(' ')[-1]) < 0:
                end_idx = int(pam_char.split(' ')[-1]) * (-1)
                pam_char = pam_char.split(' ')[0][0 : end_idx]
                pam_len = end_idx
            else:
                end_idx = int(pam_char.split(' ')[-1])
                pam_char = pam_char.split(' ')[0][end_idx * (-1):]
                pam_len = end_idx
        subprocess.run(['cp pam/' + pam + ' ' + result_dir + '/pam.txt'], shell = True)
        pam = result_dir + '/pam.txt'

    guides_file = result_dir + '/guides.txt'
    if text_guides is not None and text_guides is not '':
        save_guides_file = open(result_dir + '/guides.txt', 'w')
        text_guides = text_guides.replace('\n', 'N' * pam_len + '\n') + 'N' * pam_len    #TODO what if pam at beginning?
        save_guides_file.write(text_guides)
        save_guides_file.close()     
    else:
        decoded_guides = parse_contents(file_guides).decode('utf-8')
        save_guides_file = open(result_dir + '/guides.txt', 'w')
        save_guides_file.write(decoded_guides)
        save_guides_file.close()

    if (int(dna) == 0 and int(rna) == 0):
        index = False
        search_index = False
    max_bulges = rna
    if (int(dna) > int(rna)):
        max_bulges = dna
    if (index and (pam_char + '_' + str(max_bulges) + '_' + genome_ref + '_' + variant) in all_genomes_idx):
        index = False

    if (search_index):
        search = False
    # else:
    #     search_index = False

    if variant is 'None':
        genome_idx = pam_char + '_' + str(max_bulges) + '_' + genome_ref
    else:
        genome_idx = pam_char + '_' + str(max_bulges) + '_' + genome_ref + '_' + variant

    #Create Params.txt file
    with open(result_dir + '/Params.txt', 'w') as p:
        p.write('Genome_ref\t' + genome_enr + '\n')
        if search_index:
            p.write('Genome_idx\t' + genome_idx + '\n')
        else:
            p.write('Genome_idx\t' + 'None\n')
        p.write('Variant\t' + str(variant) + '\n')
        p.write('Pam\t' + pam_char + '\n')
        p.write('Max_bulges\t' + str(max_bulges) + '\n')
        p.write('Mismatches\t' + str(mms) + '\n')
        p.write('DNA\t' + str(dna) + '\n')
        p.write('RNA\t' + str(rna) + '\n')
        p.close()


    
    #Check if input parameters (mms, bulges, pam, guides, genome) are the same as a previous search
    all_result_dirs = [f for f in listdir('Results') if isdir(join('Results', f))]
    all_result_dirs.remove(job_id)
    for check_param_dir in all_result_dirs:
        if os.path.exists('Results/' + check_param_dir + '/Params.txt'):
            print('checkparamdir:', check_param_dir)
            if (filecmp.cmp('Results/' + check_param_dir + '/Params.txt', result_dir + '/Params.txt' )):
                search = False
                search_index = False
                subprocess.run(['ln -s $PWD/Results/' + check_param_dir + '/' + check_param_dir + '* ' + result_dir + '/'], shell = True) #TODO copy result from one directory to the current one or create simlink
                subprocess.run(['ln -s $PWD/Results/' + check_param_dir + '/*.png ' + result_dir + '/'], shell = True)
                subprocess.run(['rename \'s/' + check_param_dir + '/' + job_id + '/g\' ' + result_dir + '/*'], shell = True)
                break           #BUG manca il controllo sulle guide
    #Annotation
    if (not search and not search_index):
        annotation = False      #TODO copy result from one directory to the current one or create simlink
 
    #Generate report
    if (not enrich and not index and not search and not search_index):
        report = False          #TODO copy result from one directory to the current one or create simlink
    #TODO if human genome -> annotation = human genome. mouse -> annotation mouse etc
    subprocess.Popen(['assets/./submit_job.sh ' + 'Results/' + job_id + ' ' + str(variant) + ' ' + 'Genomes/' + genome_ref + ' ' + 'variants_genome/SNPs_genome/' + genome_enr + ' ' + 'genome_library/' + genome_idx + (
        ' ' + pam + ' ' + guides_file + ' ' + str(mms) + ' ' + str(dna) + ' ' + str(rna) + ' ' + str(enrich) + ' ' + str(index) + ' ' + str(search_index) + ' ' + str(search) + (
            ' ' + str(annotation) + ' ' + str(report))) ], shell = True)
    return '/load','?job=' + job_id

#When url changed, load new page
@app.callback(
    [Output('page-content', 'children'),
    Output('job-link', 'children')],
    [Input('url', 'pathname')],
    [State('url','href'), 
    State('url','search')
    ]
)
def changePage(path, href, search):

    if path == '/load':
        return load_page, 'http://127.0.0.1:8050/load' + search #NOTE change the url part when DNS are changed
    if path == '/result':
        return result_page, 'http://127.0.0.1:8050/load' + search
    return index_page, ''

#Check end job
@app.callback(
    [Output('view-results', 'style'),
    Output('add-variants-status', 'children'),
    Output('search-status', 'children'),
    Output('generate-report-status', 'children'),
    Output('view-results','href')],
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
                add_var_status = html.P('To do', style = {'color':'red'})
                search_status = html.P('To do', style = {'color':'red'})
                report_status = html.P('To do', style = {'color':'red'})
                current_log = log.read()
                if ('Add-variants\tDone' in current_log):
                    add_var_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                if ('Search-index\tDone' in current_log or 'Search\tDone' in current_log):
                    search_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                if ('Report\tDone' in current_log):
                    report_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                if all_done == 3:
                    return {'visibility':'visible'}, add_var_status, search_status, report_status, '/result?job=' + dir_name.split('=')[-1]
                else:
                    return {'visibility':'hidden'}, add_var_status, search_status, report_status,''
    raise PreventUpdate

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

#Callback to populate the tab, note that it's called when the result_page is loaded (dash implementation), so we do not use raise update to block this first callback
@app.callback(
    [Output('signal','children'),
    Output('result-table','page_current'),
    Output('result-table', "sort_by"), 
    Output('result-table','filter_query')],
    [Input('url', 'pathname')],
    [State('url', 'search')]
)
def populateTable(pathname, search):
    if pathname != '/result':
        raise PreventUpdate

    job_id = search.split('=')[-1]
    job_directory = 'Results/' + job_id + '/'
    #print('JOB_ID: ', job_id)
    global_store(job_id)
    return job_id, 0, [], ''

#Send the data when next or prev button is clicked on the result table
@app.callback(
    [Output('result-table', 'data'),
    Output('guide-table', 'data'),
    Output('mms-dropdown','options')],
    [Input('signal', 'children'),
     Input('result-table', "page_current"),
     Input('result-table', "page_size"),
     Input('result-table', "sort_by"),
     Input('result-table', 'filter_query')]
)
def update_table(value, page_current,page_size, sort_by, filter):
    #print('signal_children', value)
    if value is None:
        raise PreventUpdate

    
    filtering_expressions = filter.split(' && ')    
    df = global_store(value)
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
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )
    
    #Load guide table
    df_guide = pd.read_csv('Results/' + value + '/guides.txt', names = ['Guides'])

    #Load mismatches
    with open('Results/' + value + '/Params.txt') as p:
       
        mms = (next(s for s in p.read().split('\n') if 'Mismatches' in s)).split('\t')[-1]
    mms = int(mms[0])
    mms_values = [{'label':i, 'value':i} for i in range(mms + 1) ]

    return dff.iloc[
        page_current*page_size:(page_current+ 1)*page_size
    ].to_dict('records'), df_guide.to_dict('records'), mms_values

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

#Show images
@app.callback(
    [Output('barplot-img', 'src'),
    Output('link-barplot', 'href'),
    Output('radar-img','src'),
    Output('link-radar','href')],
    [Input('guide-table','selected_cells'),
    Input('mms-dropdown','value')],
    [State('guide-table', 'data'),
    State('url','search')]
)
def testcel(sel_cel, mms, all_guides, job_id):
    if sel_cel is None or mms is None:
        raise PreventUpdate
    job_id = job_id.split('=')[-1]
    barplot_img = 'summary_histogram_' + str(mms) + 'mm.png'
    try:            #NOTE serve per non generare errori se il barplot non è stato fatto
        barplot_src = 'data:image/png;base64,{}'.format(base64.b64encode(open('Results/' + job_id + '/' + barplot_img, 'rb').read()).decode())
        barplot_href = 'assets/Img/' + job_id + '/' + barplot_img
    except:
        barplot_src = ''
        barplot_href = ''
    guide = all_guides[int(sel_cel[0]['row'])]['Guides'] 
    radar_img = 'summary_single_guide_' + guide + '_' + str(mms) + 'mm.png'
    radar_src = 'data:image/png;base64,{}'.format(base64.b64encode(open('Results/' + job_id + '/' + radar_img, 'rb').read()).decode())
    radar_href = 'assets/Img/' + job_id + '/' + radar_img
    return barplot_src, barplot_href, radar_src, radar_href

#Show filename if user upload a file
@app.callback(
    [Output('uploaded-filename', 'children'),
    Output('uploaded-filename', 'style')],
    [Input('upload-guides', 'filename')]
)
def showUploadedFilename(name):
    if name is None:
        raise PreventUpdate
    return 'Uploaded file: ' + name, {'visibility':'visible'}

if __name__ == '__main__':
    app.run_server(debug=True)
    cache.clear()       #delete cache when server is closed

    #TODO se faccio l'annotazione (stessi parametri) dei targets ottenuti da enr e ref genomes, poi posso usare i loro summary counts per fare il barplot, che dipende solo dai mm e non dalle guide
    #BUG quando faccio scores, se ho dei char IUPAC nei targets, nel terminale posso vedere 150% 200% etc perche' il limite massimo e' basato su wc -l dei targets, ma possono aumentare se ho molti
    #Iupac
    #BUG emx1.txt error on loading extended_profile


    #TODO bootstrap per fare il menu, aggiungere selezione per gecko e il barplot, inserire email quando finito, magari un menu avanzate per opzioni avanzate, divisione tra genomi e genomi enr (al posto
    # della select varfile), cambiare il nome delle pam togliendo txt e mettendo nome giusto