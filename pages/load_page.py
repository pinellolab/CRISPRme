import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from app import app, current_working_directory, URL
import os
from os.path import isfile, isdir, join  # for getting directories
from os import listdir

# Check end job


@app.callback(
    [Output('view-results', 'style'),
     Output('index-status', 'children'),
     Output('search-status', 'children'),
     Output('post-process-status', 'children'),
     Output('merge-status', 'children'),
     Output('images-status', 'children'),
     Output('database-status', 'children'),
     Output('integrate-status', 'children'),
     Output('view-results', 'href'),
     Output('no-directory-error', 'children')],
    [Input('load-page-check', 'n_intervals')],
    [State('url', 'search')]
)
def refreshSearch(n, dir_name):
    '''
    Il componente Interval chiama questa funzione ogni 3 secondi. Essa controlla lo stato del lavoro e aggiorna la pagina se una parte del lavoro
    è stata fatta.
    Quando la ricerca è finita, visualizza un link per passare alla pagina dei risultati
    Se il job non esiste, ritorna un avviso di errore
    '''
    if n is None:
        raise PreventUpdate

    onlydir = [f for f in listdir(current_working_directory + 'Results')
               if isdir(join(current_working_directory + 'Results', f))]
    current_job_dir = current_working_directory + \
        'Results/' + dir_name.split('=')[-1] + '/'
    if dir_name.split('=')[-1] in onlydir:
        onlyfile = [f for f in listdir(
            current_job_dir) if isfile(join(current_job_dir, f))]
        if os.path.exists(current_job_dir + 'guides.txt'):
            with open(current_job_dir + 'guides.txt') as guides:
                n_guides = len(guides.read().strip().split('\n'))
        else:
            n_guides = -1
        if 'log.txt' in onlyfile:
            with open(current_job_dir + 'log.txt') as log:
                all_done = 0

                index_status = html.P('To do', style={'color': 'red'})
                search_status = html.P('To do', style={'color': 'red'})
                post_process_status = html.P('To do', style={'color': 'red'})
                merge_status = html.P('To do', style={'color': 'red'})
                images_status = html.P('To do', style={'color': 'red'})
                database_status = html.P('To do', style={'color': 'red'})
                integrate_status = html.P('To do', style={'color': 'red'})
                current_log = log.read()

                variant = False
                with open(current_job_dir + 'Params.txt') as f:
                    if "Ref_comp\tTrue" in f.read():
                        variant = True

                if variant:
                    if "Index-genome Variant\tEnd" in current_log:
                        index_status = html.P('Done', style={'color': 'green'})
                        all_done = all_done + 1
                    elif "Index-genome Variant\tStart" in current_log:
                        index_status = html.P(
                            'Indexing Enriched Genome...' + ' ' + 'Step [4/4]', style={'color': 'orange'})
                    elif "Index-genome Reference\tStart" in current_log:
                        index_status = html.P(
                            'Indexing Reference Genome...' + ' ' + 'Step [3/4]', style={'color': 'orange'})
                    elif "Indexing Indels\tStart" in current_log:
                        index_status = html.P(
                            'Indexing Indels Genome...' + ' ' + 'Step [2/4]', style={'color': 'orange'})
                    elif 'Add-variants\tStart' in current_log:
                        index_status = html.P(
                            'Adding variants...' + ' ' + 'Step [1/4]', style={'color': 'orange'})
                    elif 'Search Reference\tStart' in current_log:
                        index_status = html.P('Done', style={'color': 'green'})
                        all_done = all_done + 1
                else:
                    if "Index-genome Reference\tEnd" in current_log:
                        index_status = html.P('Done', style={'color': 'green'})
                        all_done = all_done + 1
                    elif "Index-genome Reference\tStart" in current_log:
                        index_status = html.P(
                            'Indexing Reference Genome...' + ' ' + 'Step [1/1]', style={'color': 'orange'})
                    elif 'Search Reference\tStart' in current_log:
                        index_status = html.P('Done', style={'color': 'green'})
                        all_done = all_done + 1

                if variant:
                    if "Search Reference\tEnd" in current_log and 'Search Variant\tEnd' in current_log and "Search INDELs\tEnd" in current_log:
                        search_status = html.P(
                            'Done', style={'color': 'green'})
                        all_done = all_done + 1
                    elif "Search Reference\tStart" in current_log or 'Search Variant\tStart' in current_log:
                        search_status = html.P(
                            'Searching...', style={'color': 'orange'})
                else:
                    if "Search Reference\tEnd" in current_log:
                        search_status = html.P(
                            'Done', style={'color': 'green'})
                        all_done = all_done + 1
                    elif "Search Reference\tStart" in current_log:
                        search_status = html.P(
                            'Searching...', style={'color': 'orange'})

                if variant:
                    if 'Post-analysis SNPs\tEnd' in current_log and 'Post-analysis INDELs\tEnd' in current_log:
                        post_process_status = html.P(
                            'Done', style={'color': 'green'})
                        all_done = all_done + 1
                    elif 'Post-analysis SNPs\tEnd' in current_log:
                        post_process_status = html.P(
                            'Post-analysis on INDELs...' + ' ' + 'Step [2/2]', style={'color': 'orange'})
                    elif 'Post-analysis SNPs\tStart' in current_log:
                        post_process_status = html.P(
                            'Post-analysis on SNPs...' + ' ' + 'Step [1/2]', style={'color': 'orange'})
                else:
                    if 'Post-analysis\tEnd' in current_log:
                        post_process_status = html.P(
                            'Done', style={'color': 'green'})
                        all_done = all_done + 1
                    elif 'Post-analysis\tStart' in current_log:
                        post_process_status = html.P(
                            'Post-analysis...' + ' ' + 'Step [1/1]', style={'color': 'orange'})

                if 'Merging Close Targets\tEnd' in current_log:
                    merge_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Merging Close Targets\tStart' in current_log:
                    merge_status = html.P(
                        'Processing...' + ' ' + 'Step [1/1]', style={'color': 'orange'})

                if 'Creating images\tEnd' in current_log:
                    images_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Creating images\tStart' in current_log:
                    images_status = html.P(
                        'Processing...' + ' ' + 'Step [1/1]', style={'color': 'orange'})

                if 'Creating database\tEnd' in current_log:
                    database_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Creating database\tStart' in current_log:
                    database_status = html.P(
                        'Processing...' + ' ' + 'Step [1/1]', style={'color': 'orange'})
                
                if 'Integrating results\tEnd' in current_log:
                    integrate_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Integrating results\tStart' in current_log:
                    integrate_status = html.P(
                        'Processing...' + ' ' + 'Step [1/1]', style={'color': 'orange'})
                '''
                if ('Search-index\tDone' in current_log and 'Search\tDone' in current_log):
                    search_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                elif os.path.exists(current_job_dir + 'output.txt'):                #Extract % of search done 
                    with open(current_job_dir + 'output.txt', 'r') as output_status:
                        line = output_status.read().strip()
                        if 'Indexing_Reference' in line:
                            search_status = html.P('Indexing Reference Genome...' + ' ' + 'Step [1/2]', style = {'color':'orange'})
                        if 'Indexing_Enriched' in line:
                            search_status = html.P('Indexing Enriched Genome...' + ' ' + 'Step [2/2]', style = {'color':'orange'})
                        if 'Search_output' in line:
                            if 'both' in line:
                                last_percent = line.rfind('%')
                                if last_percent > 0:
                                    last_percent = line[line[:last_percent].rfind(' '): last_percent]
                                    search_status_message = last_percent + '%'
                                else:
                                    search_status_message = 'Searching...'

                                steps = 'Step [1/2]'
                                if 'Search_output_ref' in line:
                                    steps = 'Step [2/2]'
                                
                            else:
                                last_percent = line.rfind('%')
                                if last_percent > 0:
                                    last_percent = line[line[:last_percent].rfind(' '): last_percent]
                                    search_status_message = last_percent + '%'
                                else:
                                    search_status_message = 'Searching...'
                                steps = ''
                            search_status = html.P(search_status_message + ' ' + steps, style = {'color':'orange'})

                if ('Report\tDone' in current_log):
                    report_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                elif os.path.exists(current_job_dir + 'output.txt'):                #Extract % of search done
                    with open(current_job_dir + 'output.txt', 'r') as output_status:
                        line = output_status.read().strip()
                        if 'Generate_report' in line:
                            if n_guides < 0:
                                report_status = html.P('Generating Report...', style = {'color':'orange'}) 
                            else:
                                status_message = round((len(line.split('\n')) - 1) / n_guides, 2)
                                report_status = html.P(str(status_message * 100) + '%', style = {'color':'orange'})
                if ('PostProcess\tDone' in current_log):
                    post_process_status = html.P('Done', style = {'color':'green'})
                    all_done = all_done + 1
                elif os.path.exists(current_job_dir + 'output.txt'):                #Extract % of search done
                    with open(current_job_dir + 'output.txt', 'r') as output_status:
                        line = output_status.read().strip()
                        if 'PostProcess_output' in line:
                            line = line.split('\n')
                            if line[-1] == 'PostProcess_output':
                                post_process_status = html.P('Processing data...', style = {'color':'orange'})    
                            else:
                                if 'Annotating...' in line:
                                    last_percent = line[-1].rfind('%')
                                    if last_percent > 0:
                                        last_percent = line[line[:last_percent].rfind(' '): last_percent]
                                        status_message = last_percent + '%'
                                    else:
                                        status_message = 'Annotating...'
                                else:
                                    status_message = line[-1]
                                post_process_status = html.P(status_message, style = {'color':'orange'})
                if all_done == 3 or 'Job\tDone' in current_log:
                    return {'visibility':'visible'},  search_status, report_status, post_process_status ,'/result?job=' + dir_name.split('=')[-1], ''
                else:
                    return {'visibility':'hidden'},  search_status, report_status, post_process_status,'', ''
                '''
                if all_done == 7 or 'Job\tDone' in current_log:
                    return {'visibility': 'visible'}, index_status, search_status, post_process_status, merge_status, images_status, database_status, integrate_status, '/result?job=' + dir_name.split('=')[-1], ''
                else:
                    return {'visibility': 'hidden'}, index_status, search_status, post_process_status, merge_status, images_status, database_status, integrate_status, '', ''
        elif 'queue.txt' in onlyfile:
            return {'visibility': 'hidden'}, html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('To do', style={'color': 'red'}), html.P('To do', style={'color': 'red'}), html.P('To do', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}),'', dbc.Alert("Job submitted. Current status: in queue", color="info")
    return {'visibility': 'hidden'},  html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}),html.P('Not available', style={'color': 'red'}),'', dbc.Alert("The selected result does not exist", color="danger")

# Load Page


def load_page():
    final_list = []
    final_list.append(
        html.Div(
            html.Div(
                html.Div(
                    [
                        html.P(
                            'Job submitted. Copy this link to view the status and the result page '),
                        html.Div(
                            html.P(
                                'link', id='job-link', style={'margin-top': '0.75rem', 'font-size': 'large'}),
                            style={'border-radius': '5px', 'border': '2px solid', 'border-color': 'blue',
                                   'width': '100%', 'display': 'inline-block', 'margin': '5px'}
                        ),
                        html.P('Results will be kept available for 3 days')
                    ],
                    style={'display': 'inline-block'}
                ),
                style={'display': 'inline-block', 'background-color':
                       'rgba(154, 208, 150, 0.39)', 'border-radius': '10px', 'border': '1px solid black', 'width': '70%'}
            ),
            style={'text-align': 'center'}
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
                                    html.Li('Indexing Genome(s)'),
                                    html.Li('Searching crRNA'),
                                    html.Li('Post processing'),
                                    html.Li('Merge Targets'),
                                    html.Li('Generating images'),
                                    html.Li('Populating database'),
                                    html.Li('Integrating results'),
                                    #html.Li('Annotating result'),
                                    #html.Li('Generating report')
                                ]
                            ),
                            style={'flex': '0 0 20%'}
                        ),
                        html.Div(
                            html.Ul(
                                [
                                    html.Li('To do', style={
                                            'color': 'red'}, id='index-status'),
                                    html.Li('To do', style={
                                            'color': 'red'}, id='search-status'),
                                    html.Li('To do', style={
                                            'color': 'red'}, id='post-process-status'),
                                    html.Li('To do', style={
                                            'color': 'red'}, id='merge-status'),
                                    html.Li('To do', style={
                                            'color': 'red'}, id='images-status'),
                                    html.Li('To do', style={
                                            'color': 'red'}, id='database-status'),
                                    html.Li('To do', style={
                                            'color': 'red'}, id='integrate-status'),
                                    #html.Li('To do', style = {'color':'red'}, id = 'annotate-result-status'),
                                    #html.Li('To do', style = {'color':'red'}, id = 'generate-report-status')
                                ],
                                style={'list-style-type': 'none'}
                            )
                        )
                    ],
                    className='flex-status'
                ),
                html.Div(
                    [
                        dcc.Link('View Results', style={
                                 'visibility': 'hidden'}, id='view-results', href=URL),
                        html.Div(id='no-directory-error')
                    ]
                )
            ],
            id='div-status-report'
        )
    )

    final_list.append(html.P('', id='done'))

    final_list.append(dcc.Interval(id='load-page-check', interval=3*1000))
    load_page = html.Div(final_list, style={'margin': '1%'})

    return load_page
