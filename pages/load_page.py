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
     Output('no-directory-error', 'children'),
     Output('button-remove-result', 'hidden')],
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
                with open(current_job_dir + '.Params.txt') as f:
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

                if 'Merging Targets\tEnd' in current_log:
                    merge_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Merging Targets\tStart' in current_log:
                    merge_status = html.P(
                        'Processing...' + ' ' + 'Step [1/1]', style={'color': 'orange'})

                if 'Annotating results\tStart' in current_log:
                    images_status = html.P(
                        'Annotating... Step[1/2]', style={'color': 'orange'})

                if 'Creating images\tEnd' in current_log:
                    images_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Creating images\tStart' in current_log:
                    images_status = html.P(
                        'Generating images...' + ' ' + 'Step [2/2]', style={'color': 'orange'})

                if 'Integrating results\tEnd' in current_log:
                    integrate_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Integrating results\tStart' in current_log:
                    integrate_status = html.P(
                        'Processing...' + ' ' + 'Step [1/1]', style={'color': 'orange'})

                if 'Creating database\tEnd' in current_log:
                    database_status = html.P('Done', style={'color': 'green'})
                    all_done = all_done + 1
                elif 'Creating database\tStart' in current_log:
                    database_status = html.P(
                        'Inserting data...' + ' ' + 'Step [1/1]', style={'color': 'orange'})

                if os.path.isfile(f"{current_job_dir}/log_error.txt") and os.path.getsize(f"{current_job_dir}/log_error.txt") > 0:
                    return {'visibility': 'hidden'},  html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), '', dbc.Alert("The selected result encountered some errors, please remove it and try to submit again.", color="danger"), False
                if all_done == 7 or 'Job\tDone' in current_log:
                    return {'visibility': 'visible'}, index_status, search_status, post_process_status, merge_status, images_status, integrate_status, database_status, URL+'/result?job=' + dir_name.split('=')[-1], '', True
                else:
                    return {'visibility': 'hidden'}, index_status, search_status, post_process_status, merge_status, images_status, integrate_status, database_status, '', '', True
        elif 'queue.txt' in onlyfile:
            return {'visibility': 'hidden'}, html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('To do', style={'color': 'red'}), html.P('To do', style={'color': 'red'}), html.P('To do', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), '', dbc.Alert("Job submitted. Current status: in queue", color="info"), True
    return {'visibility': 'hidden'},  html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), html.P('Not available', style={'color': 'red'}), '', dbc.Alert("The selected result does not exist", color="danger"), True


@app.callback(
    Output('result-deleted', 'children'),
    [Input('button-remove-result', 'n_clicks')],
    [State('url', 'search')]
)
def removeResult(n, dir_name):
    if n == 0:
        raise PreventUpdate
    if n == 1:
        current_job_dir = current_working_directory + \
            'Results/' + dir_name.split('=')[-1]
        os.system(f'rm -rf {current_job_dir}')
        return html.P('Result deleted')
    return None


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
                        # html.P('Results will be kept available for 3 days')
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
                                    html.Li('Indexing genome(s)'),
                                    html.Li('Searching spacer'),
                                    html.Li('Post processing'),
                                    html.Li('Merge targets'),
                                    html.Li(
                                        'Annotating and generating images'),
                                    html.Li('Integrating results'),
                                    html.Li('Populating database'),
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
                        html.Div(id='no-directory-error'),
                        html.Div([html.Button(
                            'Remove result', id='button-remove-result', n_clicks=0, hidden=True)]),
                        html.Div(id='result-deleted')
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
