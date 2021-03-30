
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


def helpPage():
    final_list = []
    final_list.append(
        html.Div([
            html.H3('About'),
            html.P([
                'CRISPRme  performs  predictive analysis and result assessment on population and individual specific CRISPR/Cas experiments.' +
                ' CRISPRme enumerates on- and off-target accounting simultaneously for  substitutions, DNA/RNA bulges and common genetic variants from the 1000 genomes project.'
            ]),
            html.P(['Open this ', html.A('example', href='http://crispritz.di.univr.it/result?job=Q47PXDTBC8',
                                         target='_blank'), ' to navigate the results we show in this page'])

        ])

    )

    final_list.append(html.H3('Main Page'))

    final_list.append(
        html.Div([
            html.P(
                ['In the main page of CRISPRme, users can select a wide range of options to personalize their searches. The input phase is divided into three main steps:',
                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open('assets/main_page.png', 'rb').read()).decode()), width="100%", height="auto")]
            ),
            html.Ul(
                [
                    html.Li(
                        [
                            'STEP 1: Genome and PAM selection',
                            html.Ul(
                                [
                                    html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                        open('assets/helpPage/genome.png', 'rb').read()).decode()), width='30%'),
                                    html.Li(
                                        'Genome: here you can select a genome from the ones present, some genomes have also the variant version enriched (indicated with a \'+\' symbol) with genetic variant from the 1000genome project'),
                                    html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                        open('assets/helpPage/pam.png', 'rb').read()).decode()), width='30%'),
                                    html.Li(
                                        'PAM: here you can select a Protospacer Adjacent Motif for a specific Cas protein.'),
                                ], style={'padding': '15px'}
                            )
                        ]
                    ),
                    html.Li(
                        [
                            'STEP 2: Guide insertion and threshold configuration',
                            html.Ul(
                                [
                                    html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                        open('assets/helpPage/guides.PNG', 'rb').read()).decode()), width='40%'),
                                    html.Li(
                                        'Guides: a list of crRNAs sequences, consisting in 1 or more sequences (max 1000 sequences) to search on the genome'),
                                    html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                        open('assets/helpPage/sequence.PNG', 'rb').read()).decode()), width='40%'),
                                    html.Li('Sequence: one or more genetic sequences (max 1000 characters), each sequence MUST BE separated with the header \'>name\'. The sequence can be also submitted with a ' +
                                            'chromosome range, also provided with an header. The region will be extracted from the Genome selected in STEP 1'),
                                    html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                        open('assets/helpPage/crRNA.PNG', 'rb').read()).decode()), width='20%'),
                                    html.Li(
                                        'Allowed mismatches: number of tolerated mismatches in a target'),
                                    html.Li(
                                        'Bulge DNA size: size of bubbles tolerated on the DNA sequence (can be consecutive(AA--AA) or interleaved(AA-A-AA)).'),
                                    html.Li(
                                        'Bulge RNA size: size of bubbles tolerated on the RNA sequence (can be consecutive(AA--AA) or interleaved(AA-A-AA))'),
                                    # html.Img(src = 'data:image/png;base64,{}'.format(base64.b64encode(open('assets/helpPage/crRNA.PNG', 'rb').read()).decode()), width='20%' ),
                                    html.Li(
                                        'crRNA length: available only when a genetic sequence is given as input, represents the length of the guides (without PAM) that you want to extract from the sequence.')
                                ], style={'padding': '15px'}
                            )
                        ]
                    ),
                    html.Li(
                        [
                            'STEP 3 (Advanced options): Select various comparisons',
                            html.Ul(
                                [
                                    html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                        open('assets/helpPage/advOpt.PNG', 'rb').read()).decode()), width='40%'),
                                    html.Li(
                                        'Compare your results with the GeCKO v2 library: selected by default, compares the results of your guides with the results obtained in a previous search with guides from the well-known GeCKO library.'),
                                    html.Li('Compare your results with the corresponding reference genome: selected by default when an enriched genome is chosen, compares the results with the respective reference genome to evaluate differences when variant are added.'),
                                    html.Li(
                                        'Notify me by email: if selected, let you insert an email to receive a notification when your job is terminated.'),
                                ], style={'padding': '15px'}
                            )
                        ]
                    )
                ], style={'padding': '15px'}
            )
        ])
    )

    final_list.append(
        html.P(
            ['After selecting the desired inputs, click on the Submit button to start the search']
        )
    )

    final_list.append(
        html.Div(
            [
                dbc.Alert(
                    [
                        'WARNING: If some inputs are missing, a warning popup will be displayed', html.P(),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                            open('assets/helpPage/warning.png', 'rb').read()).decode()), width='100%'),
                    ], color='warning', fade=False, style={'width': '70%'}
                )
            ]
        )
    )

    final_list.append(
        html.P(
            [
                'After the submission, the status of the search will be displayed in a new page',
                html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    open('assets/waitPage/loadPage.png', 'rb').read()).decode()), width='100%'),
                'When the job is complete, the result link will appear at the end of the status report',
                html.P(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    open('assets/waitPage/jobDone.png', 'rb').read()).decode()), width='20%'))
            ]
        )
    )

    final_list.append(html.H3('Result Page'))
    final_list.append(
        html.P(
            [
                'At the top of the page, you find a table with the list of sgRNAs used during the search phase. This table summarizes the results obtained for each input guide.',
                html.P(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                    'assets/resultPage/resultsSummary.png', 'rb').read()).decode()), width='100%')),
                html.Ul(
                    [
                        html.Li('CFD: Off-Target Cutting Frequency Determination Score, calculates how much is the affinity of the guides with the off-targets, basically tells you the likelihood of the guide to perform cut in off-target regions.'),
                        html.Li('Doench 2016: On-Target Efficacy Scoring (Azimuth 2.0), it’s a trained machine learning model that gives you the likelihood of on-target activity for the selected guide.'),
                        html.Li(
                            'On-Targets Reference: shows how many possible On-Targets the guide can have in the Reference Genome.'),
                        html.Li(['Samples in Class 0 - 0+ - 1 - 1+: shows the number of samples grouped by Sample Class:',
                                 html.Ul([
                                     html.Li(
                                         'Class 0: Samples that does not have any On-Targets'),
                                     html.Li(
                                         'Class 0+: Samples that have a subset of the Reference Genome On-Targets'),
                                     html.Li(
                                         'Class 1: Samples that have the same On-Targets as the Reference Genome'),
                                     html.Li(
                                         'Class 1+: Samples that creates at least a new On-Target, that is not present in the Reference Genome')
                                 ])
                                 ]),
                        html.Li('Off-Targets Reference (0 - n Mismatches + Bulges): shows how many possible Off-Targets the guide can have in the Reference Genome. Targets are also grouped by Mismatch + Bulge value.'),
                        html.Li('Off-Targets Enriched (0 - n Mismatches + Bulges): shows how many possible Off-Targets the guide can have in the Enriched Genome. Targets are also grouped by Mismatch + Bulge value.')
                    ], style={'padding': '15px'}
                ),
                html.P(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                    'assets/resultPage/populationDistribution.png', 'rb').read()).decode()), width='100%')),
                'The Show Target Distribution in Populations button opens a section where informations about the number of targets found in each Superpopulation (EAS, EUR, AFR, AMR, SAS) are provided by means of a barplot for each Mismatch + Bulge value. ',
                html.P('In the middle of the page there are four tabs:'),
                html.Ul(
                    [
                        html.Li([html.Span('Summary by Guide: ', style={
                                'color': 'red'}), 'This table collects all the possible On-/Off- Targets grouped by mismatch/bulge couples.']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                            'assets/resultPage/summaryByGuide.png', 'rb').read()).decode()), width='100%'),
                        html.Ul(
                            [
                                html.Li(
                                    'Bulge Type: type of bulge of the targets, can be X (no bulge), DNA or RNA'),
                                html.Li(
                                    'Bulge Size: size of the bulge present in the targets'),
                                html.Li(
                                    'Mismatches: number of mismatches present in the targets'),
                                html.Li(
                                    'Targets in Reference: number of targets found in the Reference Genome for that combination of mismatch/bulge'),
                                html.Li(
                                    'Targets in Enriched: number of targets found in the Enriched Genome for that combination of mismatch/bulge. Each of these targets is associated with at least one sample'),
                                html.Li(
                                    'PAM Creation: number of possible created PAMs due to variants addition'),
                                html.Li(
                                    'Show Targets: open a new page to display all the targets of the row of interest')
                            ]
                        ),
                        html.Li([html.Span('Summary by Sample: ', style={
                                'color': 'red'}), 'This table collects all the possible On-/Off- Targets grouped by sample name.']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                            'assets/resultPage/summaryBySamples.png', 'rb').read()).decode()), width='100%'),
                        html.Ul(
                            [
                                html.Li(
                                    'Population: population which the sample belong to'),
                                html.Li(
                                    'Super Population: continent which the sample belong to'),
                                html.Li(
                                    'Targets in Enriched: number of targets found in the Enriched Genome that are generated by that sample'),
                                html.Li(
                                    'Targets in Population: number of targets found in the Enriched Genome that are generated by all the sample of the population'),
                                html.Li(
                                    'Targets in Super Population: number of targets found in the Enriched Genome that are generated by all the populations'),
                                html.Li(
                                    'PAM Creation: number of possible created PAMs due to variants addition'),
                                html.Li(
                                    'Class: Sample Class (0 - 0+ - 1 - 1+) associated with the sample'),
                                html.Li(
                                    'Show Targets: open a new page to display all the targets of the row of interest')
                            ]
                        ),
                        html.Li([html.Span('Summary by Position: ', style={
                                'color': 'red'}), 'This table collects all the possible On-/Off- Targets grouped by position in the genome (composed by chromosome and relative position)']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                            'assets/resultPage/summaryByPosition.png', 'rb').read()).decode()), width='100%'),
                        html.Ul(
                            [
                                html.Li(
                                    'Position: chromosome relative position of the first letter of the guide'),
                                html.Li(
                                    'Best Target: best target found in that position'),
                                html.Li(
                                    'Min Mismatch: minimum number of mismatches present in the targets in that position'),
                                html.Li(
                                    'Min Bulge: minimum number of bulges present in the targets in that position'),
                                html.Li(
                                    'Bulge: number of bulges present in the targets in that position'),
                                html.Li(
                                    'Targets in Cluster by Mismatch Value: Matrix showing the distribution of the targets grouped by mismatch/bulge count'),
                            ]

                        ),
                        html.Li([html.Span('Graphical Report: ', style={
                                'color': 'red'}), 'This page shows graphics about a specific guide, including genomic annotation and motif logos. The main feature introduced is the possibility to visualize graphical reports at individual level.']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(
                            'assets/resultPage/summaryByGraphicaGecko.png', 'rb').read()).decode()), width='100%'),
                        html.Ul(
                            [
                                html.Li(
                                    'Select a Mismatch Value: generate graphics with the specified mismatch value'),
                                html.Li(
                                    'Select Individual Data: generate individual data, by selecting Super Population, Population and Sample')
                            ]
                        )
                    ], style={'padding': '15px'}
                )

            ]
        )
    )

    final_list.append(
        html.Div(
            [
                html.H3('Browser Compatibility'),
                html.Div([
                    dash_table.DataTable(
                        data=[{'OS': 'Linux', 'V': 'Ubuntu 18.04.2 LTS', 'Ch': '79.0', 'S': 'n/a', 'ME': 'n/a', 'F': '71.0'},
                              {'OS': 'MacOS', 'V': 'Mojave', 'Ch': ' 79.0',
                               'S': '13.0.4', 'ME': 'n/a', 'F': '71.0'},
                              {'OS': 'Windows', 'V': '10', 'Ch': '79.0', 'S': 'n/a', 'ME': '44.18362.449.0', 'F': '71.0'}],

                        columns=[{'id': 'OS', 'name': 'Operative System'}, {'id': 'V', 'name': 'Version'}, {'id': 'Ch', 'name': 'Chrome'},
                                 {'id': 'S', 'name': 'Safari'}, {'id': 'ME', 'name': 'Microsoft Edge'}, {
                                     'id': 'F', 'name': 'Firefox'}
                                 ],

                        style_cell={
                            'textAlign': 'center',
                            'width': '20'
                        },

                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': 'rgb(248, 248, 248)'
                            }
                        ],
                        style_header={
                            'backgroundColor': 'rgb(230, 230, 230)',
                            'fontWeight': 'bold'
                        }
                    )
                ],   style={'display': 'inline-block', 'width': '48%'})
            ]
        )
    )
    page = html.Div(final_list, style={'margin': '1%'})
    return page
