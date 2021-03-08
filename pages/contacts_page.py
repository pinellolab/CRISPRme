import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc


def contactPage():
    f = []
    f.append(
        html.Div(
            [
                html.H3('CONTACTS'),
                html.P(
                    'CRISPRme was developed by:'
                ),
                html.Ul(
                    [
                        html.Li(
                            [
                                ('Elia Dirupo, Samuele Cancellieri, Nicola Bombieri and Rosalba Giugno in the Department of Computer Science, University of Verona, Italy, InfOmics Lab ('),
                                html.A('https://infomics.github.io/InfOmics/index.html',
                                       href='https://infomics.github.io/InfOmics/index.html', target='_blank'),
                                ')'
                            ]
                        ),
                        html.Li(
                            [
                                ('Luca Pinello Molecular Pathology Unit, Center for Computational and Integrative Biology, Center for Cancer Research, Massachusetts General Hospital, Charlestown, MA, USA. Department of Pathology, Harvard Medical School, Boston, MA, USA. Broad Institute of MIT and Harvard, Cambridge, MA, USA. Pinello Lab ('),
                                html.A(
                                    'http://pinellolab.org/', href='http://pinellolab.org/', target='_blank'),
                                ')'
                            ]
                        ),
                    ], style={'padding': '15px'}
                ),
                html.P(
                    'Please send any comment or bug to:'
                ),
                html.Ul(
                    [
                        html.Li(
                            'rosalba DOT giugno AT univr DOT it'
                        ),
                        html.Li(
                            'lpinello AT jimmy DOT harvard DOT edu'
                        )
                    ], style={'padding': '15px'}
                )
            ]
        )
    )
    return f
