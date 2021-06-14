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
                                ('Samuele Cancellieri, Francesco Masillo, Elia Dirupo, Nicola Bombieri, and Rosalba Giugno,  InfOmics Lab, Department of Computer Science, University of Verona, Italy, '),
                                html.A('(https://infomics.github.io/InfOmics/index.html)',
                                       href='https://infomics.github.io/InfOmics/index.html', target='_blank')
                            ]
                        ),
                        html.Li(
                            [
                                ('Luca Pinello, Molecular Pathology Unit and Center for Cancer Research, Massachusetts General Hospital, Charlestown, MA, USA. Department of Pathology, Harvard Medical School, Boston, MA, USA. Broad Institute of MIT and Harvard, Cambridge, MA, USA. Pinello Lab '),
                                html.A(
                                    '(http://pinellolab.org/)', href='http://pinellolab.org/', target='_blank'),
                            ]
                        ),
                        html.Li(
                            [
                                ('Linda Yingqi Lin and Daniel E. Bauer, Division of Hematology/Oncology, Boston Children’s Hospital, Department of Pediatric Oncology, Dana-Farber Cancer Institute, Harvard Stem Cell Institute, Broad Institute, Department of Pediatrics, Harvard Medical School, Boston, Massachusetts 02115, USA '),
                                html.A(
                                    '(http://bauerlab.org/)', href='http://bauerlab.org/', target='_blank'),
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
                            'lpinello AT mgh DOT harvard DOT edu'
                        ),
                        html.Li(
                            'bauer AT bloodgroup DOT tch DOT harvard DOT edu'
                        ),
                    ], style={'padding': '15px'}
                ),
                html.Div(
                    [
                        html.P('Alternatively, please open an issue on GitHub: '),
                        html.A('https://github.com/pinellolab/CRISPRme/issues',
                               href='https://github.com/pinellolab/CRISPRme/issues', target='_blank'),
                    ]
                )
            ], style={'margin-left': '1%'}
        )
    )
    return f
