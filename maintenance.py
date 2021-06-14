import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

import flask
import pandas as pd
import time
import os

server = flask.Flask('app')

app = dash.Dash('app', server=server)

app.layout = html.Div(
    [
        html.P('We are sorry, CRISPRme website is under maintenance, we will be back online as soon as possible. Thanks for your patience.'),
        html.Div(['In the meantime you can check our offline version on ', html.A(
            'Github', target='_blank', href='https://github.com/pinellolab/CRISPRme')]),
    ]
)
if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=80, debug=False,
                   dev_tools_ui=False, dev_tools_props_check=False)
