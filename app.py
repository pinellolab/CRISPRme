"""Define the Flask-based server app basics """

from flask_caching import Cache

import dash_bootstrap_components as dbc

import concurrent.futures  
import flask
import dash
import sys
import os


IPADDRESS = "127.0.0.1:8080"


def __start_message() -> None:
    """ (PRIVATE) 
    Write server start message to stderr.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    sys.stderr.write("SERVER STARTED\n")
    sys.stderr.write(f"GO TO {IPADDRESS} TO USE THE WEB APP\n\n")


URL = ""  # server URL
external_stylesheets = [
    "https://codepen.io/chriddyp/pen/bWLwgP.css", dbc.themes.BOOTSTRAP
]  # CSS stylesheet used to style the website
server = flask.Flask(__name__)  # define flask app.server
app = dash.Dash(
        __name__, 
        external_stylesheets=external_stylesheets, 
        suppress_callback_exceptions=True, 
        server=server
)  # initialize server app
app_directory = os.path.dirname(os.path.realpath(__file__))  # current location
__start_message()  # print server start message
current_working_directory = os.getcwd() + '/'  # This for files TODO: remove
app.title = "CRISPRme"  # assign flask app name
# necessary if update element in a callback generated in another callback
# app.config['suppress_callback_exceptions'] = True
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True
# define filtering operators used when querying tables
operators = [
    ["ge ", ">="], 
    ["le ", "<="], 
    ["lt ", "<"], 
    ["gt ", ">"], 
    ["ne ", "!="],
    ["eq ", "="],
    ["contains "]
] 
ONLINE = False  # NOTE change to True for online version, False for offline
DISPLAY_OFFLINE = ""
DISPLAY_ONLINE = ""
if ONLINE:
    DISPLAY_OFFLINE = "none"
else:
    DISPLAY_ONLINE = "none"
# set to execute multiple jobs at time
pool_executor = concurrent.futures.ProcessPoolExecutor(max_workers=2)
# configure caching
CACHE_CONFIG = {
    # try 'filesystem' if you don't want to setup redis
    "CACHE_TYPE": "filesystem",
    "CACHE_DIR": ("Cache")  # os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()  # initialize cache
cache.init_app(app.server, config=CACHE_CONFIG)  # start web-app

