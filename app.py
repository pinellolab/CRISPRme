"""Main module for the CRISPRme web application.

This module sets up the Flask and Dash applications, configures the server, and initializes the web application. It also defines stylesheets, caching, and various application settings.

Attributes:
    WEBADDRESS (str): The web address for accessing the application.
    IPADDRESS (str): The IP address and port for local access to the web application.
    URL (str): The server URL.
    external_stylesheets (list): List of external stylesheets for the web application.
    server (Flask): The Flask server instance.
    app (Dash): The Dash application instance.
    app_directory (str): The directory of the current application.
    current_working_directory (str): The current working directory.
    operators (list): List of filtering operators used for querying tables.
    ONLINE (bool): Flag indicating if the application is online or offline.
    DISPLAY_OFFLINE (str): CSS display property for offline mode.
    DISPLAY_ONLINE (str): CSS display property for online mode.
    pool_executor (ProcessPoolExecutor): Executor for running multiple jobs 
        concurrently.
    CACHE_CONFIG (dict): Configuration settings for caching.
    cache (Cache): Cache instance for the application.
"""

from flask_caching import Cache

import dash_bootstrap_components as dbc

import concurrent.futures
import flask
import dash
import sys
import os

WEBADDRESS = "http://crisprme.di.univr.it"
IPADDRESS = "127.0.0.1:8080"


def _start_message() -> None:
    """Prints a startup message to the standard error stream.

    This function outputs a message indicating that the server has started and
    provides the URL to access the web application.

    Args:
        None

    Returns:
        None
    """

    sys.stderr.write("SERVER STARTED\n")
    sys.stderr.write(f"GO TO http://{IPADDRESS} TO USE THE WEB APP\n\n")


URL = ""  # server URL
external_stylesheets = [
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
    dbc.themes.BOOTSTRAP,
]  # CSS stylesheet used to style the website
server = flask.Flask(__name__)  # define flask app.server
app = dash.Dash(
    __name__,
    external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True,
    server=server,  # type: ignore
)  # initialize server app
app_directory = os.path.dirname(os.path.realpath(__file__))  # current location
_start_message()  # print server start message
current_working_directory = os.getcwd() + "/"  # This for files TODO: remove
app.title = "CRISPRme"  # type: ignore - assign flask app name
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
    ["contains "],
]
ONLINE = True  # NOTE change to True for online version, False for offline
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
    "CACHE_DIR": ("Cache"),  # os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()  # initialize cache
cache.init_app(app.server, config=CACHE_CONFIG)  # start web-app
