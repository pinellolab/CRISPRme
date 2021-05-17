import dash
import dash_bootstrap_components as dbc
import os
import concurrent.futures  # For workers and queue
from flask_caching import Cache  # for cache of .targets or .scores
# from index import DISPLAY_HISTORY

URL = ''
external_stylesheets = [
    'https://codepen.io/chriddyp/pen/bWLwgP.css', dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
                suppress_callback_exceptions=True)

app_location = os.path.realpath(__file__)
print('started')
app_main_directory = os.path.dirname(app_location) + '/'  # This for scripts
current_working_directory = os.getcwd() + '/'  # This for files

app.title = 'CRISPRme'
# necessary if update element in a callback generated in another callback
# app.config['suppress_callback_exceptions'] = True
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

app_location = os.path.dirname(os.path.abspath(__file__)) + '/'
operators = [['ge ', '>='],
             ['le ', '<='],
             ['lt ', '<'],
             ['gt ', '>'],
             ['ne ', '!='],
             ['eq ', '='],
             ['contains ']]  # for filtering

ONLINE = False  # NOTE change to True for online version, False for offline
DISPLAY_OFFLINE = ''
DISPLAY_ONLINE = ''
if ONLINE:
    DISPLAY_OFFLINE = 'none'
    DISPLAY_ONLINE = ''
else:
    DISPLAY_OFFLINE = ''
    DISPLAY_ONLINE = 'none'

exeggutor = concurrent.futures.ProcessPoolExecutor(max_workers=2)

CACHE_CONFIG = {
    # try 'filesystem' if you don't want to setup redis
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': ('Cache')  # os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)
