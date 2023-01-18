"""
"""

# from pages import (
#     main_page, 
#     navbar_creation, 
#     results_page, 
#     load_page, 
#     history_page, 
#     help_page, 
#     contacts_page
# )
# from app import app, URL, current_working_directory, cache
# from utils import check_directories

# from dash.dependencies import Input, Output, State
# from typing import Tuple

# import dash_core_components as dcc
# import dash_html_components as html

# import sys
# import os


# MODEFILE = ".mode_type.txt"  # running mode report file
# HOST = "0.0.0.0"  # server host
# PORTWEB = 80  # website port
# PORTLOCAL = 8080  # local server port

# # initialize the webpage
# server = app.server  # start server
# navbar = navbar_creation.navbar()  # create navigation bar on top of the webpage
# # display multipage website


# # switch between the website pages
# @app.callback(
#     [Output("page-content", "children"), Output("job-link", "children")],
#     [Input("url", "href"), Input("url", "pathname"), Input("url", "search")],
#     [State("url", "hash")]
# )
# def change_page(href: str, path: str, search: str, hash_guide: str) -> Tuple:
#     """The function switches between the selected pages.

#     ...

#     Parameters
#     ----------
#     href : str
#         URL
#     path : str
#         Current path
#     search : str
#         Current search
#     hash_guide : str
#         Guide hash

#     Returns
#     -------
#     Tuple
#     """

#     if not isinstance(href, str):
#         raise TypeError(f"Expected {str.__name__}, got {type(href).__name__}")
#     if not isinstance(path, str):
#         raise TypeError(f"Expected {str.__name__}, got {type(path).__name__}")
#     if not isinstance(search, str):
#         raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
#     if not isinstance(hash_guide, str):
#         raise TypeError(f"Expected {str.__name__}, got {type(hash_guide).__name__}")
#     if path == "/load":  # show loading page
#         return (
#                 load_page.load_page(), 
#                 os.path.join("".join(href.split("/")[:-1]), "/load", search)
#         )
#     if path == "/result":  # display results page
#         job_id = search.split("=")[-1]
#         if not hash_guide or hash_guide is None:
#             return results_page.result_page(job_id), os.path.join(URL, "load", search)
#         elif "new" in hash_guide:  # TODO: change name to guide page
#             return (
#                 results_page.guidePagev3(job_id, hash_guide.split("#")[1]),
#                 os.path.join(URL, "load", search)
#             )
#         elif "-Sample-" in hash_guide:
#             return (
#                 results_page.sample_page(job_id, hash_guide.split("#")[1]),
#                 os.path.join(URL, "load", search)
#             )
#         elif "-Pos-" in hash_guide:
#             return (
#                 results_page.cluster_page(job_id, hash_guide.split("#")[1]),
#                 os.path.join(URL, "load", search)
#             )
#         return results_page.result_page(job_id), os.path.join(URL, "load", search)

#     if path == "/user-guide":  # display manual page
#         return help_page.helpPage(), os.path.join(URL, "load", search)
#     if path == "/contacts":  # display contacts page
#         return contacts_page.contact_page(), os.path.join(URL, "load", search)
#     if path == "/history":  # display results history page
#         return history_page.history_page(), os.path.join(URL, "load", search)
#     if path == "/index":  # display main page
#         return main_page.index_page(), "/index"
#     return main_page.index_page(), "/index"


# def index():
#     """
#     """

#     # check CRISPRme directory tree consistency
#     check_directories(current_working_directory)
#     # check if debug mode is active 
#     # TODO: replace using argparse in crisprme.py
#     debug = False
#     if "--debug" in sys.argv[1:]:
#         debug = True
#     # check if local server or website
#     website = False
#     if "--website" in sys.argv[1:]:
#         website = True
#     # keep track of running mode
#     try:
#         handle = open(os.path.join(current_working_directory, MODEFILE), mode="w")
#         if website:  # 'server' mode
#             handle.write("server")
#         else:  # 'local' mode
#             handle.write("local")
#     except OSError:
#         raise OSError(f"An error occurred while writing {MODEFILE}")
#     finally:
#         handle.close()
#     if website:
#         app.run_server(
#             host=HOST, port=PORTWEB, debug=debug, dev_tools_ui=debug, dev_tools_props_check=debug
#         )
#         cache.clear()  # clear cache once server is closed
#     else:  # local web-interface running
#         app.run_server(
#             host=HOST, port=PORTLOCAL, debug=debug, dev_tools_ui=debug, dev_tools_props_check=debug
#         )
#         app.run_server(
#             host=HOST, port=PORTWEB, debug=debug, dev_tools_ui=debug, dev_tools_props_check=debug
#         )
#         cache.clear()  # clear cache once server is closed

#!/usr/bin/env python
# DASH IMPORT
from pickle import TRUE
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
# SYSTEM IMPORT
import os
import sys
# APP IMPORT
from app import app, URL, current_working_directory, cache
# PAGES IMPORT
# from pages import UpdateDict
from pages import main_page
from pages import navbar_creation
# from pages import GUImessage
# from pages import creazione_dizionari
# from pages import ChooseFiles
# from pages import annotations
# from pages import send_mail
from pages import results_page
from pages import load_page
from pages import history_page
from pages import help_page
# from pages import personalization_page
from pages import contacts_page
# from pages import genome_database


server = app.server  # I add this part here

navbar = navbar_creation.navbar()
# For multipage
app.layout = html.Div([
    navbar,
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content'),
    html.P(id='signal', style={'visibility': 'hidden'})
])


def directoryCheck():
    # function to check the main directory status, if some directory is missing, create it
    directoryList = ['Genomes', 'Results', 'Dictionaries',
                     'VCFs', 'Annotations', 'PAMs', 'samplesIDs']
    for directory in directoryList:
        if not os.path.exists(current_working_directory+directory):
            os.makedirs(current_working_directory+directory)


@app.callback(
    [Output('page-content', 'children'),
     Output('job-link', 'children')],
    [Input('url', 'href'),
     Input('url', 'pathname'),
     Input('url', 'search')],
    [State('url', 'hash')]
)
def changePage(href, path, search, hash_guide):
    '''
    Controllo della pagina da mostrare in base all'url
    '''

    # print(href)
    if path == '/load':
        # print("CHANGING TO LOAD PAGE")
        # return load_page.load_page(), href.strip().split('/')[0] + '/load' + search
        return load_page.load_page(), ''.join(href.split('/')[:-1]) + '/load' + search
    if path == '/result':
        job_id = search.split('=')[-1]
        if hash_guide is None or hash_guide == '':
            return results_page.result_page(job_id), URL + '/load' + search
        if 'new' in hash_guide:  # TODO cambiare nome alla pagina delle guide
            return results_page.guidePagev3(job_id, hash_guide.split('#')[1]), URL + '/load' + search
        if '-Sample-' in hash_guide:
            return results_page.sample_page(job_id, hash_guide.split('#')[1]), URL + '/load' + search
        if '-Pos-' in hash_guide:
            return results_page.cluster_page(job_id, hash_guide.split('#')[1]), URL + '/load' + search
        return results_page.result_page(job_id), URL + '/load' + search
    # if path == '/genome-dictionary-management':
        # return personalization_page.genomeAndDictionaryManagement(), URL + '/load' + search
    if path == '/user-guide':
        return help_page.helpPage(), URL + '/load' + search
    if path == '/contacts':
        return contacts_page.contact_page(), URL + '/load' + search
    if path == '/history':
        return history_page.history_page(), URL + '/load' + search
    # if path == '/genomes':
    #     genomes_page = html.Div(genome_database.get_genomes(
    #         current_working_directory), style={'margin': '1%'})
        # return genomes_page, URL + '/load' + search
    if path == '/index':
        return main_page.index_page(), '/index'
    return main_page.index_page(), '/index'


if __name__ == '__main__':
    # check if all necessary directories exist
    directoryCheck()
    # keep track of mode of activation (local or server)
    mode_file = open(current_working_directory+'.mode_type.txt', 'w')
    # debug mode
    debug_mode = False
    if '--debug' in sys.argv[1:]:
        debug_mode = True
    if '--website' in sys.argv[1:]:
        # write 'server' in mode_file
        mode_type = 'server'
        mode_file.write(mode_type)
        mode_file.close()
        app.run_server(host='0.0.0.0', port=80, debug=debug_mode,
                       dev_tools_ui=debug_mode, dev_tools_props_check=debug_mode)

        cache.clear()  # delete cache when server is closed
    else:
        # write 'local' in mode_file
        mode_type = 'local'
        mode_file.write(mode_type)
        mode_file.close()
        app.run_server(host='0.0.0.0', port=8080, debug=debug_mode,
                       dev_tools_ui=debug_mode, dev_tools_props_check=debug_mode)
        cache.clear()  # delete cache when server is close