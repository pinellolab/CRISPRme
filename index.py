#!/usr/bin/env python
"""Main module for the CRISPRme web application.

This module initializes the web application, sets up the server, and defines the 
layout and callbacks for page navigation. It handles the routing of different pages 
based on the URL and manages the server's running mode.

Attributes:
    MODEFILE (str): The filename for storing the running mode.
    HOST (str): The host address for the server.
    PORTWEB (int): The port number for the website.
    PORTLOCAL (int): The port number for the local server.
    server (Flask): The Flask server instance.
    navbar (html.Div): The navigation bar component of the web application.
    app.layout (html.Div): The layout structure of the Dash application.
"""


from pages import (
    main_page,
    navbar_creation,
    results_page,
    load_page,
    history_page,
    help_page,
    contacts_page,
)
from app import (
    URL,
    IPADDRESS,
    WEBADDRESS,
    ONLINE,
    app,
    current_working_directory,
    cache,
)

from dash.dependencies import Input, Output, State
from typing import Tuple

import dash_core_components as dcc
import dash_html_components as html

import sys
import os


MODEFILE = ".mode_type.txt"  # running mode report file
HOST = "0.0.0.0"  # server host
PORTWEB = 80  # website port
PORTLOCAL = 8080  # local server port
CRISPRME_DIRS = [
    "Genomes",
    "Results",
    "Dictionaries",
    "VCFs",
    "Annotations",
    "PAMs",
    "samplesIDs",
]  # crisprme directory tree

# initialize the webpage
server = app.server  # start server
navbar = navbar_creation.navbar()  # create navigation bar on top of the webpage
# display multipage website
app.layout = html.Div(
    [
        navbar,
        dcc.Location(id="url", refresh=False),
        html.Div(id="page-content"),
        html.P(id="signal", style={"visibility": "hidden"}),
    ]
)


def check_directories(basedir: str) -> None:
    """Check and create necessary directories in the specified base directory.

    This function verifies if the provided base directory exists and creates
    any missing subdirectories defined in the CRISPRME_DIRS list.

    Args:
        basedir (str): The base directory to check and create subdirectories in.

    Raises:
        TypeError: If basedir is not a string.
        FileNotFoundError: If the base directory does not exist.
    """

    if not isinstance(basedir, str):
        raise TypeError(f"Expected {str.__name__}, got {type(basedir).__name__}")
    if not os.path.exists(basedir):
        raise FileNotFoundError(f"Unable to locate {basedir}")
    for d in CRISPRME_DIRS:
        if not os.path.exists(os.path.join(basedir, d)):
            os.makedirs(os.path.join(basedir, d))


# switch between the website pages
@app.callback(
    [Output("page-content", "children"), Output("job-link", "children")],
    [Input("url", "href"), Input("url", "pathname"), Input("url", "search")],
    [State("url", "hash")],
)
def change_page(href: str, path: str, search: str, hash_guide: str) -> Tuple:
    """Handles page changes based on the current URL and its components.

    This callback function updates the content of the web application based on
    the provided URL parameters. It determines which page to display and returns
    the corresponding content and link based on the path and search parameters.

    Args:
        href (str): The full URL of the current page.
        path (str): The path component of the URL.
        search (str): The query string of the URL.
        hash_guide (str): The hash component of the URL.

    Returns:
        Tuple: A tuple containing the children for the page content and the job
            link.

    Raises:
        TypeError: If any of the input parameters are not of type str.
    """

    if not isinstance(href, str):
        raise TypeError(f"Expected {str.__name__}, got {type(href).__name__}")
    if not isinstance(path, str):
        raise TypeError(f"Expected {str.__name__}, got {type(path).__name__}")
    if not isinstance(search, str):
        raise TypeError(f"Expected {str.__name__}, got {type(search).__name__}")
    if not isinstance(hash_guide, str):
        raise TypeError(f"Expected {str.__name__}, got {type(hash_guide).__name__}")
    if path == "/load":  # show loading page
        # define url to display on load page to check on job status
        # if online show the webaddress, show the ip address otherwise
        job_loading_url = WEBADDRESS if ONLINE else IPADDRESS
        return (load_page.load_page(), f"{job_loading_url}/load{search}")
    if path == "/result":  # display results page
        job_id = search.split("=")[-1]  # recover job id from url
        if not hash_guide or hash_guide is None:
            return results_page.result_page(job_id), os.path.join(URL, "load", search)
        elif "new" in hash_guide:  # targets table tab
            return (
                results_page.guidePagev3(job_id, hash_guide.split("#")[1]),
                os.path.join(URL, "load", search),
            )
        elif "-Sample-" in hash_guide:  # sample tab
            return (
                results_page.sample_page(job_id, hash_guide.split("#")[1]),
                os.path.join(URL, "load", search),
            )
        elif "-Pos-" in hash_guide:  # genomic region tab
            return (
                results_page.cluster_page(job_id, hash_guide.split("#")[1]),
                os.path.join(URL, "load", search),
            )
        return results_page.result_page(job_id), os.path.join(URL, "load", search)
    if path == "/user-guide":  # display manual page
        return help_page.helpPage(), os.path.join(URL, "load", search)
    if path == "/contacts":  # display contacts page
        return contacts_page.contact_page(), os.path.join(URL, "load", search)
    if path == "/history":  # display results history page
        return history_page.history_page(), os.path.join(URL, "load", search)
    return main_page.index_page(), "/index"  # display main page


def index():
    """Starts the CRISPRme web application server.

    This function checks the directory structure for consistency, determines the
    running mode (local or website), and starts the Dash application server
    accordingly. It also handles the creation of a mode file to track the current
    running mode and clears the cache before starting the server.

    Args:
        None

    Returns:
        None

    Raises:
        OSError: If there is an issue writing the mode file.
    """

    # check CRISPRme directory tree consistency
    check_directories(current_working_directory)
    # TODO: replace using argparse in crisprme.py
    debug = "--debug" in sys.argv[1:]  # check if debug mode is active
    website = "--website" in sys.argv[1:]  # check if local server or website
    try:  # keep track of the running mode (debugging purposes)
        modefname = os.path.join(current_working_directory, MODEFILE)
        with open(modefname, mode="w") as outfile:
            mode = "server" if website else "local"
            outfile.write(mode)
    except IOError as e:
        raise OSError("Cannot write mode file") from e
    if website:  # online web-interface running
        app.run_server(
            host=HOST,
            port=PORTWEB,  # type: ignore
            debug=debug,
            dev_tools_ui=debug,
            dev_tools_props_check=debug,
        )
    else:  # local web-interface running
        app.run_server(
            host=HOST,
            port=PORTLOCAL,  # type: ignore
            debug=debug,
            dev_tools_ui=debug,
            dev_tools_props_check=debug,
        )
    cache.clear()  # clear cache once server is closed


if __name__ == "__main__":
    index()
