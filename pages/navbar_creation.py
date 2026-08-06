"""Create the navigation bar displayed on top of CRISPRme 
website.
"""

from app import URL

import dash_bootstrap_components as dbc
from dash import html

import os


PLOTLY_LOGO = "assets/crisprme-icon.svg"
CRISPRME_LOGO = "assets/crisprme-logo.png"  # full wordmark (correct Outfit font baked in)
DISPLAY_HISTORY = ""


def create_search_bar() -> dbc.Row:
    """Create search bar on top of CRISPRme webpage.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    dbc.Row
    """

    # start bar creation
    search_bar = dbc.Row(
        [
            dbc.Col(
                dbc.NavLink(
                    # home button
                    html.A(
                        "HOME",
                        href=os.path.join(URL, "index"),
                        target="",
                        style={"text-decoration": "none", "color": "#37474f"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "#37474f",
                        "font-size": "1.5rem",
                        "padding": "0 0.9rem",
                        "white-space": "nowrap",
                    },
                )
            ),
            dbc.Col(
                dbc.NavLink(
                    # help button (CRISPRme manual)
                    html.A(
                        "MANUAL",
                        href=os.path.join(URL, "user-guide"),
                        target="",
                        style={"text-decoration": "none", "color": "#37474f"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "#37474f",
                        "font-size": "1.5rem",
                        "padding": "0 0.9rem",
                        "white-space": "nowrap",
                    },
                )
            ),
            dbc.Col(
                dbc.NavLink(
                    # contacts button
                    html.A(
                        "CONTACTS",
                        href=os.path.join(URL, "contacts"),
                        target="",
                        style={"text-decoration": "none", "color": "#37474f"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "#37474f",
                        "font-size": "1.5rem",
                        "padding": "0 0.9rem",
                        "white-space": "nowrap",
                    },
                )
            ),
            dbc.Col(
                dbc.NavLink(
                    html.A(
                        # user history button
                        "HISTORY",
                        href=os.path.join(URL, "history"),
                        target="",
                        style={"text-decoration": "none", "color": "#37474f"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "#37474f",
                        "font-size": "1.5rem",
                        "padding": "0 0.9rem",
                        "white-space": "nowrap",
                    },
                )
            ),
            dbc.Col(
                dbc.NavLink(
                    html.A(
                        # settings / data manager button
                        "SETTINGS",
                        href=os.path.join(URL, "settings"),
                        target="",
                        style={"text-decoration": "none", "color": "#37474f"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "#37474f",
                        "font-size": "1.5rem",
                        "padding": "0 0.9rem",
                        "white-space": "nowrap",
                    },
                )
            ),
        ],
        className="g-0 ml-auto flex-nowrap mt-3 mt-md-0",
        align="center",
    )
    return search_bar


def navbar():
    """Create the navigation bar of CRISPRme website.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    dbc.Navbar
    """

    # create search bar
    search_bar = create_search_bar()
    # construct the navigation bar: a clean white bar so the full-colour CRISPRme+
    # wordmark logo reads well (rendered PNG keeps the exact Outfit font/colours)
    navbar = dbc.Navbar(
        [
            html.A(
                html.Img(src=CRISPRME_LOGO, height="52px"),
                href=os.path.join(URL, "index"),
                style={"text-decoration": "none"},
            ),
            dbc.NavbarToggler(id="navbar-toggler"),
            dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
        ],
        color="white",
        dark=False,
        style={
            "box-shadow": "0 2px 10px rgba(0, 0, 0, 0.08)",
            "border-bottom": "1px solid #eaeaea",
            "padding": "0.4rem 1.2rem",
        },
    )
    return navbar
