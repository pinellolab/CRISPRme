"""Create the navigation bar displayed on top of CRISPRme 
website.
"""

from app import URL

import dash_bootstrap_components as dbc
import dash_html_components as html

import os


PLOTLY_LOGO = "assets/favicon.png"
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
                        style={"text-decoration": "none", "color": "white"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "white",
                        "font-size": "1.5rem",
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
                        style={"text-decoration": "none", "color": "white"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "white",
                        "font-size": "1.5rem",
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
                        style={"text-decoration": "none", "color": "white"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "white",
                        "font-size": "1.5rem",
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
                        style={"text-decoration": "none", "color": "white"},
                    ),
                    active=True,
                    className="testHover",
                    style={
                        "text-decoration": "none",
                        "color": "white",
                        "font-size": "1.5rem",
                    },
                )
            ),
        ],
        no_gutters=True,
        className="ml-auto flex-nowrap mt-3 mt-md-0",
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
    # construct the navigation bar
    navbar = dbc.Navbar(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src=PLOTLY_LOGO, height="60px")),
                        dbc.Col(
                            dbc.NavbarBrand(
                                "CRISPRme",
                                className="ml-2",
                                style={"font-size": "30px"},
                            )
                        ),
                    ],
                    align="center",
                    no_gutters=True,
                ),
                href=os.path.join(URL, "index"),
            ),
            dbc.NavbarToggler(id="navbar-toggler"),
            dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
        ],
        color="dark",
        dark=True,
    )
    return navbar
