import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from app import URL
# from index import DISPLAY_HISTORY
PLOTLY_LOGO = 'assets/favicon.png'

# DISPLAY_OFFLINE = ''
DISPLAY_HISTORY=''
search_bar = dbc.Row(
    [
        #dbc.Col(dbc.Input(type="search", placeholder="Search")),
        dbc.Col(dbc.NavLink(
            html.A('HOME', href=URL + '/index', target='', style={
                   'text-decoration': 'none', 'color': 'white'}),
            active=True,
            className='testHover', style={'text-decoration': 'none', 'color': 'white', 'font-size': '1.5rem'})),
        # dbc.Col(dbc.NavLink(
        #     html.A('PERSONAL_DATA_MANAGEMENT', href=URL + '/genome-dictionary-management', target='',
        #            style={'text-decoration': 'none', 'color': 'white'}),
        #     active=True,
        #     className='testHover', style={'text-decoration': 'none', 'color': 'white', 'font-size': '1.5rem', 'display': DISPLAY_OFFLINE})),
        dbc.Col(dbc.NavLink(
            html.A('MANUAL', href=URL + '/user-guide', target='',
                   style={'text-decoration': 'none', 'color': 'white'}),
            active=True,
            className='testHover', style={'text-decoration': 'none', 'color': 'white', 'font-size': '1.5rem'})),
        dbc.Col(dbc.NavLink(
            html.A('CONTACTS', href=URL + '/contacts', target='',
                   style={'text-decoration': 'none', 'color': 'white'}),
            active=True,
            className='testHover', style={'text-decoration': 'none', 'color': 'white', 'font-size': '1.5rem'})),
        dbc.Col(dbc.NavLink(
            html.A('HISTORY', href=URL + '/history', target='',
                   style={'text-decoration': 'none', 'color': 'white'}),
            active=True,
            className='testHover', style={'text-decoration': 'none', 'color': 'white', 'font-size': '1.5rem', 'display': DISPLAY_HISTORY}))
    ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)


def Navbar():
    navbar = dbc.Navbar(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src=PLOTLY_LOGO, height="60px")),
                        dbc.Col(dbc.NavbarBrand(
                            "CRISPRme", className="ml-2", style={'font-size': '30px'}))
                    ],
                    align="center",
                    no_gutters=True,
                ),
                href=URL + '/index',
            ),
            dbc.NavbarToggler(id="navbar-toggler"),
            dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
        ],
        color="dark",
        dark=True,
    )
    return navbar
