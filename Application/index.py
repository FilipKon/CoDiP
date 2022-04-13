import dash_bootstrap_components as dbc
import dash_uploader as du
from dash import callback, Input, Output, dcc, html, Dash

from apps import overview, bindingsites, contact

# styling the sidebar
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

# padding for the page content
CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

external_stylesheets = [
    {
        "href": "https://fonts.googleapis.com/css2?"
                "family=Lato:wght@400;700&display=swap",
        "rel": "stylesheet",
    },
]
app = Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "CoDiP: Conformational Difference of Protein structures"
app._favicon = ("D:\\ConFam\\logo_codip.png")
# app = Dash(__name__, external_stylesheets=[dbc.themes.JOURNAL])
server = app.server


nav_item = dbc.Nav(
    [
        dbc.NavItem(dbc.NavLink("Overview", active="exact", href="/apps/overview")),
        dbc.NavItem(dbc.NavLink("Input", active="exact", href="/apps/input")),
        dbc.NavItem(dbc.NavLink("Global", active="exact", href="/apps/globala")),
        dbc.NavItem(dbc.NavLink("Binding sites", active="exact", href="/apps/bindingsites")),
        dbc.NavItem(dbc.NavLink("Contact", active="exact", href="/apps/contact"))
    ],
    fill=False,
    vertical=True,
)

navbar = dbc.Navbar(
    html.A(
        dbc.Row([
                dbc.Col(dbc.Col(html.Img(src=app.get_asset_url("logo_codip_new.png"), height="30px"))),
                dbc.Col(dbc.NavbarBrand(html.H1("CoDiP", className="ml-2")))])
    ),
    dbc.NavbarToggler(id="navbar-toggler"),
    dbc.Nav([nav_item], navbar=True, className="mr-auto"),
    color="dark",
    dark=True,
    className="menu")

sidebar = html.Div([
    html.H3("CAP", className='display-4'),
    html.Hr(),
    html.P("Choose type of analyzes", className="lead"),
    dbc.Nav(
        [
            dbc.NavLink("Overview", href='/apps/overview', active="exact"),
            dbc.NavLink("Input", href='/apps/input', active="exact"),
            dbc.NavLink("Global", href='/apps/globala', active="exact"),
            dbc.NavLink("Binding sites", href='/apps/bindingsites', active="exact"),
            # dbc.NavLink("Tryouts", href='/apps/raw_function_tryouts', active="exact"),
            dbc.NavLink("Contact", href='/apps/contact', active="exact")
        ],
        vertical=True,
        pills=True,
    ),
], style=SIDEBAR_STYLE)


app.layout = html.Div([
    html.Div(
        children=[navbar],
        className="menu",
    ),
    html.Div(children=[
        html.P(children="![](logo_codip_new.png)", className="header-emoji"),
        html.H1(
            children="CoDiP", className="header-title"
        ),
        html.P(
            children="Analyze protein structures"
                     " and their similarities by calculating "
                     " distances between residues and performing PCA on that",
            className="header-description",
        ),
    ],
        className="header", ),
    dcc.Location(id='url', refresh=False),
    #navbar,
    html.Div(id='page-content_newest_1yo', children=[], style=CONTENT_STYLE),
    dcc.Store(id='global_link'),
    dcc.Store(id='global_data'),
    dcc.Store(id='global_labels'),
    dcc.Store(id='binding_site_3_link'),
    dcc.Store(id='binding_site_4_link'),
    dcc.Store(id='binding_site_5_link')
])

UPLOAD_FOLDER = r"/tmp"
du.configure_upload(app, UPLOAD_FOLDER)


@callback(Output('page-content_newest_1yo', 'children'),
          [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/apps/overview':
        return overview.layout
    #elif pathname == '/apps/input':
    #    return input.layout
    #elif pathname == '/apps/globala':
    #    return globala.layout
    elif pathname == '/apps/bindingsites':
        return bindingsites.layout
    elif pathname == '/apps/contact':
        return contact.layout
    else:
        return overview.layout


if __name__ == '__main__':
    app.run_server(debug=True)
