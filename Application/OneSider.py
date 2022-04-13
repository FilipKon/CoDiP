import base64
import dash_bootstrap_components as dbc
from dash import html, Dash
from apps import overview, bindingsites, contact, howto

image_filename = 'D:\\ConFam\\logo_codip_new.png'
#image_filename = '/home/fkoniuszewski/codip/logo_codip_new.png'
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

external_stylesheets2 = [
    {
        "href": "https://fonts.googleapis.com/css2?"
                "family=Lato:wght@400;700&display=swap",
        "rel": "stylesheet",
    },
]
external_stylesheets = [dbc.themes.LUX, external_stylesheets2]
app = Dash(__name__, external_stylesheets=[dbc.themes.LUX])
#           assets_external_path="D:\\ConFam\\assets\\bootstrap.css",
#           assets_folder='D:\\ConFam\\assets\\',
#           external_stylesheets=["D:\\ConFam\\assets\\bootstrap.css", dbc.themes.LUX])
# external_stylesheets=external_stylesheets
app.title = "CoDiP: Conformational Difference of Protein structures"


# nav_item = dbc.Nav(
nav_item = [
    dbc.NavItem(dbc.NavLink("Overview", active="exact", href="#overview", external_link=True)),
    dbc.NavItem(dbc.NavLink("How to", active="exact", href="#howto", external_link=True)),
    # dbc.NavItem(dbc.NavLink("Input", active="exact", href="#input", external_link=True)),
    # dbc.NavItem(dbc.NavLink("Global", active="exact", href="#globala", external_link=True)),
    dbc.NavItem(dbc.NavLink("Binding sites", active="exact", href="#bindingsites", external_link=True)),
    dbc.NavItem(dbc.NavLink("Contact", active="exact", href="#contact", external_link=True))
]
# fill=False,
# vertical=True,
# )
"""
navbar = dbc.NavbarSimple(
    children=[
        #nav_item,
        dbc.NavItem(dbc.NavLink("Overview", href="#")),
        dbc.NavItem(dbc.NavLink("Input", active="exact", href="/apps/input")),
        dbc.NavItem(dbc.NavLink("Global", active="exact", href="/apps/globala")),
        dbc.NavItem(dbc.NavLink("Binding sites", active="exact", href="/apps/bindingsites")),
        dbc.NavItem(dbc.NavLink("Contact", active="exact", href="/apps/contact"))
    ],
    brand="CoDiP",
    brand_href="#",
    color="primary",
    dark=True,
)
nav1 = dbc.Nav(nav_item, pills=True, fill=True)
"""
navbar = dbc.Navbar(
    dbc.Container(
        [
            #html.A(
            #   children=[
                dbc.Row(
                    [
                        dbc.Col(
                            children=[
                                html.Img(id='pici', src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                         style={'height': '13%', 'width': '13%'}),
                                html.H5("CoDiP")], width=5)
                    ],
                    justify="start",
                    align="center",
                    className="g-0",
                ),
                #dbc.Row(
                #    [
                #        dbc.Col(dbc.NavbarBrand("CoDiP", className="mr-auto"), width=5)
                #    ],
                #    justify="start",
                #    align="center",
                #    className="mr-auto",
                #),
            #]),
            dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
            dbc.Row(dbc.Nav(nav_item, pills=True, className='mr-auto'),
                    className="navigations",
                    justify="end",
                    align="center")
        ],
        fluid=False
    ),
    sticky="top",
    color="dark",
    dark=True,
    style={"position": "fixed", "z-index": 9999,
           'marginRight': '0px', "width": "100%"},
    className="navbar"
)

app.layout = html.Div([
    html.Div([navbar] + [html.Br()] * 6),
    html.Div([html.P("Overview", id="overview", className="fs-1 text fw-bolder"), overview.layout] + [html.Br()] * 10),
    html.Div([html.P("How to", id="howto", className="fs-1 text fw-bolder"), howto.layout] + [html.Br()] * 10),
    # html.Div([html.P("Input", id="input", className="fs-1 text fw-bolder"), input.layout] + [html.Br()]*10),
    html.Div([html.P("Binding sites", id="bindingsites", className="fs-1 text fw-bolder"), bindingsites.layout] + [
        html.Br()] * 10),
    html.Div([html.P("Contact", id="contact", className="fs-1 text fw-bolder"), contact.layout] + [html.Br()] * 10)
], className="dash-bootstrap")

if __name__ == '__main__':
    app.run_server(debug=True, port=8022, use_reloader=False, host='127.0.0.1')
