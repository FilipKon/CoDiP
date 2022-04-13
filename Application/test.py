import base64

import dash_bootstrap_components as dbc
from dash import html, dash

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

# the style arguments for the sidebar. We use position:fixed and a fixed width
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

sidebar = html.Div(
    [
        html.H2("Sidebar", className="display-4"),
        html.Hr(),
        html.P(
            "A simple sidebar layout with navigation links for same page", className="lead"
        ),
        dbc.Nav(
            [
                dbc.NavLink("Home", href="#start", external_link=True),
                dbc.NavLink("Mid", href="#mid",external_link=True),
                dbc.NavLink("End", href="#end",external_link=True),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)

image_filename = 'D:\\ConFam\\logo_codip_new.png'
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

content = html.Div([
            html.Div([html.Img(src=app.get_asset_url('D:\\ConFam\\logo_codip.png.png'))]),
            html.Div([html.P("Hello World", id="start"),
                      html.Div([
                          html.Img(id='Picture', src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                   style={'height': '5%', 'width': '5%', "float": "left", 'margin-top': 0})
                          # html.Img(src=app.get_asset_url('Figure7.png')
                      ], style={'display': 'inline-block'}),
                      ] +
                                  [html.Br()]*50 +
                                  [html.P("Middle World", id="mid")] +
                                  [html.Br()]*50 +
                                  [html.P("Goodbye world", id="end")]
                                 ),
    ], id="page-content", style=CONTENT_STYLE)


app.layout = html.Div([sidebar, content])



if __name__ == '__main__':
    app.run_server(debug=True, port=8022,use_reloader=False, host='127.0.0.1')