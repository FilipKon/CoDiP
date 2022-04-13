from dash import html

layout = html.Div([
    html.Div([html.Br()]*2),
    html.Div([html.H3("Ernst Laboratoy - Medical University of Vienna")] + [html.Br()]*3),
    html.Div([html.H4("Filip Koniuszewski: filip.koniuszewski@meduniwien.ac.at")] + [html.Br()]),
    html.Div([html.H5("Do you have any questions or suggestions? Just write me and I'll be happy to help you")])
    #html.Div([html.H4("Margot Ernst: margot.ernst@meduniwien.ac.at")])
])
