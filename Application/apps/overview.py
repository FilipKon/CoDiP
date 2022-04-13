import base64

from dash import html

image_filename2 = 'D:\\ConFam\\Figure7.png'
#image_filename = '/home/fkoniuszewski/codip/Figure7.png'
encoded_image2 = base64.b64encode(open(image_filename2, 'rb').read())

layout = html.Div([
    html.H5(children="Conformational Differences of Protein structures for TMD binding sites of pLGICs",
            style={"textAlign": "center"}),
    html.P("Methods see Koniuszewski et al., 2022", className="text-center"),
    html.P("Overview of the three different binding sites used in the analysis", className="text-center"),
    #html.Div([html.Img(id='figure7', src='data:image/png;base64,{}'.format(encoded_image.decode()),
    #         style={'height': '20%', 'width': '20%'})], style={'textAlign': 'center'}),
    html.Div([
        html.A([
            html.Img(src='data:image/png;base64,{}'.format(encoded_image2.decode()),
                     style={'height': '20%', 'width': '20%'})], href='#bindingsites'),
    ], style={'textAlign': 'center'}),
    html.Div([html.Br()]*2),
    html.P("Please Cite: Koniuszewski F., Vogel D. F., Bampali K., Fabjan J., Seidl T., Scholze P., Schmiedhofer P., Langer \
T., Ernst M.; „Molecular mingling: Multimodal predictions of ligand promiscuity in pentameric \
ligand-gated ion channels“; Front. Mol. Biosci. doi: 10.3389/fmolb.2022.860246", className="text-center")
])