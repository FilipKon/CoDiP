from dash import html
import base64

image_filename = 'D:\\ConFam\\Explanation_3dmol.png'
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

layout = html.Div([
    html.Div([html.Br()]*2),
    html.Div([html.H2("Explanation of the binding site analysis")] + [html.Br()]*2),
    html.Div([
        html.H5("1) Choose a binding site of which you want to see the analysis of"),
        html.H5("2) If you want to display only some protein families than deselected the protein families which "
                "shouldn't be displayed above the binding site selection"),
        html.H5("3) Click on a data point in the 3D scatter plot to see the measured distances for this data point"),
        html.H5("4) If you want to compare multiple data points click on each one and the distances will be displayed "
                "in the table above with additional the min/max distances for all structures"),
        html.H5("5) Click on a distance in the table to see where the amino acids are in a 3d molecular structures"),
        html.H5("6) The green cylinder showing the distance which was clicked in the table"),
        html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                     style={'height': '30%', 'width': '30%', 'margin': '0px', 'marginBottom': '0px'}),
        html.H5("7) To find out which residues and chains are on a specific position just click on the structure")],
    style={'marginLeft': '20px'}),
    html.Div([html.Br()])
])
