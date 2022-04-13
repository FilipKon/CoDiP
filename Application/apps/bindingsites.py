import sys
import dash_bio as dashbio
import pandas as pd
import plotly.graph_objs as go
import dash_bootstrap_components as dbc
# import Bio
from Bio.PDB.PDBParser import PDBParser
from dash import callback, Input, Output, dcc, html, dash_table
from dash_bio_utils import PdbParser, create_mol3d_style
import base64
import global_variables
from PIL import Image

image_filename = 'D:\\ConFam\\Figure_legend.png'
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

if not sys.warnoptions:
    import warnings

    warnings.simplefilter("ignore")

ids = ["Max", "Min"]

tabs_styles = {
    'height': '30px'
}

tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold',
    'height': '30px'
}

tab_style_2 = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold',
    'height': '20px'
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#119DFF',
    'color': 'white',
    'padding': '6px'
}

layout = html.Div([
    html.Div([
        html.Div([
            html.Div([
                html.H5('Protein families to be display:'),
                dcc.Checklist(
                    options=[{'label': 'GABAAR', 'value': 'GABAAR'},
                             {'label': 'GlyR', 'value': 'GlyR'},
                             {'label': 'nAChR', 'value': 'nAChR'},
                             {'label': '5HT3R', 'value': '5HT3R'}],
                    value=['GABAAR', 'GlyR', 'nAChR', '5HT3R'], id='Families'),  # , inline=True
            ], style={'display': 'inline-block', 'width': '100%', 'margin-bottom': 50}),
            dcc.Tabs(id="tabs-BS", value='tab-1-BS', children=[
                dcc.Tab(label='site 3', value='site3', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='site 4', value='site4', style=tab_style, selected_style=tab_selected_style),
                dcc.Tab(label='site 5', value='site5', style=tab_style, selected_style=tab_selected_style),
            ], colors={
                "border": "black",
                "primary": "gold",
                "background": "#6e7d92"}, style=tabs_styles)], style={'width': '100%',
                                                                      'display': 'inline-block',
                                                                      'margin-top': 70})]),
    html.Div(id='tabs-content-BS'),
    html.Div(id='tbl_out'),
    html.Div(id='hover-BS-data'),
])


@callback(Output('tabs-content-BS', 'children'),
          Input('tabs-BS', 'value'),
          Input('Families', 'value'))
def render_content(tabBS, famvalues):
    filters = ['Min', 'Max']
    for item in famvalues:
        #print(['--', item])
        if item == 'GABAAR':
            famvalues.append('GABAAR_chimeras')
            filters += global_variables.GABAAR_ids
        if item == 'GlyR':
            famvalues.append('GlyR_chimeras')
            filters += global_variables.GlyR_ids
        if item == 'nAChR':
            filters += global_variables.nAChRs_ids
        if item == '5HT3R':
            filters += global_variables.HT3R_ids
    #print(filters)
    if tabBS == 'site3':
        df = pd.read_csv('data/BS3/BS3_PCA_data_new.csv', sep=';')
        # df = pd.read_csv('/home/fkoniuszewski/data/BS3/BS3_PCA_data_new.csv', sep=';')
        df_dists = pd.read_csv('data/BS3/Apos_full_data_bs3.csv', sep=',')
        # df_dists = pd.read_csv('/home/fkoniuszewski/data/BS3/Apos_full_data_bs3.csv', sep=';')
        df_dists_names = pd.read_csv('data/BS3/Apos_full_labels_bs3.csv', sep=',')
        # df_dists_names = pd.read_csv('/home/fkoniuszewski/data/BS3/Apos_full_labels_bs3.csv', sep=';')
        if len(famvalues) < 6:
            df = df[df['Family'].isin(famvalues)]
            df_dists = df_dists[df_dists['pdbid'].str.contains('|'.join(filters))]
            df_dists_names = df_dists_names[df_dists_names['pdbid'].str.contains('|'.join(filters))]
        df_dists = df_dists.drop(df_dists.columns[0], axis=1)
        df_dists_names = df_dists_names.drop(df_dists_names.columns[0], axis=1)
        Scene = dict(xaxis=dict(title='PC1', showgrid=True, gridwidth=1, gridcolor='black'),
                     yaxis=dict(title='PC2', showgrid=True, gridwidth=1, gridcolor='black'),
                     zaxis=dict(title='PC3', showgrid=True, gridwidth=1, gridcolor='black'))
        trace = go.Scatter3d(x=df['PC1'], y=df['PC2'], z=df['PC3'], mode='markers+text', hovertext=df['labels'],
                             marker_symbol=df['markers'], showlegend=False,
                             marker=dict(color=df['Colors'], size=10, line=dict(color='black', width=10)))
        data = [trace]
        layout = go.Layout(margin=dict(l=0, r=0), scene=Scene, height=800, width=800, showlegend=True)
        fig = go.Figure(data=data, layout=layout)
        return html.Div([
            html.Div([html.H3('Binding site 3')],
                     style={'textAlign': 'center', 'marginTop': 30}),
            html.Div([
                dbc.Row([
                    dbc.Col(dbc.CardImg(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                        style={"height": "80%", "width": "80%"}), width={"size": 2, "offset": 1},
                            lg=2, md=2, align="end"),
                    dbc.Col(dbc.Card(html.Div("After you click on a data point in the 3d scatter plot you will be "
                                              "able to also select residues in the molecular 3d viewer. Here the "
                                              "selected residue and chain will be displayed", id="SelectedMol3d"),
                                     style={'marginTop': '10px'}),
                            width={"size": 2}, lg=3, md=4, align="end"),
                ], justify="between"),
                dbc.Row([
                    dbc.Col(dcc.Graph(id='site3-scatter3d', figure=fig, hoverData={'points': [{'customdata': '6HIO'}]},
                                      style={'width': '49%', "margin": 0, 'display': 'inline-block'}),
                            width={"size": 3, "offset": 1}),
                    dbc.Col(html.Div(id="mol3d",
                                     style={'width': '49%', "margin": 0, 'display': 'inline-block'}),
                            width={"size": 3, "offset": 2}, align="start"),
                ])], style={'display': 'inline-block'}),
            dcc.Store(id='intermediate-table', data=df_dists.to_json(date_format='iso', orient='split')),
            dcc.Store(id='intermediate-table_nr', data=df_dists_names.to_json(date_format='iso', orient='split'))
        ])
    elif tabBS == 'site4':
        df = pd.read_csv('data/BS4/BS4_PCA_data_new.csv', sep=';')
        # df = pd.read_csv('/home/fkoniuszewski/data/BS4/BS4_PCA_data_new.csv', sep=';')
        df_dists = pd.read_csv('data/BS4/Apos_full_data_bs4.csv', sep=',')
        # df_dists = pd.read_csv('/home/fkoniuszewski/data/BS4/Apos_full_data_bs4.csv', sep=';')
        df_dists_names = pd.read_csv('data/BS4/Apos_full_labels_bs4.csv', sep=',')
        # df_dists_names = pd.read_csv('/home/fkoniuszewski/data/BS4/Apos_full_labels_bs4.csv', sep=';')
        if len(famvalues) < 6:
            df = df[df['Family'].isin(famvalues)]
            df_dists = df_dists[df_dists['pdbid'].str.contains('|'.join(filters))]
            df_dists_names = df_dists_names[df_dists_names['pdbid'].str.contains('|'.join(filters))]
        df_dists = df_dists.drop(df_dists.columns[0], axis=1)
        df_dists_names = df_dists_names.drop(df_dists_names.columns[0], axis=1)
        Scene = dict(xaxis=dict(title='PC1', showgrid=True, gridwidth=1, gridcolor='black'),
                     yaxis=dict(title='PC2', showgrid=True, gridwidth=1, gridcolor='black'),
                     zaxis=dict(title='PC3', showgrid=True, gridwidth=1, gridcolor='black'))
        trace = go.Scatter3d(x=df['PC1'], y=df['PC2'], z=df['PC3'], mode='markers+text', hovertext=df['labels'],
                             marker_symbol=df['markers'],  # symbol=df2['markers'],
                             marker=dict(color=df['Colors'], size=10, line=dict(color='black', width=10)))
        data = [trace]
        layout = go.Layout(margin=dict(l=0, r=0), scene=Scene, height=800, width=800)
        fig = go.Figure(data=data, layout=layout)
        return html.Div([
            html.Div([html.H3('Binding site 3')],
                     style={'textAlign': 'center', 'marginTop': 30}),
            html.Div([
                dbc.Row([
                    dbc.Col(dbc.CardImg(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                        style={"height": "80%", "width": "80%"}), width={"size": 2, "offset": 1},
                            lg=2, md=2, align="end"),
                    dbc.Col(dbc.Card(html.Div("After you click on a data point in the 3d scatter plot you will be "
                                              "able to also select residues in the molecular 3d viewer. Here the "
                                              "selected residue and chain will be displayed", id="SelectedMol3d"),
                                     style={'marginTop': '10px'}),
                            width={"size": 2}, lg=3, md=4, align="end"),
                ], justify="between"),
                dbc.Row([
                    dbc.Col(dcc.Graph(id='site3-scatter3d', figure=fig, hoverData={'points': [{'customdata': '6HIO'}]},
                                      style={'width': '49%', "margin": 0, 'display': 'inline-block'}),
                            width={"size": 3, "offset": 1}),
                    dbc.Col(html.Div(id="mol3d",
                                     style={'width': '49%', "margin": 0, 'display': 'inline-block'}),
                            width={"size": 3, "offset": 2}, align="start"),
                ])], style={'display': 'inline-block'}),
            dcc.Store(id='intermediate-table', data=df_dists.to_json(date_format='iso', orient='split')),
            dcc.Store(id='intermediate-table_nr', data=df_dists_names.to_json(date_format='iso', orient='split'))
        ])
    elif tabBS == 'site5':
        df = pd.read_csv('data/BS5/BS5_PCA_data_new.csv', sep=';')
        # df = pd.read_csv('/home/fkoniuszewski/data/BS5/BS5_PCA_data_new.csv', sep=';')
        df_dists = pd.read_csv('data/BS5/Apos_full_data_bs5_new.csv', sep=',')
        # df_dists = pd.read_csv('/home/fkoniuszewski/data/BS5/Apos_full_data_bs5_new.csv', sep=';')
        df_dists_names = pd.read_csv('data/BS5/Apos_full_labels_bs5_new.csv', sep=',')
        # df_dists_names = pd.read_csv('/home/fkoniuszewski/data/BS5/Apos_full_labels_bs5_new.csv', sep=';')
        if len(famvalues) < 6:
            df = df[df['Family'].isin(famvalues)]
            df_dists = df_dists[df_dists['pdbid'].str.contains('|'.join(filters))]
            df_dists_names = df_dists_names[df_dists_names['pdbid'].str.contains('|'.join(filters))]
        df_dists = df_dists.drop(df_dists.columns[0], axis=1)
        df_dists_names = df_dists_names.drop(df_dists_names.columns[0], axis=1)
        Scene = dict(xaxis=dict(title='PC1', showgrid=True, gridwidth=1, gridcolor='black'),
                     yaxis=dict(title='PC2', showgrid=True, gridwidth=1, gridcolor='black'),
                     zaxis=dict(title='PC3', showgrid=True, gridwidth=1, gridcolor='black'))
        trace = go.Scatter3d(x=df['PC1'], y=df['PC2'], z=df['PC3'], mode='markers+text', hovertext=df['labels'],
                             marker_symbol=df['markers'],  # symbol=df2['markers'],
                             marker=dict(color=df['Colors'], size=10, line=dict(color='black', width=10)))
        data = [trace]
        layout = go.Layout(margin=dict(l=0, r=0), scene=Scene, height=800, width=1400)
        fig = go.Figure(data=data, layout=layout)
        return html.Div([
            html.Div([html.H3('Binding site 5')],
                     style={'textAlign': 'center', 'marginTop': 30}),
            html.Div([
                dbc.Row([
                    dbc.Col(dbc.CardImg(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                                        style={"height": "80%", "width": "80%"}), width={"size": 2, "offset": 1},
                            lg=2, md=2, align="end"),
                    dbc.Col(dbc.Card(html.Div("After you click on a data point in the 3d scatter plot you will be "
                                              "able to also select residues in the molecular 3d viewer. Here the "
                                              "selected residue and chain will be displayed", id="SelectedMol3d"),
                                     style={'marginTop': '10px'}),
                            width={"size": 2}, lg=3, md=4, align="end"),
                ], justify="between"),
                dbc.Row([
                    dbc.Col(dcc.Graph(id='site3-scatter3d', figure=fig, hoverData={'points': [{'customdata': '6HIO'}]},
                                      style={'width': '49%', "margin": 0, 'display': 'inline-block'}),
                            width={"size": 3, "offset": 1}),
                    dbc.Col(html.Div(id="mol3d",
                                     style={'width': '49%', "margin": 0, 'display': 'inline-block'}),
                            width={"size": 3, "offset": 2}, align="start"),
                ])], style={'display': 'inline-block'}),
            dcc.Store(id='intermediate-table', data=df_dists.to_json(date_format='iso', orient='split')),
            dcc.Store(id='intermediate-table_nr', data=df_dists_names.to_json(date_format='iso', orient='split'))
        ])


# before hoverData - clickData
@callback(
    Output('hover-BS-data', 'children'),
    Input('site3-scatter3d', 'clickData'),
    Input('intermediate-table', 'data'))
def update_table(hoverData, jsonified_cleaned_data):
    id = hoverData['points'][0]['hovertext']
    if id not in ids:
        ids.append(id)
    df = pd.read_json(jsonified_cleaned_data, orient='split')
    # dff = df[df['pdbid'].isin(ids)].set_index('pdbid').T
    dff = df[df['pdbid'].isin(ids)]
    dff = dff.set_index('pdbid').T
    return html.Div([
        dash_table.DataTable(
            id='hovertbl',
            columns=[{"name": i, "id": i} for i in dff.columns.values],
            data=dff.to_dict('records'),
            style_header={
                'backgroundColor': 'rgb(50, 50, 50)',
                'color': 'white',
                'align': 'center',
                'fontWeight': 'bold',
                'textAlign': 'center'
            },
            style_data={
                'backgroundColor': 'rgb(100, 100, 100)',
                'color': 'white'
            },
            page_action="native",
            style_cell={'textAlign': 'center'},
        )
    ])


@callback(
    Output('tbl_out', 'children'),
    Input('site3-scatter3d', 'clickData'),
    Input('intermediate-table_nr', 'data'),
    Input('hovertbl', 'active_cell'))
def update_graphs(hoverData, jsonified_cleaned_data_nrs, active_cell):
    n = 0
    if active_cell == None:
        active_cell = {'row': 0, 'column': 0}
    id = hoverData['points'][0]['hovertext']
    if id not in ids:
        ids.append(id)
    df_nrs = pd.read_json(jsonified_cleaned_data_nrs, orient='split')
    dff = df_nrs[df_nrs.pdbid.isin(ids)].set_index('pdbid').T
    row = active_cell['row']
    col = active_cell['column']
    text = str(dff.iat[row, col])
    text = text.replace("het=", "")
    text = text.replace("resseq=", "")
    text = text.replace("icode= ", "")
    text = text.replace("id=", "")
    n = text.count('-')
    text = text.split('-')
    if n < 2:
        txt1 = 'Residue 1: ' + text[0]
        txt2 = 'Residue 2: ' + text[1]
        txt3 = ''
    else:
        txt1 = 'Residue 1: ' + text[1]
        txt2 = 'Residue 2: ' + text[2]
        txt3 = text[3] + text[4]

    return html.Div([
        html.H4('Selected distance was calculated by choosing'),
        html.H4(txt1),
        html.H4(txt2),
        html.H4(txt3)
    ])


@callback(
    Output("mol3d", "children"),
    Input('site3-scatter3d', 'clickData'),
    Input('intermediate-table_nr', 'data'),
    Input('hovertbl', 'active_cell')
)
def update_3dviewer(hoverData, jsonified_cleaned_data_nrs, active_cell):
    id = hoverData['points'][0]['hovertext']
    chain = 'No'
    if active_cell == None:
        active_cell = {'row': 0, 'column': 0}
    if '[' in id:
        id = id.split('[')
        id = id[0]
    elif '_' in id:
        id = id.split('_')
        chain = id[2]
        id = id[0]

    df_nrs = pd.read_json(jsonified_cleaned_data_nrs, orient='split')
    dff = df_nrs[df_nrs.pdbid.isin(ids)].set_index('pdbid').T
    row = active_cell['row']
    col = active_cell['column']
    text = str(dff.iat[row, col])
    #print(text)
    text = text.replace("het=", "")
    text = text.replace("resseq=", "")
    text = text.replace("icode= ", "")
    text = text.replace("id=", "")
    n = text.count('-')
    text = text.split('-')
    if n >= 3:
        pdbid = text[3]
        pdbid = pdbid[19:-2]
        print(["PDBID ----", pdbid])
        residue1 = text[1]
        residue2 = text[2]
        chain1 = residue1.split('>')
        chain1 = chain1[0].split(' ')
        chain1 = chain1[2]
    else:
        print(id)
        pdbid = id
        residue1 = text[0]
        residue2 = text[1]
        chain1 = residue1.split('>')
        chain1 = chain1[0].split(' ')
        chain1 = chain1[1]
    chain2 = residue2.split('>')
    chain2 = chain2[0].split(' ')
    chain2 = chain2[2]
    print(residue1)
    residue1 = [int(s) for s in residue1.split() if s.isdigit()]
    print(residue1)
    residue1 = residue1[0]
    print(residue1)
    residue2 = [int(s) for s in residue2.split() if s.isdigit()]
    residue2 = residue2[0]
    print(["residues: ", residue1, residue2])
    print(["chains: ", chain1, chain2])
    print(pdbid)
    link_strcts = 'C:\\Users\\filip\\PycharmProjects\\plgic\\Structures\\' + pdbid + '.pdb'
    # link_strcts = '/home/fkoniuszewski/codip/data/Structures/' + pdbid + '.pdb'
    parser = PdbParser(link_strcts)
    parser2 = PDBParser()
    structure = parser2.get_structure(id, link_strcts)
    first_model = structure[0]
    chain1 = first_model[chain1]
    chain2 = first_model[chain2]
    residue1 = chain1[int(residue1)]
    residue2 = chain2[int(residue2)]
    print(["Residues", residue1, residue2])
    ca1 = residue1["CA"]
    ca2 = residue2["CA"]
    coords1 = ca1.get_coord().tolist()
    coords2 = ca2.get_coord().tolist()
    #print(["CAs", coords1, coords2])
    data = parser.mol3d_data()
    styles = create_mol3d_style(
        data['atoms'], visualization_type='cartoon', color_element='chain'
    )
    return html.Div([
        html.Div([pdbid], id="pdbid", style={"color": "white"}),
        dashbio.Molecule3dViewer(
            id='dashbio-default-molecule3d',
            modelData=data,
            styles=styles,
            #selectionType='atoms',
            shapes=[
                {
                    'type': 'Cylinder',
                    'start': {'x': coords1[0], 'y': coords1[1], 'z': coords1[2]},
                    'end': {'x': coords2[0], 'y': coords2[1], 'z': coords2[2]},
                    'radius': 1.0,
                    'fromCap': 1,
                    'toCap': 2,
                    'color': 'green',
                    'opacity': 1
                }
            ]
        )
    ], style={'width': '49%', "margin": 0, 'display': 'inline-block', 'marginLeft': '250px', 'marginTop': '0px'})


@callback(
    Output('SelectedMol3d', 'children'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds'),
    Input('pdbid', 'children')
)
def show_selected_atoms(atom_ids, pdbid):
    pdbid = pdbid[0]
    link_strcts = 'C:\\Users\\filip\\PycharmProjects\\plgic\\Structures\\' + pdbid + '.pdb'
    parser = PdbParser(link_strcts)
    data = parser.mol3d_data()
    if atom_ids is None or len(atom_ids) == 0:
        return 'No residue has been selected. Click somewhere on the molecular \
        structure to select an atom.'
    text = 'PDB ID: ' + pdbid
    return [html.Div([
        html.Div(text),
        html.Div('Chain: {}'.format(data['atoms'][atm]['chain'])),
        html.Div('Residue name: {}'.format(data['atoms'][atm]['residue_name'])),
    ]) for atm in atom_ids]
