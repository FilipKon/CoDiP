from pathlib import Path

import dash_uploader as du
from dash import callback
from dash import dcc, html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

import global_variables
from apps import input_calculations

"""
    Different types of inputs
        - Load my mol2 file
        - Append my file to yours
    Create DF, DF_Names for each file uploaded:
        - Specify Names to remember
        - Save links
"""

UPLOAD_FOLDER = r"/tmp"

layout = html.Div([
    html.H4("Please upload your file for the Full/ECD/TMD analysis (mol2/csv/txt/xls)"),
    du.Upload(
        id='dash-uploader',
        text='Drag and Drop Here to upload!',
        text_completed='Uploaded: ',
        text_disabled='The uploader is disabled.',
        cancel_button=True,
        pause_button=False,
        disabled=False,
        filetypes=None,
        max_file_size=1024,
        chunk_size=1,
        default_style=None,
        upload_id=None,
        max_files=1,
    ),
    html.Div(id='callback-output'),
    html.Div([
        dcc.Input(
            id="ECD_nrs",
            type='number',
            min=0,
            max=2000,
            step=1,
            placeholder="ECD amino acids (first - till)",
            style={'width': '50%'}
        ),
        dcc.Input(
            id="TMD_nrs",
            type='number',
            min=0,
            max=2000,
            step=1,
            placeholder="TMD amino acids (ECD last+1 - end)",
            style={'width': '50%'}
        )
    ], style={'width': '100%', 'marginTop': 25, 'marginBottom': 25}),
    html.Div(id='range_sliders'),
    html.Button(children='Submit', id='submit-val', type='submit', n_clicks=0),
    html.Div(id='answer_1'),
    html.Div([html.H4("Please upload your file for the TMD binding site 3 analysis (mol2/csv/txt/xls)")],
             style={'marginBottom': 35, 'marginTop': 25}),
    du.Upload(
        id='dash-uploader_bs3',
        text='Drag and Drop Here to upload!',
        text_completed='Uploaded: ',
        text_disabled='The uploader is disabled.',
        cancel_button=True,
        pause_button=False,
        disabled=False,
        filetypes=None,
        max_file_size=1024,
        chunk_size=1,
        default_style=None,
        upload_id=None,
        max_files=1,
    ),
    html.Div(id='callback-output_3'),
    html.Div([html.H4("Please upload your file for the TMD binding site 4 analysis (mol2/csv/txt/xls)")],
             style={'marginBottom': 35, 'marginTop': 25}),
    du.Upload(
        id='dash-uploader_bs4',
        text='Drag and Drop Here to upload!',
        text_completed='Uploaded: ',
        text_disabled='The uploader is disabled.',
        cancel_button=True,
        pause_button=False,
        disabled=False,
        filetypes=None,
        max_file_size=1024,
        chunk_size=1,
        default_style=None,
        upload_id=None,
        max_files=1,
    ),
    html.Div(id='callback-output_4'),
    html.Div([html.H4("Please upload your file for the TMD binding site 5 analysis (mol2/csv/txt/xls)")],
             style={'marginBottom': 35, 'marginTop': 25}),
    du.Upload(
        id='dash-uploader_bs5',
        text='Drag and Drop Here to upload!',
        text_completed='Uploaded: ',
        text_disabled='The uploader is disabled.',
        cancel_button=True,
        pause_button=False,
        disabled=False,
        filetypes=None,
        max_file_size=1024,
        chunk_size=1,
        default_style=None,
        upload_id=None,
        max_files=1,
    ),
    html.Div(id='callback-output_5')
])


@callback(
    [Output('callback-output', 'children'),
     Output('global_link', 'data')],
    [Input('dash-uploader', 'isCompleted')],
    [State('dash-uploader', 'fileNames'),
     State('dash-uploader', 'upload_id')],
)
def display_files(isCompleted, fileNames, upload_id):
    if not isCompleted:
        return
    print(fileNames)
    print(isCompleted)
    print(upload_id)
    filepath = 'None'
    if fileNames is not None:
        out = []
        for filename in fileNames:
            file = Path(UPLOAD_FOLDER) / filename
            out.append(file)
        filepath = '/tmp/' + upload_id + '/' + filename
        # datas = open(filepath, 'r')
        # datas = datas.readlines()
        # print(datas)
        return html.Ul([html.Li(str(x)) for x in out]), filepath
    return html.Ul(html.Li("No Files Uploaded Yet!")), filepath


@callback(
    [Output('answer_1', 'children')],
    [Input('submit-val', 'n_clicks'),
     Input('global_link', 'data')],
    State('ECD_nrs', 'value'),
    State('TMD_nrs', 'value')
)
def update_output(n_clicks, global_link, ECD_value, TMD_value):
    if n_clicks is None:
        raise PreventUpdate
    if n_clicks:
        if global_link is None:
            return html.Div([html.H2('Please upload a file first.')])
        if 'mol' in global_link:
            data = open(global_link, 'r')
            data = data.readlines()
            # print(data)
            input_calculations.read_seq(data, ECD_value, TMD_value)
            if '/' in str(global_link):
                text = global_link.split('/')
                print(text)
                # text = text[len(text)]
                global_variables.names_df.append(text)
            else:
                global_variables.names_df.append(global_link)
            return html.Div([html.H2('Analysis done. Please visit the global section to view results.')])
    return html.Div([html.H2('...')])


"""
@callback(
    [Output('callback-output_3', 'children')],
    [Input('dash-uploader_bs3', 'isCompleted')],
    [State('dash-uploader_bs3', 'fileNames'),
     State('dash-uploader_bs3', 'upload_id')],
)
def display_files(isCompleted, fileNames, upload_id):
    if not isCompleted:
        return
    print(fileNames)
    print(isCompleted)
    print(upload_id)
    if fileNames is not None:
        out = []
        for filename in fileNames:
            file = Path(UPLOAD_FOLDER) / filename
            out.append(file)
        filepath = '/tmp/' + upload_id + '/' + filename
        # datas = open(filepath, 'r')
        # datas = datas.readlines()
        # print(datas)
        return html.Ul([html.Li(str(x)) for x in out]), filepath
    return html.Ul(html.Li("No Files Uploaded Yet!"))
@callback(
    [Output('callback-output_4', 'children')],
    [Input('dash-uploader_bs4', 'isCompleted')],
    [State('dash-uploader_bs4', 'fileNames'),
     State('dash-uploader_bs4', 'upload_id')],
)
def display_files(isCompleted, fileNames, upload_id):
    if not isCompleted:
        return
    print(fileNames)
    print(isCompleted)
    print(upload_id)
    if fileNames is not None:
        out = []
        for filename in fileNames:
            file = Path(UPLOAD_FOLDER) / filename
            out.append(file)
        filepath = '/tmp/' + upload_id + '/' + filename
        # datas = open(filepath, 'r')
        # datas = datas.readlines()
        # print(datas)
        return html.Ul([html.Li(str(x)) for x in out]), filepath
    return html.Ul(html.Li("No Files Uploaded Yet!"))
@callback(
    [Output('callback-output_5', 'children')],
    [Input('dash-uploader_bs5', 'isCompleted')],
    [State('dash-uploader_bs5', 'fileNames'),
     State('dash-uploader_bs5', 'upload_id')],
)
def display_files(isCompleted, fileNames, upload_id):
    if fileNames is None:
        raise PreventUpdate
    if not isCompleted:
        return
    if fileNames is not None:
        out = []
        for filename in fileNames:
            file = Path(UPLOAD_FOLDER) / filename
            out.append(file)
        filepath = '/tmp/' + upload_id + '/' + filename
        datas = open(filepath, 'r')
        datas = datas.readlines()
        #df, df_nr input_calculations(datas)
        return html.Ul([html.Li(str(x)) for x in out]), filepath
    return html.Ul(html.Li("No Files Uploaded Yet!"))
"""
