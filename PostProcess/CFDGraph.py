'''
Crea area graph per distribuzione del CFD. Per ogni valore di CFD (0-100) sull'asse X, indica sull'asse y il numero di target con quel
valore di CFD, distinti in REF e VAR. In input prende il file .CDFGraph.txt
TODO Modificare lo script in modo da aggiungere un dash Input, in modo che l'utente possa scrivere un sample da visualizzare, e un Button
per confermare la selezione di quel sample. Quando l'utente seleziona un sample, lo script deve caricare il file CFDGraph, come un dataframe,
e tenere solo la colonna del valore REF e quella del sample selezionato. Se non è stato selezionato alcun sample, visualizzare colonna REF e VAR.
Se il sample non esiste, far apparire un messaggio di errore (semplicemente un div il cui children è l'output della callback del filtering, se non
c'è il sample far apparire 'Il sample non esiste', altrimenti ritorna '').
TODO Potrebbe essere interessante fare in modo che se l'utente clicca sul grafico, prendere il valore di X (ovvero il CFD) e poi fare un grep sul
file dei risultati e mostrare i target con quel valore di CFD. La funzione display_selected_data tenta di fare la prima parte, ma non funziona con
questa tipologia di grafico, vedere se si può risolvere in qualche modo
NOTE    https://plotly.com/python/filled-area-plots/ 
        https://dash.plotly.com/dash-core-components/graph
        https://dash.plotly.com/interactive-graphing
NOTE 127.0.0.1:8050 per aprire la pagina (o 0.0.0.0:8050)
'''


import plotly.graph_objects as go # or plotly.express as px

# fig.add_trace( ... )
# fig.update_layout( ... )

import pandas as pd
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import sys

def createGraph(cfd_distribution, showLog):
    fig = go.Figure() # or any Plotly Express function e.g. px.bar(...)
    fig.add_trace(go.Scatter(x=list(range(101)), y=list(cfd_distribution['ref']), fill='tozeroy', name = 'Tagets in Reference',#fillcolor = 'yellow',
                    #mode='none' # override default markers+lines
                        ))
    fig.add_trace(go.Scatter(x=list(range(101)), y=list(cfd_distribution['var']), fill='tozeroy',name = 'Targets in Enriched',
                        #mode= 'none'
                        ))
    if showLog:
        fig.update_layout(yaxis_type="log", xaxis_title="CFD Value",
        yaxis_title="Number of Targets (log scale)", hovermode = 'x')
    else:
        fig.update_layout( xaxis_title="CFD Value",
        yaxis_title="Number of Targets",hovermode = 'x')
    return fig

import dash
import dash_core_components as dcc
import dash_html_components as html

def CFDGraph(CFD_total_path):
    if CFD_total_path is None:
        return ''
    final_list = []
    final_list.append(html.P('Number of targets found in the Reference and Enriched Genome for each CFD Score value (0-100)'))
    # final_list.append(
    #     dcc.Checklist(
    #         options = [
    #             {'label':'Log value for Number of Targets','value':'logval', 'disabled':False}
    #         ], 
    #         value = ['logval'],
    #         id = 'checklist-advanced',
            
    #     ),
    # )
    cfd_distribution = pd.read_csv(CFD_total_path, sep = '\t')
    final_list.append(
        html.Div(
            dcc.Graph(figure=createGraph(cfd_distribution,True), id = 'CFD-graph-id'),
            id = 'div-CFD-graph'
        )
    )
    final_list.append(html.Div(id = 'selected-data'))


    final_list.append(html.Br())
    # import plotly.express as px
    # df = px.data.gapminder()
    # fig = px.area(df, x="year", y="pop", color="continent",
    # 	      line_group="country")

    # final_list.append(dcc.Graph(figure = fig, id = 'graph-id2'))
    return final_list

# @app.callback(
#     Output('div-graph','children'),
#     [Input('checklist-advanced', 'value')]
# )
# def changeLogGraph(val):
#     if val is None:
#         raise PreventUpdate
#     cfd_distribution = pd.read_csv('/home/francesco/Progetto_Giugno/Progetto/script/total.CFDGraph.txt'
# , sep = '\t')
#     if 'logval' in val:
#         return dcc.Graph(figure=createGraph(cfd_distribution,True))
#     return dcc.Graph(figure=createGraph(cfd_distribution,False))

# @app.callback(
#     Output('selected-data', 'children'),
#     [Input('CFD-graph-id', 'clickData')])
# def display_selected_data(selectedData):
#     print(selectedData) #NOTE non funziona
#     return ''

if __name__ == '__main__':
    CFDGraph(None)