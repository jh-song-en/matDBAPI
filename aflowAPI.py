import requests
import pandas as pd
import json, sys, os
from urllib.request import urlopen
from itertools import combinations
import plotly.graph_objects as go
import math


def get_data_from_database(composition, database):
    if database == 'af':
        compo = ",".join(composition)
        SERVER = "http://aflow.org"
        API = "/API/aflux/v1.0/?"
        MATCHBOOK = f"species({compo}),nspecies({len(composition)}),enthalpy_formation_atom,composition,species,stoichiometry,Bravais_lattice_relax"
        DIRECTIVES = "$paging(0)"
        SUMMONS = MATCHBOOK + "," + DIRECTIVES
        response = json.loads(urlopen(SERVER + API + SUMMONS).read().decode("utf-8"))
    if database == 'mp':
        data = {
            'criteria': {
                'elements': {'$all': composition},
                'nelements': len(composition),
            },
            'properties': [
                'pretty_formula',
                'material_id',
                'formula',
                'formation_energy_per_atom',
            ]
        }
        r = requests.post('https://materialsproject.org/rest/v2/query',
                          headers={'X-API-KEY': 'P2tLmbTfRg9uMJDqzXV'},
                          data={k: json.dumps(v) for k, v in data.items()})
        response_content = r.json()
        response = response_content['response']

    return response

def get_every_combination(composition_list):
    total_combination = []

    total_number = len(composition_list)
    if total_number > 1:
        for i in range(1,len(composition_list)):
            total_combination += combinations(composition_list, i)

    total_combination.append(tuple(composition_list))
    return total_combination


def get_data_by_composition(composition_list, database):
    composition_combination = get_every_combination(composition_list)
    if database == "af":
        dictionary_list = []
        for composition in composition_combination:
            dictionary_list += get_data_from_database(composition, database)
        df1 = list_to_dataframe(dictionary_list)
        species = df1["species"].to_list()
        stoichiometry = df1["stoichiometry"].tolist()

        a = [dict(zip(species[i], stoichiometry[i])) for i in range(len(species))]
        print(a)
        df2 = list_to_dataframe(a, column=composition_list)
        for composition in composition_list:
            df2[composition] = df2[composition].fillna(0)
        merged_df = pd.concat([df1,df2],axis=1)
        return merged_df
    elif database == "mp":
        dictionary_list = []
        for composition in composition_combination:
            dictionary_list += get_data_from_database(composition, database)
        dataframe_keys = list(dictionary_list[0].keys()) + composition_list
        for dictionary in dictionary_list:
            keys = dictionary['formula'].keys()
            values = dictionary['formula'].values()
            atom_sum = sum(values)
            values = [value / atom_sum for value in values]
            new_formula = dict(zip(keys, values))
            dictionary.update(new_formula)

        df = list_to_dataframe(dictionary_list,column = dataframe_keys)
        for composition in composition_list:
            df[composition] = df[composition].fillna(0)
        return df

def list_to_dataframe(dictionary_list, column = None):
    if column == None:
        column = list(dictionary_list[0].keys())
    len_column = len(column)
    full_list = []
    for dictionary in dictionary_list:
        dictlist = len_column * [None]
        for key in list(dictionary.keys()):
            dictlist[column.index(key)] = dictionary[key]
        full_list.append(dictlist)
    df = pd.DataFrame(full_list, columns=column)
    return df


def show_phase_diagram(df, composition, formation_e_range, database):
    if database == 'af':
        customdata = ["compound"] + composition + ["enthalpy_formation_atom", 'auid']
        title = "phase diagram"+"("+",".join(composition)+")"+"-Aflow"
    if database == 'mp':
        customdata = ["pretty_formula"] + composition + ["formation_energy_per_atom", 'material_id']
        title = "phase diagram"+"("+",".join(composition)+")"+"-Materials Project"
        pass
    fig = go.Figure()
    if len(composition) == 2:


        fig.add_trace(go.Scatter(
                                x = df[composition[1]],
                                y = df[customdata[3]],
                                mode="markers",
                                marker = dict(
                                       size = 4,
                                       color = df[customdata[3]],
                                       colorscale = 'Viridis',
                                       cmin=formation_e_range[0],
                                       cmax=formation_e_range[1],
                                       showscale=True

                                   ),
                                customdata=df[customdata],
                                hovertemplate='<b>%{customdata[0]}</b>(%{customdata[4]})<br>'+
                                              composition[0] + ':%{customdata[1]}<br>'+
                                              composition[1] + ':%{customdata[2]}<br>'+
                                              customdata[3]+ ':%{customdata[3]}<br>'
                                )

        )

        fig.update_layout(title_text=title)
        fig.update_xaxes(title_text=customdata[2] + "(at%)",
                         range=[0,1])
        fig.update_yaxes(title_text=customdata[3] + "(eV)",
                         range=formation_e_range)
    if len(composition) == 3:
        sqrt3_over_2 = 0.5*math.sqrt(3)
        df['x'] = df[composition[1]] + 0.5 * df[composition[2]]
        df['y'] = sqrt3_over_2 * df[composition[2]]

        fig = go.Figure()

        fig.add_trace(go.Scatter3d(
            x=[0, 1, 0.5,0],
            y=[0, 0, 0.5*math.sqrt(3),0],
            z=[0,0,0,0],
            mode="lines+text",
            name="Lines, Markers and Text",
            text=composition,
            marker= dict(color = "black"),
            textposition="top center",
            textfont=dict(
                size=40
            ),
            hoverinfo='none'
        ))
        fig.add_trace(go.Scatter3d(x = df['x'], y = df['y'], z = df[customdata[4]],
                                   mode="markers",
                                   marker = dict(
                                       size = 4,
                                       color = df[customdata[4]],
                                       colorscale = 'Viridis',
                                       cmin=formation_e_range[0],
                                       cmax=formation_e_range[1],
                                       showscale=True
                                   ),
                                   customdata=df[customdata],
                                   hovertemplate='<b>%{customdata[0]}</b>(%{customdata[5]})<br>'+
                                                 composition[0] + ':%{customdata[1]}<br>'+
                                                 composition[1] + ':%{customdata[2]}<br>'+
                                                 composition[2] + ':%{customdata[3]}<br>'+
                                                 customdata[4] + ':%{customdata[4]}<br>'

                                   ))
        fig.update_scenes(xaxis_visible=False, yaxis_visible=False,zaxis_visible=True )
        fig.update_layout(title_text=title,
                          xaxis_showgrid=False,
                          yaxis_showgrid=False,
                          showlegend=False,
                          scene=dict(
                                      xaxis = dict(nticks=4, range=[0,1],),
                                      yaxis = dict(nticks=4, range=[0,1],),
                                      zaxis = dict(nticks=4, range=formation_e_range,title = customdata[4] + "(eV)"))
                          )
        df.pop("x")
        df.pop("y")
    fig.show()

if __name__ == "__main__":
    composition = ["Ni", "Fe"]
    df_mp = get_data_by_composition(composition, 'mp')
    print(df_mp)


