import requests
import pandas as pd
import numpy as np
import json, sys, os
from urllib.request import urlopen
from itertools import combinations
import plotly.graph_objects as go
import math
from scipy.spatial import ConvexHull


class APIRester():
    def __init__(self, composition, database, get_hull = False):
        self.composition = composition
        self.database = database
        self.df = self.get_data_by_composition()
        self.get_hull = get_hull
        if self.database == 'af':
            self.formation_column_name = 'enthalpy_formation_atom'
            self.id_name = 'auid'
        elif self.database == 'mp':
            self.formation_column_name = "formation_energy_per_atom"
            self.id_name = 'material_id'
        elif self.database == 'oq':
            self.formation_column_name = "delta_e"
            self.id_name = "entry_id"
        else:
            pass
        if database != "oq":
            self.np_hull, self.point_list = self.get_convex_elements()

    def get_convex_elements(self):

        columns = [self.id_name] + self.composition + [self.formation_column_name]
        return_df = self.df[columns]
        for compo in self.composition:
            return_df = return_df[return_df[compo] != 1]
        return_df = return_df[return_df[self.formation_column_name] < 0]
        if len(self.composition)==2:
            pure = [['-', 1, 0, 0, '-'], ['-', 0, 1, 0, '-']]
            return_df['original_index'] = return_df.index
            for row in pure:
                new_row = dict(zip(return_df.columns, row))

                return_df = return_df.append(new_row, ignore_index=True)
            return_df.reset_index(drop=True, inplace=True)
            return_df.dropna()
            numum = np.array(return_df[[self.composition[1], self.formation_column_name]])
            hull = ConvexHull(
                numum)  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
            convex_list = [[0, 0]]
            for index in hull.vertices:
                if numum[index][1] < 0:
                    convex_list.append(list(numum[index]))
            convex_list.append([1, 0])
            convex_array = np.array(convex_list)

            return convex_array, None
        if len(self.composition) == 3:
            pure = [['-', 1, 0, 0, 0, '-'], ['-', 0, 1, 0, 0, '-'], ['-', 0, 0, 1, 0, '-']]
            return_df['original_index'] = return_df.index
            for row in pure:
                new_row = dict(zip(return_df.columns, row))

                return_df = return_df.append(new_row, ignore_index=True)

            sqrt3_over_2 = 0.5 * math.sqrt(3)
            return_df['x'] = return_df[self.composition[1]] + 0.5 * return_df[self.composition[2]]
            return_df['y'] = sqrt3_over_2 * return_df[self.composition[2]]

            return_df.reset_index(drop=True, inplace=True)
            return_df.dropna()
            np_convex = np.array(return_df[['x', 'y', self.formation_column_name]])
            hull = ConvexHull(np_convex)
            vertices = list(hull.vertices)
            wall_lists_in_convex_hull = []
            for compo in self.composition:
                wall_list = list(return_df[return_df[compo] == 0].index)
                intersection = list(set(vertices).intersection(wall_list))
                wall_lists_in_convex_hull.append(intersection)
            wall_lists_in_convex_hull.append(vertices[-3:])
            true_simplices = []
            for simplix in list(hull.simplices):
                simp = list(simplix)

                put = True
                for compo_2_list in wall_lists_in_convex_hull:
                    intersection = list(set((simp)).intersection(compo_2_list))

                    if len(intersection) >= 3:
                        put = False
                if put:
                    true_simplices.append(list(simp))
            final_convex_3d_points = np.array(true_simplices)
            return np_convex, final_convex_3d_points



    def get_data_from_database(self, selected_composition):
        if self.database == 'af':
            compo = ",".join(selected_composition)
            SERVER = "http://aflow.org"
            API = "/API/aflux/v1.0/?"
            MATCHBOOK = f"species({compo}),nspecies({len(selected_composition)}),enthalpy_formation_atom,composition,species,stoichiometry,Bravais_lattice_relax"
            DIRECTIVES = "$paging(0)"
            SUMMONS = MATCHBOOK + "," + DIRECTIVES
            response = json.loads(urlopen(SERVER + API + SUMMONS).read().decode("utf-8"))
        elif self.database == 'mp':
            data = {
                'criteria': {
                    'elements': {'$all': selected_composition},
                    'nelements': len(selected_composition),
                },
                'properties': [
                    'pretty_formula',
                    'material_id',
                    'formula',
                    'formation_energy_per_atom',
                    'e_above_hull',
                ]
            }
            r = requests.post('https://materialsproject.org/rest/v2/query',
                              headers={'X-API-KEY': 'P2tLmbTfRg9uMJDqzXV'},
                              data={k: json.dumps(v) for k, v in data.items()})
            response_content = r.json()
            response = response_content['response']
        elif self.database == 'oq':
            compo = ",".join(selected_composition)
            column_list = ["name", "entry_id", "spacegroup", "ntypes", "band_gap", "delta_e", "stability"]
            query = "http://oqmd.org/oqmdapi/formationenergy?fields=" + ','.join(
                column_list) + "&filter=element_set=("+compo+") AND ntypes=" + str(len(selected_composition))
            reg = requests.get(query).json()
            response = reg["data"]
        return response

    def get_every_combination(self, composition_list):
        total_combination = []

        total_number = len(composition_list)
        if total_number > 1:
            for i in range(1,len(composition_list)):
                total_combination += combinations(composition_list, i)

        total_combination.append(tuple(composition_list))
        return total_combination


    def get_data_by_composition(self):
        composition_combination = self.get_every_combination(self.composition)
        if self.database == "af":
            dictionary_list = []
            for selected_composition in composition_combination:
                dictionary_list += self.get_data_from_database(selected_composition)
            df1 = self.list_to_dataframe(dictionary_list)
            species = df1["species"].to_list()
            stoichiometry = df1["stoichiometry"].tolist()

            a = [dict(zip(species[i], stoichiometry[i])) for i in range(len(species))]
            df2 = self.list_to_dataframe(a, column=self.composition)
            for atom in self.composition:
                df2[atom] = df2[atom].fillna(0)
            merged_df = pd.concat([df1,df2],axis=1)
            return merged_df
        elif self.database == "mp":
            dictionary_list = []
            for selected_composition in composition_combination:
                dictionary_list += self.get_data_from_database(selected_composition)
            dataframe_keys = list(dictionary_list[0].keys()) + self.composition
            for dictionary in dictionary_list:
                keys = dictionary['formula'].keys()
                values = dictionary['formula'].values()
                atom_sum = sum(values)
                values = [value / atom_sum for value in values]
                new_formula = dict(zip(keys, values))
                dictionary.update(new_formula)

            df = self.list_to_dataframe(dictionary_list,column = dataframe_keys)
            for atom in self.composition:
                df[atom] = df[atom].fillna(0)
            return df
        elif self.database == "oq":
            dictionary_list = []
            for selected_composition in composition_combination:
                dictionary_list += self.get_data_from_database(selected_composition)
            dataframe_keys = list(dictionary_list[0].keys())
            df = self.list_to_dataframe(dictionary_list, column=dataframe_keys)
            return df

    def list_to_dataframe(self, dictionary_list, column = None):
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


    def show_phase_diagram(self, formation_e_range, show_hull = True, show_figure = True):
        composition = self.composition
        if self.database == 'af':
            customdata = ["compound"] + composition + ["enthalpy_formation_atom", 'auid']
            title = "phase diagram"+"("+",".join(composition)+")"+"-Aflow"
        if self.database == 'mp':
            customdata = ["pretty_formula"] + composition + ["formation_energy_per_atom", 'material_id']
            title = "phase diagram"+"("+",".join(composition)+")"+"-Materials Project"
            pass
        fig = go.Figure()
        if len(composition) == 2:


            fig.add_trace(go.Scatter(
                                    x = self.df[composition[1]],
                                    y = self.df[customdata[3]],
                                    mode="markers",
                                    marker = dict(
                                           size = 4,
                                           color = self.df[customdata[3]],
                                           colorscale = 'Viridis',
                                           cmin=formation_e_range[0],
                                           cmax=formation_e_range[1],
                                           showscale=True

                                       ),
                                    customdata=self.df[customdata],
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
            if show_hull:
                fig.add_trace(go.Scatter(
                    x=self.np_hull[:, 0],
                    y=self.np_hull[:, 1],
                    mode="lines",
                    line=dict(
                        color='red',
                        width=2

                    ),
                    hoverinfo='none'
                ))
        if len(composition) == 3:
            sqrt3_over_2 = 0.5*math.sqrt(3)
            self.df['x'] = self.df[composition[1]] + 0.5 * self.df[composition[2]]
            self.df['y'] = sqrt3_over_2 * self.df[composition[2]]

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
            fig.add_trace(go.Scatter3d(x = self.df['x'], y = self.df['y'], z = self.df[customdata[4]],
                                       mode="markers",
                                       marker = dict(
                                           size = 4,
                                           color = self.df[customdata[4]],
                                           colorscale = 'Viridis',
                                           cmin=formation_e_range[0],
                                           cmax=formation_e_range[1],
                                           showscale=True
                                       ),
                                       customdata=self.df[customdata],
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
            self.df.pop("x")
            self.df.pop("y")
            if show_hull:
                for s in self.point_list:
                    s = np.append(s, s[0])
                    fig.add_trace(go.Scatter3d(
                        x=self.np_hull[s, 0], y=self.np_hull[s, 1], z=self.np_hull[s, 2],

                        mode="lines",
                        line=dict(
                            color='red',
                            width=2

                        ),
                        hoverinfo='none'
                ))
        if show_figure:
            fig.show()
        return fig



if __name__ == "__main__":
    composition = ["Ni", "Fe"]
    Rester = APIRester(composition, 'oq')
    print(Rester.df)


