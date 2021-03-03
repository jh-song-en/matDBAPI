from pymatgen import MPRester
import requests
import pandas as pd
import plotly.graph_objects as go
import numpy as np


api_key = "P2tLmbTfRg9uMJDqzXV"

def get_mp_id_list_by_composition(composition_list, api_key = api_key):
    with MPRester("P2tLmbTfRg9uMJDqzXV") as m:
        mp_entries = m.get_entries_in_chemsys(composition_list)
    mp_id_list = []

    return mp_id_list

def get_XRD_prediction_by_id(material_id, rad_source = "Cu"):
    rv = requests.get(
        'https://materialsproject.org/materials/{material_id}/xrd?symbol={rad_source}'.format(
            material_id=material_id,
            rad_source=rad_source,
        ))
    try:
        data = rv.json()
        df = pd.DataFrame(data["pattern"], columns=data["meta"])
        return True, df
    except:
        return False, "unknown material id or radiation source"
class somesomesome():
    def __init__(self, api_key):
        self.api_key = api_key
        self.mp_id_list = []
        self.xrd_dictionary = {}
    def get_entries_by_composition(self, composition_list):
        with MPRester(self.api_key) as m:
            self.mp_entries = m.get_entries_in_chemsys(composition_list)
        self.mp_id_list = []
        for entry in self.mp_entries:
            self.mp_id_list.append(entry.entry_id)
        self.xrd_dictionary = {}
        for mp_id in self.mp_id_list:
            success, df = get_XRD_prediction_by_id(mp_id)
            if success:
                self.xrd_dictionary[mp_id] = df

    def show_XRD_graph(self):
        if self.xrd_dictionary != {}:
            fig = go.Figure()
            mp_id_list = list(self.xrd_dictionary.keys())
            for id in mp_id_list:
                dff = self.xrd_dictionary[id]
                x = dff["two_theta"].to_numpy()
                y = dff["amplitude"].to_numpy()
                fig.add_trace(go.Bar(visible=False, x=x, y=y, width=2, opacity=0.8,
                                     marker=dict(color='red',
                                                 line=dict(width=0))))
            fig.data[0].visible = True
            steps = []
            for i in range(len(fig.data)):
                step = dict(
                    method="update",
                    args=[{"visible": [False] * len(fig.data)}],  # layout attribute
                    label=mp_id_list[i]
                )
                step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
                steps.append(step)
            sliders = [dict(
                active=0,
                steps=steps
            )]
            fig.update_yaxes(scaleanchor="x", scaleratio=1)
            fig.update_layout(
                xaxis_title="2Theta",
                yaxis_title="Amplitude",
                template='plotly_white',
                sliders=sliders,
                xaxis_range=[0, 100]
            )
            fig.show()


if __name__ == "__main__":
    hehe = somesomesome("P2tLmbTfRg9uMJDqzXV")
    hehe.get_entries_by_composition(["Ni", "Fe"])
    hehe.show_XRD_graph()

