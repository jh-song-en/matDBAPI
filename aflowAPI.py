import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json, sys, os
from urllib.request import urlopen
from itertools import combinations

def get_data_from_aflow(composition):
    compo = ",".join(composition)
    SERVER = "http://aflow.org"
    API = "/API/aflux/v1.0/?"
    MATCHBOOK = f"species({compo}),nspecies({len(composition)}),enthalpy_formation_atom,composition,species,stoichiometry,Bravais_lattice_relax"
    DIRECTIVES = "$paging(0)"
    SUMMONS = MATCHBOOK + "," + DIRECTIVES
    response = json.loads(urlopen(SERVER + API + SUMMONS).read().decode("utf-8"))
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
    if database == "af":
        composition_combination = get_every_combination(composition_list)
        dictionary_list = []
        for composition in composition_combination:
            dictionary_list += get_data_from_aflow(composition)
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
        pass

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


if __name__ == "__main__":
    composition = ["Ni", "Fe"]
    data = get_data_by_composition(composition, 'af')
    print(data)


