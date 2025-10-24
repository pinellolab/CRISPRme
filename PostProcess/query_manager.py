"""
This file contains two main functions to deal with querying of the database
in all the specific cases we want.
The two functions are called by the main app file and take in input the parameters 
that are chosen by the user. 
If the parameters are linked to a sort query without shold only the parameters needed are in input and passed to noshold().
Otherwise, we have all the parameters including sholds and the other stuff we need for querying. shold()
The parameters in input are modified with a prefix and a suffix to be coherent with the database structure (column names).
So, if the name of columns in the input file are modified, also these parameters inside here need to be modified.
In output the two functions return a dict that is then included in the table in the main app file.
"""

import dash
import plotly
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_table
import pandas as pd
import sqlite3

guide_column = "Spacer+PAM"


def shold(
    target,
    n_clicks,
    page_current,
    page_size,
    radio_order,
    orderdrop,
    sholddrop,
    maxdrop,
    asc1,
    url_job,
    guide,
    current_working_directory,
):
    # query with threshold
    # print('entro shold')
    if asc1 == None:
        asc1 = "DESC"

    url = url_job[5:]

    path = current_working_directory + "/Results/" + url + "/." + url + ".db"

    conn = sqlite3.connect(path)
    c = conn.cursor()
    param = [
        guide,
        page_size,
        page_current * page_size,
    ]

    if target == "CFD":
        radio_order = str(radio_order) + "_(highest_CFD)"
        orderdrop = str(orderdrop) + "_(highest_CFD)"
    elif target == "fewest":
        radio_order = str(radio_order) + "_(fewest_mm+b)"
        orderdrop = str(orderdrop) + "_(fewest_mm+b)"
    elif target == "CRISTA":
        if "CFD" in radio_order:
            radio_order = str(radio_order).replace("CFD", "CRISTA")
        try:
            if "CFD" in orderdrop:
                orderdrop = str(orderdrop).replace("CFD", "CRISTA")
        except:
            pass
        radio_order = str(radio_order) + "_(highest_CRISTA)"
        orderdrop = str(orderdrop) + "_(highest_CRISTA)"

    # print(radio_order, orderdrop)

    if orderdrop != "None_(highest_CFD)":  # ordinamento doppio
        if maxdrop == None:
            maxdrop = 1000
        # if maxdrop == None:  # min
        #     df = pd.read_sql_query("SELECT * FROM final_table WHERE \"{}\"=? AND {}>={} ORDER BY \"{}\" {},\"{}\" {} LIMIT ? OFFSET ?".format(
        #         guide_column, radio_order, sholddrop, radio_order, asc1, orderdrop, asc1), conn, params=param)
        # else:  # min e max
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? AND "{}" BETWEEN {} AND {} ORDER BY "{}" {},"{}" {}  LIMIT ? OFFSET ?'.format(
                guide_column,
                radio_order,
                sholddrop,
                maxdrop,
                radio_order,
                asc1,
                orderdrop,
                asc1,
            ),
            conn,
            params=param,
        )
    else:  # ordinamento singolo
        if maxdrop == None:  # min
            maxdrop = 1000
            # df = pd.read_sql_query("SELECT * FROM final_table WHERE \"{}\"=? AND \"{}\">={} ORDER BY \"{}\" {} LIMIT ? OFFSET ?".format(
            #     guide_column, radio_order, sholddrop, radio_order, asc1), conn, params=param)  # ok
        # else:  # min e max
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? AND "{}" BETWEEN {} AND {} ORDER BY "{}" {} LIMIT ? OFFSET ?'.format(
                guide_column, radio_order, sholddrop, maxdrop, radio_order, asc1
            ),
            conn,
            params=param,
        )  # ok

    conn.commit()
    conn.close()
    data = df
    return data


def noshold(
    target,
    n_clicks,
    page_current,
    page_size,
    radio_order,
    orderdrop,
    asc1,
    url_job,
    guide,
    current_working_directory,
):
    # print(asc1)
    # print('entro noshold')
    if asc1 == None:
        asc1 = "DESC"
    url = url_job[5:]
    path = current_working_directory + "/Results/" + url + "/." + url + ".db"

    conn = sqlite3.connect(path)
    c = conn.cursor()
    param = [
        guide,
        page_size,
        page_current * page_size,
    ]
    # print("noshold")

    if target == "CFD":
        radio_order = str(radio_order) + "_(highest_CFD)"
        orderdrop = str(orderdrop) + "_(highest_CFD)"
    elif target == "fewest":
        radio_order = str(radio_order) + "_(fewest_mm+b)"
        orderdrop = str(orderdrop) + "_(fewest_mm+b)"
    elif target == "CRISTA":
        if "CFD" in radio_order:
            radio_order = str(radio_order).replace("CFD", "CRISTA")
        try:
            if "CFD" in orderdrop:
                orderdrop = str(orderdrop).replace("CFD", "CRISTA")
        except:
            pass
        radio_order = str(radio_order) + "_(highest_CRISTA)"
        orderdrop = str(orderdrop) + "_(highest_CRISTA)"

    # print(radio_order, orderdrop)

    if orderdrop != "None_(highest_CFD)":
        # query con multiordinamento
        # print('double ordering')
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? ORDER BY "{}" {}, "{}" {} LIMIT ? OFFSET ?'.format(
                guide_column, radio_order, asc1, orderdrop, asc1
            ),
            conn,
            params=param,
        )
    else:
        # query con ordinamento singolo
        # print('single ordering')
        df = pd.read_sql_query(
            'SELECT * FROM final_table WHERE "{}"=? ORDER BY "{}" {} LIMIT ? OFFSET ?'.format(
                guide_column, radio_order, asc1
            ),
            conn,
            params=param,
        )

    # query senza soglia
    conn.commit()
    conn.close()
    data = df
    return data
