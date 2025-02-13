# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:48:20 2021

@author: franc
"""

import sqlite3
import time
import sys
import pandas as pd

# Creation of the connection to a new database and table

"""
This file contains the step to create the database from the input file sg1617.total.scored.txt
There's a simple preprocessing for the column names and after that the creation of the table with the names read inside a list.
So if you need to add a column you have to add a column element "{} TYPE" inside the string q, where type can be in this case TEXT,INTEGER,NUMERIC.
Then the rows of the file are read with list and the for loop and added in the insert() function. Again if you add some column you need to add a "?".
Finally there's the creation of indexes for the fast visualization of the results of the queries. There are really specific then if you want to add something 
you need to make couples or triples of indexes for covering the columns you need to order by query. If you change columns names you also have update names of parameters of indexes
For example  c.execute("CREATE INDEX IF NOT EXISTS rgt ON final_table(Real_Guide_1,Total_1)") create an index to speed up the order of columns by Total_1. Always include Real_Guide_1 
to permit the use of multiple guides. Same type of indexes are created in the case of multiple order using in this case "triples".

"""


def dict_pd_dtypes_to_sql_types(pd_dtype):
    if pd_dtype == "O":
        return "TEXT"
    elif pd_dtype == "int64":
        return "INTEGER"
    elif pd_dtype == "float64":
        return "NUMERIC"
    else:
        return "TEXT"


# fileIn = sys.argv[1]
# fileOut = sys.argv[2]

fileIn = sys.argv[1]
fileOut = sys.argv[2]

conn = sqlite3.connect(f"{fileOut}.db")
c = conn.cursor()
# tot_lines = 0
# with open(f"{fileIn}", "r") as f:
#     header = f.readline().split()
#     for line in f:
#         tot_lines += 1
db_schema = []
try:
    df = pd.read_csv(fileIn, sep="\t", index_col=False, na_filter=False, nrows=1000)
    types = [dict_pd_dtypes_to_sql_types(pd_dtype) for pd_dtype in df.dtypes]

    for i, col in enumerate(df.columns):
        db_schema.append(f'"{col}" {types[i]}')

    db_schema = ", ".join(db_schema)
except:
    with open(fileIn) as f_in:
        cols = f_in.readline().strip().split("\t")
        for col in cols:
            db_schema.append(f'"{col}" TEXT')
        db_schema = ", ".join(db_schema)
# print(db_schema)

q = f"CREATE TABLE IF NOT EXISTS final_table ({db_schema})"

c.execute(q)


def insert(data):
    print("Inserting chunk of data")
    question_marks = ["?"] * len(data[0])
    question_marks = ",".join(question_marks)
    c.executemany(f"INSERT INTO final_table VALUES ({question_marks})", data)
    conn.commit()


with open(f"{fileIn}", "r") as f:
    a = []
    data = []
    y = 0
    start_time = time.time()
    for line in f:
        if y == 0:
            line.split()
            y = y + 1
        else:
            g = line.strip().split("\t")
            data.append(g)
            if len(data) == 500000:  # stop to small enough chunks of data
                insert(data)
                data = []
    if len(data) != 0:  # insert last chunk of data
        insert(data)
    # create indexes
    print("Now creating indexes")

    guide_column = "Spacer+PAM"
    mm_column = "Mismatches_(highest_CFD)"
    blg_column = "Bulges_(highest_CFD)"
    total_column = "Mismatches+bulges_(highest_CFD)"
    cfd_column = "CFD_score_(highest_CFD)"
    risk_cfd_column = "CFD_risk_score_(highest_CFD)"
    pos_column = "Start_coordinate_(highest_CFD)"
    samples_column = "Variant_samples_(highest_CFD)"
    chrom_column = "Chromosome"
    bulge_t_column = "Bulge_type_(highest_CFD)"

    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_mm ON final_table("{guide_column}","{mm_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_blg ON final_table("{guide_column}","{blg_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_tot ON final_table("{guide_column}","{total_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_cfd ON final_table("{guide_column}","{cfd_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_risk ON final_table("{guide_column}","{risk_cfd_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_chrom_pos ON final_table("{guide_column}","{chrom_column}","{pos_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_samples ON final_table("{guide_column}","{samples_column}")'
    )

    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_mm_blg_blgt ON final_table("{guide_column}","{mm_column}","{blg_column}","{bulge_t_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_mm_tot ON final_table("{guide_column}","{mm_column}","{total_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_mm_cfd ON final_table("{guide_column}","{mm_column}","{cfd_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_blg_tot ON final_table("{guide_column}","{blg_column}","{total_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_blg_cfd ON final_table("{guide_column}","{blg_column}","{cfd_column}")'
    )
    c.execute(
        f'CREATE INDEX IF NOT EXISTS g_tot_cfd ON final_table("{guide_column}","{total_column}","{cfd_column}")'
    )

    conn.commit()
    conn.close()
