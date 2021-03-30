import sqlite3
import time
import sys
#Creation of the connection to a new database and table

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

fileIn = sys.argv[1]
fileOut = sys.argv[2]

conn = sqlite3.connect(f'{fileOut}.db')
c = conn.cursor()
tot_lines = 0
with open(f"{fileIn}", "r") as f:
    first = f.readline().split()
    for line in f:
        tot_lines += 1

a=int(len(first)/2)
l1 = first[0:a]
i=0
li=[]
for i in l1 :
  i=i.replace('#','')
  i=i+'_1'
  li.append(i)

l2 = first[a:len(first)]
i=0
for i in l2 :
  i=i.replace('#','')
  i=i+'_2'
  li.append(i)


q = """CREATE TABLE IF NOT EXISTS final_table (
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  INTEGER,
  {}  INTEGER,
  {}  TEXT,
  {}  INTEGER,
  {}  INTEGER,
  {}  INTEGER,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  NUMERIC,
  {}  NUMERIC,
  {}  NUMERIC,
  {}  NUMERIC,
  
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  INTEGER,
  {}  INTEGER,
  {}  TEXT,
  {}  INTEGER,
  {}  INTEGER,
  {}  INTEGER,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  TEXT,
  {}  NUMERIC,
  {}  NUMERIC,
  {}  NUMERIC,
  {}  NUMERIC
)""".format(*li);

c.execute(q)

def insert(data):
  s=time.time()
  print("Inserting chunk of data")
  c.executemany('INSERT INTO final_table VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', data)
  conn.commit()

with open(f"{fileIn}", "r") as f:
    a=[]
    data=[]
    y=0
    start_time = time.time()
    for line in f:
       if y==0:
         line.split()
         y=y+1
       else:
         g=line.split('\t')
         data.append(g)
         if len(data)==tot_lines//10:
           insert(data)
           data=[]
    insert(data)
    # create indexes
    c.execute("CREATE INDEX IF NOT EXISTS g ON final_table(Real_Guide_1)")

    c.execute("CREATE INDEX IF NOT EXISTS rgm ON final_table(Real_Guide_1,Mismatches_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgc ON final_table(Real_Guide_1,CFD_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgt ON final_table(Real_Guide_1,Total_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbs ON final_table(Real_Guide_1,Bulge_Size_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgcrs ON final_table(Real_Guide_1,Highest_CFD_Risk_Score_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rghcar ON final_table(Real_Guide_1,Highest_CFD_Absolute_Risk_Score_1)")

    c.execute("CREATE INDEX IF NOT EXISTS rgm2 ON final_table(Real_Guide_1,MMBLG_Mismatches_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgc2 ON final_table(Real_Guide_1,MMBLG_CFD_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgt2 ON final_table(Real_Guide_1,MMBLG_Total_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbs2 ON final_table(Real_Guide_1,MMBLG_Bulge_Size_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgcrs2 ON final_table(Real_Guide_1,MMBLG_CFD_Risk_Score_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rghcar2 ON final_table(Real_Guide_1,MMBLG_CFD_Absolute_Risk_Score_2)")

    c.execute("CREATE INDEX IF NOT EXISTS rgmcfd ON final_table(Real_Guide_1,Mismatches_1,CFD_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgmbul ON final_table(Real_Guide_1,Mismatches_1,Bulge_Size_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgmtot ON final_table(Real_Guide_1,Mismatches_1,Total_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgcmi ON final_table(Real_Guide_1,CFD_1,Mismatches_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgcbulg ON final_table(Real_Guide_1,CFD_1,Bulge_Size_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgctot ON final_table(Real_Guide_1,CFD_1,Total_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgtmis ON final_table(Real_Guide_1,Total_1,Mismatches_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgtcfd ON final_table(Real_Guide_1,Total_1,CFD_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgtbulg ON final_table(Real_Guide_1,Total_1,Bulge_Size_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbsmism ON final_table(Real_Guide_1,Bulge_Size_1,Mismatches_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbscfd ON final_table(Real_Guide_1,Bulge_Size_1,CFD_1)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbstot ON final_table(Real_Guide_1,Bulge_Size_1,Total_1)")

    c.execute("CREATE INDEX IF NOT EXISTS rgmcfd2 ON final_table(Real_Guide_1,MMBLG_Mismatches_2,MMBLG_CFD_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgmbul2 ON final_table(Real_Guide_1,MMBLG_Mismatches_2,MMBLG_Bulge_Size_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgmtot2 ON final_table(Real_Guide_1,MMBLG_Mismatches_2,MMBLG_Total_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgcmi2 ON final_table(Real_Guide_1,MMBLG_CFD_2,MMBLG_Mismatches_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgcbulg2 ON final_table(Real_Guide_1,MMBLG_CFD_2,MMBLG_Bulge_Size_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgctot2 ON final_table(Real_Guide_1,MMBLG_CFD_2,MMBLG_Total_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgtmis2 ON final_table(Real_Guide_1,MMBLG_Total_2,MMBLG_Mismatches_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgtcfd2 ON final_table(Real_Guide_1,MMBLG_Total_2,MMBLG_CFD_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgtbulg2 ON final_table(Real_Guide_1,MMBLG_Total_2,MMBLG_Bulge_Size_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbsmism2 ON final_table(Real_Guide_1,MMBLG_Bulge_Size_2,MMBLG_Mismatches_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbscfd2 ON final_table(Real_Guide_1,MMBLG_Bulge_Size_2,MMBLG_CFD_2)")
    c.execute("CREATE INDEX IF NOT EXISTS rgbstot2 ON final_table(Real_Guide_1,MMBLG_Bulge_Size_2,MMBLG_Total_2)")

    conn.commit()
    conn.close()


