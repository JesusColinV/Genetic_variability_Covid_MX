# %% [markdown]
# 
# Created on Sat Sep 24 17:40:00 2022
# 
# @author: JESUS ALEJANDRO COLIN VILCHIS
# 
# # target:
#   ## Analysis of the gisaid data, until August 31, 2022
# 

# %%
# Pandas
import pandas as pd

# Python
import warnings
import os

from concurrent.futures import ThreadPoolExecutor

# Plotly
import plotly.graph_objects as go
import plotly.express as px
import datetime

# Developed
from utils import Files as f
from utils import Counter as c
from utils import counter, autoadjust, add_p_value_annotation
from utils import get_voc

# %% [markdown]
# Cleaning of library alerts

# %%
warnings.filterwarnings('ignore')

# %% [markdown]
# Loading, preprocessing, normalization and data cleaning.

# %%
files = os.listdir("files")

# Reading
original_file = pd.read_csv("files/gisaid_hcov-19_2022_08_17_04.tsv",sep="\t")
dictionaries = f.load_dict(dir = "json")

# Preprocessing and standardization
completed_file = f.complete(original_file)
normalized_file = f.normalize(completed_file , dictionaries)

# Tables are filtered according to age
df_under18 = f.filter_age(normalized_file,"under18")
df_over18 = f.filter_age(normalized_file,"over18")


# %% [markdown]
# The size of the recordsets according to their age are as follows.

# %%
n_under18 = len(df_under18)
n_over18 = len(df_over18)

print(
    "Under_18: \t", n_under18 ,"\n", 
    "Over_18: \t", n_over18 ,"\n", 
    "Total: \t", len(normalized_file)
    )

# %% [markdown]
# Counting mutations, and calculating percentages

# %%
if "table_mutation_count.csv" not in files:
    # Mutations grouped by protein of occurrence are counted
    df_count_all = counter(normalized_file,dictionaries["proteins"])
    df_count_under18 = counter(df_under18,dictionaries["proteins"])
    df_count_over18 = counter(df_over18,dictionaries["proteins"])


    # column name change
    mk_table = lambda df:pd.concat(
        [df[k].rename(columns = {k:'mutation'}) 
         for k in df.keys()] 
        ).set_index('mutation')
    
    df_mutations_total = mk_table(df_count_all)
    df_mutations_under18 = mk_table(df_count_under18)
    df_mutations_over18 = mk_table(df_count_over18)

    # column name change
    rname = lambda df,txt="" : df.rename(columns = {i:i+txt if i!= 'mutation'else i for i in df.columns})
    df_m_total= rname(df_mutations_total)
    df_m_under18= rname(df_mutations_under18,"_under18")
    df_m_over18= rname(df_mutations_over18,"_over18")

    # Join the tables
    df_mutations = df_m_total.join(df_m_under18).join(df_m_over18)
    
    # Percentage calculation
    df_mutations['percentage_under18'] = [(df_mutations['count_under18'].iloc[i]/n_under18)*100 for i in range(len(df_mutations))]
    df_mutations['percentage_over18'] = [(df_mutations['count_over18'].iloc[i]/n_over18)*100  for i in range(len(df_mutations))]
    df_mutations = df_mutations[['change','position','count','count_under18','percentage_under18','count_over18','percentage_over18','full']].fillna(0)

    # save the results
    df_mutations.to_csv('files/table_mutation_count.csv')
    
else:
    
    #save the results
    df_mutations = pd.read_csv("files/table_mutation_count.csv")


# %%
if "tabla_mutations_variant.csv" not in files:
    
    
    with ThreadPoolExecutor(max_workers = 3) as executor:
        f1 = executor.submit(c.properties_mutations,normalized_file,df_mutations['change'].tolist())
        f2 = executor.submit(c.properties_mutations,df_under18,df_mutations['change'].tolist())
        f3 = executor.submit(c.properties_mutations,df_over18,df_mutations['change'].tolist())

    table_mutations_variant = f1.result()
    table_mutations_variant_under = f2.result() 
    table_mutations_variant_over18 = f3.result() 

    # Dictionaries are converted into dataframes
    table_mutations_variant = pd.DataFrame(table_mutations_variant)
    table_mutations_variant_under18 = pd.DataFrame(table_mutations_variant_under)
    table_mutations_variant_over18 = pd.DataFrame(table_mutations_variant_over18)

    # Join the tables
    table_mv = pd.concat(
        [
            table_mutations_variant.transpose(), 
            table_mutations_variant_under18.transpose(), 
            table_mutations_variant_over18.transpose()], 
        axis=1
        )
    
    #save the results
    table_mv.to_csv("files/tabla_mutations_variant.csv")
    
else:
    
    #save the results
    table_mv = pd.read_csv("files/tabla_mutations_variant.csv")

# %% [markdown]
# # Graphics

# %% [markdown]
# From the age distributions for each variant, a box plot is generated.

# %%
fig = go.Figure()

# For each variant
for vt in normalized_file["variant_type"].unique():
    
    fig.add_trace(
        go.Box(
            y=normalized_file[
                normalized_file["variant_type"] == vt]['age'],
            name=vt,
            boxpoints='outliers'
        ))
    
fig = add_p_value_annotation(fig, [[0,1], [0,2], [0,3], [0,4], [0,5]])

fig = fig.update_layout(
    autosize=False,
    width=1200,
    showlegend=True,
    template="plotly_white"
    )

fig.show()

# %% [markdown]
# Cleaning week no data uploaded
# * This is for display only

# %%


voc = get_voc(normalized_file)
voc = pd.DataFrame(voc)
voc.iloc[106] = voc.iloc[105]
voc = voc.iloc[8:]

# %% [markdown]
# 100 Percent Stacked Area Chart from variant counts

# %%
df = voc.divide(voc.sum(axis=1), axis=0)
df.to_csv("files/variant_semana.csv")
fig = px.area(
    df,
    line_shape="spline",
    template="plotly_white"
    )
fig.show()


# %% [markdown]
# Cleaning week no data uploaded
# * This is for display only

# %%
from utils import get_clado
clado = get_clado(normalized_file)
clado = pd.DataFrame(clado)
clado.iloc[106] = clado.iloc[105]#[104:107]
clado = clado.iloc[8:]
clado.to_csv("files/clado_semana.csv")

# %% [markdown]
# 100 Percent Stacked Area Chart from clade counts

# %%
fig = px.area(
    clado.divide(
        clado.sum(axis=1), 
        axis=0),
    line_shape="spline",
    template="plotly_white"
    )
fig.show()

# %%
from utils import count_representative_state, count_representative_country

writer =  pd.ExcelWriter('files/estados.xlsx')
count_representative_state(
  writer,
  normalized_file,
  dictionaries
  )


count_representative_country(
    writer,
    normalized_file,
)

writer.save()

# %%

week_state = dict()

for i in dictionaries["unique_states_types"].values():
  week_state[i] = pd.read_excel(
    open('files/estados.xlsx', 'rb'), 
    header = 1 ,
    sheet_name=str(i))[['week','predominant']]

# %%
from utils import get_states

states = get_states(normalized_file, dictionaries)
df_states = pd.DataFrame(states,columns=['week_cont', 'state_key' ,'variant_type'])

# %%
df_states.groupby('variant_type').size()

# %%
fig = px.scatter(
    df_states.query('state_key != "Extra"'), 
    y='state_key', 
    x= 'week_cont', 
    color='variant_type', 
    size_max=60,
    template='plotly_white', 
    title="variante predominante por Estado",
    color_discrete_sequence=px.colors.qualitative.Vivid
    )

fig.update_layout(yaxis={
                    'categoryorder': 'array', 
                    'categoryarray': [
                        'Baja California','Baja California Sur',
                        'Sonora','Chihuahua','Coahuila',
                        'Nuevo Leon','Colima','Chiapas',
                        'Tamaulipas','Aguascalientes','Hidalgo',
                        'Guanajuato','Durango','Jalisco',
                        'Sinaloa','Ciudad de Mexico',
                        'Estado de Mexico','Michoacan',
                        'Morelos','Nayarit','Oaxaca',
                        'Puebla','Queretaro','San Luis Potosi',
                        'Guerero','Tabasco','Tlaxcala',
                        'Veracruz','Zacatecas','Campeche',
                        'Quintana Roo','Yucatan'
                        ]},
    width=1000,
    height=1000)

fig.show()

# %%
df_states.to_csv("files/estados_por_puntos.csv")

# %%
df_states.groupby(['state_key','variant_type']).size().to_csv("files/estados_cuenta.csv")

# %%


# %%
path = "files/220901COVID19MEXICO.csv"
df_sinav = pd.read_csv(path)
df_sinav = df_sinav[df_sinav['FECHA_INGRESO'].notna()]
df_sinav = df_sinav.query('CLASIFICACION_FINAL == 3 and EDAD < 18 ')
df_sinav.head()

# %%
df_sinav['Year'] = [datetime.datetime.strptime(i,'%Y-%m-%d').year for i in df_sinav['FECHA_INGRESO']]
df_sinav['date'] = [datetime.datetime.strptime(i,'%Y-%m-%d').isocalendar() for i in df_sinav['FECHA_INGRESO']]

delay = {
    2020:0,
    2021:53,
    2022:106
    }

df_sinav["week"] = [df_sinav.date.iloc[i][1]  for i in range(len(df_sinav))]
df_sinav["week_cont"] = [df_sinav.date.iloc[i][1] + delay[df_sinav.date.iloc[i][0]] for i in range(len(df_sinav))]

# %%
len(df_sinav)

# %%
df_sinav['TIPO_PACIENTE_AMP'] = [
    df_sinav['TIPO_PACIENTE'].iloc[i] if df_sinav['FECHA_DEF'].iloc[i] == '9999-99-99' else 3  
    for i in range(len(df_sinav))
    ]
df_sinav['UCI_INTUBADO'] = [1 if df_sinav['UCI'].iloc[i] == 1 or df_sinav['INTUBADO'].iloc[i] == 1 else 2 
                            for i in range(len(df_sinav))]
df_sinav['EMBARAZO_bool'] = [1 if df_sinav['EMBARAZO'].iloc[i] == 1  else 2 for i in range(len(df_sinav))]

# %%
#
#
#table = list()
#for i in range(len(df_sinav)):
#  try:
#    val = normalized_file.query(
#        f" state_key == { df_sinav['ENTIDAD_RES'].iloc[i] } and week_cont == { df_sinav['week_cont'].iloc[i] }"
#        ).groupby(
#          ['variant_type']
#            ).size().transform(lambda x: x/x.sum())
#    table.append(val[val>0.5].keys()[0])
#  except IndexError:
#    table.append('None')
#  except Exception as ex:
#    print(ex)
#
#df_sinav['Predominant'] = table

# %%
def calculate_sinav(df_sinav):
    table = list()
    for i in range(len(df_sinav)):
        try:
            val = normalized_file.query(
                f" state_key == { df_sinav['ENTIDAD_RES'].iloc[i] } and week_cont == { df_sinav['week_cont'].iloc[i] }"
                ).groupby(
                ['variant_type']
                    ).size().transform(lambda x: x/x.sum())
            table.append(val[val>0.5].keys()[0])
        except IndexError:
            table.append('None')
        except Exception as ex:
            table.append('None')
            print(ex)
    return table

# %%
df_2 = df_sinav 
df_sinav = df_sinav[1:]
n = len(df_sinav)/5
n

# %%
with ThreadPoolExecutor(max_workers=10) as executor:
    futures = list()
    futures.append(executor.submit(calculate_sinav,df_sinav.iloc[:int(n)]))
    futures.append(executor.submit(calculate_sinav,df_sinav.iloc[int(n):int(2*n)]))
    futures.append(executor.submit(calculate_sinav,df_sinav.iloc[int(2*n):int(3*n)]))
    futures.append(executor.submit(calculate_sinav,df_sinav.iloc[int(3*n):int(4*n)]))
    futures.append(executor.submit(calculate_sinav,df_sinav.iloc[int(4*n):]))

# %%
r = list()
for i in futures:
    r = r + i.result()
df_sinav['Predominant'] = r


# %%
df_sinav.to_csv("files/sinav_predominant.csv")


