# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 17:40:00 2022

@author: JESUS ALEJANDRO COLIN VILCHIS
"""
import pandas as pd
import datetime
from concurrent.futures import ThreadPoolExecutor
import os
import re
import json
import logging
import numpy as np
from scipy import stats

class Files:
    
    @staticmethod
    def load_files(path:str) -> dict:
        df = pd.read_csv(path,sep="\t")
        return df
    
    @staticmethod
    def load_dict(dir:str) -> dict:
        """_summary_
            Cargamos los catalogos con la logica del analisis
        Args:
            None

        Returns:
            dictionaries: catalogos con la logica del analisis
        """
        dictionaries = dict()
        #dictionaries[f'{age_quinquennia}'] = age_quinquennia
        files = os.listdir(dir)
        for file in files:
            with open(f"{dir}/{file}", 'r',encoding="utf-8") as f:
                dictionaries[file[:-5]] = json.load(f)
        return dictionaries
    
    @staticmethod
    def complete(df:pd.DataFrame, *args, **kwargs) -> pd.DataFrame:
        """_summary_
            Cargamos los catalogos con la logica del analisis
        Args:
            df (DataFrame): el archivo a orginal a analizar
    
        Returns:
            dictionaries: catalogos con la logica del analisis
        """
        df_cleaned = df[df['Collection date'].str.len()>7]
        del df
        df_cleaned['date'] = [datetime.datetime.strptime(i,'%Y-%m-%d').isocalendar() for i in df_cleaned['Collection date']]
        with ThreadPoolExecutor(max_workers = 4) as executor:
            f1 = executor.submit(m_year,df_cleaned)
            f2 = executor.submit(m_week,df_cleaned)
            f3 = executor.submit(m_week_cont,df_cleaned)
            f4 = executor.submit(m_state,df_cleaned)
        df_cleaned['Year'] = f1.result()
        df_cleaned["week"] = f2.result() 
        df_cleaned["week_cont"] = f3.result() 
        df_cleaned['state'] = f4.result()
        return df_cleaned
    
    @staticmethod
    def normalize(df:pd.DataFrame,dictionaries:dict) -> pd.DataFrame:
        #dic_normalizer = {category : None for category in ("state","group_age","group_patient_status","age")}
        with ThreadPoolExecutor(max_workers = 4) as executor:
            f1 = executor.submit(linage_normalized,df['Lineage'],dictionaries["variant_types"])
            f2 = executor.submit(age_normalized,df[['Patient age','Gender']],dictionaries["age_unification"])
            f3 = executor.submit(state_group_patient_status,df['Patient status'],dictionaries["patient_status"])
            f4 = executor.submit(state_normalized,df['state'],dictionaries["states_types"],dictionaries["unique_states_types"])
        df['variant_type'] = f1.result()           
        df[['age','group_age']] = f2.result()          
        df['group_status'] =  f3.result()              
        df[['state_key','region_key']] = f4.result()
        #df = self.filter_age(df,*args,**kwargs)  
        return df
    
    @staticmethod
    def filter_age(df,age):
            if bool(len(age)):
              if age == "under18":
                df = df[df['age']<18]
              elif age == "over18":
                df = df[df['age']>=18]
            return df



class Counter:

    
    @staticmethod
    def get_mutations(df,parts_protein:dict,protein:tuple = ("Spike_","E_","N_","M_","NSP")) -> dict:
        search_amino = lambda segment,sustitution : [aa for user in sustitution 
                                                     for aa in user.replace(r'(','').replace(r')','').split(',')
                                                        if bool(re.search(segment,aa))]
        position = lambda val : int(re.findall("\d+",val.split('_')[1])[0])
        table = lambda val,count,aa: [(val[i],count[i],val[i].split('_')[0],position(val[i]), p_p(val[i],aa)) for i in range(len(val))]
        table2 = lambda val,count: [(val[i],count[i],val[i].split('_')[0],position(val[i])) for i in range(len(val))]
        v_c = lambda val: np.unique(val, return_counts=True)
        # devuelve la porcion de proteina a la que pertenece , dada una mutación y si tipo de segmento
        p_p = lambda val,aa:  [pp[0] for pp in parts_protein[aa] if position(val) in range(pp[1],pp[2]+1)][0]
        
        with ThreadPoolExecutor(max_workers=len(protein)) as executor:
            futures = list()
            unique = list()
            futures.append(executor.submit(search_amino,protein[0],df['AA Substitutions']))
            futures.append(executor.submit(search_amino,protein[1],df['AA Substitutions']))
            futures.append(executor.submit(search_amino,protein[2],df['AA Substitutions']))
            futures.append(executor.submit(search_amino,protein[3],df['AA Substitutions']))
            futures.append(executor.submit(search_amino,protein[4],df['AA Substitutions']))
            values = [f.result() for f in futures]
            for i,aa in enumerate(protein[:-1]):
                v,c = v_c(values[i])
                unique.append(executor.submit(table,v,c,aa[:-1]))
            v,c = v_c(values[-1])
            unique.append(executor.submit(table2,v,c))
        amino = {protein[i]:u.result() for i,u in enumerate(unique)}
        return amino
    
    @staticmethod
    def get_table(amino,parts_protein:dict,txt:str = ""):
        #, #Guruprasad L. Human SARS CoV-2 spike protein mutations. Proteins. 2021 May;89(5):569-576. doi: 10.1002/prot.26042. Epub 2021 Jan 17. PMID: 33423311; PMCID: PMC8014176
        futures = list()
        create_df = lambda table, segment: pd.DataFrame(
                                                table,
                                                columns=[segment,"count","full","change","position"]
                                                ) if segment != "NSP" else pd.DataFrame(
                                                    table, 
                                                    columns=[segment,"count","full","change"])
        with ThreadPoolExecutor(max_workers=len(amino.keys())) as executor:
            for k,v in amino.items():
                futures.append(executor.submit(create_df,v,k))
            dfs = {k:futures[i].result() for i,k in enumerate(amino.keys())}
        save_files(dfs)
        return dfs
    
    @staticmethod
    def properties_mutations(df,aminos:list):
        df_amino_changes = dict()
        for change in set(aminos):
            df_amino_changes[change] = {i:0 for i in df['variant_type'].unique()} 
            for i in range(len(df)): # recorre todos los registros existentes  
                if re.search(str(change),df['AA Substitutions'].iloc[i]) is not None:
                    df_amino_changes[change][df['variant_type'].iloc[i]] +=1
        return df_amino_changes
    
def save_files(files):
    for k,v in files.items():
        v.to_csv(f"files/{k}.csv")
        #print(k)

def add_p_value_annotation(fig, array_columns, subplot=None, _format=dict(interline=0.07, text_height=1.07, color='black')):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)
    
    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare 
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
    # Specify in what y_range to plot for each pair of columns
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        y_range[i] = [1.01+i*_format['interline'], 1.02+i*_format['interline']]

    # Get values from figure
    fig_dict = fig.to_dict()

    # Get indices if working with subplots
    if subplot:
        if subplot == 1:
            subplot_str = ''
        else:
            subplot_str =str(subplot)
        indices = [] #Change the box index to the indices of the data for that subplot
        for index, data in enumerate(fig_dict['data']):
            #print(index, data['xaxis'], 'x' + subplot_str)
            if data['xaxis'] == 'x' + subplot_str:
                indices = np.append(indices, index)
        indices = [int(i) for i in indices]
        print((indices))
    else:
        subplot_str = ''

    # Print the p-values
    for index, column_pair in enumerate(array_columns):
        if subplot:
            data_pair = [indices[column_pair[0]], indices[column_pair[1]]]
        else:
            data_pair = column_pair

        # Mare sure it is selecting the data and subplot you want
        #print('0:', fig_dict['data'][data_pair[0]]['name'], fig_dict['data'][data_pair[0]]['xaxis'])
        #print('1:', fig_dict['data'][data_pair[1]]['name'], fig_dict['data'][data_pair[1]]['xaxis'])

        # Get the p-value
        pvalue = stats.ttest_ind(
            fig_dict['data'][data_pair[0]]['y'],
            fig_dict['data'][data_pair[1]]['y'],
            equal_var=False,
        )[1]
        if pvalue >= 0.05:
            symbol = 'ns'
        elif pvalue >= 0.01: 
            symbol = '*'
        elif pvalue >= 0.001:
            symbol = '**'
        else:
            symbol = '***'
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][0], 
            x1=column_pair[0], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        # Horizontal line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[0], y0=y_range[index][1], 
            x1=column_pair[1], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        # Vertical line
        fig.add_shape(type="line",
            xref="x"+subplot_str, yref="y"+subplot_str+" domain",
            x0=column_pair[1], y0=y_range[index][0], 
            x1=column_pair[1], y1=y_range[index][1],
            line=dict(color=_format['color'], width=2,)
        )
        ## add text at the correct x, y coordinates
        ## for bars, there is a direct mapping from the bar number to 0, 1, 2...
        fig.add_annotation(dict(font=dict(color=_format['color'],size=14),
            x=(column_pair[0] + column_pair[1])/2,
            y=y_range[index][1]*_format['text_height'],
            showarrow=False,
            text=symbol,
            textangle=0,
            xref="x"+subplot_str,
            yref="y"+subplot_str+" domain"
        ))
    return fig

   
    
def m_year(df_cleaned:pd.DataFrame) -> list:
    return [datetime.datetime.strptime(i,'%Y-%m-%d').year for i in df_cleaned['Collection date']]


def m_week(df_cleaned:pd.DataFrame) -> list:
    return [df_cleaned.date.iloc[i][1]  for i in range(len(df_cleaned))]


def m_week_cont(df_cleaned:pd.DataFrame) -> list:
    delay = {2020:0,2021:53,2022:106}
    return [df_cleaned.date.iloc[i][1] + delay[df_cleaned.date.iloc[i][0]] for i in range(len(df_cleaned))]


def m_state(df_cleaned:pd.DataFrame) -> list:
    return [i.split('/')[2]for i in df_cleaned['Location']]



def linage_normalized(df:pd.DataFrame, dic:dict) -> list:
    variant = list()
    for _df in df:
        if re.search('AY',_df):
            variant.append('Delta')
        elif re.search('BA',_df):
            variant.append('Omicron')
        else:
            for k,v in dic.items():
                if _df in v:
                    variant.append(k)
                    break
            else:
                logging.error("No se identifica la clasificación de: "+ _df)
    return variant


def get_quinquenios(val):
    age_quinquennia = {
        -1:[-1],
        1: [0],
        2:range(1,5),
        3:range(5,10),
        4:range(10,15),
        5:range(15,20),
        6:range(20,30),
        7:range(30,40),
        8:range(40,50),
        9:range(50,60),
        10:range(60,150)
    } 
    for k,v in age_quinquennia.items():
        if val in v:
            return (val,k)


def age_normalized(df,dic)-> list:
    """_summary_
        noramliza edad y edge group
    Args:
        df (_type_): _description_
        dic (_type_): _description_

    Returns:
        list: _description_
    """
    variant = list()
    for i in range(len(df)):  
        try:
            val = int(df['Patient age'].iloc[i])
            variant.append(get_quinquenios(val)) 
        except:
            if df['Patient age'].iloc[i].lower() in dic.keys():
                val = dic[df['Patient age'].iloc[i].lower()]
                variant.append(get_quinquenios(val))
            elif df['Patient age'].iloc[i] in ['Hospitalized', 'Ambulatory','Male','Female']:
                variant.append(get_quinquenios(int(df['Gender'].iloc[i])))
            elif bool(re.search('A',df['Patient age'].iloc[i])):
                val = re.search('A',df['Patient age'].iloc[i])
                variant.append(get_quinquenios(int(df['Patient age'].iloc[i][:val.span()[0]])))
            elif bool(re.search('onths',df['Patient age'].iloc[i])):# referencia a meses
                variant.append(get_quinquenios(0))
            elif bool(re.search('\.',df['Patient age'].iloc[i])):
                val = re.search('\.',df['Patient age'].iloc[i])
                variant.append(get_quinquenios(int(df['Patient age'].iloc[i][:val.span()[0]])))
            else:
                variant.append(get_quinquenios(-1))
                logging.error("No se identificó la edad:" + df['Patient age'].iloc[i])
    return variant


def state_normalized(df,dic,dic_k) -> list:
    """_summary_
        clave y valor real del estado unificado
    Args:
        df (DataFrame): _description_
        dic (dict): _description_
    Returns:
        list<tuple>: clave y valor de estados 
    """   
    evaluate = lambda _df: (_df,dic_k[str(_df)]) if str(_df) in dic_k.keys() else (99,"Extra")
    return [evaluate(dic[_df]) for _df in df] 


def state_group_patient_status(df:pd.DataFrame, dic:dict) -> list:
    variant = list()
    for _df in df:
        try:
            variant.append(int(dic[_df.lower()]))
        except:
            if re.search('home',_df.lower()):
                variant.append(1)
            else:
                logging.error("No se identifica el estado:" + _df)
    return variant


    
def counter(df,dic):    
        amino_unique = Counter.get_mutations(df,dic)
        return Counter.get_table(amino_unique,dic)

def autoadjust(df, sheetname, writer):
    df.to_excel(writer, sheet_name=sheetname, index=False, startrow=1)
    for column in df.columns:
        column_length =  len(str(column))+1
        col_idx = df.columns.get_loc(column)
        writer.sheets[sheetname]

def get_voc(normalized_file):
    voc = {voc:[] for voc in normalized_file['variant_type'].unique()}
    df_voc = normalized_file.groupby(['week_cont','variant_type'],as_index=False).size()
    existe=list(df_voc['variant_type'].unique())
    for i in range(int(df_voc['week_cont'].max())): # for each week
        for j in range(len(df_voc)): # for table
            if i == df_voc.iloc[j]['week_cont']:
                for mivoc in df_voc['variant_type'].unique(): 
                    if mivoc == df_voc.iloc[j]['variant_type']:
                        voc[mivoc].append(df_voc.iloc[j]['size'])
                        existe.remove(df_voc.iloc[j]['variant_type'])

        for mivoc in existe: # we fill in the empty spaces
            voc[mivoc].append(0)
            
        existe=list(df_voc['variant_type'].unique())
    return voc

def get_clado(normalized_file):
    clado = {clado:[] for clado in normalized_file['Clade'].unique()}
    df_clado = normalized_file.groupby(['week_cont','Clade'],as_index=False).size()

    existe=list(df_clado.Clade.unique())
    for i in range(int(df_clado['week_cont'].max())): # for each week
        for j in range(len(df_clado)): # for table
            if i == df_clado.iloc[j]['week_cont']:
                for miclado in df_clado.Clade.unique(): # for clades
                    if miclado == df_clado.iloc[j]['Clade']:
                        clado[miclado].append(df_clado.iloc[j]['size'])
                        existe.remove(df_clado.iloc[j]['Clade'])

        for miclado in existe: # we fill in the empty spaces
            clado[miclado].append(0)
        existe=list(df_clado.Clade.unique())
    return clado

def get_states(normalized_file, dictionaries):
    table = list()
    for week in range(normalized_file['week_cont'].min(),normalized_file['week_cont'].max()+1):
        for s_k in normalized_file['state_key'].unique():
            try:
                val = normalized_file.query(f' state_key == {s_k} and week_cont == {week}').groupby(['variant_type']).size().transform(lambda x: x/x.sum())
                table.append([week,dictionaries["unique_states_types"][str(s_k)],val[val>0.5].keys()[0]])
            except IndexError:
                table.append([week,dictionaries["unique_states_types"][str(s_k)],'None'])
            except Exception as ex:
                print(ex,week,s_k)
    return table
    