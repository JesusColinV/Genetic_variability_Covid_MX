# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 17:40:00 2022

@author: JESUS ALEJANDRO COLIN VILCHIS
"""

# Pandas
import pandas as pd

# Numpy
import numpy as np

# Python
import datetime
from concurrent.futures import ThreadPoolExecutor
import os

# Regular expressions
import re

# JSON
import json

# Logging
import logging

# SciPy
from scipy import stats

class Files:
    """
    File management class

    """
        
    @staticmethod
    def load_dict( dir : str ) -> dict:
        """ We load the metadata prepared in json format

        Args:
            dir (str): json file location
    
        Returns:
            dictionaries (dict): metadata
        """
        
        dictionaries = dict()
        files = os.listdir(dir)
        
        
        for file in files:
            
            with open(f"{dir}/{file}", 'r',encoding="utf-8") as f:
                
                dictionaries[file[:-5]] = json.load(f)
                
                
        return dictionaries
    
    @staticmethod
    def complete( df : pd.DataFrame) -> pd.DataFrame:
        """ Generate the columns with the data to be processed
        
        Args:
            df (DataFrame): dataframe with the original data
    
        Returns:
            df (DataFrame): dataframe with preprocessed data
        """
        
        # dates are cleared that comply with the format
        df_cleaned = df[df['Collection date'].str.len()>7]
        
        del df
        
        # change of string data to datetime, with iso calendar format
        df_cleaned['date'] = [
            datetime.datetime.strptime(i,'%Y-%m-%d').isocalendar() 
            for i in df_cleaned['Collection date']
            ]
        
        # Parallel processes
        with ThreadPoolExecutor(max_workers = 4) as executor:
        
            f1 = executor.submit(m_year,df_cleaned)
            f2 = executor.submit(m_week,df_cleaned)
            f3 = executor.submit(m_week_cont,df_cleaned)
            f4 = executor.submit(m_state,df_cleaned)
            
        # Generating time features
        df_cleaned['Year'] = f1.result()
        df_cleaned["week"] = f2.result() 
        df_cleaned["week_cont"] = f3.result() 
        df_cleaned['state'] = f4.result()
        
        
        return df_cleaned
    
    @staticmethod
    def normalize( df : pd.DataFrame, dictionaries : dict ) -> pd.DataFrame:
        """ From the cleaned data, we normalize columns that will be used for analysis.

        Args:
            df (pd.DataFrame): cleaned data
            dictionaries (dict): metadata catalog

        Returns:
            pd.DataFrame: standardized data
        """
        
        # Parallel processes        
        with ThreadPoolExecutor(max_workers = 4) as executor:
            
            f1 = executor.submit(linage_normalized,df['Lineage'],dictionaries["variant_types"])
            f2 = executor.submit(age_normalized,df[['Patient age','Gender']],dictionaries["age_unification"])
            f3 = executor.submit(state_group_patient_status,df['Patient status'],dictionaries["patient_status"])
            f4 = executor.submit(state_normalized,df['state'],dictionaries["states_types"],dictionaries["unique_states_types"])
        
        # Assignment to new columns
        df['variant_type'] = f1.result()           
        df[['age','group_age']] = f2.result()          
        df['group_status'] =  f3.result()              
        df[['state_key','region_key']] = f4.result()

        return df
    
    @staticmethod
    def filter_age(df: pd.DataFrame, age:str) -> pd.DataFrame:
        """ Data frame filtering from an age trigger

        Args:
            df (pd.DataFrame): Complete table
            age (str): trigger

        Returns:
            pd.DataFrame: Filtered table
        """
        # Filter application conditional
        if bool(len(age)):
            
            if age == "under18":
                df = df.query('age < 18')
                
            elif age == "over18":
                df = df.query('age >= 18')
                
        return df



class Counter:

    
    @staticmethod
    def get_mutations( 
        df : pd.DataFrame,
        parts_protein : dict,
        protein : tuple = ("Spike_","E_","N_","M_","NSP")
        ) -> dict:
        """ Mutation count by protein

        Args:
            df (pd.DataFrame): the standardized table of GISAID records
            parts_protein (dict): dictionary with proteins and their segments
            protein (tuple, optional): proteins of occurrence of mutations. Defaults to ("Spike_","E_","N_","M_","NSP").

        Returns:
            dict: table of each protein with the columns of the mutation count and its characteristics
        """
        # returns the amino acids that belong to the protein delivered as a parameter
        search_amino = lambda segment,aa_sustitutions : [
            aa for record in aa_sustitutions 
                for aa in record.replace(r'(','').replace(r')','').split(',')
                    if bool(re.search(segment,aa))
                    ]
        
        # Position of the amino acid mutation
        position = lambda val : int(re.findall("\d+",val.split('_')[1])[0])
        
        # from each amino acid in protein, the characteristics are counted and obtained
        # e.g. ['Spike_', 'count', 'full', 'change', 'position']
        table_protein = lambda val,count,aa: [
            (val[i],count[i],val[i].split('_')[0],
            position(val[i]), p_p(val[i],aa)) 
            for i in range(len(val))]
        
        # from each amino acid in nonstructural protein, the characteristics are counted and obtained
        # e.g. ['NSP', 'count', 'full', 'change']
        table_nsp = lambda val,count: [
            (val[i],count[i],val[i].split('_')[0],position(val[i])) 
            for i in range(len(val))]
        
        # returns unique mutations and their count
        v_c = lambda val: np.unique(val, return_counts=True)
        
        # returns the portion of protein to which it belongs, given a mutation and its segment name
        p_p = lambda val,aa:  [
            pp[0] 
            for pp in parts_protein[aa] 
            if position(val) in range(pp[1],pp[2]+1)
            ][0]
        
        
        with ThreadPoolExecutor(max_workers=len(protein)) as executor:
            
            # parallel processes for each type of protein
            futures = [
                executor.submit(
                            search_amino,
                            protein[i],
                            df['AA Substitutions']
                        ) for i in range(5)]
            
            values = [f.result() for f in futures]
            
            unique = list()
            
            for i, aa in enumerate(protein[:-1]):
                v,c = v_c(values[i])
                unique.append(executor.submit(table_protein,v,c,aa[:-1]))
                
            v,c = v_c(values[-1])
            unique.append(executor.submit(table_nsp,v,c))
            
        amino = {
            protein[i]:u.result() 
            for i,u in enumerate(unique)
            }
        
        return amino
    
    @staticmethod
    def dic2df(amino : dict) -> dict:
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
                
            table = {
                k:futures[i].result() 
                for i,k in enumerate(amino.keys())
                }

        return table
    
    @staticmethod
    def properties_mutations(
        df : pd.DataFrame,
        aminos : list
        ) -> dict:
        """_summary_

        Args:
            df (pd.DataFrame): _description_
            aminos (list): _description_

        Returns:
            dict: _description_
        """
        
        df_amino_changes = dict()
        
        # For each amino acid
        for change in set(aminos):
            
            
            df_amino_changes[change] = {i:0 for i in df['variant_type'].unique()} 
            
            # For each of the records
            for i in range(len(df)): 
                 
                 
                if re.search(str(change),df['AA Substitutions'].iloc[i]) is not None:
                    
                    df_amino_changes[change][df['variant_type'].iloc[i]] +=1
                    
                    
        return df_amino_changes
    

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
    error=[]
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
                error.append(_df)
                #logging.error("No se identifica la clasificación de: "+ _df)
    print(set(error))
    return variant


def get_quinquenios(val:int):
    """ Given an age number recognizes to which five-year period it belongs

    Args:
        val (int): Age

    Returns:
        val (int): Age
        k (int): five-year period in which it was classified
    """
    
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
    
    # goes through the five-year period to determine which classification it belongs to
    for k,v in age_quinquennia.items():
        if val in v:
            return (val,k)


def age_normalized(df : pd.DataFrame, dic : dict)-> list:
    """ Normalizes the age of patients and recognizes which group they belong to

    Args:
        df (pd.DataFrame): data from the gisaid database
        dic (dict): values to unify the age to a single format

    Returns:
        list: _description_
    """
    
    variant = list()
    
    
    for i in range(len(df)):  
        
        try:
            
            val = int(df['Patient age'].iloc[i])
            variant.append(get_quinquenios(val)) 
            
        except:
            
            # Adapts the age data, previously recognizing in which fields it can be found, 
            # due to a probelma in the storage
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
    """ Unify the names in the status column
    
    Args:
        df (DataFrame): _description_
        dic (dict): _description_
    Returns:
        list<tuple>: clave y valor de estados 
    """  
    
    #for each State value is assigned to a tuple the state and its key
    evaluate = lambda _df: (
        _df,dic_k[str(_df)]
        ) if str(_df) in dic_k.keys() else (99,"Extra")
    
    # The list of states with unified names is returned
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


    
def counter( df : pd.DataFrame, dic : dict):
    """ Generate the table with the count of amino acids, for each group (over or under 18 years)

    Args:
        df (pd.DataFrame): filtered group table
        dic (dict): metadata, with dictionaries 

    Returns:
        table (dict): table with mutation count and the rest of the characteristics
    """
    amino_unique = Counter.get_mutations(df,dic)
    table = Counter.dic2df(amino_unique)
    return table

def autoadjust( df : pd.DataFrame, sheetname : str, writer : pd.ExcelWriter ):
    """ formatting for Excel generation

    Args:
        df (pd.DataFrame): _description_
        sheetname (str): _description_
        writer (pd.ExcelWriter): _description_
    """
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
    
def count_representative_state(
    writer:pd.ExcelWriter,
    normalized_file:pd.DataFrame,
    dictionaries:dict ,
    table:dict = dict(), 
    period:dict = dict()
    ):

    for state_key in normalized_file['state_key'].unique():
        my_key = dictionaries["unique_states_types"][str(state_key)]
        df_2 = normalized_file.query(f'state_key == {state_key}')[['week_cont','variant_type']].groupby(['week_cont','variant_type']).size(
                        ).reset_index(name='patients')
        table[my_key]=dict()
        period[my_key]={'Alpha':[], 'Delta':[], 'Gamma':[], 'Omicron':[],'Otro':[],'Ninguno':[]}
        
        counting_representative(table, period, my_key, df_2, writer)
        
def count_representative_country(
    writer:pd.ExcelWriter,
    normalized_file:pd.DataFrame,
    table:dict = dict(), 
    period:dict = dict(),
    my_key:str = 'Pais'
):   
    df_2 = normalized_file[['week_cont','variant_type']].groupby(['week_cont','variant_type']).size(
                    ).reset_index(name='patients')
    table[my_key]=dict()
    period[my_key]={'Alpha':[], 'Delta':[], 'Gamma':[], 'Omicron':[],'Otro':[],'Ninguno':[]}
    counting_representative(table, period, my_key, df_2, writer)

def counting_representative(table, period, my_key, df_2, writer):
    # empty data/ asigna un valor 0  a todas las posiciones
        for week_num in range(df_2['week_cont'].max()): # for each week
            table[my_key][week_num+1]=dict()
            for variant_type in ['Alpha', 'Delta', 'Gamma', 'Omicron','Otro']:
                table[my_key][week_num+1][variant_type]=0
                
        # completa data
        for j in range(len(df_2)): # for each week
            table[my_key][df_2['week_cont'].iloc[j]][df_2['variant_type'].iloc[j]]=df_2['patients'].iloc[j]          
        
        # calculate values
        for week_num in range(df_2['week_cont'].max()): # for each week
            table[my_key][week_num+1]['week'] =week_num+1
            for variant_type in ['Alpha', 'Delta', 'Gamma', 'Omicron','Otro']:
                if not bool(week_num): 
                    if bool(table[my_key][week_num+1][variant_type]):
                        table[my_key][week_num+1][variant_type+'_Growth_rate']=100
                    else:
                        table[my_key][week_num+1][variant_type+'_Growth_rate']=0
            
                elif bool(table[my_key][week_num][variant_type]) and bool(table[my_key][week_num+1][variant_type]):
                    table[my_key][week_num+1][variant_type+'_Growth_rate'] = (table[my_key][week_num+1][variant_type]/table[my_key][week_num][variant_type])*100
                else:
                    if not bool(table[my_key][week_num][variant_type]) and bool(table[my_key][week_num+1][variant_type]):
                        table[my_key][week_num+1][variant_type+'_Growth_rate']=100
                    elif bool(table[my_key][week_num][variant_type]) and not bool(table[my_key][week_num+1][variant_type]):
                        table[my_key][week_num+1][variant_type+'_Growth_rate']=0
                    else:
                        table[my_key][week_num+1][variant_type+'_Growth_rate']=0
            total , max , aux = 0,0,0
            
            for variant_type in ['Alpha', 'Delta', 'Gamma', 'Omicron','Otro']:
            # create parameters
                max = table[my_key][week_num+1][variant_type]
                total= total + max
                if max > aux:
                    aux = max
                    table[my_key][week_num+1]['predominant']=variant_type
            
            if aux == 0:
                table[my_key][week_num+1]['predominant']='Ninguno'
                
            table[my_key][week_num+1]['total'] = total
            period[my_key][table[my_key][week_num+1]['predominant']].append(week_num+1)
            
            for variant_type in ['Alpha', 'Delta', 'Gamma', 'Omicron','Otro']:
                if table[my_key][week_num+1][variant_type] == 0 or total == 0:
                    table[my_key][week_num+1]['%_'+variant_type] = 0
                else:
                    table[my_key][week_num+1]['%_'+variant_type] = (table[my_key][week_num+1][variant_type]/total)*100
        autoadjust(pd.DataFrame(table[my_key]).transpose(),my_key,writer)

def fromisocalendar(y,w,d):
    date = datetime.datetime.strptime( "%04dW%02d-%d"%(y,w-1,d), "%YW%W-%w")
    return date.date().strftime("%d/%m/%Y")
