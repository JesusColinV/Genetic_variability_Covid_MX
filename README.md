Methodology

For the development of this code, the following Python modules were used: Pandas, Warnings, Os, Concurrent, Plotly, Datetime. Starting from the information of the Gisaid database and some pre-built json files with base knowledge to classifier information, we proceeded to complete and normalize the main table: 

Complete

A) Cleaning the collection date: ignoring string length greater than 7 and changing the string data type to datetime. 
B) Generate the time characteristics: reading and assigning the corresponding year, calendar week, continuous week, corresponding state of each record
a.	The state is read from the location feature by doing a division by "/" for the third position.

Normalize

A) Normalized/standardized lineage: using the json variant_types file, where each clade was assigned to a variant studied.
B) Age normalized: needed due to an error in the database, where some patient age records were assigned to the Gender characteristic and to fix the problem the generalized logic was written in the json file age_unification, this returns both the age variable and the age group.
C) Patient group status: using the json file patient_status to classify each patient into 3 categories according to the severity of the patient's recorded condition in the patient status function.
D) Normalized geographic status: using the states_types json file to unify all record names in the status function, which returns the state_key and region_key variables.

Finally, the main table is filtered according to the analysis, filtering where under 18 and equal and over 18 and also the original normalized table. 

From the maximum and minimum date characteristic, the support and end recorded for the analysis record are obtained and each of the three tables mentioned above are counted.

The second part of the code is to count and measure the percentage of each mutation variant and display the 10 most frequent by counting and grouping each variant.

A) Count the variant: launch the get_mutation function on the Counter class, the value of the unique amino acid for each of the Spike, E, N, N, M, and NSP type proteins is obtained.
a.	The Search_amino function allows to return each of the amino acid values from the "AA Substitutions" function of the main table. 
b.	The Position function allows to return the position value of the mutation.
c.	The Table_protein function, given the current dictionary as "val", updates the count of each amino acid and has the same function for table_nsp that applies to NSP mutations. 
d.	V_C function counts for each unique mutation. 
e.	P_P function returns the portion of the protein to which it belongs means it is a classifier. 
f.	Dic2df converts the dictionary variable into a data frame to be displayed as a table. 

B) Clustering is performed by repeating the process for each of the 3 tables, < 18, > 18 and normalized into a single table with new name.

C) The percentage count is performed by dividing the counted result by the total of each record and multiplied by 100, and for nan values it is replaced by 0 using the "fillna" function.

Finally we can identify the segment and part of the protein by grouping the "full" function.

The table to identify the most representative variant by state first by preprocessing using the count_representative_state function that obtains parameters such as growth for each week and the percentage of each variant by state, and second using the get_states function where the unique_state_types feature is grouped and when the value is greater than 50% it is considered as the most representative variant for these weeks, otherwise the value is None, it means that there is no representative variant in these states and these week. 
