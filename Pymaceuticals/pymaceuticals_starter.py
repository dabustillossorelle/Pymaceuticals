#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Dependencies and Setup
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import sem

# Hide warning messages in notebook
import warnings
warnings.filterwarnings('ignore')

# File to Load (Remember to Change These)
mouse_drug_data_to_load = "data/mouse_drug_data.csv"
clinical_trial_data_to_load = "data/clinicaltrial_data.csv"

# Read the Mouse and Drug Data and the Clinical Trial Data
mouse_drug_df=pd.read_csv(mouse_drug_data_to_load)
clinical_trial_df=pd.read_csv(clinical_trial_data_to_load)

# Combine the data into a single dataset
combined_pymaceuticals_df=pd.merge(mouse_drug_df, clinical_trial_df, how='outer', on='Mouse ID')

# Display the data table for preview
combined_pymaceuticals_df.head()


# ## Tumor Response to Treatment

# In[2]:


pymaceuticals_grpby_df=combined_pymaceuticals_df.groupby(['Drug', 'Timepoint'])
print(pymaceuticals_grpby_df)


# In[3]:


pymaceuticals_mean=pymaceuticals_grpby_df["Tumor Volume (mm3)"].mean()
pymaceuticals_mean


# In[4]:


pymaceuticals_means_df=pd.DataFrame(pymaceuticals_mean).reset_index()
pymaceuticals_means_df.head()


# In[5]:


pymaceuticals_sem=pymaceuticals_grpby_df["Tumor Volume (mm3)"].sem()
pymaceuticals_sem


# In[6]:


pymaceuticals_sem_df=pd.DataFrame(pymaceuticals_sem).reset_index()
pymaceuticals_sem_df.head()


# In[7]:


Restruct_Data=pymaceuticals_means_df.pivot(index='Timepoint', columns='Drug', values='Tumor Volume (mm3)')
Restruct_Data.head()


# In[8]:


fig, ax=plt.subplots()
Capomulin_error = pymaceuticals_sem_df.loc[pymaceuticals_sem_df["Drug"] == "Capomulin", "Tumor Volume (mm3)"]
Infubinol_error = pymaceuticals_sem_df.loc[pymaceuticals_sem_df["Drug"] == "Infubinol", "Tumor Volume (mm3)"]
Ketapril_error = pymaceuticals_sem_df.loc[pymaceuticals_sem_df["Drug"] == "Ketapril", "Tumor Volume (mm3)"]
Placebo_error = pymaceuticals_sem_df.loc[pymaceuticals_sem_df["Drug"] == "Placebo", "Tumor Volume (mm3)"]
Time = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
ax.errorbar(Time, Restruct_Data["Capomulin"] , yerr= Capomulin_error, label= "Capomulin", marker= "o", color="b")
ax.errorbar(Time, Restruct_Data["Infubinol"] , yerr= Infubinol_error, label= "Infubinol", marker= "^", color="r")
ax.errorbar(Time, Restruct_Data["Ketapril"] , yerr= Ketapril_error, label= "Ketapril", marker= "D", color="g")
ax.errorbar(Time, Restruct_Data["Placebo"] , yerr= Placebo_error , label= "Placebo", marker= "s", color="y")
plt.legend()
plt.title("Mean Tumor Response to Treatment ")
ax.set_xlabel("Timepoint")
ax.set_ylabel("Tumor Volume (mm3)")
plt.show()


# ## Metastatic Response to Treatment

# In[9]:


# Store the Mean Met. Site Data Grouped by Drug and Timepoint 

# Convert to DataFrame

# Preview DataFrame
pymaceuticals_grpby_df=combined_pymaceuticals_df.groupby(['Drug', 'Timepoint'])
print(pymaceuticals_grpby_df)


# In[10]:


pymaceuticals_mean2=pymaceuticals_grpby_df["Metastatic Sites"].mean()
pymaceuticals_mean2


# In[11]:


pymaceuticals_means2_df=pd.DataFrame(pymaceuticals_mean2).reset_index()
pymaceuticals_means2_df.head()


# In[12]:


pymaceuticals_sem2=pymaceuticals_grpby_df["Metastatic Sites"].sem()
pymaceuticals_sem2


# In[13]:


pymaceuticals_sem2_df=pd.DataFrame(pymaceuticals_sem2).reset_index()
pymaceuticals_sem2_df.head()


# In[14]:


# Minor Data Munging to Re-Format the Data Frames

# Preview that Reformatting worked
Restruct_Data2=pymaceuticals_means2_df.pivot(index='Timepoint', columns='Drug', values='Metastatic Sites')
Restruct_Data2.head()


# In[15]:


fig, ax=plt.subplots()
Capomulin_error = pymaceuticals_sem2_df.loc[pymaceuticals_sem2_df["Drug"] == "Capomulin", "Metastatic Sites"]
Infubinol_error = pymaceuticals_sem2_df.loc[pymaceuticals_sem2_df["Drug"] == "Infubinol", "Metastatic Sites"]
Ketapril_error = pymaceuticals_sem2_df.loc[pymaceuticals_sem2_df["Drug"] == "Ketapril", "Metastatic Sites"]
Placebo_error = pymaceuticals_sem2_df.loc[pymaceuticals_sem2_df["Drug"] == "Placebo", "Metastatic Sites"]
Time = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
ax.errorbar(Time, Restruct_Data2["Capomulin"] , yerr= Capomulin_error, label= "Capomulin", marker= "o", color="b")
ax.errorbar(Time, Restruct_Data2["Infubinol"] , yerr= Infubinol_error, label= "Infubinol", marker= "^", color="r")
ax.errorbar(Time, Restruct_Data2["Ketapril"] , yerr= Ketapril_error, label= "Ketapril", marker= "D", color="g")
ax.errorbar(Time, Restruct_Data2["Placebo"] , yerr= Placebo_error , label= "Placebo", marker= "s", color="y")
plt.legend()
plt.title("Mean Number of Metastatic Sites")
ax.set_xlabel("Timepoint")
ax.set_ylabel("Metastatic Sites")
plt.show()


# ## Survival Rates

# In[16]:


# Store the Count of Mice Grouped by Drug and Timepoint (W can pass any metric)

# Convert to DataFrame

# Preview DataFrame
pymaceuticals_grpby_df=combined_pymaceuticals_df.groupby(['Drug', 'Timepoint'])
print(pymaceuticals_grpby_df)


# In[17]:


pymaceuticals_counts=pymaceuticals_grpby_df["Mouse ID"].count()
pymaceuticals_counts


# In[18]:


pymaceuticals_counts_df=pd.DataFrame(pymaceuticals_counts).reset_index()
pymaceuticals_counts_df.head()


# In[19]:


pymaceuticals_counts_df=pymaceuticals_counts_df.rename(columns={"Mouse ID":"Mouse Counts"})
pymaceuticals_counts_df.head()


# In[20]:


# Minor Data Munging to Re-Format the Data Frames

# Preview the Data Frame
Restruct_Data3=pymaceuticals_counts_df.pivot(index='Timepoint', columns='Drug', values='Mouse Counts')
Restruct_Data3.head()


# In[22]:


# Generate the Plot (Accounting for percentages)

# Save the Figure

# Show the Figure
Time = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
plt.plot(Time, (Restruct_Data3["Capomulin"]/25)*100 ,label= "Capomulin", marker= "o", color="b")
plt.plot(Time, (Restruct_Data3["Infubinol"]/25)*100 , label= "Infubinol", marker= "^", color="r")
plt.plot(Time, (Restruct_Data3["Ketapril"]/25)*100 , label= "Ketapril", marker= "D", color="g")
plt.plot(Time, (Restruct_Data3["Placebo"]/25)*100 , label= "Placebo", marker= "s", color="y")
plt.legend()
plt.title("Survival Rates")
ax.set_xlabel("Timepoint")
ax.set_ylabel("Survival Rate (%)")
plt.show()


# ## Summary Bar Graph

# In[23]:


# Calculate the percent changes for each drug

# Display the data to confirm
Start_Tumor_Vol= 45
pymaceuticals_per_chng=((Restruct_Data.loc[45, :]-Start_Tumor_Vol)/Start_Tumor_Vol)*100
pymaceuticals_per_chng.head()


# In[25]:


pymaceuticals_perchng_df=pd.DataFrame(pymaceuticals_per_chng).reset_index()
pymaceuticals_perchng_df.head()


# In[28]:


pymaceuticals_perchng_df=pymaceuticals_perchng_df.rename(columns={45:"Percent Change"})
pymaceuticals_perchng_df.head()


# In[33]:


pymaceuticals_perchng_df=pd.DataFrame(pymaceuticals_perchng_df).round(decimals=2)
pymaceuticals_perchng_df.head()


# In[34]:


# Store all Relevant Percent Changes into a Tuple

Drugs_perchng=list(pymaceuticals_perchng_df[['Drug', 'Percent Change']].itertuples(index=False, name=None))
Drugs_perchng


# In[51]:


# Show the Figure
x_labels = [val[0] for val in Drugs_perchng]
y_labels = [val[1] for val in Drugs_perchng]
plt.figure(figsize=(12, 6))
ax=pd.Series(y_labels).plot(kind='bar')
ax.set_xticklabels(x_labels)

rects=ax.patches

for rect, label in zip(rects, y_labels):
    height=rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2, height-4, label, ha='center', va='bottom')

