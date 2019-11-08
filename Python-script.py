import pandas as pd
pd.options.display.float_format = '{:,.2f}'.format
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score,train_test_split
from sklearn import svm
from sklearn import preprocessing
%matplotlib inline

# Creating whole info Table ##
# All tag 2 uniq_rec (input_path)
all_tag_2_uniq_record = pd.read_csv('/scratch1/gao022/0_IPK_pangenome_dome/uni_record_2_seqs/hashes_2000_pars0/seqs_bino_ids.csv',
             names=['tag','uniq_record'])

# UAMT, uniquely Aligned Morex Tag (input_path)
UAMT = pd.read_table('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1497977_dome_UAMT',
             names=['tag', 'real_chr', 'real_pos'])

# Creating whole tabel df: including UAMT result in
df = pd.merge(all_tag_2_uniq_record, UAMT, on = 'tag', how = 'outer')

#############################################
##Sum machine learning feature (input_path)##
#############################################
names =['Sum', 'uniq_record', 'Total_Sig_Num', 'Best_SNP', 'Best_value', 'Sig_Num_chr1','Sig_Num_chr2','Sig_Num_chr3',
        'Sig_Num_chr4','Sig_Num_chr5','Sig_Num_chr6','Sig_Num_chr7','Sig_Num_chr8','Ratio2ed_odd','RatioM_odd','G_width','best_value',
       'sec_value', 'mid_value', 'Ratio2ed', 'RatioM']
uniq_record_Sum_result = pd.read_table('/scratch1/gao022/0_IPK_pangenome_dome/bino_result_2000_pars0_machine_feature/all_dome_feature_sum',
                                            names = names)
del uniq_record_Sum_result['Sum']

###################################
## Add TaxaCount into Sum result ##
###################################
# 1.creat TaxaCount (input_path)
uniq_record_2_Taxa = pd.read_csv('/scratch1/gao022/0_IPK_pangenome_dome/uni_record_2_seqs/hashes_2000_pars0/unique_records_ids.csv',names=['uniq_record', 'ids'])
uniq_record_2_Taxa['uniq_record'] = uniq_record_2_Taxa['uniq_record'].astype(str) + '_bino'
TaxaCount = []
for i in range(len(uniq_record_2_Taxa['ids'])):
    TaxaCount.append( uniq_record_2_Taxa['ids'][i].count('ERX') )
uniq_record_2_TaxaCount = uniq_record_2_Taxa[['uniq_record']]
uniq_record_2_TaxaCount['TaxaCount'] = TaxaCount

# 2.Merge TaxaCount to uniq_record_Sum_result
uniq_record_Sum_result = pd.merge(uniq_record_Sum_result, uniq_record_2_TaxaCount, on = 'uniq_record', how = 'outer')

###################################################
## Creating whole tabel df: including sum result ##
###################################################
df = pd.merge(df, uniq_record_Sum_result, on = 'uniq_record', how = 'outer')

###########################################################################################
## Creating whole tabel df: adding trans_pos; RP_distance and calculating other features ##
###########################################################################################

predicted_position = df["Best_SNP"].str.split("_", n=1, expand = True)

df["pre_chr"] = predicted_position[0]
df["pre_pos"] = predicted_position[1]

#replace str in df.columns
df["pre_chr"] = df["pre_chr"].str.replace('S','')
df["real_chr"] = df["real_chr"].str.replace('chr','')
df["real_chr"] = df["real_chr"].str.replace('H','')
df["real_chr"] = df["real_chr"].str.replace('Un','8')

# Creat Trans_Real_Pos
df['Trans_Real_Pos'] = df.real_chr.astype('float64') * 1000000000 + df.real_pos.astype('float64')

# Creat Trans_Pre_Pos
df['Trans_Pre_Pos'] = df.pre_chr.astype('float64') * 1000000000 + df.pre_pos.astype('float64')

# Creat PR_distance
df['RP_distance'] = abs(df.Trans_Pre_Pos - df.Trans_Real_Pos)

# factorize RP_distance
df['Y'] =  np.where(df['RP_distance'] <= 1000, 1, 
                   (np.where(df['RP_distance'] <= 100000, 2, 
                            np.where(df['RP_distance'] <= 10000000, 3, 4
                            ))))

# calculate the Sig_Num on best_chr AND ad it back to df
df_sig_num =  df[['pre_chr', 'Sig_Num_chr1','Sig_Num_chr2', 'Sig_Num_chr3','Sig_Num_chr4', 'Sig_Num_chr5', 
                  'Sig_Num_chr6', 'Sig_Num_chr7','Sig_Num_chr8']]
list = df_sig_num['pre_chr'].astype('int64')
df['Sig_Num_BestChr'] = df_sig_num.values[range(len(list)),list]

# creat (Sig_Num_BestChr / Total_Sig_Num)%
df['Sig_Num_BC_percent'] = df['Sig_Num_BestChr'] / df['Total_Sig_Num']

# Add in "RatioGSN" = "Gwidth / Sig_Num_BestChr" ##Not important
df['RatioGSN'] = df['G_width'] / df['Sig_Num_BestChr']
df['RootGSN'] = df['G_width'] ** (1 / df['Sig_Num_BestChr'])


######################################################
###### Creating machine_learing dataset from df ######
######################################################
df_ML = df[['tag',
            'uniq_record',
            'Best_value',
            'RatioM_odd', 
            'Ratio2ed_odd',
            'G_width',
            'Total_Sig_Num', 
            'Sig_Num_BestChr', 
            'Sig_Num_BC_percent',
            'RatioGSN',
            'RootGSN',
            'TaxaCount',
            'Y',
            'RP_distance']].dropna().drop_duplicates().reset_index(drop = True)

# UAMT description
print("UAMT number:",len(df_ML))
print("UAMT number of class 1 (<1Kb):",len(df_ML[df_ML['Y']==1]))
print("UAMT number of class 2 (1Kb~100Kb):",len(df_ML[df_ML['Y']==2]))
print("UAMT number of class 3 (100Kb~10Mb):",len(df_ML[df_ML['Y']==3]))

#Add random col
df_ML['random'] = np.random.uniform(0, 1, len(df_ML))



#########################################################################
# Train model for tuning parameter and plot Mulitple classes ROC curves #
#########################################################################
from itertools import cycle
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier

# Get data ready (take n% of df_ML)
X = df_ML[df_ML['random']>=0.8].iloc[:,2:12]
y = df_ML[df_ML['random']>=0.8].iloc[:,12]

# Binarize the output
y = label_binarize(y, classes=[1, 2, 3, 4])
n_classes = y.shape[1]

# Shuffle and split training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.25, random_state=0)

# Learn to predict each class against the other
classifier = OneVsRestClassifier(RandomForestClassifier(random_state=0, n_estimators= 90, min_samples_leaf= 2, 
                                                        n_jobs = -1, max_features = 6))
y_score = classifier.fit(X_train, y_train).predict_proba(X_test)

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

# Plot all ROC curves
plt.figure(figsize=(10,6.6))
plt.plot(fpr["micro"], tpr["micro"], label='micro-average (AUC = {0:0.2f})' ''.format(roc_auc["micro"]),
         color='k', linestyle='-', linewidth=1.5)

colors = cycle(['r', 'm', 'darkorange', 'b'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=1.5,label='Class {0} (AUC = {1:0.2f})' ''.format(i+1, roc_auc[i]))

plt.plot([0, 1], [0, 1], '--', lw=2, color ='dimgray')
plt.xlim([-0.02, 1.02])
plt.ylim([-0.02, 1.02])
plt.xlabel('False Positive Rate', fontsize = 20, labelpad=0)
plt.ylabel('True Positive Rate', fontsize = 20, labelpad=0)
plt.title('ROC Curves (domesticated dataset)',fontsize = 22, y = 1.04)
plt.legend(loc="lower right", prop={'size': 17}, frameon=False)

ax = plt.gca()
ax.tick_params(direction = 'out', labelsize = 20, width = 1, top=False, right=False)

plt.savefig('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/roc_wild.png', format = "png")

# train classifier_final using whole df_ML
clf_final = RandomForestClassifier(random_state=0, n_estimators= 90, min_samples_leaf= 2, n_jobs = -1, max_features = 6)

# Whole training set
train_all = df_ML.iloc[:,2:12]

# Training the classifier
clf_final.fit(train_all, df_ML['Y'])

# Plot feature ranking
features = df_ML.columns[2:12]
importances = clf_final.feature_importances_
y_pos=np.arange(len(importances))

plt.figure(figsize=(11,6.6))

plt.barh(y_pos, importances, align = 'center', height =0.7, lw = 0, color = 'b')
plt.ylim([-0.7,9.7])
ax = plt.gca()
ax.grid(axis = 'x', lw =1, ls = ':', c = 'k')
ax.set_yticks(y_pos)
ax.set_yticklabels(features)
ax.tick_params(direction = 'out', labelsize = 15, width = 1, top=False, right=False)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Importance', fontdict={'size':20})
ax.set_title('Importance of features', fontdict={'size':22},y=1.04)
plt.tight_layout()
plt.savefig('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/feature_ranking.png',format='png')

################################################################################################
# Applying the trained Classifier to the test and tidy up and visualize final predicted result #
################################################################################################

preds = clf_final.predict(df[features])

# Final result df_F
df_F = df[['tag','uniq_record','Trans_Pre_Pos']]
df_F['class'] = preds

# Merge sam position wiht df_F class <=3
df_F1= df_F[df_F['class'] <=3].reset_index(drop=True)

# sam_position info
sam_pos = pd.read_table('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/14M_dome_tag_2_Morex.sam_pos',
                        names = ['tag','sam_chr','sam_pos','multi_hit'])

# Meger df_F1 and sam_pos to df_final_correct df_FC
df_FC = pd.merge(df_F1, sam_pos, on ='tag', how = 'inner')

#replace str in df.columns
df_FC["sam_chr"] = df_FC["sam_chr"].str.replace('chr','')
df_FC["sam_chr"] = df_FC["sam_chr"].str.replace('H','')
df_FC["sam_chr"] = df_FC["sam_chr"].str.replace('Un','8')
df_FC['Trans_Sam_Pos'] = df_FC.sam_chr.astype('float64') * 1000000000 + df_FC.sam_pos.astype('float64')

# Sorting PAV-I (absent on Morex)
print("PAV-I total:", len(df_FC[df_FC['sam_pos'] == 0]))
print("PAV-I class 1:", len(df_FC[(df_FC['sam_pos'] == 0) & (df_FC['class'] == 1)]))
print("PAV-I class 2:", len(df_FC[(df_FC['sam_pos'] == 0) & (df_FC['class'] == 2)]))
print("PAV-I class 3:", len(df_FC[(df_FC['sam_pos'] == 0) & (df_FC['class'] == 3)]))

# Output PAV-I data for distribution plot
n, bins, patches = plt.hist(df_FC[df_FC['sam_pos'] == 0]['Trans_Pre_Pos'], 1000)
hist = pd.DataFrame({"bins": bins})
hist['frequency'] = np.append([0],n)
hist.to_csv('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/PVA_I_hist.csv',index =False)

# Sorting PAV-II (Over100Mb)
cond = abs(df_FC['Trans_Pre_Pos'] - df_FC['Trans_Sam_Pos']) >= 100000000
print("PAV-II (Over100Mb) total:", len(df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) ]))
print("PAV-II (Over100Mb) class 1:", len(df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) & (df_FC['class'] == 1)]))
print("PAV-II (Over100Mb) class 2:", len(df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) & (df_FC['class'] == 2)]))
print("PAV-II (Over100Mb) class 3:", len(df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) & (df_FC['class'] == 3)]))

# Reduce PAVII result for plotting
PAVII = df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) ][['tag','uniq_record','class','Trans_Pre_Pos','Trans_Sam_Pos']]

PAVII['Occ_No'] =  PAVII.groupby("uniq_record").cumcount()+1

PAVII_1 = PAVII[PAVII['Occ_No']==1]
PAVII_1['rand'] = np.random.uniform(0, 1, len(PAVII_1))

PAVII_2 = PAVII_1[PAVII_1['rand']<0.04]


# Output PAV-II data for circos plot and validation
PAVII_2.to_csv('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/PAV_II.csv', index = False)


# pre_pos-sam_pos matched tag: to show the capability of ML model
cond_1 = abs(df_FC['Trans_Pre_Pos'] - df_FC['Trans_Sam_Pos']) < 100000000
x = df_FC[cond_1]['Trans_Sam_Pos']
y = df_FC[cond_1]['Trans_Sam_Pos']

plt.figure(figsize=(9,9))
plt.scatter(x-900000000, y-900000000, s = 1, c="K")
plt.grid(axis = 'both', lw =1, ls = ':')

plt.xlim([1, 8e9])
plt.ylim([1, 8e9])
plt.xlabel('Physical position', fontsize = 20, labelpad=0)
plt.ylabel('GWAS-ML predicted position', fontsize = 20, labelpad=0)
plt.title('Non-PAV tags (class 1-3) filtered by ML model',fontsize = 22, y = 1.04)
ax = plt.gca()
ax.tick_params(direction = 'out', labelsize = 20, width = 1, top=False, right=False)
plt.savefig('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/Non-PAV.png', format = "png")

# Creat pie chart for class 1-3 tags filtered by ML-model
labels = 'Non-PAV', 'PAV-I Class2', 'PAV-I Class1', 'PAV-I Class3', 'PAV-II Class1', 'PAV-II Class2', 'PAV-II Class3'
sizes = [921709, 21668, 102615, 122397, 61720, 14199, 136157 ]
explode = (0,0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
colors = ['gray', 'red','red', 'red', 'b', 'b', 'b']

plt.figure(figsize=(20,10))
plt.pie(sizes, explode=explode,labels=labels, colors=colors, autopct='%1.1f%%',shadow=False, startangle=57.5,
       textprops = {'fontsize': 30, 'color': 'k'},pctdistance=1.17, labeldistance=1.38,radius=0.2,
       wedgeprops = {'linewidth': 1,"edgecolor":"k"})
# plt.title('444k Class1-3 tags filtered by ML model',fontsize = 22, y = 1.04)
plt.axis('equal')
plt.savefig('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/pie_chart.png', format='png')

# Creat pie chart for class 1-3 tags filtered by ML-model
labels = 'Non-PAV', 'PAV-I Class2', 'PAV-I Class1', 'PAV-I Class3', 'PAV-II Class1', 'PAV-II Class2', 'PAV-II Class3'
sizes = [921709, 21668, 102615, 122397, 61720, 14199, 136157 ]
explode = (0,0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
colors = ['gray', 'red','red', 'red', 'b', 'b', 'b']

plt.figure(figsize=(10,10))
plt.pie(sizes, explode=explode, colors=colors,shadow=False, startangle=57.5,
       textprops = {'fontsize': 30, 'color': 'k'},pctdistance=1.17, labeldistance=1.38,radius=0.2,
       wedgeprops = {'linewidth': 0,"edgecolor":"k"})
# plt.title('444k Class1-3 tags filtered by ML model',fontsize = 22, y = 1.04)
plt.axis('equal')
plt.savefig('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/1_plot/pie_chart_2.png', format='png')


# Output PAVI tag
df_FC[(df_FC['sam_pos'] == 0)][['tag','Trans_Pre_Pos']].to_csv('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/2_anchor_pantranscriptome/PAVI.csv', 
                                               index = False)

# Output PAVII tag
df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) & (df_FC['class'] <= 3)][['tag','Trans_Pre_Pos']].to_csv('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/2_anchor_pantranscriptome/PAVII.csv', 
                                               index = False)

tag_2_SNP = df[['tag','Best_SNP']]

# Output PAVI tag
PAVI_result = pd.merge(df_FC[(df_FC['sam_pos'] == 0)][['tag','Trans_Pre_Pos']], tag_2_SNP, on = 'tag', how = 'inner')

PAVI_result.to_csv('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/2_anchor_pantranscriptome/PAVI.csv',
                   index = False)

# Output PAVII tag
PAVII_result = pd.merge(df_FC[cond & (df_FC['Trans_Sam_Pos'] != 0) & (df_FC['class'] <= 3)][['tag','Trans_Pre_Pos']], tag_2_SNP, on = 'tag', how = 'inner')

PAVII_result.to_csv('/scratch1/gao022/0_IPK_pangenome_dome/1_MACHINE_LEARNING_pars0/2_anchor_pantranscriptome/PAVII.csv', 
                    index = False)

plt.figure(figsize=(20,30))
_ = plt.hist(PAVII_result['Trans_Pre_Pos'],1000)