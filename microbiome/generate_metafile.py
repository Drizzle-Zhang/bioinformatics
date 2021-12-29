import pandas as pd
import re


file_in = '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
df_meta = pd.read_csv(file_in, sep='\t', encoding='gbk')
df_meta = df_meta.iloc[:, :3]
df_meta.columns = ['Sample', 'origin_group', 'Description']
df_meta = df_meta.dropna()
df_meta['Time'] = df_meta['origin_group'].apply(lambda x: x.split(' ')[1][0])
df_meta['Group'] = df_meta['origin_group'].apply(lambda x: x.split(' ')[0])
df_meta = df_meta.replace({'Group': 'Group'}, value='Treat')
df_meta['Dose'] = df_meta['origin_group'].apply(lambda x: x.split(' ')[1][1])
df_meta = df_meta.replace({'Dose': '4'}, value='0')
df_meta['Gender'] = df_meta['origin_group'].apply(
    lambda x: x.split(' ')[1][3:-1])
pattern_id = re.compile(r'[0-9]+')
df_meta['ID_Mouse'] = df_meta['Description'].apply(
    lambda x: pattern_id.match(x).group())
file_out = '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df_meta.to_csv(file_out, sep='\t', index=None)
