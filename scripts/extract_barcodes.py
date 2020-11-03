#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd


def find_amplicon_pos(df):
	"""
	Find amplicon position from a list of matches.
	Returns a tuple (start, stop).
	
	Will return the first occurence of a match on the + strand followed by a match on the minus strand.
	"""
	dfs = df.sort_values(by='start')
	start = None
	for index, row in dfs.iterrows():
		if row["strand"] == 'plus' and not start: 
			start = row['start']
		elif row['strand'] == 'minus' and start:
			end = row['start'] # blast reports end as the 3' position of the primer on the reverse strand
			return (start, end)
			

def main(blastfile, reportout):
	df = pd.read_csv(blastfile, 
						sep="\t", 
						header=0)

	# empty df to store amplicon informations
	dfout = pd.DataFrame(columns = ['seqid', 'taxid', 'start', 'end', 'length'])

	for s in set(df.seqid):
		sub_df = df.loc[df['seqid'] == s]
		pos = find_amplicon_pos(sub_df)
		if pos:
			df_s = pd.DataFrame([{'seqid' : s,
									'taxid' : list(sub_df['taxid'])[0],
									'start' : pos[0],
									'end' : pos[1],
									'length' : int(pos[1]) - int(pos[0])
									}])
			dfout = dfout.append(df_s, ignore_index=True)
	dfout.to_csv(reportout, sep='\t', header=False, index=False) 


if __name__ == '__main__':
	main(snakemake.input[0], snakemake.output[0])