#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from taxidTools.taxidTools import Taxdump

	
def main(taxid_file, parent, output, rankedlineage_dmp, nodes_dmp, failed):
	txd = Taxdump(rankedlineage_dmp, nodes_dmp)
	
	with open(taxid_file, "r") as fin:
		db_entries = set(fin.read().splitlines()[1:])
	
	with open(output, "w") as fout, open(failed, 'w') as ffailed:
		for taxid in db_entries:
			try:
				if txd.isDescendantOf(str(taxid).strip(), str(parent).strip()):
					fout.write(taxid + "\n")
				else:
					pass
			except KeyError:
				#print("WARNING: taxid %s missing from Taxonomy reference, it will be ignored" % taxid)
				ffailed.write(taxid + "\n")
				
if __name__ == '__main__':
	main(snakemake.input[0], 
			snakemake.params["taxid"], 
			snakemake.output["mask"], 
			snakemake.params["lineage"], 
			snakemake.params["nodes"],
			snakemake.output["failed"])