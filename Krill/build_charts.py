import subprocess
import os


def build_charts(path):
	input_file_path = os.path.join(path,'DBs_BGCs_with_Hits.tsv')
	output_dir = os.path.join(path,'Charts')
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	subprocess.run(
		f'Rscript Krill_visualizations.R {input_file_path} {output_dir}',
		shell=True,
		executable='/bin/bash',
		stdout=subprocess.DEVNULL,
    	stderr=subprocess.STDOUT
	)