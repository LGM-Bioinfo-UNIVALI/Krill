import subprocess
import os


def build_charts(path):
	project_dir = os.path.dirname(os.path.abspath(__file__))
	input_file_path = os.path.join(path,'DBsReportOutput/DBs_BGCs_with_Hits.tsv')
	output_dir = os.path.join(path,'DBsReportOutput/Charts')
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	subprocess.run(
		f'Rscript {project_dir}/Krill_visualizations.R {input_file_path} {output_dir}',
		shell=True,
		executable='/bin/bash',
		stdout=subprocess.DEVNULL,
    	stderr=subprocess.STDOUT
	)


if __name__ == '__main__':
	build_charts("/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/NCBI_PROJECTS/AtlanticoSul")