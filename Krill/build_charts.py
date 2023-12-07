import subprocess
import os


def build_charts(path):
	project_dir = os.path.dirname(os.path.abspath(__file__))
	input_file_path = os.path.join(path,'DBsReportOutput/DBs_BGCs_with_Hits.tsv')
	output_dir_hybrids = os.path.join(path,'DBsReportOutput/Charts/WithHybrids')
	output_dir_no_hybrids = os.path.join(path,'DBsReportOutput/Charts/WithoutHybrids')
	
	if not os.path.exists(output_dir_hybrids):
		os.makedirs(output_dir_hybrids)

	if not os.path.exists(output_dir_no_hybrids):
		os.makedirs(output_dir_no_hybrids)

	subprocess.run(
		f'Rscript {project_dir}/Krill_visualizations_with_hybrids.R {input_file_path} {output_dir_hybrids}',
		shell=True,
		executable='/bin/bash',
		stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT
	)

	subprocess.run(
		f'Rscript {project_dir}/Krill_visualizations_without_hybrids.R {input_file_path} {output_dir_no_hybrids}',
		shell=True,
		executable='/bin/bash',
		stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT
	)


if __name__ == '__main__':
	build_charts("/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/NCBI_PROJECTS/AtlanticoSul")