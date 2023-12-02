import subprocess, os, tempfile
import pandas as pd


def df2xlsx(path, sheet_name, df):
    df = df.copy().reset_index()
    writer = pd.ExcelWriter(path, engine='xlsxwriter')

    df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)
    worksheet = writer.sheets[sheet_name]

    (max_row, max_col) = df.shape
    column_settings = [{'header': column} for column in df.columns]
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
    worksheet.set_column(0, max_col - 1, 15)

    writer._save()


def run(path, ext, threads):
    os.chdir(path)

    orffinder = os.path.join(
                            os.path.dirname(__file__),
                            'ORF_Finder','ORFfinder'
                            )

    bash="""cp {orffinder} {path} && \
    chmod u+x ORFfinder && \
    echo "FILE,NT,ORFS" && \
    mkdir -p ORFs && \
    parallel -j{threads} 'count=$(grep -v ">" {{}} | tr -d '\\n' | wc -c) && file=$(echo "{{.}}") && ./ORFfinder -in {{}} -g 11 -s 1 -n true -out ORFs/"{{}}" && ORFS=$(grep -c ">" ORFs/"{{}}") && echo -n "$file," && echo -n "$count," && echo $ORFS' ::: *.{ext} && \
    rm ORFfinder""".format(orffinder=orffinder,path=path,ext=ext,threads=threads)

    tmp = tempfile.NamedTemporaryFile(suffix='.sh',mode='w+',delete=False)
    tmp.write(bash)
    tmp.seek(0)
    tmp.close()
    
    os.chmod(tmp.name,0o777)

    output, err = subprocess.Popen([tmp.name],
                           shell=True,
                           executable='/bin/bash',
                           universal_newlines=True,
                           text=True,
                           stdout=subprocess.PIPE).communicate()

    os.remove(tmp.name)

    unformatted_tsv = [x.split(',') for x in output.splitlines()]
    df = pd.DataFrame(unformatted_tsv[1:], columns=unformatted_tsv[0])
    df.to_csv(os.path.join(path,'ReportOutput/DNABases_and_ORFs_count.tsv'),sep='\t',index=False)
    df2xlsx(os.path.join(path,'ReportOutput/DNABases_and_ORFs_count.xlsx'), 'DNABases and ORFs count', df)