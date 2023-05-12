import subprocess, os, tempfile
import pandas as pd

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
    df.to_csv(os.path.join(path,'DNABases_and_ORFs_count.tsv'),sep='\t',index=False)