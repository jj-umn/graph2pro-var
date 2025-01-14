script=$(basename $(readlink -nf $0))
parfile=$1
withhead=$2
if [ -z $parfile ]; then
  echo "Error: parameter file not given"
  exit 1
fi
if [ ! -f $parfile ]; then
   echo "Error: $parfile NOT found"
   exit 1
fi
if [ -z $withhead ]; then
  withhead=true
fi
#get parameters from input parfile
IFS=$'\n'; set -f; par=($(<$parfile))
n=${#par[@]}
fdr=0.01
for apar in "${par[@]}"; do
   IFS='=' read -r -a item <<< "$apar"
   case "${item[0]}" in
     id) exp_d=${item[1]};;
     ms) mgf=${item[1]};;
     fdr) fdr=${item[1]};;
   esac
done
if [ -z $exp_d ]; then
  echo "Error: exp_d not found"
  exit 1
fi
if [ -z $mgf ]; then
  echo "Error: mgf not found"
  exit 1
fi

spec_count=`grep "BEGIN IONS" $mgf -c`
if [ -z $exp_d ]; then
   echo "error: prefix of the filenames including the folder is needed"
   exit 1
fi
contig_info=`getUniquePeptides_files.py -i ${exp_d}.fgs.tsv.$fdr.tsv -b`
other_info="NA NA NA NA NA NA"
if [ -f ${exp_d}.uncordant.tsv.$fdr.tsv ]; then
   other_info=`getUniquePeptides_files.py -i ${exp_d}.graph2pro.tsv.$fdr.tsv -v ${exp_d}.var2pep.tsv.$fdr.tsv -b -o ${exp_d}.final-peptide.txt`
fi
if $withhead; then
   echo "spectra-count contig-peptide-all contig-peptide-unique graph2pro-peptide-all graph2pro-peptide-unique var2peponly-peptide-all var2peponly-peptide-unique graph2pro&var2pep-peptide-all graph2pro&var2pep-peptide-unique"
fi

echo "$exp_d $spec_count $contig_info $other_info"
