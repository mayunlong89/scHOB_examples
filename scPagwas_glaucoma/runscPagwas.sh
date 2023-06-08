#/bin/sh
sh_dir=$(cd $(dirname $0); cd ../; pwd)

show_usage="args:[-d|--dir]\
                 [-f|--vcf]\
                 [-n|--num]\
                 [-t|--trait]"
dir=""                 
trait=""
seuDir=""
gwasDes=""
GETOPT_ARGS=`getopt -o d:f:n:t:s:g: -al dir:,vcf:,num:,trait:,seuDir:,gwasDes: -- "$@"`
while [ -n "$1" ]
do
  case "$1" in 
    -d|--dir) dir=$2; shift 2 ;;
    -g|--gwasDes) gwasDes=$2; shift 2 ;;
    -t|--trait) trait=$2; shift 2 ;;
    -s|--seuDir) seuDir=$2; shift 2 ;;
    --) break ;;
    *)echo $show_usage; break ;;
  esac
done

cd $dir
results_dir="./results/scPagwas_"$trait"/"
if [ ! -d $results_dir ]; then
  mkdir -pv $results_dir
fi

if [ ! -e $results_dir$trait"scPagwas.rds" ]; then
    echo 'run scPagwas', $trait
    Rscript $sh_dir/lib/scPagwas.R --dir $dir --seuDir $seuDir --gwas $gwasDes --trait $trait
  else
    echo $trait' file has existed, continue!'
fi

#End
