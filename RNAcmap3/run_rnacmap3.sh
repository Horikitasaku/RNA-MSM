#!/bin/bash

# set -euxo pipefail

if [ $# -eq 0 ]; then
    echo "No parameters given, try \"$0 -h\" for usage introductions"
    exit 0
fi

PARSED_OPTIONS=$(getopt -n "$0"  -o hi:n:b:c:d:rg \
    --long "help,reuse_cm1,ignore_vol_err,dca_method:,blast_db:,infernal_db:,ncpu:,input:" -- "$@")

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ];
then
  exit 1
fi
 
# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

while true; do
    case "$1" in
        -h|--help)
            echo "usage: $0 -h|--help
usage: $0 [-i|--input INPUT] [-n|--ncpu NCPU ] [-b|--blast_db BLAST_DB] [-c|--infernal_db INFERNAL_DB]
          [-d|--dca_method DCA_METHOD] [-r|--reuse_cm1] [-g|--ignore_vol_err]
   
Parameters:
    INPUT: input file, a nucleotide sequence in fasta format, default: ./seq.fasta
           outputs and intermediate files will be written to the directory where the input file locates

    NCPU: number of cpus to use, optional, default: 24. NCPU does not support cross-node tasks 

    BLAST_DB: absolute path to the location of .nal file (without .nal suffix) of blast database,
              accepts a single path or a quoted, whitespace-separated list of paths,
              optional, default: $PWD/db/blast/MARS/MARS
              
    INFERNAL_DB: absolute path to the directory of fasta database for infernal,
                 accepts a single path, optional, default: $PWD/db/infernal/MARS
                 fasta files in the directory must end with '.fasta'

    DCA_METHOD: method to run dca analysis, accept: gremlin, mfdca, plmc, plmdca
                optional, default: gremlin

Options:
    -r  --reuse_cm1
        By default, RNAcmap3 always regenerates cm_msa1 file from the BLAST result.
        Use this option to use existing cm_msa1 file, take effect at rerun only

    -g  --ignore_vol_err
        By default, RNAcmap3 exits immediately after an infernal search Error.
        Use this option to ignore infernal search error on individual database volumes and 
        continue search on the next volume.
    "
            shift
            exit 0;;
        -i|--input)
            input="$2"
            shift 2;;
        -n|--ncpu)
            ncpu="$2"
            shift 2;;
        -b|--blast_db)
            blast_db="$2"
            shift 2;;
        -c|--infernal_db)
            infernal_db="$2"
            shift 2;;
        -d|--dca_method)
            dca="$2"
            shift 2;;
        -r|--reuse_cm1)
            reuse_cm1=T
            shift 1;;
        -g|--ignore_vol_err)
            ignore_vol_err=T
            shift 1;;
        --)
            shift
            break;;
        *)
            echo "Unknown option $1"
            exit 1
    esac
done

start=`date +%s`

sel_dca=${dca:-gremlin}

echo "will use dca method $sel_dca"

input0=${input:-./seq.fasta}
input="$(cd "$(dirname "$input0")"; pwd)/$(basename "$input0")"
input_dir=$(dirname $input)
seq_id=$(basename $(basename $input) | cut -d. -f1)
program_dir=$(dirname $(readlink -f $0))

path_blastn_database=(${blast_db:-$PWD/db/blast/MARS/MARS})
inf_db=${infernal_db:-$PWD/db/infernal/MARS}
path_infernal_database=($(find $inf_db -name "*.fasta"|sort))
np=${ncpu:-24}  # number of threads
use_cm1=${reuse_cm1:-F}
ignore_volError=${ignore_vol_err:-F}

mkdir -p $input_dir/${seq_id}_features #&& mkdir -p $input_dir/${seq_id}_outputs
echo ">"$seq_id > $input_dir/${seq_id}_features/$seq_id.fasta
awk -i inplace '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $input 
tail -n1 $input >> $input_dir/${seq_id}_features/$seq_id.fasta
echo "" >> $input_dir/${seq_id}_features/$seq_id.fasta

feature_dir=$input_dir/${seq_id}_features


###### check if aligned homologous sequences file already exists ############
if [ -f $feature_dir/$seq_id.a2m_msa2 ];	then
        echo ""
        echo "==========================================================================="
        echo "    MSA file $feature_dir/$seq_id.a2m_msa2 from Infernal Pipeline already  "
        echo "    exists for query sequence $feature_dir/$seq_id.fasta.                  "
        echo "                                                                           "
        echo "    Delete existing $feature_dir/$seq_id.a2m_msa2 if want to generate new  "
        echo "    alignment file                                                         "
        echo "==========================================================================="
    	echo ""
else

    #################### check if blastn alignment file ready exists ######################
    if [ -f $feature_dir/$seq_id.bla_msa1 ];       then
	    echo ""
	    echo "============================================================================"
	    echo "    MSA-1 file $feature_dir/$seq_id.bla_msa1 from Infernal Pipeline already "
	    echo "    exists for query sequence $feature_dir/$seq_id.fasta.                   "
	    echo "                                                                            "
	    echo "    Delete existing $feature_dir/$seq_id.bla_msa1 if want to generate new   "
	    echo "    alignment file.                                                         "
	    echo "============================================================================"
		echo ""
    else
        echo ""
        echo "==========================================================================================================================="
        echo "      Running BLASTN for first round of homologous sequence search for query sequence $feature_dir/$seq_id.fasta.          "
        echo "      May take 5 mins to few hours depending on sequence length and no. of homologous sequences in database.               "
        echo "==========================================================================================================================="
        echo ""
        blastn -db "${path_blastn_database[*]}" -query $feature_dir/$seq_id.fasta -out $feature_dir/$seq_id.bla_msa1 -evalue 0.001 -num_descriptions 1 -num_threads $np -line_length 1000 -num_alignments 50000
    fi
			
	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      First round of MSA-1 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "=============================================================================="
        echo "        Error occured while formatting the nt database.                       "
        echo "=============================================================================="
        echo ""
        exit 1
    fi

    if [ "$use_cm1" == T ] && [ -s $feature_dir/$seq_id.cm_msa1 ]; then
        echo ""
        echo "==========================================================================="
        echo "    Use existing $feature_dir/$seq_id.cm_msa1"
        echo "==========================================================================="
    	echo ""
    else
	    ######## reformat the output ################
        echo ""
        echo "=================================================================================================="
        echo "         Converting $feature_dir/$seq_id.bla_msa1 from BLASTN to $feature_dir/$seq_id.sto_msa1.   "
        echo "=================================================================================================="
        echo ""
	    $program_dir/utils/parse_blastn_local.pl $feature_dir/$seq_id.bla_msa1 $feature_dir/$seq_id.fasta $feature_dir/$seq_id.aln_msa1
	    sed -i 's/\s.*$//' $feature_dir/$seq_id.aln_msa1
	    sed -i "s/$seq_id/$seq_id E=0.0/g" $feature_dir/$seq_id.aln_msa1

	    $program_dir/utils/seqkit rmdup -n $feature_dir/$seq_id.aln_msa1 > $feature_dir/temp.aln_msa1
	    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/temp.aln_msa1 | sed '/^$/d' > $feature_dir/$seq_id.aln_msa1

	    $program_dir/utils/reformat.pl fas sto $feature_dir/$seq_id.aln_msa1 $feature_dir/$seq_id.sto_msa1
        ## Klark note: variate duplicate seqnames in sto alignments
        mv -f $feature_dir/$seq_id.sto_msa1 $feature_dir/$seq_id.sto_msa1.tmp
        python $program_dir/utils/unique_seqname.py $feature_dir/$seq_id.sto_msa1.tmp >$feature_dir/$seq_id.sto_msa1 && rm $feature_dir/$seq_id.sto_msa1.tmp
        ## Klark done

	    if [ $? -eq 0 ]; then
	        echo ""
	        echo "=========================================="
            echo "      Converison completed successfully.  "
	        echo "=========================================="
	        echo ""
	    else
            echo ""
            echo "==================================================================================================="
            echo "   Error occured while Converting $feature_dir/$seq_id.bla_msa1 to $feature_dir/$seq_id.sto_msa1   "
            echo "===================================================================================================="
            echo ""
            exit 1
        fi

	    ######## predict secondary structure from RNAfold ################
        echo ""
        echo "==============================================================================================================================="
        echo "       Predicting Consensus Secondary Structure (CSS) of query sequence $feature_dir/$seq_id.fasta using RNAfold predictor.   "
        echo "==============================================================================================================================="
        echo ""

	    RNAfold $feature_dir/$seq_id.fasta | awk '{print $1}' | tail -n +3 > $feature_dir/$seq_id.db

	    ################ reformat ss with according to gaps in reference sequence of .sto file from blastn ################
	    for i in `awk '{print $2}' $feature_dir/$seq_id.sto_msa1 | head -n5 | tail -n1 | grep -b -o - | sed 's/..$//'`; do sed -i "s/./&-/$i" $feature_dir/$seq_id.db; done

	    #########  add reformated ss from last step to .sto file of blastn ##############
	    head -n -1 $feature_dir/$seq_id.sto_msa1 > $feature_dir/temp.sto
	    echo "#=GC SS_cons                     "`cat $feature_dir/$seq_id.db` > $feature_dir/temp.txt
	    cat $feature_dir/temp.sto $feature_dir/temp.txt > $feature_dir/$seq_id.sto_msa1
	    echo "//" >> $feature_dir/$seq_id.sto_msa1

	    if [ $? -eq 0 ]; then
	        echo ""
	        echo "=================================================================="
            echo "      Consensus Secondary Structure (CSS) generated successfully. "
	        echo "=================================================================="
	        echo ""
	    else
            echo ""
            echo "================================================================================"
            echo "             Error occured while generating structure from RNAfold.             "
            echo "================================================================================"
            echo ""
            exit 1
        fi

	    ######## run infernal ################
        echo ""
        echo "========================================================================================================================"
        echo "      Building Covariance Model from BLASTN alignment (with SS from RNAfold) from $feature_dir/$seq_id.sto_msa1 file.  "
        echo "========================================================================================================================"
        echo ""
	    cmbuild --hand -F $feature_dir/$seq_id.cm_msa1 $feature_dir/$seq_id.sto_msa1

	    if [ $? -eq 0 ]; then
	        echo ""
	        echo "================================================================================="
            echo "    Covariance Model (CM) built successfully from $feature_dir/$seq_id.sto_msa1. "
	        echo "================================================================================="
	        echo ""
	    else
            echo ""
            echo "==============================================================================================="
            echo "     Error occured while building Covariance Model (CM) from cmbuild.           "
            echo "==============================================================================================="
            echo ""
            exit 1
        fi

        echo ""
        echo "======================================================================="
        echo "       Calibrating the Covariance Model $feature_dir/$seq_id.cm_msa1.  "
        echo "======================================================================="
        echo ""
	    cmcalibrate --cpu $np $feature_dir/$seq_id.cm_msa1

	    if [ $? -eq 0 ]; then
	        echo ""
	        echo "==========================================================="
            echo "    CM calibrated $feature_dir/$seq_id.cm_msa1 successfully.    "
	        echo "==========================================================="
	        echo ""
	    else
            echo ""
            echo "================================================================================"
            echo "     Error occured while calibrating $feature_dir/$seq_id.cm_msa1.              "
            echo "================================================================================"
            echo ""
            exit 1
        fi
    fi

    echo ""
    echo "==========================================================================================================================="
    echo "        Second round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm_msa1.    "
    echo "                 May take 15 mins to few hours for this step.                                                              "
    echo "==========================================================================================================================="
    echo ""
    ipart=0
    rm -f $feature_dir/$seq_id.msalist_msa2
    for p1 in "${path_infernal_database[@]}"; do
        if [ -s $feature_dir/$seq_id.${ipart}.msa_msa2 ]; then
            # echo "find non-empty and existing cmsearch ${ipart}.msa2, skipping"
            echo "$feature_dir/$seq_id.${ipart}.msa_msa2" >>"$feature_dir/$seq_id.msalist_msa2"
            let ipart+=1
            continue
        fi
        echo "cmsearch on database $p1 ..."
	    cmsearch -o $feature_dir/$seq_id.${ipart}.out_msa2 -A $feature_dir/$seq_id.${ipart}.msa_msa2 --cpu $np --incE 10.0 $feature_dir/$seq_id.cm_msa1 "$p1"
        if [ $? -ne 0 ] && [ "$ignore_volError" == F ]; then
             echo ""
             echo "========================================================================================="
             echo "     Error occured during the second round search using CM $feature_dir/$seq_id.cm_msa1."
             echo "========================================================================================="
             echo ""
             exit 1
        fi
        if [ -s "$feature_dir/$seq_id.${ipart}.msa_msa2" ]; then
            echo "$feature_dir/$seq_id.${ipart}.msa_msa2" >>"$feature_dir/$seq_id.msalist_msa2"
        fi
        let ipart+=1
    done
    if [ -s "$feature_dir/$seq_id.msalist_msa2" ]; then
        esl-alimerge -o "$feature_dir/$seq_id.msa_msa2" --list "$feature_dir/$seq_id.msalist_msa2"
    else  # all .msa_msa2 are empty
        cp "$feature_dir/$seq_id.0.msa_msa2" "$feature_dir/$seq_id.msa_msa2"
    fi

	# cmsearch -o $feature_dir/$seq_id.out_msa2 -A $feature_dir/$seq_id.msa_msa2 --cpu 72 --incE 10.0 $feature_dir/$seq_id.cm_msa1 $path_infernal_database

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      Second round of MSA search (MSA-2) completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================="
        echo "     Error occured during the second round search using CM $feature_dir/$seq_id.cm_msa1. "
        echo "========================================================================================="
        echo ""
        exit 1
    fi

	######### reformat the alignment without gaps and dashes  ###############
    echo ""
    echo "============================================================================"
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa_msa2   "
    echo "============================================================================"
    echo ""

	##### check if .msa_msa2 is not empty  #########
	if [[ -s $feature_dir/$seq_id.msa_msa2 ]] 
	  then 
		esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa_msa2 > $feature_dir/temp.a2m_msa2 || exit 1
		sed -i 's/\s.*$//' $feature_dir/temp.a2m_msa2  || exit 1 # remove everything after space
		sed -i "s/$seq_id/$seq_id E=0.0/g" $feature_dir/temp.a2m_msa2 || exit 1  # 
	else 
	  cat $feature_dir/$seq_id.fasta > $feature_dir/temp.a2m_msa2
	  cat $feature_dir/$seq_id.fasta >> $feature_dir/temp.a2m_msa2
	  sed -i '$ s/.$/./' $feature_dir/temp.a2m_msa2
	fi

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================"
        echo "   Reformatted the $feature_dir/$seq_id.msa_msa2 successfully.  "
	    echo "================================================================"
	    echo ""
	else
        echo ""
        echo "============================================================================================="
        echo "     Error occured during the refomatting the alignment file $feature_dir/$seq_id.msa_msa2.  "
        echo "============================================================================================="
        echo ""
        exit 1
    fi

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s $feature_dir/temp.a2m_msa2 > $feature_dir/temp1.a2m_msa2 || exit 1
	$program_dir/utils/seqkit rmdup -n $feature_dir/temp1.a2m_msa2 > $feature_dir/temp.a2m_msa2 || exit 1
	cp $feature_dir/temp.a2m_msa2 $feature_dir/$seq_id.rmdup.a2m_msa2

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================="
        echo "   Duplicate sequences removed successfully.   "
	    echo "==============================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the removel of duplicates from MSA-2.  "
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.rmdup.a2m_msa2 | sed '/^$/d' > $feature_dir/temp.a2m_msa2 
	############# add query sequence at the top of MSA file and consider top 50000 RNAs  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m_msa2 | sed '/^[[:space:]]*$/d' > $feature_dir/temp1.a2m_msa2
    head -n100000 $feature_dir/temp1.a2m_msa2 > $feature_dir/$seq_id.a2m_msa2 
	

fi


############### check Neff-value of $feature_dir/$seq_id.a2m_msa2  alignment  ############
neff=`$program_dir/GREMLIN_CPP/gremlin_cpp -only_neff -alphabet rna -i $feature_dir/$seq_id.a2m_msa2 | grep 'NEFF' | awk '{ print $3}'`
thres=50.0
echo $neff

if (( $(echo "$neff > $thres" |bc -l) )) || [[ ! -s $feature_dir/$seq_id.msa_msa2 ]]; then
	
	echo ""
	echo "=============================================================="
	echo "Neff-value greater than 50 or No-hit found by Infernal search " 
	echo "               Evalute DCA from MSA-2                         "
	echo "=============================================================="
	echo ""

	############### run dca predictors ############
	if [[ $sel_dca = "gremlin" ]]; then

		$program_dir/GREMLIN_CPP/gremlin_cpp -alphabet rna -i $feature_dir/$seq_id.a2m_msa2 -o $feature_dir/$seq_id.dca_gremlin &> $feature_dir/$seq_id.log_gremlin

	elif [[ $sel_dca = "plmc" ]]; then

		$program_dir/plmc/bin/plmc -c $feature_dir/$seq_id.dca_plmc -a -.ACGUNX -le 20 -lh 0.01 -m 50 $feature_dir/$seq_id.a2m_msa2 &> $feature_dir/$seq_id.log_plmc

	elif [[ $sel_dca = "mfdca" ]]; then

		mfdca compute_fn rna $feature_dir/$seq_id.a2m_msa2 --apc --pseudocount 0.5 --verbose &> temp.log

	elif [[ $sel_dca = "plmdca" ]]; then

		plmdca compute_fn rna $feature_dir/$seq_id.a2m_msa2 --max_iterations 500 --num_threads $np --apc --verbose &> temp.log
		
	fi

else

	echo ""
	echo "=============================="
	echo "  Neff-value less than 50     "
	echo "  Going for of MSA-3 search   "
	echo "=============================="
	echo ""

	$program_dir/utils/reformat.pl a2m sto $feature_dir/$seq_id.a2m_msa2 $feature_dir/$seq_id.sto_msa2
    ## Klark note: variate duplicate seqnames in sto alignments
    mv -f $feature_dir/$seq_id.sto_msa2 $feature_dir/$seq_id.sto_msa2.tmp
    python $program_dir/utils/unique_seqname.py $feature_dir/$seq_id.sto_msa2.tmp >$feature_dir/$seq_id.sto_msa2
    ## Klark done
	sed -i 's/#=GF DE/#=GF DE                          E=0.0/g' $feature_dir/$seq_id.sto_msa2  # missing line at the top with E=0.0

	RNAfold $feature_dir/$seq_id.fasta | awk '{print $1}' | tail -n +3 > $feature_dir/$seq_id.db

	################ reformat ss with according to gaps in reference sequence of .sto file from blastn ################
	for i in `awk '{print $2}' $feature_dir/$seq_id.sto_msa2 | head -n5 | tail -n1 | grep -b -o - | sed 's/..$//'`; do sed -i "s/./&-/$i" $feature_dir/$seq_id.db; done

	#########  add reformated ss from last step to .sto file of blastn ##############
	head -n -1 $feature_dir/$seq_id.sto_msa2 > $feature_dir/temp.sto
	echo "#=GC SS_cons                     "`cat $feature_dir/$seq_id.db` > $feature_dir/temp.txt
	cat $feature_dir/temp.sto $feature_dir/temp.txt > $feature_dir/$seq_id.sto_msa2
	echo "//" >> $feature_dir/$seq_id.sto_msa2


	######## run infernal round-2 ################
    echo ""
    echo "==============================================================================================================================="
    echo "      Building Covariance Model from infernal MSA-2 alignment (with SS from RNAfold) from $feature_dir/$seq_id.sto_msa2 file.  "
    echo "==============================================================================================================================="
    echo ""
	cmbuild --hand -F $feature_dir/$seq_id.cm_msa2 $feature_dir/$seq_id.sto_msa2

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================================="
        echo "    Covariance Model (CM) built successfully from $feature_dir/$seq_id.sto_msa2. "
	    echo "================================================================================="
	    echo ""
	else
        echo ""
        echo "==============================================================================================="
        echo "     Error occured while building Covariance Model (CM) from cmbuild.           "
        echo "==============================================================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "======================================================================="
    echo "       Calibrating the Covariance Model $feature_dir/$seq_id.cm_msa2.  "
    echo "======================================================================="
    echo ""
	cmcalibrate --cpu $np $feature_dir/$seq_id.cm_msa2

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "    CM calibrated $feature_dir/$seq_id.cm_msa2 successfully.    "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "================================================================================"
        echo "     Error occured while calibrating $feature_dir/$seq_id.cm_msa2.              "
        echo "================================================================================"
        echo ""
        exit 1
    fi

    echo ""
    echo "==========================================================================================================================="
    echo "        Third round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm_msa2.    "
    echo "                 May take 15 mins to few hours for this step.                                                              "
    echo "==========================================================================================================================="
    echo ""
    ipart=0
    rm -f $feature_dir/$seq_id.msalist_msa3
    for p1 in "${path_infernal_database[@]}"; do
        if [ -s "$feature_dir/$seq_id.${ipart}.msa_msa3" ]; then
            # echo "find non-empty and existing cmsearch ${ipart}.msa3, skipping"
            echo "$feature_dir/$seq_id.${ipart}.msa_msa3" >>"$feature_dir/$seq_id.msalist_msa3"
            let ipart+=1
            continue
        fi
        echo "cmsearch on database $p1 ..."
	    cmsearch -o $feature_dir/$seq_id.${ipart}.out_msa3 -A $feature_dir/$seq_id.${ipart}.msa_msa3 --cpu $np --incE 10.0 $feature_dir/$seq_id.cm_msa2 "$p1"
	    if [ $? -ne 0 ] && [ "$ignore_volError" == F ]; then
            echo ""
            echo "========================================================================================="
            echo "     Error occured during the third round search using CM $feature_dir/$seq_id.cm_msa2. "
            echo "========================================================================================="
            echo ""
            exit 1
        fi
        if [ -s "$feature_dir/$seq_id.${ipart}.msa_msa3" ]; then
            echo "$feature_dir/$seq_id.${ipart}.msa_msa3" >>"$feature_dir/$seq_id.msalist_msa3"
        fi
        let ipart+=1
    done
    if [ -s "$feature_dir/$seq_id.msalist_msa3" ]; then
        esl-alimerge -o "$feature_dir/$seq_id.msa_msa3" --list "$feature_dir/$seq_id.msalist_msa3"
    else  # all .msa_msa3 are empty
        cp "$feature_dir/$seq_id.0.msa_msa3" "$feature_dir/$seq_id.msa_msa3"
    fi

# 	cmsearch -o $feature_dir/$seq_id.out_msa3 -A $feature_dir/$seq_id.msa_msa3 --cpu 16 --incE 10.0 $feature_dir/$seq_id.cm_msa2 $path_infernal_database

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      Third round of MSA search (MSA-3) completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================="
        echo "     Error occured during the third round search using CM $feature_dir/$seq_id.cm_msa2. "
        echo "========================================================================================="
        echo ""
        exit 1
    fi

	######### reformat the alignment without gaps and dashes  ###############
    echo ""
    echo "============================================================================"
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa_msa3   "
    echo "============================================================================"
    echo ""

	##### check if .msa_msa3 is not empty  #########
	if [[ -s $feature_dir/$seq_id.msa_msa3 ]] 
	  then 
		esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa_msa3 > $feature_dir/temp.a2m_msa3 || exit 1
	else 
	  cat $feature_dir/$seq_id.fasta > $feature_dir/temp.a2m_msa3 || exit 1
	  cat $feature_dir/$seq_id.fasta >> $feature_dir/temp.a2m_msa3 || exit 1
	  sed -i '$ s/.$/./' $feature_dir/temp.a2m_msa3 || exit 1
	fi

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "================================================================"
        echo "   Reformatted the $feature_dir/$seq_id.msa_msa3 successfully.  "
	    echo "================================================================"
	    echo ""
	else
        echo ""
        echo "============================================================================================="
        echo "     Error occured during the refomatting the alignment file $feature_dir/$seq_id.msa_msa3.  "
        echo "============================================================================================="
        echo ""
        exit 1
    fi

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s $feature_dir/temp.a2m_msa3 > $feature_dir/$seq_id.rmdup.a2m_msa3

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================="
        echo "   Duplicate sequences removed successfully.   "
	    echo "==============================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the removel of duplicates from MSA-3.  "
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.rmdup.a2m_msa3 | sed '/^$/d' > $feature_dir/temp.a2m_msa3 || exit 1
	############# add query sequence at the top of MSA file and consider top 50000 RNAs  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m_msa3 | sed '/^[[:space:]]*$/d' > $feature_dir/temp1.a2m_msa3 || exit 1
    head -n100000 $feature_dir/temp1.a2m_msa3 > $feature_dir/$seq_id.a2m_msa3 || exit 1
#	sed -i '/^[[:space:]]*$/d' $feature_dir/$seq_id.a2m_msa3 


	############### run dca predictors ############
	if [[ $sel_dca = "gremlin" ]]; then

		$program_dir/GREMLIN_CPP/gremlin_cpp -alphabet rna -i $feature_dir/$seq_id.a2m_msa3 -o $feature_dir/$seq_id.dca_gremlin &> $feature_dir/$seq_id.log_gremlin

	elif [[ $sel_dca = "plmc" ]]; then

		$program_dir/plmc/bin/plmc -c $feature_dir/$seq_id.dca_plmc -a -.ACGUNX -le 20 -lh 0.01 -m 50 $feature_dir/$seq_id.a2m_msa3 &> $feature_dir/$seq_id.log_plmc

	elif [[ $sel_dca = "mfdca" ]]; then

		mfdca compute_fn rna $feature_dir/$seq_id.a2m_msa3 --apc --pseudocount 0.5 --verbose &> temp.log

	elif [[ $sel_dca = "plmdca" ]]; then

		plmdca compute_fn rna $feature_dir/$seq_id.a2m_msa3 --max_iterations 500 --num_threads $np --apc --verbose &> temp.log
		
	fi
fi

end=`date +%s`

runtime=$((end-start))

echo -e "\ncomputation time = "$runtime" seconds"

