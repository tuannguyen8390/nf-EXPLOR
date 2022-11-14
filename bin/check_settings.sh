# check if params are both true

#####################
nanofilt=`grep params.enable_nanofilt ../nextflow.config | gawkp.scr - 3` 
filtlong=`grep params.enable_filtlong  ../nextflow.config | gawkp.scr - 3` 

minimap2=`grep params.enable_minimap2 ../nextflow.config | gawkp.scr - 3` 
ngmlr=`grep params.enable_ngmlr ../nextflow.config | gawkp.scr - 3` 
winnowmap2=`grep params.enable_winnowmap2 ../nextflow.config | gawkp.scr - 3` 

#####################
echo "Checking if params are OK for analysis"
echo "######################################"
echo "QC       ||| Nanofilt + FLRMC: $nanofilt , filtlong: $filtlong"
echo "Mapping  ||| Minimap: $minimap2 , NGMLR: $ngmlr, Winnowmap: $winnowmap2"
echo "######################################"
echo ""
echo "Results  "
echo ""
echo "######################################"
if [[ $nanofilt == "false" ]] && [[ $filtlong == "false" ]]; then
        echo "ERROR QC PARAMS: Please select at least one QC method"      
fi

if [[ "$nanofilt" != "false" && "$nanofilt" != "false" ]] 
then
    echo "ERROR QC PARAMS: Please select only one mapper"
fi

#####################


if [[ $minimap2 == "false" ]] && [[ $ngmlr == "false" ]] && [[ $winnowmap2 == "false" ]]; then
        echo "ERROR MAPPING PARAMS: Please select at least one mapper"      
fi

if [[ "$minimap2" != "false" && "$ngmlr" != "false" ]] || [[ "$minimap2" != "false" && "$winnowmap2" != "false" ]] || [[ "$ngmlr" != "false" && "$winnowmap2" != "false" ]] 
then
    echo "ERROR MAPPING PARAMS: Please select only one mapper"
fi
echo "######################################"

