# conda activate alevin

################################
## create index for GTF files ##
################################

# export PATH="/bi/group/reik/ricard/data/software:$PATH"
# gtf_transposable_elements="/bi/group/reik/ricard/data/mm10_sequence/transposable_elements/gtf/mm10_rmsk_TE.gtf"
# job 25 1 "TEtranscripts_indexer --afile $gtf_transposable_elements --itype TE"
# gtf_genes="/bi/group/reik/ricard/data/ensembl/mouse/v87/GTF/all/Mus_musculus.GRCm38.87.chr.gtf"
# job 25 1 "TEtranscripts_indexer --afile $gtf_genes --itype gene"

##############
## Settings ##
##############

gtf_transposable_elements="/bi/group/reik/ricard/data/mm10_sequence/transposable_elements/gtf/mm10_rmsk_TE.gtf.ind"
gtf_genes="/bi/group/reik/ricard/data/ensembl/mouse/v87/GTF/all/Mus_musculus.GRCm38.87.chr.gtf.ind"

samples=( "SRR10055443" "SRR10055444" "SRR10055445" "SRR10055446" "SRR10055447" "SRR10055448" "SRR10055449" "SRR10055450" "SRR10055451" "SRR10055452" "SRR10055453" "SRR10055454" "SRR10055455" "SRR10055456" "SRR10055457" "SRR10055458" "SRR10055459" "SRR10055460" "SRR10055461" "SRR10055462" "SRR10055463" "SRR10055464" "SRR10055465" "SRR10055466" "SRR10055467" "SRR10055468" "SRR10055469" "SRR10055470" "SRR10055471" "SRR10055472" "SRR10055473" "SRR10055474" "SRR10055475" "SRR10055476" "SRR10055477" "SRR10055478" "SRR10055479" "SRR10055480" "SRR10055481" "SRR10055482" "SRR10055483" "SRR10055484" "SRR10055485" "SRR10055486" "SRR10055487" "SRR10055488" "SRR10055489" "SRR10055490" "SRR10055491" "SRR10055492" "SRR10055493" "SRR10055494" "SRR10055495" "SRR10055496" "SRR10055497" "SRR10055498" "SRR10055499" "SRR10055500" "SRR10055501" "SRR10055502" "SRR10055503" "SRR10055504" "SRR10055505" "SRR10055506" "SRR10055507" "SRR10055508" "SRR10055509" "SRR10055510" "SRR10055511" "SRR10055512" "SRR10055513" "SRR10055514" "SRR10055515" "SRR10055516" "SRR10055517" "SRR10055518" "SRR10055519" "SRR10055520" "SRR10055521" "SRR10055522" "SRR10055523" "SRR10055524" "SRR10055525" "SRR10055526" "SRR10055527" "SRR10055528" "SRR10055529" "SRR10055530" "SRR10055531" "SRR10055532" "SRR10055533" "SRR10055534" "SRR10055535" "SRR10055536" "SRR10055537" "SRR10055538" "SRR10055539" "SRR10055540" "SRR10055541" "SRR10055542" "SRR10055543" "SRR10055544" "SRR10055545" "SRR10055546" "SRR10055547" "SRR10055548" "SRR10055549" "SRR10055550" "SRR10055551" "SRR10055552" "SRR10055553" "SRR10055554" "SRR10055555" "SRR10055556" "SRR10055557" "SRR10055558" "SRR10055559" "SRR10055560" "SRR10055561" "SRR10055562" "SRR10055563" "SRR10055564" "SRR10055565" "SRR10055566" "SRR10055567" "SRR10055568" "SRR10055569" "SRR10055570" "SRR10055571" "SRR10055572" "SRR10055573" "SRR10055574" "SRR10055575" "SRR10055576" "SRR10055577" "SRR10055578" "SRR10055579" "SRR10055580" "SRR10055581" "SRR10055582" "SRR10055583" "SRR10055584" "SRR10055585" "SRR10055586" "SRR10055587" "SRR10055588" "SRR10055589" "SRR10055590" "SRR10055591" "SRR10055592" "SRR10055593" "SRR10055594" "SRR10055595" "SRR10055596" "SRR10055597" "SRR10055598" "SRR10055599" "SRR10055600" "SRR10055601" "SRR10055602" "SRR10055603" "SRR10055604" "SRR10055605" "SRR10055606" "SRR10055607" "SRR10055608" "SRR10055609" "SRR10055610" "SRR10055611" "SRR10055612" "SRR10055613" "SRR10055614" "SRR10055615" "SRR10055616" "SRR10055617" "SRR10055618" "SRR10055619" "SRR10055620" "SRR10055621" "SRR10055622" "SRR10055623" "SRR10055624" "SRR10055625" "SRR10055626" "SRR10055627" "SRR10055628" "SRR10055629" "SRR10055630" "SRR10055631" "SRR10055632" "SRR10055633" "SRR10055634" "SRR10055635" "SRR10055636" "SRR10055637" "SRR10055638" "SRR10055639" "SRR10055640" "SRR10055641" "SRR10055642" "SRR10055643" "SRR10055644" "SRR10055645" "SRR10055646" "SRR10055647" "SRR10055648" "SRR10055649" "SRR10055650" "SRR10055651" "SRR10055652" "SRR10055653" "SRR10055654" "SRR10055655" "SRR10055656" "SRR10055657" "SRR10055658" "SRR10055659" "SRR10055660" "SRR10055661" "SRR10055662" "SRR10055663" "SRR10055664" "SRR10055665" "SRR10055666" "SRR10055667" "SRR10055668" "SRR10055669" "SRR10055670" "SRR10055671" "SRR10055672" "SRR10055673" "SRR10055674" "SRR10055675" "SRR10055676" "SRR10055677" "SRR10055678" "SRR10055679" "SRR10055680" "SRR10055681" "SRR10055682" "SRR10055683" "SRR10055684" "SRR10055685" "SRR10055686" "SRR10055687" "SRR10055688" "SRR10055689" "SRR10055690" "SRR10055691" "SRR10055692" "SRR10055693" "SRR10055694" "SRR10055695" "SRR10055696" "SRR10055697" "SRR10055698" "SRR10055699" "SRR10055700" "SRR10055701" "SRR10055702" "SRR10055703" "SRR10055704" "SRR10055705" "SRR10055706" "SRR10055707" "SRR10055708" "SRR10055709" "SRR10055710" "SRR10055711" "SRR10055712" "SRR10055713" "SRR10055714" "SRR10055715" "SRR10055716" "SRR10055717" "SRR10055718" "SRR10055719" "SRR10055720" "SRR10055721" "SRR10055722" "SRR10055723" "SRR10055724" "SRR10055725" "SRR10055726" "SRR10055727" "SRR10055728" "SRR10055729" "SRR10055730" "SRR10055731" "SRR10055732" "SRR10055733" "SRR10055734" "SRR10055735" "SRR10055736" "SRR10055737" "SRR10055738" "SRR10055739" "SRR10055740" "SRR10055741" "SRR10055742" "SRR10055743" "SRR10055744" "SRR10055745" "SRR10055746" "SRR10055747" "SRR10055748" "SRR10055749" "SRR10055750" "SRR10055751" "SRR10055752" "SRR10055753" "SRR10055754" "SRR10055755" "SRR10055756" "SRR10055757" "SRR10055758" "SRR10055759" "SRR10055760" "SRR10055761" "SRR10055762" "SRR10055763" "SRR10055764" "SRR10055765" "SRR10055766" "SRR10055767" "SRR10055768" "SRR10055769" "SRR10055770" "SRR10055771" "SRR10055772" "SRR10055773" "SRR10055774" "SRR10055775" "SRR10055776" "SRR10055777" "SRR10055778" "SRR10055779" "SRR10055780" "SRR10055781" "SRR10055782" "SRR10055783" "SRR10055784" "SRR10055785" "SRR10055786" "SRR10055787" "SRR10055788" "SRR10055789" "SRR10055790" "SRR10055791" "SRR10055792" "SRR10055793" "SRR10055794" "SRR10055795" "SRR10055796" "SRR10055797" "SRR10055798" "SRR10055799" "SRR10055800" "SRR10055801" "SRR10055802" "SRR10055803" "SRR10055804" "SRR10055805" "SRR10055806" "SRR10055807" "SRR10055808" "SRR10055809" "SRR10055810" "SRR10055811" "SRR10055812" "SRR10055813" "SRR10055814" "SRR10055815" "SRR10055816" "SRR10055817" "SRR10055818" "SRR10055819" "SRR10055820" "SRR10055821" "SRR10055822" "SRR10055823" "SRR10055824" "SRR10055825" "SRR10055826" "SRR10055827" "SRR10055828" "SRR10055829" "SRR10055830" "SRR10055831" "SRR10055832" "SRR10055833" "SRR10055834" "SRR10055835" "SRR10055836" "SRR10055837" "SRR10055838" "SRR10055839" "SRR10055840" "SRR10055841" "SRR10055842" "SRR10055843" )

################################
## create index for bam files ##
################################

# cd /bi/group/reik/ricard/data/Wang2021_mouse_preimplantation/star

# i="SRR1041728"
# for i in "${samples[@]}"; do
# 	if ! [ -f "${i}_Aligned.sortedByCoord.out.bam.bai" ]; then
# 		echo $i
# 	    job 10 1 "samtools index ${i}_Aligned.sortedByCoord.out.bam"
# 	fi
# done

#############
## TEcount ##
#############

bam_dir="/bi/group/reik/ricard/data/Wang2021_mouse_preimplantation/star"
outdir="/bi/group/reik/ricard/data/Wang2021_mouse_preimplantation/TEcounts"

mkdir $outdir; cd $outdir

# i="SRR1041728"
for i in "${samples[@]}"; do

	if ! [ -f "${i}.cntTable.gz" ]; then
		echo $i
		file="${bam_dir}/${i}_Aligned.sortedByCoord.out.bam"
		cmd="TEcount --GTF $gtf_genes --TE $gtf_transposable_elements --sortByPos --format BAM --mode multi -b $file --project $i --outdir $outdir"
	    job 40 1 "$cmd"
	fi
done


# SAMPLES THAT FAILED
# SRR10055692
# SRR10055443