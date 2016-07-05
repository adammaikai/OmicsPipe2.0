Oncogenomics Report for patient ID <code class="knitr inline">SRR1027184</code>
========================================================

<!-- dataTables -->
<script src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.js"></script>
<link rel="stylesheet" href="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/css/jquery.dataTables.css"/>
<script src="http://next.datatables.net/release-datatables/examples/resources/bootstrap/3/dataTables.bootstrap.js"></script>
<link rel="stylesheet" href="http://next.datatables.net/release-datatables/examples/resources/bootstrap/3/dataTables.bootstrap.css"/>

<!-- TableTools -->
<script src="http://cdn.aldu.net/jquery.datatables/extras/tabletools/2.0.3/js/TableTools.js"></script>
<link rel="stylesheet" href="http://cdn.aldu.net/jquery.datatables/extras/tabletools/2.0.3/css/TableTools.css"/>

<!-- ZeroClipboard -->
<script src="http://cdn.aldu.net/jquery.datatables/extras/tabletools/2.0.3/js/ZeroClipboard.js"></script>

<!-- ColVis -->
<script src="http://cdn.eventcore.com/js/plugins/dataTables.ColVis-1.07.js"></script>

<!-- Scientific notation sorting -->
<script src="../media/src/scientific.js"></script>

<!-- read in external core code -->


<!-- read in all optional modules from the module folder -->






# Quality Control & Preprocessing

## Raw data QC

The quality of the RNA-seq reads was assessed using FASTQC <code class="knitr inline">(, 2014)</code>. The table below lists basic QC statistics. Click on the provided links to open up the full FastQC report.


<table id="fastqc_table", class="stdTable">
<tbody>
  <tr>
   <td align="left">  </td>
   <td align="left"> Pair_1 </td>
   <td align="left"> Pair_2 </td>
  </tr>
  <tr>
   <td align="left"> file </td>
   <td align="left"> SRR1027184_1.fastq </td>
   <td align="left"> SRR1027184_2.fastq </td>
  </tr>
  <tr>
   <td align="left"> basic_statistics </td>
   <td align="left"> pass </td>
   <td align="left"> pass </td>
  </tr>
  <tr>
   <td align="left"> file_type </td>
   <td align="left"> Conventional base calls </td>
   <td align="left"> Conventional base calls </td>
  </tr>
  <tr>
   <td align="left"> encoding </td>
   <td align="left"> Illumina 1.5 </td>
   <td align="left"> Illumina 1.5 </td>
  </tr>
  <tr>
   <td align="left"> total_seq </td>
   <td align="left"> 68888018 </td>
   <td align="left"> 68888018 </td>
  </tr>
  <tr>
   <td align="left"> filtered_seq </td>
   <td align="left"> 0 </td>
   <td align="left"> 0 </td>
  </tr>
  <tr>
   <td align="left"> seq_length </td>
   <td align="left"> 60 </td>
   <td align="left"> 60 </td>
  </tr>
  <tr>
   <td align="left"> percent_gc </td>
   <td align="left"> 52 </td>
   <td align="left"> 53 </td>
  </tr>
  <tr>
   <td align="left"> link_to_report </td>
   <td align="left"> <a href =  '/gpfs/group/su/kfisch/Cancer/Kumar/QC/SRR1027184_1_fastqc/fastqc_report.html' target="_blank" > FastQC report </a> </td>
   <td align="left"> <a href =  '/gpfs/group/su/kfisch/Cancer/Kumar/QC/SRR1027184_2_fastqc/fastqc_report.html' target="_blank" > FastQC report </a> </td>
  </tr>
</tbody>
</table>


******

## Alignment

Reads were aligned to the human genome (hg19) using the STAR aligner <code class="knitr inline">(Dobin, Davis, Schlesinger, Drenkow, Zaleski, Jha, Batut, Chaisson, and Gingeras, 2012)</code>. Resulting .SAM files were sorted and converted to .BAM files using SAMtools software <code class="knitr inline">(Li, Handsaker, Wysoker, Fennell, Ruan, Homer, Marth, Abecasis, and Durbin, 2009)</code>. The table below lists the alignment statistics as provided by the STAR output:

<div style="width:40%;padding:0 10pt 0 0;float:left;">
<table id="star_table", class="stdTable">
<tbody>
  <tr>
   <td align="left"> Metric </td>
   <td align="left"> Value </td>
  </tr>
  <tr>
   <td align="left"> Started job on  </td>
   <td align="left"> Jun 29 18:08:04 </td>
  </tr>
  <tr>
   <td align="left"> Started mapping on  </td>
   <td align="left"> Jun 29 18:20:30 </td>
  </tr>
  <tr>
   <td align="left"> Finished on  </td>
   <td align="left"> Jun 29 18:40:34 </td>
  </tr>
  <tr>
   <td align="left"> Mapping speed, Million of reads per hour  </td>
   <td align="left"> 205.98 </td>
  </tr>
  <tr>
   <td align="left"> Number of input reads  </td>
   <td align="left"> 68888018 </td>
  </tr>
  <tr>
   <td align="left"> Average input read length  </td>
   <td align="left"> 120 </td>
  </tr>
  <tr>
   <td align="left"> UNIQUE READS: </td>
   <td align="left">  </td>
  </tr>
  <tr>
   <td align="left"> Uniquely mapped reads number  </td>
   <td align="left"> 45311287 </td>
  </tr>
  <tr>
   <td align="left"> Uniquely mapped reads %  </td>
   <td align="left"> 65.78% </td>
  </tr>
  <tr>
   <td align="left"> Average mapped length  </td>
   <td align="left"> 117.14 </td>
  </tr>
  <tr>
   <td align="left"> Number of splices: Total  </td>
   <td align="left"> 2616281 </td>
  </tr>
  <tr>
   <td align="left"> Number of splices: Annotated (sjdb)  </td>
   <td align="left"> 0 </td>
  </tr>
  <tr>
   <td align="left"> Number of splices: GT/AG  </td>
   <td align="left"> 2596502 </td>
  </tr>
  <tr>
   <td align="left"> Number of splices: GC/AG  </td>
   <td align="left"> 18053 </td>
  </tr>
  <tr>
   <td align="left"> Number of splices: AT/AC  </td>
   <td align="left"> 1726 </td>
  </tr>
  <tr>
   <td align="left"> Number of splices: Non-canonical  </td>
   <td align="left"> 0 </td>
  </tr>
  <tr>
   <td align="left"> Mismatch rate per base, %  </td>
   <td align="left"> 1.26% </td>
  </tr>
  <tr>
   <td align="left"> Deletion rate per base  </td>
   <td align="left"> 0.01% </td>
  </tr>
  <tr>
   <td align="left"> Deletion average length  </td>
   <td align="left"> 1.89 </td>
  </tr>
  <tr>
   <td align="left"> Insertion rate per base  </td>
   <td align="left"> 0.01% </td>
  </tr>
  <tr>
   <td align="left"> Insertion average length  </td>
   <td align="left"> 1.65 </td>
  </tr>
  <tr>
   <td align="left"> MULTI-MAPPING READS: </td>
   <td align="left">  </td>
  </tr>
  <tr>
   <td align="left"> Number of reads mapped to multiple loci  </td>
   <td align="left"> 3615953 </td>
  </tr>
  <tr>
   <td align="left"> % of reads mapped to multiple loci  </td>
   <td align="left"> 5.25% </td>
  </tr>
  <tr>
   <td align="left"> Number of reads mapped to too many loci  </td>
   <td align="left"> 56282 </td>
  </tr>
  <tr>
   <td align="left"> % of reads mapped to too many loci  </td>
   <td align="left"> 0.08% </td>
  </tr>
  <tr>
   <td align="left"> UNMAPPED READS: </td>
   <td align="left">  </td>
  </tr>
  <tr>
   <td align="left"> % of reads unmapped: too many mismatches  </td>
   <td align="left"> 0.00% </td>
  </tr>
  <tr>
   <td align="left"> % of reads unmapped: too short  </td>
   <td align="left"> 28.83% </td>
  </tr>
  <tr>
   <td align="left"> % of reads unmapped: other  </td>
   <td align="left"> 0.07% </td>
  </tr>
</tbody>
</table>

</div>

<div style="width:60%;padding:0 10pt 0 0;float:right;">
<iframe src='
figure/starQCPie.html
' scrolling='no' seamless class='rChart 
nvd3
 '
id=iframe-
chartfe937a9efbb
></iframe>
<style>iframe.rChart{ width: 100%; height: 400px;}</style>

</div>
<div style="clear:both;"></div>

## Insert size distribution

Insert size distribution is calculated using Picard tools [CollectInsertSizeMetrics](http://picard.sourceforge.net/command-line-overview.shtml#CollectInsertSizeMetrics) function <code class="knitr inline">(, 2014)</code>. The table lists statistics as provided by this tool:

<div style="width:40%;padding:0 10pt 0 0;float:left;">
<table id="insertSize_table", class="stdTable">
<tbody>
  <tr>
   <td align="left"> Metric </td>
   <td align="left"> Value </td>
  </tr>
  <tr>
   <td align="left"> MEDIAN_INSERT_SIZE </td>
   <td align="left"> 185 </td>
  </tr>
  <tr>
   <td align="left"> MEDIAN_ABSOLUTE_DEVIATION </td>
   <td align="left"> 29 </td>
  </tr>
  <tr>
   <td align="left"> MIN_INSERT_SIZE </td>
   <td align="left"> 20 </td>
  </tr>
  <tr>
   <td align="left"> MAX_INSERT_SIZE </td>
   <td align="left"> 2342127 </td>
  </tr>
  <tr>
   <td align="left"> MEAN_INSERT_SIZE </td>
   <td align="left"> 173.657428 </td>
  </tr>
  <tr>
   <td align="left"> STANDARD_DEVIATION </td>
   <td align="left"> 51.598197 </td>
  </tr>
  <tr>
   <td align="left"> READ_PAIRS </td>
   <td align="left"> 48927240 </td>
  </tr>
  <tr>
   <td align="left"> PAIR_ORIENTATION </td>
   <td align="left"> FR </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_10_PERCENT </td>
   <td align="left"> 9 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_20_PERCENT </td>
   <td align="left"> 19 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_30_PERCENT </td>
   <td align="left"> 29 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_40_PERCENT </td>
   <td align="left"> 41 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_50_PERCENT </td>
   <td align="left"> 59 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_60_PERCENT </td>
   <td align="left"> 83 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_70_PERCENT </td>
   <td align="left"> 115 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_80_PERCENT </td>
   <td align="left"> 167 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_90_PERCENT </td>
   <td align="left"> 257 </td>
  </tr>
  <tr>
   <td align="left"> WIDTH_OF_99_PERCENT </td>
   <td align="left"> 25671 </td>
  </tr>
</tbody>
</table>

</div>

<div style="width:60%;padding:0 10pt 0 0;float:right;">
<img src='/gpfs/group/su/kfisch/Cancer/Kumar/QC/SRR1027184/insertSizeHist.png'>
</div>
<div style="clear:both;"></div>

## RNA sequencing metrics

Metrics are assessed using Picard tools [CollectRnaSeqMetrics](http://picard.sourceforge.net/command-line-overview.shtml#CollectRnaSeqMetrics) function <code class="knitr inline">(, 2014)</code>. Metrics about the alignment of RNA to various functional classes of loci in the genome are collected: coding, intronic, UTR, intergenic, ribosomal. The table below lists the metrics provided by this tool:

<div style="width:40%;padding:0 10pt 0 0;float:left;">
<table id="rnaseqmetrics_table", class="stdTable">
<tbody>
  <tr>
   <td align="left"> Metric </td>
   <td align="left"> Value </td>
  </tr>
  <tr>
   <td align="left"> PF_BASES </td>
   <td align="left"> 58712688000 </td>
  </tr>
  <tr>
   <td align="left"> PF_ALIGNED_BASES </td>
   <td align="left"> 57260965640 </td>
  </tr>
  <tr>
   <td align="left"> RIBOSOMAL_BASES </td>
   <td align="left">             NA </td>
  </tr>
  <tr>
   <td align="left"> CODING_BASES </td>
   <td align="left">  7011864570 </td>
  </tr>
  <tr>
   <td align="left"> UTR_BASES </td>
   <td align="left">  4802658840 </td>
  </tr>
  <tr>
   <td align="left"> INTRONIC_BASES </td>
   <td align="left"> 23449223220 </td>
  </tr>
  <tr>
   <td align="left"> INTERGENIC_BASES </td>
   <td align="left"> 21997219010 </td>
  </tr>
  <tr>
   <td align="left"> IGNORED_READS </td>
   <td align="left">          00 </td>
  </tr>
  <tr>
   <td align="left"> CORRECT_STRAND_READS </td>
   <td align="left">          00 </td>
  </tr>
  <tr>
   <td align="left"> INCORRECT_STRAND_READS </td>
   <td align="left">          00 </td>
  </tr>
  <tr>
   <td align="left"> PCT_RIBOSOMAL_BASES </td>
   <td align="left">             NA </td>
  </tr>
  <tr>
   <td align="left"> PCT_CODING_BASES </td>
   <td align="left">          0.122 </td>
  </tr>
  <tr>
   <td align="left"> PCT_UTR_BASES </td>
   <td align="left">          0.084 </td>
  </tr>
  <tr>
   <td align="left"> PCT_INTRONIC_BASES </td>
   <td align="left">          0.410 </td>
  </tr>
  <tr>
   <td align="left"> PCT_INTERGENIC_BASES </td>
   <td align="left">          0.384 </td>
  </tr>
  <tr>
   <td align="left"> PCT_MRNA_BASES </td>
   <td align="left">          0.206 </td>
  </tr>
  <tr>
   <td align="left"> PCT_USABLE_BASES </td>
   <td align="left">          0.201 </td>
  </tr>
  <tr>
   <td align="left"> PCT_CORRECT_STRAND_READS </td>
   <td align="left">          00 </td>
  </tr>
  <tr>
   <td align="left"> MEDIAN_CV_COVERAGE </td>
   <td align="left">          0.534 </td>
  </tr>
  <tr>
   <td align="left"> MEDIAN_5PRIME_BIAS </td>
   <td align="left">          0.146 </td>
  </tr>
  <tr>
   <td align="left"> MEDIAN_3PRIME_BIAS </td>
   <td align="left">          0.310 </td>
  </tr>
  <tr>
   <td align="left"> MEDIAN_5PRIME_TO_3PRIME_BIAS </td>
   <td align="left">          0.537 </td>
  </tr>
  <tr>
   <td align="left"> SAMPLE </td>
   <td align="left">             NA </td>
  </tr>
  <tr>
   <td align="left"> LIBRARY </td>
   <td align="left">             NA </td>
  </tr>
  <tr>
   <td align="left"> READ_GROUP </td>
   <td align="left">             NA </td>
  </tr>
</tbody>
</table>

</div>

<div style="width:60%;padding:0 10pt 0 0;float:right;">
<iframe src='
figure/RNASeqMetricsPie.html
' scrolling='no' seamless class='rChart 
nvd3
 '
id=iframe-
chartfe920276d5e
></iframe>
<style>iframe.rChart{ width: 100%; height: 400px;}</style>

</div>
<div style="clear:both;"></div>

******

## Gene expression quantification 

Gene expression quantification was done using [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) function within the python HTSeq analysis package <code class="knitr inline">(, 2014)</code>, which counts all reads overlapping known exons.

<div style="width:40%;padding:0 10pt 0 0;float:left;">
<table id="htseq_table", class="stdTable">
<tbody>
  <tr>
   <td align="left"> Metric </td>
   <td align="left"> Value </td>
  </tr>
  <tr>
   <td align="left"> __no_feature </td>
   <td align="left"> 35025875 </td>
  </tr>
  <tr>
   <td align="left"> __ambiguous </td>
   <td align="left"> 135432 </td>
  </tr>
  <tr>
   <td align="left"> __too_low_aQual </td>
   <td align="left"> 0 </td>
  </tr>
  <tr>
   <td align="left"> __not_aligned </td>
   <td align="left"> 0 </td>
  </tr>
  <tr>
   <td align="left"> __alignment_not_unique </td>
   <td align="left"> 10174185 </td>
  </tr>
</tbody>
</table>

</div>

<div style="width:60%;padding:0 10pt 0 0;float:right;">
<iframe src='
figure/htseqQCPie.html
' scrolling='no' seamless class='rChart 
nvd3
 '
id=iframe-
chartfe92e130f7f
></iframe>
<style>iframe.rChart{ width: 100%; height: 400px;}</style>

</div>
<div style="clear:both;"></div>

## Gene expression reliability



Gene expresssion values are grouped into "reliable expressed genes" and "not so reliable expressed genes" based on a background distribution of <code class="knitr inline">136</code> genes which have a median expression value of 0 in the <code class="knitr inline">1057</code> TCGA reference cohort:

* **reliable expressed** Number of genes reliable expressed in the sample
* **marginal expressed** Number of genes that are marginally (between 0.01 & 0.02) reliable expressed in the sample
* **not reliable expressed** Number of genes that are not reliable expressed in the sample



<iframe src='
figure/reliableGEPPlot.html
' scrolling='no' seamless class='rChart 
highcharts
 '
id=iframe-
chartfe97ac665ca
></iframe>
<style>iframe.rChart{ width: 100%; height: 400px;}</style>


******

Expression density plot

<div class="chunk" id="reliableGEPDensityPlot"><div class="rcode"><div class="rimage center"><img src="figure/reliableGEPDensityPlot.png" title="plot of chunk reliableGEPDensityPlot" alt="plot of chunk reliableGEPDensityPlot" class="plot" /></div>
</div></div>

******

# Analysis

## Differential Gene Expression

Differential gene expression is calculated based on a model using the negative binomial distribution as implemented in the DESeq2 package <code class="knitr inline">(Anders and Huber, 2010)</code>. 

<!-- run the differential expression on the patient samples -->


<iframe src='
figure/diffExpPlot.html
' scrolling='no' seamless class='rChart 
highcharts
 '
id=iframe-
chartfe94816c7d1
></iframe>
<style>iframe.rChart{ width: 100%; height: 400px;}</style>


******

## Gene Expression Table

<!-- create gene expression table including results from diff. expression calling if available -->
* **Gene** HUGO Symbol for the gene 
* **Diff_expressed** [upregulated/downregulated] above/below median expression value within the normal control +/- three standard deviations 
* **Rank** gene rank within differentially expressed genes
* **perRank** percentile rank within the distribution of all genes within the patient sample
* **Patient** Deseq scaled and log2 transformed read count value of the patient sample
* **Normal** Deseq scaled and log2 transformed median read count value of the normal reference cohort (n=<code class="knitr inline">4</code>)
* **Normal_sd** standard deviation for each scaled read count value of the normal reference cohort
* **d_normal** delta between patient sample read count value and median read count value of the normal reference cohort
* **Ref** Deseq scaled and log2 transformed median read count value of the TCGA reference cohort (n=<code class="knitr inline">947</code>)
* **perRankRef** percentile rank within the distribution of all genes within the reference cohort
* **Ref_sd** standard deviation for each scaled read count value of the reference cohort
* **d_ref** delta between patient sample read count value and median read count value of the reference cohort
* **List_Gene** genes within the supplied gene list



<table id="gep_table">
    <thead>
        <tr>
            <th>Gene</th>
            <th>Diff_expressed</th>
            <th>Reliable_expressed</th>
            <th>Rank</th>
            <th>perRank</th>
            <th>Patient</th>
            <th>Normal</th>
            <th>Normal_sd</th>
            <th>d_normal</th>
            <th>Ref</th>
            <th>perRankRef</th>
            <th>Ref_sd</th>
            <th>d_ref</th>
            <th>Pvalue</th>
            <th>List_Gene</th>
            <th>GeneWiki</th>
        </tr>
    </thead>
    <tbody></tbody>
</table>

<script type="text/javascript" charset="utf-8">
    $(document).ready(function() {
    var table = $('#gep_table').dataTable({
      "sDom": 'T<"clear">lfrtip',
      "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
            ]
      },
      "aoColumns": [ 
      null,
      null,
      null,
      {"bVisible": false},
		  {"bVisible": false},
      null,
      null,
      {"bVisible": false},
      null,
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"sType": "scientific"},
      {"bVisible": false},
      null
		], 
    "iDisplayLength": 10,
    "bPaginate": true,
    "bProcessing": true,
    "sAjaxSource": 'json/gepall.txt'
    //"sScrollXInner": "110%",
    });
    });    
</script>

******

## Prediction of Receptor Status



* **Receptor** Receptor type
* **Status** predicted receptor status for each receptor type
* **Prediction_probability** prediction probability for each receptor type

<table id="receptor_table">
 <thead>
  <tr>
   <th align="left"> Receptor </th>
   <th align="left"> Status </th>
   <th align="right"> Prediction_probability </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> ER </td>
   <td align="left"> Negative </td>
   <td align="right"> 0.97 </td>
  </tr>
  <tr>
   <td align="left"> PR </td>
   <td align="left"> Negative </td>
   <td align="right"> 0.87 </td>
  </tr>
  <tr>
   <td align="left"> HER2 </td>
   <td align="left"> Positive </td>
   <td align="right"> 0.73 </td>
  </tr>
</tbody>
</table>


<script type="text/javascript" charset="utf-8">
    $(document).ready( function () {
    $('#receptor_table').dataTable({
          "sDom": 'T<"clear">lfrtip',
          "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
          }
    });
    });
</script>

******

## Prediction of molecular subtype

* **Subtype** prediction of one of the subtypes [Luminal A / Luminal B / Basal / Her2]
* **Prediction_probability** prediction probability for each receptor type

<table id="mol_table">
 <thead>
  <tr>
   <th align="left"> Subtype </th>
   <th align="right"> Prediction_probability </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> Basal </td>
   <td align="right"> 1 </td>
  </tr>
</tbody>
</table>


<script type="text/javascript" charset="utf-8">
    $(document).ready( function () {
    $('#mol_table').dataTable({
          "sDom": 'T<"clear">lfrtip',
          "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
          }
    });
    });
</script>

******

## Pathways





### KEGG

<div class="error"><pre class="knitr r">## Error: could not open file 'hsa05322.SRR1027184.png'
</pre></div>


<script type="text/javascript" charset="utf-8">
    $(document).ready( function () {
    $('#kegg_table').dataTable({
          "sDom": 'T<"clear">lfrtip',
          "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
          },
      "aoColumns": [ 
      null,
			{"bVisible": false},
      {"bVisible": false},
		  {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"sType": "scientific"},
      {"bVisible": false},
      null,
      null
		]
    });
    });
</script>

******

### Biocarta



<table id="biocarta_table">
 <thead>
  <tr>
   <th align="left"> Name </th>
   <th align="right"> pSize </th>
   <th align="right"> NDE </th>
   <th align="right"> pNDE </th>
   <th align="right"> tA </th>
   <th align="right"> pPERT </th>
   <th align="right"> pG </th>
   <th align="right"> pGFdr </th>
   <th align="right"> pGFWER </th>
   <th align="left"> Status </th>
   <th align="left"> Image </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> pertussis toxin-insensitive ccr5 signaling in macrophage </td>
   <td align="right"> 9 </td>
   <td align="right"> 4 </td>
   <td align="right"> 0.0616 </td>
   <td align="right"> -16.75 </td>
   <td align="right"> 0 </td>
   <td align="right"> 0 </td>
   <td align="right"> 8e-04 </td>
   <td align="right"> 8e-04 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/BiocartaPlots1.png target="_blank"> <img height=50 width=50 src=figure/BiocartaPlots1.png></a> </td>
  </tr>
</tbody>
</table>


<script type="text/javascript" charset="utf-8">
    $(document).ready( function () {
    $('#biocarta_table').dataTable({
          "sDom": 'T<"clear">lfrtip',
          "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
          },
      "aoColumns": [ 
      null,
  		{"bVisible": false},
      {"bVisible": false},
		  {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"sType": "scientific"},
      {"bVisible": false},
      null,
      null
		]        
    });
    });
</script>

******

### NCI



<table id="nci_table">
 <thead>
  <tr>
   <th align="left"> Name </th>
   <th align="right"> pSize </th>
   <th align="right"> NDE </th>
   <th align="right"> pNDE </th>
   <th align="right"> tA </th>
   <th align="right"> pPERT </th>
   <th align="right"> pG </th>
   <th align="right"> pGFdr </th>
   <th align="right"> pGFWER </th>
   <th align="left"> Status </th>
   <th align="left"> Image </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> Aurora B signaling </td>
   <td align="right"> 39 </td>
   <td align="right"> 19 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 65.037 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/NCIPlots1.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots1.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Signaling events mediated by HDAC Class III </td>
   <td align="right"> 35 </td>
   <td align="right"> 19 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -1.726 </td>
   <td align="right"> 0.679 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0009 </td>
   <td align="right"> 0.0017 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/NCIPlots2.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots2.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Circadian rhythm pathway </td>
   <td align="right"> 13 </td>
   <td align="right"> 9 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> -9.236 </td>
   <td align="right"> 0.173 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0059 </td>
   <td align="right"> 0.0177 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/NCIPlots3.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots3.png></a> </td>
  </tr>
  <tr>
   <td align="left"> E2F transcription factor network </td>
   <td align="right"> 63 </td>
   <td align="right"> 25 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 3.443 </td>
   <td align="right"> 0.399 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 0.0063 </td>
   <td align="right"> 0.0251 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/NCIPlots4.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots4.png></a> </td>
  </tr>
  <tr>
   <td align="left"> BARD1 signaling events </td>
   <td align="right"> 29 </td>
   <td align="right"> 14 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 11.581 </td>
   <td align="right"> 0.141 </td>
   <td align="right"> 0.0003 </td>
   <td align="right"> 0.0072 </td>
   <td align="right"> 0.0361 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/NCIPlots5.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots5.png></a> </td>
  </tr>
  <tr>
   <td align="left"> p73 transcription factor network </td>
   <td align="right"> 66 </td>
   <td align="right"> 26 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 1.000 </td>
   <td align="right"> 0.0004 </td>
   <td align="right"> 0.0081 </td>
   <td align="right"> 0.0486 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/NCIPlots6.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots6.png></a> </td>
  </tr>
  <tr>
   <td align="left"> FOXM1 transcription factor network </td>
   <td align="right"> 36 </td>
   <td align="right"> 17 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 1.000 </td>
   <td align="right"> 0.0006 </td>
   <td align="right"> 0.0101 </td>
   <td align="right"> 0.0709 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/NCIPlots7.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots7.png></a> </td>
  </tr>
  <tr>
   <td align="left"> p53 pathway </td>
   <td align="right"> 163 </td>
   <td align="right"> 49 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 1.887 </td>
   <td align="right"> 0.923 </td>
   <td align="right"> 0.0011 </td>
   <td align="right"> 0.0160 </td>
   <td align="right"> 0.1278 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/NCIPlots8.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots8.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Regulation of p38-alpha and p38-beta </td>
   <td align="right"> 127 </td>
   <td align="right"> 37 </td>
   <td align="right"> 0.0014 </td>
   <td align="right"> -17.646 </td>
   <td align="right"> 0.217 </td>
   <td align="right"> 0.0028 </td>
   <td align="right"> 0.0366 </td>
   <td align="right"> 0.3295 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/NCIPlots9.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots9.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Validated transcriptional targets of deltaNp63 isoforms </td>
   <td align="right"> 53 </td>
   <td align="right"> 20 </td>
   <td align="right"> 0.0005 </td>
   <td align="right"> 1.327 </td>
   <td align="right"> 0.741 </td>
   <td align="right"> 0.0036 </td>
   <td align="right"> 0.0424 </td>
   <td align="right"> 0.4245 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/NCIPlots10.png target="_blank"> <img height=50 width=50 src=figure/NCIPlots10.png></a> </td>
  </tr>
</tbody>
</table>


<script type="text/javascript" charset="utf-8">
    $(document).ready( function () {
    $('#nci_table').dataTable({
          "sDom": 'T<"clear">lfrtip',
          "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
          },
      "aoColumns": [ 
      null,
  		{"bVisible": false},
      {"bVisible": false},
		  {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"sType": "scientific"},
      {"bVisible": false},
      null,
      null
		]        
    });
    });
</script>

******

### Reactome



<table id="reactome_table">
 <thead>
  <tr>
   <th align="left"> Name </th>
   <th align="right"> pSize </th>
   <th align="right"> NDE </th>
   <th align="right"> pNDE </th>
   <th align="right"> tA </th>
   <th align="right"> pPERT </th>
   <th align="right"> pG </th>
   <th align="right"> pGFdr </th>
   <th align="right"> pGFWER </th>
   <th align="left"> Status </th>
   <th align="left"> Image </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> RNA Polymerase I Promoter Opening </td>
   <td align="right"> 28 </td>
   <td align="right"> 23 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -136.421 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots1.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots1.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Cell Cycle, Mitotic </td>
   <td align="right"> 292 </td>
   <td align="right"> 112 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 715.417 </td>
   <td align="right"> 0.328 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots2.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots2.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Eukaryotic Translation Elongation </td>
   <td align="right"> 83 </td>
   <td align="right"> 41 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -110.097 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots3.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots3.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Peptide chain elongation </td>
   <td align="right"> 80 </td>
   <td align="right"> 40 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -106.140 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots4.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots4.png></a> </td>
  </tr>
  <tr>
   <td align="left"> DNA Replication </td>
   <td align="right"> 182 </td>
   <td align="right"> 77 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -608.598 </td>
   <td align="right"> 0.505 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots5.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots5.png></a> </td>
  </tr>
  <tr>
   <td align="left"> RNA Polymerase I Transcription </td>
   <td align="right"> 53 </td>
   <td align="right"> 29 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -139.543 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots6.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots6.png></a> </td>
  </tr>
  <tr>
   <td align="left"> RNA Polymerase I Promoter Clearance </td>
   <td align="right"> 51 </td>
   <td align="right"> 28 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -139.751 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots7.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots7.png></a> </td>
  </tr>
  <tr>
   <td align="left"> M Phase </td>
   <td align="right"> 82 </td>
   <td align="right"> 44 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 3.744 </td>
   <td align="right"> 0.113 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots8.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots8.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Mitotic M-M/G1 phases </td>
   <td align="right"> 162 </td>
   <td align="right"> 68 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 26.345 </td>
   <td align="right"> 0.104 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots9.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots9.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Meiotic Recombination </td>
   <td align="right"> 38 </td>
   <td align="right"> 22 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 73.706 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots10.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots10.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Mitotic Prometaphase </td>
   <td align="right"> 78 </td>
   <td align="right"> 41 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 3.744 </td>
   <td align="right"> 0.120 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots11.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots11.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Translation </td>
   <td align="right"> 114 </td>
   <td align="right"> 43 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -106.936 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots12.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots12.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Chromosome Maintenance </td>
   <td align="right"> 67 </td>
   <td align="right"> 36 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 4.491 </td>
   <td align="right"> 0.645 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots13.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots13.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Regulation of gene expression in beta cells </td>
   <td align="right"> 86 </td>
   <td align="right"> 42 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -1.418 </td>
   <td align="right"> 0.509 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots14.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots14.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Regulation of beta-cell development </td>
   <td align="right"> 89 </td>
   <td align="right"> 42 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -1.418 </td>
   <td align="right"> 0.547 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots15.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots15.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Activation of the pre-replicative complex </td>
   <td align="right"> 29 </td>
   <td align="right"> 19 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 28.518 </td>
   <td align="right"> 0.012 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots16.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots16.png></a> </td>
  </tr>
  <tr>
   <td align="left"> RNA Polymerase I, RNA Polymerase III, and Mitochondrial Transcription </td>
   <td align="right"> 89 </td>
   <td align="right"> 31 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> -143.207 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots17.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots17.png></a> </td>
  </tr>
  <tr>
   <td align="left"> G2/M Checkpoints </td>
   <td align="right"> 39 </td>
   <td align="right"> 24 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 4.519 </td>
   <td align="right"> 0.547 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots18.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots18.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Influenza Infection </td>
   <td align="right"> 132 </td>
   <td align="right"> 51 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 7.843 </td>
   <td align="right"> 0.396 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0001 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots19.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots19.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Insulin Synthesis and Processing </td>
   <td align="right"> 113 </td>
   <td align="right"> 45 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -10.410 </td>
   <td align="right"> 0.197 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0001 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots20.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots20.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Influenza Life Cycle </td>
   <td align="right"> 128 </td>
   <td align="right"> 48 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 8.114 </td>
   <td align="right"> 0.383 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0004 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots21.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots21.png></a> </td>
  </tr>
  <tr>
   <td align="left"> GTP hydrolysis and joining of the 60S ribosomal subunit </td>
   <td align="right"> 100 </td>
   <td align="right"> 41 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 1.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0005 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots22.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots22.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Activation of ATR in response to replication stress </td>
   <td align="right"> 34 </td>
   <td align="right"> 20 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 4.118 </td>
   <td align="right"> 0.530 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0005 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots23.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots23.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Telomere Maintenance </td>
   <td align="right"> 47 </td>
   <td align="right"> 24 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 4.848 </td>
   <td align="right"> 0.577 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0011 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots24.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots24.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Cap-dependent Translation Initiation </td>
   <td align="right"> 107 </td>
   <td align="right"> 41 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -2.889 </td>
   <td align="right"> 0.729 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0028 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots25.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots25.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Eukaryotic Translation Initiation </td>
   <td align="right"> 107 </td>
   <td align="right"> 41 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -2.294 </td>
   <td align="right"> 0.777 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0030 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots26.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots26.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Gene Expression </td>
   <td align="right"> 363 </td>
   <td align="right"> 75 </td>
   <td align="right"> 0.1052 </td>
   <td align="right"> -87.168 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0035 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots27.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots27.png></a> </td>
  </tr>
  <tr>
   <td align="left"> G2/M Transition </td>
   <td align="right"> 80 </td>
   <td align="right"> 30 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 36.418 </td>
   <td align="right"> 0.024 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 0.0046 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots28.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots28.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Transcription </td>
   <td align="right"> 173 </td>
   <td align="right"> 37 </td>
   <td align="right"> 0.1444 </td>
   <td align="right"> -145.009 </td>
   <td align="right"> 0.000 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 0.0047 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots29.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots29.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Fanconi Anemia pathway </td>
   <td align="right"> 22 </td>
   <td align="right"> 13 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 27.124 </td>
   <td align="right"> 0.075 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0003 </td>
   <td align="right"> 0.0093 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots30.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots30.png></a> </td>
  </tr>
  <tr>
   <td align="left"> DNA strand elongation </td>
   <td align="right"> 30 </td>
   <td align="right"> 17 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> -1.841 </td>
   <td align="right"> 0.750 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0003 </td>
   <td align="right"> 0.0107 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots31.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots31.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Centrosome maturation </td>
   <td align="right"> 68 </td>
   <td align="right"> 25 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 10.377 </td>
   <td align="right"> 0.013 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 0.0005 </td>
   <td align="right"> 0.0148 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots32.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots32.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Recruitment of mitotic centrosome proteins and complexes </td>
   <td align="right"> 68 </td>
   <td align="right"> 25 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 6.570 </td>
   <td align="right"> 0.038 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0012 </td>
   <td align="right"> 0.0400 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots33.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots33.png></a> </td>
  </tr>
  <tr>
   <td align="left"> DNA replication initiation </td>
   <td align="right"> 6 </td>
   <td align="right"> 6 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 7.168 </td>
   <td align="right"> 0.236 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0013 </td>
   <td align="right"> 0.0442 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots34.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots34.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Telomere C-strand synthesis initiation </td>
   <td align="right"> 6 </td>
   <td align="right"> 6 </td>
   <td align="right"> 0.0000 </td>
   <td align="right"> 6.864 </td>
   <td align="right"> 0.253 </td>
   <td align="right"> 0.0001 </td>
   <td align="right"> 0.0013 </td>
   <td align="right"> 0.0471 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots35.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots35.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Loss of proteins required for interphase microtubule organization from the centrosome </td>
   <td align="right"> 60 </td>
   <td align="right"> 23 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 3.744 </td>
   <td align="right"> 0.088 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 0.0021 </td>
   <td align="right"> 0.0760 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots36.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots36.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Loss of Nlp from mitotic centrosomes </td>
   <td align="right"> 60 </td>
   <td align="right"> 23 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 3.744 </td>
   <td align="right"> 0.097 </td>
   <td align="right"> 0.0002 </td>
   <td align="right"> 0.0022 </td>
   <td align="right"> 0.0831 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots37.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots37.png></a> </td>
  </tr>
  <tr>
   <td align="left"> DNA Repair </td>
   <td align="right"> 102 </td>
   <td align="right"> 30 </td>
   <td align="right"> 0.0032 </td>
   <td align="right"> 43.489 </td>
   <td align="right"> 0.009 </td>
   <td align="right"> 0.0003 </td>
   <td align="right"> 0.0038 </td>
   <td align="right"> 0.1442 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots38.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots38.png></a> </td>
  </tr>
  <tr>
   <td align="left"> G1/S Transition </td>
   <td align="right"> 100 </td>
   <td align="right"> 32 </td>
   <td align="right"> 0.0005 </td>
   <td align="right"> 35.261 </td>
   <td align="right"> 0.112 </td>
   <td align="right"> 0.0006 </td>
   <td align="right"> 0.0067 </td>
   <td align="right"> 0.2630 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots39.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots39.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Phosphorylation of Emi1 </td>
   <td align="right"> 6 </td>
   <td align="right"> 5 </td>
   <td align="right"> 0.0010 </td>
   <td align="right"> 8.991 </td>
   <td align="right"> 0.062 </td>
   <td align="right"> 0.0006 </td>
   <td align="right"> 0.0069 </td>
   <td align="right"> 0.2775 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots40.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots40.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Mitotic G1-G1/S phases </td>
   <td align="right"> 111 </td>
   <td align="right"> 34 </td>
   <td align="right"> 0.0008 </td>
   <td align="right"> 38.523 </td>
   <td align="right"> 0.116 </td>
   <td align="right"> 0.0010 </td>
   <td align="right"> 0.0104 </td>
   <td align="right"> 0.4263 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots41.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots41.png></a> </td>
  </tr>
  <tr>
   <td align="left"> M/G1 Transition </td>
   <td align="right"> 80 </td>
   <td align="right"> 24 </td>
   <td align="right"> 0.0060 </td>
   <td align="right"> 28.393 </td>
   <td align="right"> 0.027 </td>
   <td align="right"> 0.0016 </td>
   <td align="right"> 0.0164 </td>
   <td align="right"> 0.6868 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots42.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots42.png></a> </td>
  </tr>
  <tr>
   <td align="left"> DNA Replication Pre-Initiation </td>
   <td align="right"> 80 </td>
   <td align="right"> 24 </td>
   <td align="right"> 0.0060 </td>
   <td align="right"> 27.681 </td>
   <td align="right"> 0.037 </td>
   <td align="right"> 0.0021 </td>
   <td align="right"> 0.0212 </td>
   <td align="right"> 0.9106 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots43.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots43.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Telomere C-strand (Lagging Strand) Synthesis </td>
   <td align="right"> 22 </td>
   <td align="right"> 11 </td>
   <td align="right"> 0.0006 </td>
   <td align="right"> 6.167 </td>
   <td align="right"> 0.444 </td>
   <td align="right"> 0.0026 </td>
   <td align="right"> 0.0253 </td>
   <td align="right"> 1.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots44.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots44.png></a> </td>
  </tr>
  <tr>
   <td align="left"> EGFR downregulation </td>
   <td align="right"> 23 </td>
   <td align="right"> 5 </td>
   <td align="right"> 0.4012 </td>
   <td align="right"> -31.712 </td>
   <td align="right"> 0.001 </td>
   <td align="right"> 0.0035 </td>
   <td align="right"> 0.0341 </td>
   <td align="right"> 1.0000 </td>
   <td align="left"> Inhibited </td>
   <td align="left"> <a href = figure/ReactomePlots45.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots45.png></a> </td>
  </tr>
  <tr>
   <td align="left"> Cyclin A/B1 associated events during G2/M transition </td>
   <td align="right"> 14 </td>
   <td align="right"> 7 </td>
   <td align="right"> 0.0064 </td>
   <td align="right"> 26.128 </td>
   <td align="right"> 0.099 </td>
   <td align="right"> 0.0053 </td>
   <td align="right"> 0.0499 </td>
   <td align="right"> 1.0000 </td>
   <td align="left"> Activated </td>
   <td align="left"> <a href = figure/ReactomePlots46.png target="_blank"> <img height=50 width=50 src=figure/ReactomePlots46.png></a> </td>
  </tr>
</tbody>
</table>


<script type="text/javascript" charset="utf-8">
    $(document).ready( function () {
    $('#reactome_table').dataTable({
          "sDom": 'T<"clear">lfrtip',
          "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
          },
      "aoColumns": [ 
      null,
  		{"bVisible": false},
      {"bVisible": false},
		  {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"sType": "scientific"},
      {"bVisible": false},
      null,
      null
		]          
    });
    });
</script>

******

## Fusion Gene Candidates

Fusion gene candidated are assessed using [Fusioncatcher](https://code.google.com/p/fusioncatcher/) <code class="knitr inline">(Edgren, Murumagi, Kangaspeska, Nicorici, Hongisto, Kleivi, Rye, Nyberg, Wolf, Borresen-Dale, and Kallioniemi, 2011)</code>. The oncogenicc potential of the fusion candidates is predicted using [OncoFuse](http://www.unav.es/genetica/oncofuse.html) <code class="knitr inline">(Shugay, de
Mendibil, Vizmanos, and Novo, 2013)</code>.

[1] "No FusionCatcher results available"


<script type="text/javascript" charset="utf-8">
    $(document).ready(function() {
    var stdTable1 = $('#fusion_table').dataTable({
    "sDom": 'T<"clear">lfrtip',
    "oTableTools": {
          "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
          "aButtons": [
                "copy",
                "csv"
                ]
    },
    "aoColumns": [ 
  		null,
			null,
      null,
		  {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      {"bVisible": false},
      null,
      null      
		], 
    //"iDisplayLength": 10,
    //"bPaginate": true,
    //"bProcessing": true,
    //"sScrollXInner": "110%",
    });
    
    //var tableId = 'fusion_table';
    //$('<div style="width: 100%; overflow: auto"></div>').append($('#' + tableId)).insertAfter($('#' + tableId + '_wrapper div').first())
    });
</script>

******

## SNPs

* **CHROM** The mutation’s chromosome.
* **POS** The mutation’s start position.
* **rsID**
* **Gene** The HUGO symbol of the gene.
* **Intogen_Driver**
* **SnpEff_Effect**
* **ClinVar_Rating**
* **ClinVar_ClinicalSignificance**
* **COSMIC_histology**
* **COSMIC_reference**
* **PharmGbk_Reaction**
* **PharmGbk_Drugs**
* **DrugBank_drug**
* **DrugBank_reaction**
* **DrugBank_reference**

<table id="SNPs_table">
 <thead>
  <tr>
   <th align="left"> CHROM </th>
   <th align="left"> POS </th>
   <th align="left"> rsID </th>
   <th align="left"> Gene </th>
   <th align="left"> Diff_Exp </th>
   <th align="left"> Intogen_driver </th>
   <th align="right"> CADD_PHRED </th>
   <th align="left"> SnpEff_Effect </th>
   <th align="left"> ClinVar_Rating </th>
   <th align="left"> ClinVar_ClinicalSignificance </th>
   <th align="left"> COSMIC_histology </th>
   <th align="left"> COSMIC_reference </th>
   <th align="left"> PharmGkb_Reaction </th>
   <th align="left"> PharmGkb_Drugs </th>
   <th align="left"> DrugBank_drug </th>
   <th align="left"> DrugBank_reaction </th>
   <th align="left"> DrugBank_reference </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 110465998 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   .    ' target="_blank" > .   .     </a> </td>
   <td align="left"> CSF1 </td>
   <td align="left"> Yes(CSF1) </td>
   <td align="left"> No </td>
   <td align="right"> 15.020 </td>
   <td align="left"> CSF1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGg/cAg, R252Q;   CSF1: MODIFIER, INTRON, cGg/cAg, R252Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> oesophagus carcinoma adenocarcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=23525077' target="_blank" > 23525077 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 113231661 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MOV10 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.820 </td>
   <td align="left"> MOV10: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGc/cAc, R81H;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 119683231 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=139548132' target="_blank" > rs139548132 </a> </td>
   <td align="left"> WARS2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 17.020 </td>
   <td align="left"> WARS2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgg/Ggg, W13G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 145109583 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2794053' target="_blank" > rs2794053 </a> </td>
   <td align="left"> SEC22B; NBPF9 </td>
   <td align="left"> Yes(NBPF9) </td>
   <td align="left"> No </td>
   <td align="right"> 16.160 </td>
   <td align="left"> SEC22B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCg/aAg, T82K;   NBPF9: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 145112414 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2590131' target="_blank" > rs2590131 </a> </td>
   <td align="left"> SEC22B; NBPF9 </td>
   <td align="left"> Yes(NBPF9) </td>
   <td align="left"> No </td>
   <td align="right"> 17.570 </td>
   <td align="left"> SEC22B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgt/Cgt, C130R;   NBPF9: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 145115810 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2655557' target="_blank" > rs2655557 </a> </td>
   <td align="left"> SEC22B; NBPF9 </td>
   <td align="left"> Yes(NBPF9) </td>
   <td align="left"> No </td>
   <td align="right"> 17.300 </td>
   <td align="left"> SEC22B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAt/cGt, H190R;   NBPF9: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 150933093 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SETDB1; CERS2 </td>
   <td align="left"> Yes(SETDB1) </td>
   <td align="left"> No </td>
   <td align="right"> 23.000 </td>
   <td align="left"> SETDB1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gAg/gTg, E852V;   CERS2: MODIFIER, DOWNSTREAM, 4556;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 156449435 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MEF2D </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 18.680 </td>
   <td align="left"> MEF2D: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cca/Tca, P184S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 157497650 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   .    ' target="_blank" > .   .     </a> </td>
   <td align="left"> FCRL5 </td>
   <td align="left"> Yes(FCRL5) </td>
   <td align="left"> No </td>
   <td align="right"> 0.856 </td>
   <td align="left"> FCRL5: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ccc/Tcc, P573S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> skin trunk malignant_melanoma nodular </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=22197931' target="_blank" > 22197931 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 161068548 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> KLHDC9; PFDN2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.210 </td>
   <td align="left"> KLHDC9: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ggc/Cgc, G75R;   KLHDC9: MODIFIER, EXON, Ggc/Cgc, G75R;   PFDN2: MODIFIER, DOWNSTREAM, 1798;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 197111733 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ASPM </td>
   <td align="left"> Yes(ASPM) </td>
   <td align="left"> No </td>
   <td align="right"> 17.290 </td>
   <td align="left"> ASPM: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cCa/cTa, P550L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 200960181 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> KIF21B </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.770 </td>
   <td align="left"> KIF21B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tcg/Ccg, S851P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 201868719 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> LMOD1 </td>
   <td align="left"> Yes(LMOD1) </td>
   <td align="left"> No </td>
   <td align="right"> 22.000 </td>
   <td align="left"> LMOD1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, caG/caC, Q474H;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 205028322 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> CNTN2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 19.690 </td>
   <td align="left"> CNTN2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tca/Cca, S200P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 22198822 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=145185113' target="_blank" > rs145185113 </a> </td>
   <td align="left"> HSPG2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 21.400 </td>
   <td align="left"> HSPG2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aac/Gac, N1360D;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 233802526 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> KCNK1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 5.808 </td>
   <td align="left"> KCNK1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atc/Gtc, I181V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 33065947 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=704886' target="_blank" > rs704886 </a> </td>
   <td align="left"> ZBTB8A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 7.825 </td>
   <td align="left"> ZBTB8A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gGt/gCt, G418A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 40533266 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=11207440' target="_blank" > rs11207440 </a> </td>
   <td align="left"> CAP1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 9.400 </td>
   <td align="left"> CAP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgt/Ggt, C229G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 40533287 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=6665926' target="_blank" > rs6665926 </a> </td>
   <td align="left"> CAP1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 7.939 </td>
   <td align="left"> CAP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgc/Ggc, C236G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 51873879 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=150935900' target="_blank" > rs150935900 </a> </td>
   <td align="left"> EPS15 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 15.190 </td>
   <td align="left"> EPS15: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, agT/agG, S153R, S467R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 70641538 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> LRRC40 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 15.820 </td>
   <td align="left"> LRRC40: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cTa/cCa, L311P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 93160880 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=143611208' target="_blank" > rs143611208 </a> </td>
   <td align="left"> EVI5 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 17.180 </td>
   <td align="left"> EVI5: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTc/aCc, I343T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 1 </td>
   <td align="left"> 97544603 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> DPYD </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.110 </td>
   <td align="left"> DPYD: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aaa/Gaa, K1003E;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 121355952 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TIAL1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 19.060 </td>
   <td align="left"> TIAL1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGg/cTg, R9L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 126136300 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> NKX1-2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 19.540 </td>
   <td align="left"> NKX1-2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ttc/Ctc, F211L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 13699338 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=199968440' target="_blank" > rs199968440 </a> </td>
   <td align="left"> FRMD4A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 20.400 </td>
   <td align="left"> FRMD4A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Acc/Ccc, T751P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> liver carcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=   ' target="_blank" >     </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 28824522 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> WAC; WAC-AS1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 22.500 </td>
   <td align="left"> WAC: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aGt/aTt, S37I;   WAC: MODIFIER, EXON, aGt/aTt, S37I;   WAC-AS1: MODIFIER, UPSTREAM, 3239;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 31137611 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=149636280' target="_blank" > rs149636280 </a> </td>
   <td align="left"> ZNF438 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 20.100 </td>
   <td align="left"> ZNF438: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtg/Atg, V526M, V565M, V575M;   ZNF438: MODIFIER, EXON, Gtg/Atg, V526M, V565M, V575M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 46999151 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=3127820' target="_blank" > rs3127820 </a> </td>
   <td align="left"> GPRIN2 </td>
   <td align="left"> Yes(GPRIN2) </td>
   <td align="left"> No </td>
   <td align="right"> 1.236 </td>
   <td align="left"> GPRIN2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgg/Cgg, W91R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 5789309 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> FAM208B </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 3.663 </td>
   <td align="left"> FAM208B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gcc/Tcc, A1309S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 10 </td>
   <td align="left"> 70730032 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> DDX21 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.130 </td>
   <td align="left"> DDX21: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Caa/Aaa, Q370K, Q438K;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 108163417 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ATM </td>
   <td align="left"> Yes(ATM) </td>
   <td align="left"> Yes(ATM)   </td>
   <td align="right"> 9.608 </td>
   <td align="left"> ATM: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAg/cGg, Q1503R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 10822029 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> EIF4G2; SNORD97 </td>
   <td align="left"> Yes(SNORD97) </td>
   <td align="left"> No </td>
   <td align="right"> 1.094 </td>
   <td align="left"> EIF4G2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atc/Gtc, I566V, I604V;   SNORD97: MODIFIER, DOWNSTREAM, 985;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 118773050 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=148999907' target="_blank" > rs148999907 </a> </td>
   <td align="left"> BCL9L </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 12.240 </td>
   <td align="left"> BCL9L: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ccc/Tcc, P468S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 45937267 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=10742772' target="_blank" > rs10742772 </a> </td>
   <td align="left"> PEX16 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 6.552 </td>
   <td align="left"> PEX16: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtc/Atc, V116I;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 57512055 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> BTBD18; C11orf31; TMX2; TMX2-CTNND1 </td>
   <td align="left"> Yes(BTBD18) </td>
   <td align="left"> No </td>
   <td align="right"> 14.880 </td>
   <td align="left"> BTBD18: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gag/Aag, E564K;   C11orf31: MODIFIER, DOWNSTREAM, 1172;   TMX2: MODIFIER, DOWNSTREAM, 3610;   TMX2-CTNND1: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 62400116 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=138726440' target="_blank" > rs138726440 </a> </td>
   <td align="left"> GANAB </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 7.701 </td>
   <td align="left"> GANAB: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aAc/aGc, N306S, N328S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 67258391 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4930199' target="_blank" > rs4930199 </a> </td>
   <td align="left"> AIP; PITPNM1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 1.754 </td>
   <td align="left"> AIP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAg/cGg, Q307R;   PITPNM1: MODIFIER, DOWNSTREAM, 848;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 67957518 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2512606' target="_blank" > rs2512606 </a> </td>
   <td align="left"> SUV420H1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 8.504 </td>
   <td align="left"> SUV420H1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTc/aAc, I9N;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 73076531 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=141322412' target="_blank" > rs141322412 </a> </td>
   <td align="left"> ARHGEF17 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 22.700 </td>
   <td align="left"> ARHGEF17: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtg/Atg, V1883M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 74617391 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> XRRA1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 2.616 </td>
   <td align="left"> XRRA1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAt/cGt, H291R, H58R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 11 </td>
   <td align="left"> 74904362 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=1621378' target="_blank" > rs1621378 </a> </td>
   <td align="left"> SLCO2B1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.013 </td>
   <td align="left"> SLCO2B1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTc/aCc, I248T, I370T, I392T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 12 </td>
   <td align="left"> 11183642 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=73049067' target="_blank" > rs73049067 </a> </td>
   <td align="left"> TAS2R31; PRH1-PRR4 </td>
   <td align="left"> Yes(TAS2R31) </td>
   <td align="left"> No </td>
   <td align="right"> 12.300 </td>
   <td align="left"> TAS2R31: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cTt/cCt, L98P;   PRH1-PRR4: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 12 </td>
   <td align="left"> 49443911 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MLL2 </td>
   <td align="left"> No </td>
   <td align="left"> Yes(MLL2)   </td>
   <td align="right"> 11.690 </td>
   <td align="left"> MLL2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ctc/Ttc, L1154F;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 12 </td>
   <td align="left"> 51065113 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> DIP2B </td>
   <td align="left"> Yes(DIP2B) </td>
   <td align="left"> No </td>
   <td align="right"> 19.990 </td>
   <td align="left"> DIP2B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCg/aTg, T191M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 13 </td>
   <td align="left"> 24443512 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=149856612' target="_blank" > rs149856612 </a> </td>
   <td align="left"> MIPEP </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.520 </td>
   <td align="left"> MIPEP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Agc/Cgc, S288R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 13 </td>
   <td align="left"> 45147617 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TSC22D1; TSC22D1-AS1 </td>
   <td align="left"> Yes(TSC22D1) </td>
   <td align="left"> No </td>
   <td align="right"> 6.486 </td>
   <td align="left"> TSC22D1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aAt/aGt, N865S;   TSC22D1: MODIFIER, INTRAGENIC, aAt/aGt, N865S;   TSC22D1: MODIFIER, INTRON, aAt/aGt, N865S;   TSC22D1-AS1: MODIFIER, UPSTREAM, 2415;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 13 </td>
   <td align="left"> 53035925 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=7335867' target="_blank" > rs7335867 </a> </td>
   <td align="left"> CKAP2 </td>
   <td align="left"> Yes(CKAP2) </td>
   <td align="left"> No </td>
   <td align="right"> 0.000 </td>
   <td align="left"> CKAP2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ata/Gta, I322V, I323V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 13 </td>
   <td align="left"> 53421242 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> PCDH8 </td>
   <td align="left"> Yes(PCDH8) </td>
   <td align="left"> No </td>
   <td align="right"> 15.010 </td>
   <td align="left"> PCDH8: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gcg/Acg, A444T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 104493079 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=144686761' target="_blank" > rs144686761 </a> </td>
   <td align="left"> TDRD9 </td>
   <td align="left"> Yes(TDRD9) </td>
   <td align="left"> No </td>
   <td align="right"> 15.660 </td>
   <td align="left"> TDRD9: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgg/Cgg, W1029R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 38724041 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> CLEC14A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 15.040 </td>
   <td align="left"> CLEC14A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tCt/tTt, S396F;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 45700437 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MIS18BP1 </td>
   <td align="left"> Yes(MIS18BP1) </td>
   <td align="left"> No </td>
   <td align="right"> 0.003 </td>
   <td align="left"> MIS18BP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ata/Gta, I501V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 55236828 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SAMD4A </td>
   <td align="left"> Yes(SAMD4A) </td>
   <td align="left"> No </td>
   <td align="right"> 22.700 </td>
   <td align="left"> SAMD4A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aca/Gca, T126A, T447A, T535A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 75936090 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> JDP2 </td>
   <td align="left"> Yes(JDP2) </td>
   <td align="left"> No </td>
   <td align="right"> 27.400 </td>
   <td align="left"> JDP2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAc/cGc, H135R, H146R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 90755147 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=150219025' target="_blank" > rs150219025 </a> </td>
   <td align="left"> NRDE2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 4.245 </td>
   <td align="left"> NRDE2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ccc/Gcc, P858A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 14 </td>
   <td align="left"> 94847415 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=6647' target="_blank" > rs6647 </a> </td>
   <td align="left"> SERPINA1 </td>
   <td align="left"> Yes(SERPINA1) </td>
   <td align="left"> No </td>
   <td align="right"> 7.506 </td>
   <td align="left"> SERPINA1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTg/gCg, V237A;  </td>
   <td align="left"> 2 </td>
   <td align="left"> Pathogenic;other </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 39881529 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> THBS1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.041 </td>
   <td align="left"> THBS1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtc/Ctc, V634L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 43299402 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=142285781' target="_blank" > rs142285781 </a> </td>
   <td align="left"> UBR1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 24.300 </td>
   <td align="left"> UBR1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCg/aTg, T1097M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> large_intestine caecum carcinoma adenocarcinoma Stage:I </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=   ' target="_blank" >     </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 48443699 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2470103' target="_blank" > rs2470103 </a> </td>
   <td align="left"> MYEF2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 8.460 </td>
   <td align="left"> MYEF2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAa/cGa, Q426R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 48807637 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4775765' target="_blank" > rs4775765 </a> </td>
   <td align="left"> FBN1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 8.028 </td>
   <td align="left"> FBN1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tGc/tAc, C472Y;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 49048567 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=201342438' target="_blank" > rs201342438 </a> </td>
   <td align="left"> CEP152 </td>
   <td align="left"> Yes(CEP152) </td>
   <td align="left"> No </td>
   <td align="right"> 19.230 </td>
   <td align="left"> CEP152: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgg/Cgg, W960R;  </td>
   <td align="left"> 1 </td>
   <td align="left"> Uncertain significance </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 56386577 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=7170589' target="_blank" > rs7170589 </a> </td>
   <td align="left"> RFX7 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 2.277 </td>
   <td align="left"> RFX7: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aat/Cat, N1117H;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 64980951 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=182312681' target="_blank" > rs182312681 </a> </td>
   <td align="left"> OAZ2; ZNF609 </td>
   <td align="left"> Yes(OAZ2) </td>
   <td align="left"> No </td>
   <td align="right"> 6.357 </td>
   <td align="left"> OAZ2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gat/Cat, D175H;   ZNF609: MODIFIER, DOWNSTREAM, 2685;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 75664472 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=142724257' target="_blank" > rs142724257 </a> </td>
   <td align="left"> SIN3A; MAN2C1 </td>
   <td align="left"> Yes(MAN2C1) </td>
   <td align="left"> No </td>
   <td align="right"> 16.680 </td>
   <td align="left"> SIN3A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cgt/Tgt, R1224C;   MAN2C1: MODIFIER, UPSTREAM, 3504;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 15 </td>
   <td align="left"> 93588336 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4238485' target="_blank" > rs4238485 </a> </td>
   <td align="left"> RGMA </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 1.101 </td>
   <td align="left"> RGMA: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gaT/gaG, D399E, D415E, D423E;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 20809122 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=76921366' target="_blank" > rs76921366 </a> </td>
   <td align="left"> ERI2; ACSM3 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 3.249 </td>
   <td align="left"> ERI2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCa/aAa, T667K;   ACSM3: MODIFIER, DOWNSTREAM, 643;   ERI2: MODIFIER, INTRON, aCa/aAa, T667K;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 2213989 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TRAF7 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 17.530 </td>
   <td align="left"> TRAF7: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gAc/gTc, D23V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 24580865 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> RBBP6 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 18.730 </td>
   <td align="left"> RBBP6: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gag/Cag, E918Q, E952Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 2815196 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SRRM2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.850 </td>
   <td align="left"> SRRM2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aGa/aCa, R1556T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 29818703 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MAZ; KIF22; PRRT2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 17.470 </td>
   <td align="left"> MAZ: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, ttC/ttG, F176L, F199L;   KIF22: MODIFIER, DOWNSTREAM, 1997;   MAZ: MODIFIER, INTRON, ttC/ttG, F176L, F199L;   PRRT2: MODIFIER, UPSTREAM, 4706;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 29853318 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=146411255' target="_blank" > rs146411255 </a> </td>
   <td align="left"> MVP </td>
   <td align="left"> Yes(MVP) </td>
   <td align="left"> No </td>
   <td align="right"> 21.200 </td>
   <td align="left"> MVP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cgt/Tgt, R507C;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 30750116 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=149217909' target="_blank" > rs149217909 </a> </td>
   <td align="left"> SRCAP; LOC100862671 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.303 </td>
   <td align="left"> SRCAP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ctt/Ttt, L2919F;   LOC100862671: MODIFIER, UPSTREAM, 1396;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 3433570 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZSCAN32 </td>
   <td align="left"> Yes(ZSCAN32) </td>
   <td align="left"> No </td>
   <td align="right"> 13.780 </td>
   <td align="left"> ZSCAN32: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cCa/cTa, P247L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 3714393 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TRAP1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 9.371 </td>
   <td align="left"> TRAP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tCa/tTa, S431L, S484L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 3724401 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TRAP1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 10.550 </td>
   <td align="left"> TRAP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCc/aTc, T275I, T328I;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 427463 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TMEM8A; LOC100134368 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 4.439 </td>
   <td align="left"> TMEM8A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gCc/gAc, A141D;   LOC100134368: MODIFIER, UPSTREAM, 4778;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 4861709 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2085329' target="_blank" > rs2085329 </a> </td>
   <td align="left"> GLYR1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 4.601 </td>
   <td align="left"> GLYR1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, caC/caG, H459Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 58030654 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZNF319; USB1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.740 </td>
   <td align="left"> ZNF319: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cgg/Tgg, R506W;   USB1: MODIFIER, UPSTREAM, 4623;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 16 </td>
   <td align="left"> 731848 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> STUB1; JMJD8; RHBDL1; WDR24 </td>
   <td align="left"> Yes(STUB1) </td>
   <td align="left"> No </td>
   <td align="right"> 15.490 </td>
   <td align="left"> STUB1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cgg/Tgg, R194W;   JMJD8: MODIFIER, UTR_3_PRIME, 946;   RHBDL1: MODIFIER, DOWNSTREAM, 3581;   WDR24: MODIFIER, DOWNSTREAM, 2854;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 10608663 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ADPRM </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 9.274 </td>
   <td align="left"> ADPRM: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gaG/gaC, E140D;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 17398822 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> RASD1; MED9 </td>
   <td align="left"> Yes(RASD1) </td>
   <td align="left"> No </td>
   <td align="right"> 21.700 </td>
   <td align="left"> RASD1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtg/Atg, 20, V155M;   MED9: MODIFIER, DOWNSTREAM, 2288;   RASD1: MODIFIER, UTR_3_PRIME, Gtg/Atg, 20, V155M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 3514028 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=224496' target="_blank" > rs224496 </a> </td>
   <td align="left"> SHPK; TRPV1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.007 </td>
   <td align="left"> SHPK: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gaC/gaG, D421E;   TRPV1: MODIFIER, UPSTREAM, 1323;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 46929982 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=17849804' target="_blank" > rs17849804 </a> </td>
   <td align="left"> CALCOCO2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.670 </td>
   <td align="left"> CALCOCO2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Act/Gct, T201A, T231A, T273A, T294A, T297A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 62517603 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=200923114' target="_blank" > rs200923114 </a> </td>
   <td align="left"> CEP95 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 2.841 </td>
   <td align="left"> CEP95: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cct/Tct, P225S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 71193103 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> COG1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 15.050 </td>
   <td align="left"> COG1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gcc/Acc, A209T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 17 </td>
   <td align="left"> 7417086 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4344809' target="_blank" > rs4344809 </a> </td>
   <td align="left"> POLR2A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 3.541 </td>
   <td align="left"> POLR2A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gcc/Acc, A1835T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 18 </td>
   <td align="left"> 28934581 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   .    ' target="_blank" > .   .     </a> </td>
   <td align="left"> DSG1 </td>
   <td align="left"> Yes(DSG1) </td>
   <td align="left"> No </td>
   <td align="right"> 8.633 </td>
   <td align="left"> DSG1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gta/Ata, V808I;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> lung carcinoma adenocarcinoma Stage:IV </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=22980975' target="_blank" > 22980975 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 10461586 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=55886939' target="_blank" > rs55886939 </a> </td>
   <td align="left"> TYK2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.890 </td>
   <td align="left"> TYK2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gAg/gGg, E1163G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 12298712 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZNF136 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 18.380 </td>
   <td align="left"> ZNF136: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aag/Gag, K507E;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 13065096 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=145414355' target="_blank" > rs145414355 </a> </td>
   <td align="left"> GADD45GIP1; RAD23A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 6.978 </td>
   <td align="left"> GADD45GIP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aag/Gag, K199E;   RAD23A: MODIFIER, DOWNSTREAM, 639;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 21992459 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZNF43 </td>
   <td align="left"> Yes(ZNF43) </td>
   <td align="left"> No </td>
   <td align="right"> 2.582 </td>
   <td align="left"> ZNF43: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAc/cGc, H121R, H127R, H136R, H62R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 38230396 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZNF573 </td>
   <td align="left"> Yes(ZNF573) </td>
   <td align="left"> No </td>
   <td align="right"> 3.901 </td>
   <td align="left"> ZNF573: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tGt/tAt, C244Y, C274Y, C330Y, C332Y;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 4652428 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TNFAIP8L1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.830 </td>
   <td align="left"> TNFAIP8L1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gag/Aag, E183K;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 19 </td>
   <td align="left"> 8645786 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=7252299' target="_blank" > rs7252299 </a> </td>
   <td align="left"> ADAMTS10; MYO1F </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 6.689 </td>
   <td align="left"> ADAMTS10: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, caT/caG, H1101Q;   MYO1F: MODIFIER, UPSTREAM, 3455;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 10074020 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TAF1B </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 18.270 </td>
   <td align="left"> TAF1B: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, caT/caA, H558Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 166900434 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   .    ' target="_blank" > .   .     </a> </td>
   <td align="left"> SCN1A </td>
   <td align="left"> Yes(SCN1A) </td>
   <td align="left"> No </td>
   <td align="right"> 16.050 </td>
   <td align="left"> SCN1A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, agC/agA, S596R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> kidney(2x) carcinoma(2x) clear_cell_renal_cell_carcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=   ' target="_blank" >     </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 179437054 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TTN; MIR548N; TTN-AS1 </td>
   <td align="left"> No </td>
   <td align="left"> Yes(TTN; MIR548N; TTN-AS1)    </td>
   <td align="right"> 14.040 </td>
   <td align="left"> TTN: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gGg/gTg, G15537V, G15662V, G15729V, G22034V, G22961V, G24602V;   MIR548N: MODIFIER, INTRON;   TTN-AS1: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 179474590 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TTN; MIR548N; TTN-AS1 </td>
   <td align="left"> No </td>
   <td align="left"> Yes(TTN; MIR548N; TTN-AS1)    </td>
   <td align="right"> 12.620 </td>
   <td align="left"> TTN: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTt/aAt, I14619N, I15546N, I17187N, I8122N, I8247N, I8314N;   MIR548N: MODIFIER, INTRON;   TTN-AS1: MODIFIER, DOWNSTREAM, 4458;   TTN-AS1: MODIFIER, INTRAGENIC, 4458;   TTN-AS1: MODIFIER, INTRON, 4458;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 242066796 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=56033464' target="_blank" > rs56033464 </a> </td>
   <td align="left"> PASK </td>
   <td align="left"> Yes(PASK) </td>
   <td align="left"> No </td>
   <td align="right"> 9.964 </td>
   <td align="left"> PASK: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Act/Gct, T477A, T512A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 26203678 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=1465878' target="_blank" > rs1465878 </a> </td>
   <td align="left"> KIF3C </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.737 </td>
   <td align="left"> KIF3C: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAg/cGg, Q370R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 32667182 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=60197615' target="_blank" > rs60197615 </a> </td>
   <td align="left"> BIRC6 </td>
   <td align="left"> Yes(BIRC6) </td>
   <td align="left"> No </td>
   <td align="right"> 10.320 </td>
   <td align="left"> BIRC6: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtt/Ctt, V1332L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 33585796 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4422143' target="_blank" > rs4422143 </a> </td>
   <td align="left"> LTBP1 </td>
   <td align="left"> No </td>
   <td align="left"> Yes(LTBP1)   </td>
   <td align="right"> 1.642 </td>
   <td align="left"> LTBP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTg/gCg, V1010A, V1052A, V1378A, V957A, V999A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 37586801 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> QPCT </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 23.600 </td>
   <td align="left"> QPCT: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ggg/Tgg, G116W;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 2 </td>
   <td align="left"> 39964026 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> THUMPD2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 1.197 </td>
   <td align="left"> THUMPD2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gCg/gTg, A454V;   THUMPD2: MODIFIER, EXON, gCg/gTg, A454V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 20 </td>
   <td align="left"> 45839444 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=142922966' target="_blank" > rs142922966 </a> </td>
   <td align="left"> ZMYND8 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 10.170 </td>
   <td align="left"> ZMYND8: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ctc/Ttc, L1124F, L1149F, L1177F;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 20 </td>
   <td align="left"> 57415780 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> GNAS; GNAS-AS1 </td>
   <td align="left"> Yes(GNAS-AS1) </td>
   <td align="left"> Yes(GNAS; GNAS-AS1)    </td>
   <td align="right"> 12.340 </td>
   <td align="left"> GNAS: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tcg/Acg, S207T;   GNAS: MODIFIER, INTRAGENIC, Tcg/Acg, S207T;   GNAS-AS1: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 20 </td>
   <td align="left"> 60888251 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> LAMA5; ADRM1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 22.000 </td>
   <td align="left"> LAMA5: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ggc/Agc, G2950S;   ADRM1: MODIFIER, DOWNSTREAM, 4333;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 21 </td>
   <td align="left"> 27348269 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> APP </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 23.400 </td>
   <td align="left"> APP: LOW, SPLICE_SITE_REGION, Cag/Aag, Q302K, Q323K, Q358K, Q377K, Q409K, Q414K, Q433K;   APP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cag/Aag, Q302K, Q323K, Q358K, Q377K, Q409K, Q414K, Q433K;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 21 </td>
   <td align="left"> 34923594 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SON </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.560 </td>
   <td align="left"> SON: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCg/aTg, T686M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 21 </td>
   <td align="left"> 34948697 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=199930883' target="_blank" > rs199930883 </a> </td>
   <td align="left"> SON; DONSON </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 2.112 </td>
   <td align="left"> SON: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, agA/agC, R2416S;   DONSON: MODIFIER, DOWNSTREAM, 1162;   SON: MODIFIER, INTRAGENIC, agA/agC, R2416S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> central_nervous_system haematopoietic_and_lymphoid_tissue large_intestine lung brain colon carcinoma(2x) glioma lymphoid_neoplasm adenocarcinoma hairy_cell_leukaemia </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=23856246(4x)   ' target="_blank" > 23856246(4x)    </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 21 </td>
   <td align="left"> 47768999 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> PCNT </td>
   <td align="left"> Yes(PCNT) </td>
   <td align="left"> No </td>
   <td align="right"> 11.750 </td>
   <td align="left"> PCNT: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAt/cGt, H369R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 22 </td>
   <td align="left"> 22300380 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=139013152' target="_blank" > rs139013152 </a> </td>
   <td align="left"> PPM1F </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.013 </td>
   <td align="left"> PPM1F: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aGt/aCt, S14T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 22 </td>
   <td align="left"> 25425282 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=116253946' target="_blank" > rs116253946 </a> </td>
   <td align="left"> KIAA1671 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 13.320 </td>
   <td align="left"> KIAA1671: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aAg/aGg, K439R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 22 </td>
   <td align="left"> 38897266 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=150253069' target="_blank" > rs150253069 </a> </td>
   <td align="left"> DDX17 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.040 </td>
   <td align="left"> DDX17: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ggt/Agt, G103S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 22 </td>
   <td align="left"> 42071121 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=200179892' target="_blank" > rs200179892 </a> </td>
   <td align="left"> NHP2L1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 20.100 </td>
   <td align="left"> NHP2L1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAc/cCc, H68P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> liver carcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=   ' target="_blank" >     </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 22 </td>
   <td align="left"> 45599808 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=138154810' target="_blank" > rs138154810 </a> </td>
   <td align="left"> KIAA0930; MIR1249 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.150 </td>
   <td align="left"> KIAA0930: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aGt/aAt, S192N, S197N;   MIR1249: MODIFIER, UPSTREAM, 2908;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 22 </td>
   <td align="left"> 50660065 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TUBGCP6; SELO </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.890 </td>
   <td align="left"> TUBGCP6: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cCg/cTg, P908L;   SELO: MODIFIER, DOWNSTREAM, 4020;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 13612381 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> FBLN2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 4.934 </td>
   <td align="left"> FBLN2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ggt/Agt, G176S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 150290165 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=200096005' target="_blank" > rs200096005 </a> </td>
   <td align="left"> EIF2A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 23.500 </td>
   <td align="left"> EIF2A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTt/gGt, V411G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> autonomic_ganglia(2x) neuroblastoma(2x)   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=23334666(2x)   ' target="_blank" > 23334666(2x)    </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 156396253 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TIPARP; TIPARP-AS1 </td>
   <td align="left"> Yes(TIPARP) </td>
   <td align="left"> No </td>
   <td align="right"> 19.130 </td>
   <td align="left"> TIPARP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gAt/gGt, D256G;   TIPARP-AS1: MODIFIER, UPSTREAM, 2751;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 180685929 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> FXR1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 32.000 </td>
   <td align="left"> FXR1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGa/cAa, R345Q, R430Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 183957200 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> VWA5B2; ALG3; MIR1224 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 12.840 </td>
   <td align="left"> VWA5B2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gag/Aag, E741K;   ALG3: MODIFIER, DOWNSTREAM, 2917;   MIR1224: MODIFIER, UPSTREAM, 1993;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 186281868 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TBCCD1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.083 </td>
   <td align="left"> TBCCD1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTg/aCg, M84T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 40503520 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=200018880' target="_blank" > rs200018880 </a> </td>
   <td align="left"> RPL14 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.019 </td>
   <td align="left"> RPL14: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Act/Gct, T149A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 52398958 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> DNAH1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.320 </td>
   <td align="left"> DNAH1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gAg/gGg, E1814G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 57131871 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> IL17RD </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 7.557 </td>
   <td align="left"> IL17RD: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, caC/caG, H620Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 58141791 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> FLNB </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.060 </td>
   <td align="left"> FLNB: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atg/Gtg, M2269V, M2282V, M2293V, M2324V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 63982092 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ATXN7 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 23.300 </td>
   <td align="left"> ATXN7: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTc/gCc, V720A, V865A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 75786080 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=144586950' target="_blank" > rs144586950 </a> </td>
   <td align="left"> ZNF717; MIR4273 </td>
   <td align="left"> Yes(ZNF717) </td>
   <td align="left"> No </td>
   <td align="right"> 0.103 </td>
   <td align="left"> ZNF717: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gaC/gaA, D898E;   MIR4273: MODIFIER, UPSTREAM, 1351;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 75786753 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=3009020' target="_blank" > rs3009020 </a> </td>
   <td align="left"> ZNF717; MIR4273 </td>
   <td align="left"> Yes(ZNF717) </td>
   <td align="left"> No </td>
   <td align="right"> 5.123 </td>
   <td align="left"> ZNF717: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGt/cAt, R674H;   MIR4273: MODIFIER, UPSTREAM, 678;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 75787809 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=142370851' target="_blank" > rs142370851 </a> </td>
   <td align="left"> ZNF717; MIR4273 </td>
   <td align="left"> Yes(ZNF717) </td>
   <td align="left"> No </td>
   <td align="right"> 4.274 </td>
   <td align="left"> ZNF717: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tAt/tGt, Y322C;   MIR4273: MODIFIER, DOWNSTREAM, 295;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 75788010 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=151311432' target="_blank" > rs151311432 </a> </td>
   <td align="left"> ZNF717; MIR4273 </td>
   <td align="left"> Yes(ZNF717) </td>
   <td align="left"> No </td>
   <td align="right"> 3.725 </td>
   <td align="left"> ZNF717: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTg/gGg, V255G;   MIR4273: MODIFIER, DOWNSTREAM, 496;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 3 </td>
   <td align="left"> 75790444 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=141084845' target="_blank" > rs141084845 </a> </td>
   <td align="left"> ZNF717; MIR4273 </td>
   <td align="left"> Yes(ZNF717) </td>
   <td align="left"> No </td>
   <td align="right"> 3.974 </td>
   <td align="left"> ZNF717: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cCa/cTa, P87L;   MIR4273: MODIFIER, DOWNSTREAM, 2930;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 1165163 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SPON2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.990 </td>
   <td align="left"> SPON2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gAg/gGg, E111G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 170913122 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MFAP3L </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.650 </td>
   <td align="left"> MFAP3L: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tcc/Ccc, S110P, S213P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 37863193 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=76138483' target="_blank" > rs76138483 </a> </td>
   <td align="left"> PGM2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 17.980 </td>
   <td align="left"> PGM2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tTc/tCc, F600S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 38799710 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4833095' target="_blank" > rs4833095 </a> </td>
   <td align="left"> TLR1 </td>
   <td align="left"> Yes(TLR1) </td>
   <td align="left"> No </td>
   <td align="right"> 1.852 </td>
   <td align="left"> TLR1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aAt/aGt, N248S;  </td>
   <td align="left"> 1 </td>
   <td align="left"> risk factor </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 84377287 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=199771286' target="_blank" > rs199771286 </a> </td>
   <td align="left"> MRPS18C; FAM175A; HELQ </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 3.632 </td>
   <td align="left"> MRPS18C: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, ttG/ttC, L19F;   FAM175A: MODIFIER, DOWNSTREAM, 4807;   HELQ: MODIFIER, UPSTREAM, 262;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 87277 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2006764' target="_blank" > rs2006764 </a> </td>
   <td align="left"> ZNF595; ZNF718 </td>
   <td align="left"> Yes(ZNF718) </td>
   <td align="left"> No </td>
   <td align="right"> 2.646 </td>
   <td align="left"> ZNF595: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gGt/gCt, G627A;   ZNF718: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 4 </td>
   <td align="left"> 95539267 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=966845' target="_blank" > rs966845 </a> </td>
   <td align="left"> PDLIM5 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.014 </td>
   <td align="left"> PDLIM5: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gct/Act, A20T, A223T, A236T, A242T, A345T, A374T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 138729374 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> PROB1; MZB1; SPATA24 </td>
   <td align="left"> Yes(MZB1) </td>
   <td align="left"> No </td>
   <td align="right"> 11.070 </td>
   <td align="left"> PROB1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gCt/gAt, A466D;   MZB1: MODIFIER, UPSTREAM, 3769;   SPATA24: MODIFIER, DOWNSTREAM, 3082;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 140750044 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=6860609' target="_blank" > rs6860609 </a> </td>
   <td align="left"> PCDHGB3; PCDHGA1; PCDHGA2; PCDHGA3; PCDHGA4; PCDHGA5; PCDHGA6; PCDHGB1; PCDHGB2 </td>
   <td align="left"> Yes(PCDHGB3) </td>
   <td align="left"> No </td>
   <td align="right"> 0.001 </td>
   <td align="left"> PCDHGB3: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTt/gCt, V28A;   PCDHGA1: MODIFIER, INTRAGENIC;   PCDHGA1: MODIFIER, INTRON;   PCDHGA2: MODIFIER, INTRON;   PCDHGA3: MODIFIER, INTRAGENIC;   PCDHGA3: MODIFIER, INTRON;   PCDHGA4: MODIFIER, INTRON;   PCDHGA5: MODIFIER, DOWNSTREAM, 3705;   PCDHGA5: MODIFIER, INTRON, 3705;   PCDHGA5: MODIFIER, INTRAGENIC, 3705;   PCDHGA6: MODIFIER, UPSTREAM, 3607;   PCDHGB1: MODIFIER, INTRON;   PCDHGB1: MODIFIER, INTRAGENIC;   PCDHGB2: MODIFIER, INTRAGENIC;   PCDHGB2: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 159520708 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> PWWP2A </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 18.570 </td>
   <td align="left"> PWWP2A: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Ccc/Tcc, P317S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 176778602 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=143985896' target="_blank" > rs143985896 </a> </td>
   <td align="left"> LMAN2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 15.760 </td>
   <td align="left"> LMAN2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tGc/tAc, C16Y;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 176931051 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> DOK3 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 12.790 </td>
   <td align="left"> DOK3: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aAg/aTg, K475M;   DOK3: MODIFIER, INTRON, aAg/aTg, K475M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 35861068 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=1494558' target="_blank" > rs1494558 </a> </td>
   <td align="left"> IL7R </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.016 </td>
   <td align="left"> IL7R: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTc/aCc, I66T;  </td>
   <td align="left"> 1 </td>
   <td align="left"> Pathogenic </td>
   <td align="left"> stomach carcinoma adenocarcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=   ' target="_blank" >     </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 36182129 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SKP2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.790 </td>
   <td align="left"> SKP2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cTa/cCa, L210P, L424P;   SKP2: MODIFIER, INTRON, cTa/cCa, L210P, L424P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 37205506 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> C5orf42 </td>
   <td align="left"> Yes(C5orf42) </td>
   <td align="left"> No </td>
   <td align="right"> 16.930 </td>
   <td align="left"> C5orf42: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTt/aAt, I1067N;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 42719239 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=6180' target="_blank" > rs6180 </a> </td>
   <td align="left"> GHR </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 14.260 </td>
   <td align="left"> GHR: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atc/Ctc, 1041, 672, I522L, I544L, I551L;   GHR: MODIFIER, UTR_3_PRIME, Atc/Ctc, 1041, 672, I522L, I544L, I551L;  </td>
   <td align="left"> 1 </td>
   <td align="left"> risk factor </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 68728426 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MARVELD2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 10.280 </td>
   <td align="left"> MARVELD2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atg/Gtg, M407V, M419V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 73981270 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=820878' target="_blank" > rs820878 </a> </td>
   <td align="left"> HEXB </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 3.141 </td>
   <td align="left"> HEXB: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tTg/tCg, L62S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> breast carcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=   ' target="_blank" >     </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 78936682 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> PAPD4 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 13.050 </td>
   <td align="left"> PAPD4: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, caG/caT, Q158H;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 5 </td>
   <td align="left"> 98115642 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=76606616' target="_blank" > rs76606616 </a> </td>
   <td align="left"> RGMB </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.110 </td>
   <td align="left"> RGMB: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, ttG/ttT, L206F;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> large_intestine carcinoma adenocarcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=21892161' target="_blank" > 21892161 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 12122645 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=6900196' target="_blank" > rs6900196 </a> </td>
   <td align="left"> HIVEP1 </td>
   <td align="left"> No </td>
   <td align="left"> Yes(HIVEP1)   </td>
   <td align="right"> 0.863 </td>
   <td align="left"> HIVEP1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Act/Gct, T873A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 135511443 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MYB </td>
   <td align="left"> No </td>
   <td align="left"> Yes(MYB)    </td>
   <td align="right"> 22.400 </td>
   <td align="left"> MYB: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cTg/cCg, L162P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 144416537 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=199548273' target="_blank" > rs199548273 </a> </td>
   <td align="left"> SF3B5 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 24.700 </td>
   <td align="left"> SF3B5: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTg/gGg, V33G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 158924338 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=34395018' target="_blank" > rs34395018 </a> </td>
   <td align="left"> TULP4 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 5.298 </td>
   <td align="left"> TULP4: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cca/Tca, P1215S;   TULP4: MODIFIER, INTRON, Cca/Tca, P1215S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 26234929 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=75255708' target="_blank" > rs75255708 </a> </td>
   <td align="left"> HIST1H1D </td>
   <td align="left"> Yes(HIST1H1D) </td>
   <td align="left"> No </td>
   <td align="right"> 20.400 </td>
   <td align="left"> HIST1H1D: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aAc/aCc, N78T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> large_intestine carcinoma adenocarcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=21892161' target="_blank" > 21892161 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 32712999 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=75289558' target="_blank" > rs75289558 </a> </td>
   <td align="left"> HLA-DQA2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 5.147 </td>
   <td align="left"> HLA-DQA2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCc/aGc, T49S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 33048640 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=1042140' target="_blank" > rs1042140 </a> </td>
   <td align="left"> HLA-DPB1; HLA-DPA1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 10.030 </td>
   <td align="left"> HLA-DPB1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aag/Gag, K98E;   HLA-DPA1: MODIFIER, UPSTREAM, 85;  </td>
   <td align="left"> 1 </td>
   <td align="left"> risk factor </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 43173848 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   .    ' target="_blank" > .   .     </a> </td>
   <td align="left"> CUL9 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 20.800 </td>
   <td align="left"> CUL9: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atg/Gtg, M1633V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> large_intestine colon carcinoma adenocarcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=22895193' target="_blank" > 22895193 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 46672943 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=1051931' target="_blank" > rs1051931 </a> </td>
   <td align="left"> PLA2G7; TDRD6 </td>
   <td align="left"> Yes(PLA2G7) </td>
   <td align="left"> No </td>
   <td align="right"> 15.750 </td>
   <td align="left"> PLA2G7: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gTa/gCa, V379A;   TDRD6: MODIFIER, DOWNSTREAM, 887;  </td>
   <td align="left"> 1 </td>
   <td align="left"> risk factor </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 57398201 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=62398997' target="_blank" > rs62398997 </a> </td>
   <td align="left"> PRIM2 </td>
   <td align="left"> Yes(PRIM2) </td>
   <td align="left"> No </td>
   <td align="right"> 16.080 </td>
   <td align="left"> PRIM2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgt/Cgt, C302R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 57398264 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=62398998' target="_blank" > rs62398998 </a> </td>
   <td align="left"> PRIM2 </td>
   <td align="left"> Yes(PRIM2) </td>
   <td align="left"> No </td>
   <td align="right"> 13.970 </td>
   <td align="left"> PRIM2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aca/Gca, T323A;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 57398270 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=62398999' target="_blank" > rs62398999 </a> </td>
   <td align="left"> PRIM2 </td>
   <td align="left"> Yes(PRIM2) </td>
   <td align="left"> No </td>
   <td align="right"> 15.420 </td>
   <td align="left"> PRIM2: HIGH, NONSENSE, STOP_GAINED, Cag/Tag, Q325* ;      </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 57512565 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4294008' target="_blank" > rs4294008 </a> </td>
   <td align="left"> PRIM2 </td>
   <td align="left"> Yes(PRIM2) </td>
   <td align="left"> No </td>
   <td align="right"> 0.057 </td>
   <td align="left"> PRIM2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tct/Cct, S465P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 7182296 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> RREB1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 35.000 </td>
   <td align="left"> RREB1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGg/cAg, R51Q;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 74171665 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MTO1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 13.700 </td>
   <td align="left"> MTO1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gcg/Acg, A30T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 84925066 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=147915749' target="_blank" > rs147915749 </a> </td>
   <td align="left"> KIAA1009 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 15.390 </td>
   <td align="left"> KIAA1009: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTt/aCt, I146T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 87968677 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=41273279' target="_blank" > rs41273279 </a> </td>
   <td align="left"> ZNF292 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 9.969 </td>
   <td align="left"> ZNF292: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, tCt/tTt, S1777F;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 6 </td>
   <td align="left"> 99365371 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=143154211' target="_blank" > rs143154211 </a> </td>
   <td align="left"> FBXL4 </td>
   <td align="left"> Yes(FBXL4) </td>
   <td align="left"> No </td>
   <td align="right"> 5.206 </td>
   <td align="left"> FBXL4: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTa/aCa, I246T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 111503593 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   .    ' target="_blank" > .   .     </a> </td>
   <td align="left"> DOCK4 </td>
   <td align="left"> Yes(DOCK4) </td>
   <td align="left"> No </td>
   <td align="right"> 19.140 </td>
   <td align="left"> DOCK4: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtg/Atg, V770M;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> haematopoietic_and_lymphoid_tissue large_intestine(2x) colon(2x) carcinoma(2x) haematopoietic_neoplasm acute_myeloid_leukaemia_associated_with_MDS adenocarcinoma </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=21909114  22895193(2x)   ' target="_blank" > 21909114  22895193(2x)    </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 151884449 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=80039782' target="_blank" > rs80039782 </a> </td>
   <td align="left"> MLL3 </td>
   <td align="left"> No </td>
   <td align="left"> Yes(MLL3)   </td>
   <td align="right"> 14.460 </td>
   <td align="left"> MLL3: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Acg/Ccg, T1636P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 36396890 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> KIAA0895 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 6.451 </td>
   <td align="left"> KIAA0895: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCt/aTt, T112I, T12I, T150I, T163I;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 3658807 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SDK1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 20.700 </td>
   <td align="left"> SDK1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgg/Cgg, W132R;   SDK1: MODIFIER, INTRAGENIC, Tgg/Cgg, W132R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 56136260 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4245575' target="_blank" > rs4245575 </a> </td>
   <td align="left"> SUMF2; CCT6A </td>
   <td align="left"> Yes(CCT6A) </td>
   <td align="left"> No </td>
   <td align="right"> 0.127 </td>
   <td align="left"> SUMF2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gaC/gaA, D70E;   CCT6A: MODIFIER, DOWNSTREAM, 4578;   SUMF2: MODIFIER, INTRON, gaC/gaA, D70E;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 82582423 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> PCLO </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 0.133 </td>
   <td align="left"> PCLO: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtt/Ttt, V2616F;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 91714911 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=1063242' target="_blank" > rs1063242 </a> </td>
   <td align="left"> AKAP9 </td>
   <td align="left"> Yes(AKAP9) </td>
   <td align="left"> Yes(AKAP9)   </td>
   <td align="right"> 4.014 </td>
   <td align="left"> AKAP9: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cct/Tct, P2971S, P2979S;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 7 </td>
   <td align="left"> 99023166 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ATP5J2-PTCD1; PTCD1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 26.400 </td>
   <td align="left"> ATP5J2-PTCD1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGg/cTg, R379L;   PTCD1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGg/cTg, R330L;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 121457750 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MTBP; MRPL13 </td>
   <td align="left"> Yes(MTBP) </td>
   <td align="left"> No </td>
   <td align="right"> 20.700 </td>
   <td align="left"> MTBP: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Aaa/Gaa, K14E;   MRPL13: MODIFIER, UPSTREAM, 103;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 136555000 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> KHDRBS3 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 28.200 </td>
   <td align="left"> KHDRBS3: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aGa/aAa, R104K;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 142178240 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=138519822' target="_blank" > rs138519822 </a> </td>
   <td align="left"> DENND3 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.880 </td>
   <td align="left"> DENND3: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gat/Aat, D551N;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 144940706 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=112377501' target="_blank" > rs112377501 </a> </td>
   <td align="left"> EPPK1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 16.790 </td>
   <td align="left"> EPPK1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cGc/cAc, R2239H;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> autonomic_ganglia breast central_nervous_system haematopoietic_and_lymphoid_tissue(3x) meninges pancreas prostate brain carcinoma(3x) haematopoietic_neoplasm(3x) meningioma neuroblastoma primitive_neuroectodermal_tumour-medulloblastoma acinar_carcinoma Stage:1B acute_myeloid_leukaemia adenocarcinoma ductal_carcinoma Metastatic site:bone large_cell Grade:Some Grade data are given in publication </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid=22832583  23334666  23348505  23563269  24293293' target="_blank" > 22832583  23334666  23348505  23563269  24293293 </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 144946599 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> EPPK1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 5.816 </td>
   <td align="left"> EPPK1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gcc/Acc, A275T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 20110798 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> LZTS1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 18.780 </td>
   <td align="left"> LZTS1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, cAg/cCg, Q215P;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 8 </td>
   <td align="left"> 56699007 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TGS1 </td>
   <td align="left"> Yes(TGS1) </td>
   <td align="left"> No </td>
   <td align="right"> 8.939 </td>
   <td align="left"> TGS1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Gtt/Att, V184I;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 9 </td>
   <td align="left"> 124632837 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> TTLL11 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 27.300 </td>
   <td align="left"> TTLL11: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Cgg/Ggg, R648G;   TTLL11: MODIFIER, INTRAGENIC, Cgg/Ggg, R648G;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 9 </td>
   <td align="left"> 131186498 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=144617260' target="_blank" > rs144617260 </a> </td>
   <td align="left"> CERCAM </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 22.000 </td>
   <td align="left"> CERCAM: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Atg/Gtg, M170V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 9 </td>
   <td align="left"> 139990693 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> MAN1B1 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 5.923 </td>
   <td align="left"> MAN1B1: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aCa/aTa, T157I;   MAN1B1: MODIFIER, EXON, aCa/aTa, T157I;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 9 </td>
   <td align="left"> 2073321 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> SMARCA2 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 35.000 </td>
   <td align="left"> SMARCA2: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gCc/gTc, A619V;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> 9 </td>
   <td align="left"> 37440852 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZBTB5; GRHPR </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 10.680 </td>
   <td align="left"> ZBTB5: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gGa/gAa, G566E;   GRHPR: MODIFIER, DOWNSTREAM, 3866;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> X </td>
   <td align="left"> 134426392 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?.   ' target="_blank" > .    </a> </td>
   <td align="left"> ZNF75D </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 10.480 </td>
   <td align="left"> ZNF75D: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, gCc/gAc, A140D;   ZNF75D: MODIFIER, INTRON, gCc/gAc, A140D;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> X </td>
   <td align="left"> 30872602 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=5927629' target="_blank" > rs5927629 </a> </td>
   <td align="left"> TAB3 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 11.280 </td>
   <td align="left"> TAB3: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Tgg/Cgg, W394R;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> X </td>
   <td align="left"> 54209387 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=2495783' target="_blank" > rs2495783 </a> </td>
   <td align="left"> FAM120C </td>
   <td align="left"> Yes(FAM120C) </td>
   <td align="left"> No </td>
   <td align="right"> 0.019 </td>
   <td align="left"> FAM120C: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, aTt/aCt, I82T;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
  <tr>
   <td align="left"> Y </td>
   <td align="left"> 21154466 </td>
   <td align="left"> <a href =  'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=10465460' target="_blank" > rs10465460 </a> </td>
   <td align="left"> CD24; TTTY14 </td>
   <td align="left"> No </td>
   <td align="left"> No </td>
   <td align="right"> 2.235 </td>
   <td align="left"> CD24: MODERATE, MISSENSE, NON_SYNONYMOUS_CODING, Act/Tct, T44S;   TTTY14: MODIFIER, INTRON;  </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left"> <a href =  'http://cancer.sanger.ac.uk/cosmic/study/overview?pmid= ' target="_blank" >   </a> </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
   <td align="left">   </td>
  </tr>
</tbody>
</table>


<script type="text/javascript" charset="utf-8">
    $(document).ready(function() {
    var dogtable = $('#SNPs_table').dataTable({
        "sDom": 'T<"clear">lfrtip',
        "oTableTools": {
            "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
            "aButtons": [
                "copy",
                "csv"
            ]
        },
        "aoColumns": [ 
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null,
          null
       ]         
      })
      });
</script>

******

# Drug matching based on overexpressed genes

<div class="warning"><pre class="knitr r">## Warning: cannot open file
## '/gpfs/group/su/kfisch/Cancer/Kumar/results/reports/SRR1027184/dog.csv':
## Permission denied
</pre></div>
<div class="error"><pre class="knitr r">## Error: cannot open the connection
</pre></div>


<script type="text/javascript" charset="utf-8">
    $(document).ready(function() {
    var dogtable = $('#dog_table').dataTable({
        "sDom": 'T<"clear">lfrtip',
        "oTableTools": {
            "sSwfPath": "../media/swf/copy_csv_xls_pdf.swf",
            "aButtons": [
                "copy",
                "csv"
            ]
        }
      })
      });
</script>

******

# References 

[1] _Babraham Bioinformatics - FastQC A Quality Control tool for
High Throughput Sequence Data_. <URL:
http://www.bioinformatics.babraham.ac.uk/projects/fastqc>. 2014.
<URL: http://www.bioinformatics.babraham.ac.uk/projects/fastqc>.

[2] _HTSeq: Analysing high-throughput sequencing data with Python â
HTSeq v0.5.4p2 documentation_. <URL:
http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html>.
2014. <URL:
http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html>.

[3] _Picard_. <URL: http://picard.sourceforge.net/>. 2014. <URL:
http://picard.sourceforge.net/>.

[4] S. Anders and W. Huber. "Differential expression analysis for
sequence count data". In: _Genome Biol_ 11.10 (2010), p. R106.
DOI: 10.1186/gb-2010-11-10-r106. <URL:
http://dx.doi.org/10.1186/gb-2010-11-10-r106>.

[5] A. Dobin, C. A. Davis, F. Schlesinger, et al.
"$\lbrace$STAR$\rbrace$: ultrafast universal
$\lbrace$RNA$\rbrace$-seq aligner". In: _Bioinformatics_ 29.1
(Oct. 2012), pp. 15-21. DOI: 10.1093/bioinformatics/bts635. <URL:
http://dx.doi.org/10.1093/bioinformatics/bts635>.

[6] H. Edgren, A. Murumagi, S. Kangaspeska, et al. "Identification
of fusion genes in breast cancer by paired-end
$\lbrace$RNA$\rbrace$-sequencing". In: _Genome Biol_ 12.1 (2011),
p. R6. DOI: 10.1186/gb-2011-12-1-r6. <URL:
http://dx.doi.org/10.1186/gb-2011-12-1-r6>.

[7] H. Li, B. Handsaker, A. Wysoker, et al. "The Sequence
Alignment/Map format and $\lbrace$SAMtools$\rbrace$". In:
_Bioinformatics_ 25.16 (Jun. 2009), pp. 2078-2079. DOI:
10.1093/bioinformatics/btp352. <URL:
http://dx.doi.org/10.1093/bioinformatics/btp352>.

[8] M. Shugay, I. O. de Mendibil, J. L. Vizmanos, et al.
"Oncofuse: a computational framework for the prediction of the
oncogenic potential of gene fusions". In: _Bioinformatics_ 29.20
(Aug. 2013), pp. 2539-2546. DOI: 10.1093/bioinformatics/btt445.
<URL: http://dx.doi.org/10.1093/bioinformatics/btt445>.



# Session Info

<div class="chunk" id="unnamed-chunk-1"><div class="rcode"><div class="output"><pre class="knitr r">## R version 3.0.1 (2013-05-16)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=C                 LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] tcltk     splines   parallel  methods   stats     graphics  grDevices
##  [8] utils     datasets  base     
## 
## other attached packages:
##  [1] stringr_0.6.2           igraph_0.7.1           
##  [3] ReactomePA_1.6.1        pathview_1.2.4         
##  [5] graphite_1.8.1          SPIA_2.14.0            
##  [7] KEGGgraph_1.20.0        graph_1.40.1           
##  [9] pamr_1.54.1             survival_2.37-7        
## [11] RJSONIO_1.3-0           DESeq2_1.2.10          
## [13] RcppArmadillo_0.4.320.0 Rcpp_0.11.2            
## [15] GenomicRanges_1.14.4    XVector_0.2.0          
## [17] IRanges_1.20.7          cluster_1.15.2         
## [19] RColorBrewer_1.0-5      edgeR_3.4.2            
## [21] limma_3.18.13           knitcitations_1.0-1    
## [23] KEGGREST_1.2.2          pander_0.3.8           
## [25] xtable_1.7-3            XML_3.98-1.1           
## [27] annotate_1.40.1         gdata_2.13.3           
## [29] dplyr_0.1.3.0.99        plyr_1.8.1             
## [31] data.table_1.9.2        org.Hs.eg.db_2.10.1    
## [33] RSQLite_0.11.4          DBI_0.2-7              
## [35] AnnotationDbi_1.24.0    DESeq_1.14.0           
## [37] lattice_0.20-29         locfit_1.5-9.1         
## [39] Biobase_2.22.0          BiocGenerics_0.8.0     
## [41] rCharts_0.4.2           yaml_2.1.13            
## [43] knitrBootstrap_0.9.0    knitr_1.6              
## 
## loaded via a namespace (and not attached):
##  [1] assertthat_0.1     bibtex_0.3-6       Biostrings_2.30.1 
##  [4] colorspace_1.2-4   digest_0.6.4       DO.db_2.7         
##  [7] DOSE_2.0.0         evaluate_0.5.5     formatR_0.10      
## [10] genefilter_1.44.0  geneplotter_1.40.0 ggplot2_1.0.0     
## [13] GO.db_2.10.1       GOSemSim_1.20.3    grid_3.0.1        
## [16] gtable_0.1.2       gtools_3.4.1       httr_0.4          
## [19] lubridate_1.3.3    markdown_0.7.2     MASS_7.3-33       
## [22] memoise_0.2.1      munsell_0.4.2      png_0.1-7         
## [25] proto_0.3-10       qvalue_1.36.0      RCurl_1.95-4.3    
## [28] reactome.db_1.46.1 RefManageR_0.8.3   reshape2_1.4      
## [31] Rgraphviz_2.8.1    scales_0.2.4       stats4_3.0.1      
## [34] tools_3.0.1        whisker_0.3-2
</pre></div>
<div class="output"><pre class="knitr r">## [1] "Thu Sep  4 15:57:40 2014"
</pre></div>
</div></div>
