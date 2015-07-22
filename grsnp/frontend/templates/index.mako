<body bgcolor="#A8D5FF" style="width: 1000px">

	<div class="topbar">
		<div class="topbar-inner">
			<div class="container-fluid">
				<a class="brand" href="./" style="font-family:Futura,'Helvetica Neue', Arial">
					<img width="200px" src="static/images/GRLogo.png" />
				</a>
				<ul class="pull-right">
					<li><a href="./overview">Overview</a></li>
					<li><a href="./news">News</a></li>
					<li><a href="./demo">Demo</a></li>
					<!-- <li><a href="./cite">How to Cite</a></li>
					<li><a href="http://sourceforge.net/projects/genomerunner/">GenomeRunner on SourceForge</a></li>
					<li><a href="./roadmap">Roadmap</a></li> -->
					<li><a href="https://mdozmorov.github.io/grdocs/" target="_blank">Help</a></li>
<!-- <div>
					<li><a href="./google">Google Group</a></li>
					<li><img width="30px" src="static/new-icon.jpg" alt="New: GenomeRunnerSNP Google Groups" /></li>
				</div> -->
			</ul>
			<!-- <img height="40px" src="static/images/logo-reversed-small.jpg" align="right" /> -->
		</div>
	</div>
</div>
<form id="frmQuery" name="frmQuery" action="query" method="post" enctype="multipart/form-data">
	<div id="content">
		<div class="well" style="margin-top: -15px; padding: 0px">			
			<h3>Select Database Version:</h3>
			${database_versions}		
			<h3>GenomeRunner: Functional interpretation of SNPs within regulatory context</h3>
			<p>
				<span style="font-size: 16px;"><span style="font-family:arial,helvetica,sans-serif;">GenomeRunner is a tool for annotation and enrichment analysis of the SNP sets by considering SNPs co-localization with functional/regulatory genome annotation data. GenomeRunner is particularly useful for interpretation of the collective regulatory impact of SNPs in non-protein coding regions. An example of GenomeRunner&#39;s results can be found in the analysis of Sjogren&#39;s syndrome GWAS (<em><a href="http://www.nature.com/ng/journal/v45/n11/full/ng.2792.html" target="_blank">Nature Genetics</a></em>) where it identified RFX5 transcription factor binding site as strongly enriched with the disease-associated SNPs.</span></p>
				<p>
					<br />
					<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">GenomeRunner calculates <a href="http://mdozmorov.github.io/grdocs/hypergeom4/enrichment.html">the enrichment p-values</a> by evaluating whether a SNP set co-localizes with regulatory datasets more often that could happen by chance. For three or more SNP sets, GenomeRunner performs <a href="https://mdozmorov.github.io/grdocs/hypergeom4/episimilarity.html"> a regulatory similarity analysis</a> by correlating SNP set-specific regulatory enrichment profiles. The downloadable results are visualized as interactive heatmaps and tables <a href="result?id=example">(Example)</a>.</span></span></p>
					<p>
					</p>
					</div>
					<div class="well">
						<div style="float:right;margin-top: 10px;">
							<b style="font-size:150%;">Organism:</b>
							${paths.org_as_html(id="org")} 
							<img class="helptooltip" title="Select organism-specific genome assembly" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</div> 
						<h3>1. Select sets of SNPs of interest <img class="helptooltip" title="Upload files with rsIDs or genomic coordinates in BED format of SNPs of interest. At least 5 SNPs per set are required. Multiple file upload supported. Note: Avoid special characters and extra dots in file names. Do not use SNP IDs other than rsIDs." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</h3>
						<div id="div_upload_fois"	style="float: left;margin-right: 13px;margin-top: 14px;"	>
							<h4 style="float:left;">Files:</h4><input type="file" id="inputbedfile" style="margin:5px" name="bed_file" multiple="multiple"/>
							<a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">What should my data look like?</a>
						</div>
						<div id="div_demo_fois">
							<h4 style="font-size: 110%;float: left;margin-top: 15px;margin-right: 9px;;">Demo SNPs file sets: </h4>
							<div class="btn-group" id="btngroup_demo_fois" data-toggle="buttons-radio" data-toggle-name="demo_fois">
								<input type="hidden" style="margin-top: -3px;" name="demo_fois"/>
								${demo_snps}
							</div>
						</div>
						<div id="accfoi" class="accordion" style="padding-bottom: 1em;list-style:none;margin-top: 20px">        
							<h3  id="accordionheader"><a href="#" style="font-size:120%">Paste data in .BED format</a></h3>
							<div>          
								<table>
									<tr>
										<td>
											<textarea id="inputbeddata" rows=5 cols=30 style="margin:10px" name="bed_data" wrap="off" disabled>
											</textarea>
										</td>
										<td>
											<label stlye="float:left">
												Paste a list of rsIDs or tab-separated genomic coordinates (recommended) of SNPs of interest in .BED format (no headers). Pasting the data, or submitting one set of SNPs, restricts the analysis to the enrichment results only. To get the enrichment and epigenomic similarity heatmaps, upload multiple sets of SNPs as separate files using 'Choose Files' button.
											</label>
											<p style="font-size: 120%">Use <a href="https://www.ncbi.nlm.nih.gov/projects/SNP/dbSNP.cgi?list=rslist">the NCBI conversion tool</a> to convert rsIDs of SNPs into genomic coordinates.</p>
										</td>
									</tr>
								</table>
							</div>
						</div>
					</div>	
					<div class="well">
						<h3 style="float:left">2. Define the background: ${default_background}<img class="helptooltip" title="By default, all common SNPs are used as a 'universe' to calculate the probability of SNPs in a set to be enriched with regulatory dataset" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</h3>
						<!-- <input type="checkbox" style="font-size:120%;margin-top:1em" name="run_random">Run randomization test</input> -->
						<b style=" font-size:120%; margin-left: 10px"></b>
						<div class="accordion" style="padding-bottom: 1em;list-style:none;margin-top:2em">        
							<h3  id="accordionheader"><a href="#" style="font-size:120%">Upload custom background</a></h3>
							<div style="height=100px">       
								Bed File:
								<input type="file" id="inputbackgroundfile" style="margin:5px" name="background_file" />		
								<div id="accback" class="accordion" style="padding-bottom: 1em;list-style:none;margin-top: 20px">        
									<h3  id="accordionheader"><a href="#" style="font-size:120%">Paste data in .BED format</a></h3>							        
									<table style="margin-bottom:0px; padding-bottom:0px">
										<tr>
											<td>
												<textarea id="inputbackgrounddata" rows=5 cols=30 style="margin:10px" name="background_data" wrap="off" disabled>
												</textarea>
											</td>
											<td>
												<label stlye="float:left">
													The 'background', or 'universe' of SNPs assessed in a study is critical for the correct p-value calculation. A set of SNPs of interest should be a subset of the background, else the p-values may be wrong.<br><br>
													Default background selection &#40;all common SNPs from the latest organism-specific database&#41; is suitable when a genome-wide study was performed. When a microarray was used for SNPs profiling, it is advisable to use all SNPs on that array as a background.
												</label>
											</td>
										</tr>
									</table>
								</div>
							</div>
						</div>
					</div>
					<div class="well">
						<h3 style="float: left;margin-right: 10px;margin-top: 1px;">3. Select regulatory datasets:</h3>
						${custom_gfs}
						<br>

						<div id="accordGFS" class="accordion" style="padding-bottom: 1em;list-style:none; margin-top: 3em;" onClick="renderCheckBoxTree()"> 
							<h3  id="accordionheader"><a id='gfselheader' href="#" style="font-size:120%; height: 100%">Choose regulatory datasets</a></h3>
							<div >	
								<ul>
									<li id="list-bedbackground">					
										Bed Files (.bed or .gz):
										<input type="file" id="inputgenomicfeaturefile" style="margin:5px" name="genomicfeature_file" multiple="multiple"/>
										<a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">What should my data  look like?</a>
									</li>				
								</ul>
								<div id="grfdroplist" style="display: table;">	
									<div style="display: table-row;">
										<div id="divCheckBox" style="width: 70%; margin:10px; display: table-cell; verticle-align: top; visibility: hidden">
											<label >Select (categories of) regulatory datasets: </label>
											<div>
												<a class="btn" style="margin-top: 9px;"onClick="$('#jstree_gfs').jstree('open_all');">Expand All</a>
												<a class="btn" style="margin-top: 9px;" onClick="$('#jstree_gfs').jstree('close_all');">Collapse All</a>
												<a class="btn" style="margin-top: 9px;" onClick="$('#jstree_gfs').jstree('check_all');">Select All</a>
												<a class="btn" style="margin-top: 9px" onClick="$('#jstree_gfs').jstree('uncheck_all');">Deselect all</a>
												<a class="btn" style="margin-top: 9px;" id="descriptions">Track Descriptions</a>
											</div>
											<label >Search genomic features</label>
											<input id="txt_gfs_search" class='input' type="text"></input>
											<a class="btn" style="margin-top: 9px;" id="treeSelect" onClick="treeviewSelectSearchedClick()">Select</a>

											<div id="jstree_gfs"></div>
										</div>
										<div style="margin:10px; display: table-cell;  verticle-align: top; padding-left:10px">
											<label style="width: 100%; margin: 5px;">Regulatory datasets (aka functional/regulatory/epigenomic data) are mirrored from and organized according to the UCSC genome database scheme. If you know names of the tracks you want to run enrichment analyses with, start typing their names in the search box. Or simply search for keywords, like 'H3K4me1', to see which tracks are available. Use checkboxes in the collapsible TreeView to select categories of regulatory elements</label><br>
											<input type="checkbox" style="font-size:120%;margin-top:1em" name="run_annot">Run annotation analysis</input><img class="helptooltip" title="Annotate each SNP in each set for the number of overlaps with the selected regulatory elements." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
										</div>

									</div>
								</div>			 

							</div>
						</div>
					</div>	

					<div class="well">
						<table width="100%">
							<tr>
								<th width=30% style="vertical-align:bottom">
									P-Value Adjustment:	
									<select name="padjust" id="padjust">
										<option value="None">none</option>
										<option value="bonferroni">bonferroni</option>
										<option value="holm">holm</option>
										<option value="hochberg">hochberg</option>
										<option value="hommel">hommel</option>
										<option value="BH">BH</option>
										<option value="BY">BY</option>
										<option value="fdr" selected>fdr</option>
									</select>
									<img class="helptooltip" title="Sets multiple tests correction method when testing a set of SNPs against multiple regulatory datasets" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								</th>
								<th width="30%" style="vertical-align:bottom">
									Percent score threshold: ${pct_scores}
									<img class="helptooltip" title="Increasing this number filters out more low-level signal in the regulatory datasets. If a regulatory dataset does not have a score, this setting is ignored" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								</th>
								<th width="30%" style="vertical-align:bottom">
									Strand selection: 
									<select name="strand">
										<option value="both" selected>Both</option>
										<option value="plus">Plus</option>
										<option value="minus">Minus</option>
									</select>
									<img class="helptooltip" title="Sets whether or not to use strand-specific regulatory datasets, if available. If a regulatory dataset does not have a strand, this setting is ignored" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								</th>
							</tr>
							<tr>

								<td><input type="checkbox" id="disclaimer" checked="checked" >I certify that I understand that GenomeRunner is for research purposes only.</input></td>
								<td id="td_submit" style="width:90px">
									<button id="btnSubmit" class="btn btn-primary" onclick="submit_job()" type="submit" >Submit job</button>
									<img class="helptooltip" title="Submits the job for enrichment/regulatory similarity analyses" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								</td>
								<td id="td_submit" style="width:170px">
									<h3 id="upmessage" style="visibility:hidden;margin-left: -94px;margin-top: 3px;">Uploading files. Please do not refresh the page.</h3>
								</td>	
							</tr>
						</table>
						<br>
						<div class="ui-state-highlight ui-corner-all">
							<p><span class="ui-icon ui-icon-info" style="float: left; margin-right: .3em;">
							</span>GenomeRunner works best in Chrome (Windows) or Firefox. Mac users, use Safari.</p>
						</div>

					</div>
				</div>
			</form>
			<p style="text-align: center;">
				<span style="font-family:arial,helvetica,sans-serif;">You are the&nbsp;<img alt="stats counter" border="0" src="http://www.easycounter.com/counter.php?mdozmorov" />&nbsp;visitor</span>
			</p>
		</div>		
	</body>
