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
					<li><a href="./help">Help</a></li>
<!-- <div>
					<li><a href="./google">Google Group</a></li>
					<li><img width="30px" src="static/new-icon.jpg" alt="New: GenomeRunnerSNP Google Groups" /></li>
				</div> -->
			</ul>
			<img height="40px" src="static/images/logo-reversed-small.jpg" align="right" />
		</div>
	</div>
</div>
<form name="frmQuery" action="query" method="post" enctype="multipart/form-data">
	<div id="content">
		<div class="well" style="margin-top: 40px; padding: 0px">
			This website is free and open to all users and there is no login requirement. <br>
					<p>
			<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">News: GenomeRunner will be presented on a software demonstration session at <a href="http://agbt.org/about.html" style="font-family: arial, helvetica, sans-serif; font-size: 14px;" target="_blank">AGBT 2014</a>, February 12-15.</span></span></p>
			<h3>GenomeRunner Web: Functional interpretation of SNPs within epigenomic context</h3>
			<p>
				<span style="font-size: 16px;"><span style="font-family:arial,helvetica,sans-serif;">GenomeRunner Web is a tool for functional interpretation of sets of SNPs&nbsp;</span></span><span style="font-family: arial, helvetica, sans-serif; font-size: 16px;">by considering their co-localization with functional/regulatory genome annotation data (epigenomic elements). It is particularly useful for the interpretation of functional roles of SNPs in non-protein coding regions, and rare variants. An example of GenomeRunner&#39;s results can be found in the analysis of Sjogren&#39;s syndrome GWAS (<em><a href="http://www.nature.com/ng/journal/v45/n11/full/ng.2792.html" target="_blank">Nature Genetics</a></em>), where it identified RFX5 transcription factor binding site as strongly associated with the disease&#39; SNPs.</span></p>
				<p>
					<br />
					<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">As an <a href="result?id=example">output</a>, GenomeRunner Web calculates <a href="help#enrichment">enrichment p-values</a>&nbsp;by evaluating whether a&nbsp;set of SNPs co-localizes with regulatory elements more often that could happen by chance. For three or more sets of SNPs, GenomeRunner Web performs <a href="help#episimilarity">&#39;epigenomic similarity&#39; analysis</a>&nbsp;by correlating set-specific profiles of enrichment p-values. Downloadable results are visualized as interactive heatmaps and tables.</span></span></p>
					<p>
						&nbsp;</p>
					</div>
					<div class="well">
						<div style="float:right;margin-top: 10px;">
							<b style="font-size:150%;">Organism:</b>
							${paths.org_as_html(id="org")} 
							<img class="helptooltip" title="Select organism-specific genome assembly" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</div> 
						<h3>1. Select sets of SNPs of interest <img class="helptooltip" title="Upload .BED formatted files with genomic coordinates of SNPs of interest. Multiple file upload supported." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</h3>
						<div id="div_upload_fois"	style="float: left;margin-right: 13px;margin-top: 14px;"	>
							<h4 style="float:left;">Bed Files:</h4><input type="file" id="inputbedfile" style="margin:5px" name="bed_file" multiple="multiple"/>
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
												Paste tab-separated genomic coordinates of SNPs of interest in .BED format. Pasting the data, or submitting one set of SNPs, restricts the analysis to the enrichment results only. To get the enrichment and epigenomic similarity heatmaps, upload multiple sets of SNPs as separate files using 'Choose Files' button.
											</label>
										</td>
									</tr>
								</table>
							</div>
						</div>
					</div>	
					<div class="well">
						<h3 style="float:left">2. Define the background: ${default_background}<img class="helptooltip" title="By default, all common SNPs are used as a 'population' to calculate the probability of SNPs in a set to be enriched with an epigenomic element" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
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
													The background, or population selection is critical for the correct p-value calculation. The background represents all SNPs being evaluated in a user's study. A set of SNPs of interest should be a subset of the background. If that is not the case, a warning is generated.<br><br>
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
						<h3 style="float: left;margin-right: 10px;margin-top: 1px;">3. Select genome annotation features:</h3>
						${custom_gfs}
						<br>

						<div id="accordGFS" class="accordion" style="padding-bottom: 1em;list-style:none; margin-top: 3em;" onClick="renderCheckBoxTree()"> 
							<h3  id="accordionheader"><a id='gfselheader' href="#" style="font-size:120%; height: 100%">Choose genome annotation features</a></h3>
							<div >	
								<ul>
									<li id="list-bedbackground">					
										Bed Files (.bed or .gz):
										<input type="file" id="inputgenomicfeaturefile" style="margin:5px" name="genomicfeature_file" multiple="multiple"/>
										<a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">What should my data  look like?</a>
									</li>				
								</ul>
								<div id="grfdroplist" style="display: table;">	
									<label>Enter genome annotation names:</label>
									<ol style="left: 9px;position: relative; padding-top:10px;padding-bottom:10px;">        
										<li id="grf-list" class="input-text">
											<input style="float: left;" type="text" value="" id="gfs" class"grf"/>
										</li>
										<div id="grfs-auto">
											<div class="default">Type the names of the genome annotation features to run</div>
											<ul id="feed">
											</ul>
										</div>
									</ol>	
									<div style="display: table-row;">
										<div id="divCheckBox" style="width: 70%; margin:10px; display: table-cell; verticle-align: top; visibility: hidden">
											<label >Select genome annotation features: </label>
											<div>
												<a class="btn" style="margin-top: 9px;"onClick="changeCheckedStateTreeView(true)">Expand All</a>
												<a class="btn" style="margin-top: 9px;" onClick="changeCheckedStateTreeView(false)">Collapse All</a>
												<a class="btn" style="margin-top: 9px;" onClick="treeviewCheckAll()">Select All</a>
												<a class="btn" style="margin-top: 9px" onClick="treeviewUncheckAll()">Deselect all</a>
											</div>
											<div id="treeview-outer" style="padding-top:10px;">
												<div id="treeview-inner" onClick="viewBoxClick()">

												</div>
											</div>
										</div>
										<div style="margin:10px; display: table-cell;  verticle-align: top; padding-left:10px">
											<label style="width: 100%; margin: 5px;">Genome annotation data (functional/regulatory/epigenomic data) are mirrored from and organized according to the UCSC genome database scheme. If you know names of the tracks you want to run enrichment analyses with, start typing their names in the search box. Or simply search for keywords, like 'H3K4me1', to see which tracks are available. Use checkboxes in the collapsible TreeView to select groups of epigenomic elements</label><br>
											<!-- <input type="checkbox" style="font-size:120%;margin-top:1em" name="run_annot">Run annotation analysis</input> -->
										</div>

									</div>
								</div>			 

							</div>
						</div>
					</div>	

					<div class="well">
						<table width="100%">
							<tr>

								<td id="td_submit" style="width:90px">
									<button class="btn btn-primary" type="submit" >Submit job</button>
									<img class="helptooltip" title="Submits the job for enrichment/epigenomic association analyses" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								</td>
								<td id="td_submit" style="width:170px">

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
		</div>
	</body>