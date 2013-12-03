<h3>
			<span style="font-family:arial,helvetica,sans-serif;">How GenomeRunner Web can help me?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">GenomeRunner Web helps to interpret the regulatory effect of SNPs by identifying functional elements most statistically significantly co-localized with the SNPs of interest (see <a href="#enrichment">Enrichment analysis</a>). These regulatory elements may provide understanding which mechanisms (and in which cell type) may be affected by the SNPs of interest.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">If one analyzes three or more sets of SNPs, GenomeRunner Web visualizes their epigenomic similarity (see <a href="#episimilarity">Epigenomic Similarity analysis</a>). This information may be used to classify the sets of SNPs by their effect upon regulatory landscape.</span></div>
		<div>
			&nbsp;</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">How my input data should look like?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Genomic coordinates of SNPs of interest in <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1" target="_blank">.BED format</a>. One can upload a file, or copy-paste tab-separated coordinates. Multiple file upload/analysis is possible (e.g., individual-specific sets of SNPs).</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Several pre-defined sets of SNPs are available.&nbsp;</span><span style="font-family: arial, helvetica, sans-serif;">&nbsp;For&nbsp;</span><em style="font-family: arial, helvetica, sans-serif;">Homo Sapiens</em><span style="font-family: arial, helvetica, sans-serif;">&nbsp;these include:</span></div>
		<div>
			&nbsp;</div>
		<div>
			<table border="1" cellpadding="1" cellspacing="1" style="width: 1000px;">
				<tbody>
					<tr>
						<td style="text-align: center;">
							<strong><span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">Pre-defined sets of SNPs</span></span></strong></td>
						<td style="text-align: center;">
							<strong><span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">When to use</span></span></strong></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">gwasDemo - several disease-associated SNP sets from <a href="http://www.genome.gov/26525384" target="_blank">GWAScatalog</a>.</span></span></td>
						<td>
							<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">Used for <a href="./demo">Demo</a> run</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">gwasRand - randomly selected sets of SNPs from&nbsp;&nbsp;<a href="http://www.genome.gov/26525384" target="_blank">GWAScatalog</a>.</span></span></td>
						<td>
							<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">Use them for <a href="./demo">Demo</a> run, to ensure random SNPs do not show significant associations</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">snp138Rand - randomly selected sets of SNPs from the <a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=snp138&amp;hgTracksConfigPage=configure" target="_blank">snp138 </a>database.</span></span></td>
						<td>
							<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">Use for &#39;sanity check&#39;, with snp138+ background and any (category of) genome annotation elements, to ensure random SNPs do not show significant associations</span></span></td>
					</tr>
				</tbody>
			</table>
		</div>
		<div>
			&nbsp;</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">What is &ldquo;background&rdquo;?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">The background is the &ldquo;universe&rdquo; of all SNPs assessed in a study, from which SNPs of interest came from. Several pre-defined backgrounds are provided,&nbsp;</span><span style="font-family: arial, helvetica, sans-serif;">&nbsp;for&nbsp;</span><em style="font-family: arial, helvetica, sans-serif;">Homo Sapiens</em><span style="font-family: arial, helvetica, sans-serif;">&nbsp;these include:</span></div>
		<div>
			&nbsp;</div>
		<div>
			<table border="1" cellpadding="1" cellspacing="1" style="width: 1000px;">
				<tbody>
					<tr>
						<td style="text-align: center;">
							<strong><span style="font-family:arial,helvetica,sans-serif;">Pre-defined background</span></strong></td>
						<td style="text-align: center;">
							<strong><span style="font-family:arial,helvetica,sans-serif;">When to use</span></strong></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">snp138+ (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=snp138&amp;hgTracksConfigPage=configure" target="_blank">All&nbsp;Simple Nucleotide Polymorphisms (dbSNP 138)</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">For SNPs from whole-genome GWA studies</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">snp138Common+ (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=snp138Common&amp;hgTracksConfigPage=configure" target="_blank">Simple Nucleotide Polymorphisms (dbSNP 138) Found in &gt;= 1% of Samples</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">For SNPs from studies where rare variants were ignored</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">gwascatalog+ (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=gwasCatalog&amp;hgTracksConfigPage=configure" target="_blank">NHGRI Catalog of Published Genome-Wide Association Studies</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">For demo testing, to observe regulatory associations of disease-specific sets of SNPs, as compared with randomly selected from all gwascatalog</span></td>
					</tr>
				</tbody>
			</table>
		</div>
		<div>
			&nbsp;</div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">For a GWAS, the background is likely to be all SNPs (snp138+). For a study using microarrays, the background would contain coordinates of all SNPs on the microarray - upload or copy/paste them.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;"><strong>Important:</strong> SNPs of interest should be a subset of the background SNPs. If some SNPs of interest do not overlap the background, a non-critical error is issued. Create a custom background that includes all SNPs of interest.</span></div>
		<div>
			<div>
				<span style="font-family: arial, helvetica, sans-serif;"><strong>Important:</strong>&nbsp;Ensure the <em>end</em> coordinates of the SNPs of interest and the background SNPs are <em>start+1</em>. This is necessary for correct calculations of overlap statistics.</span></div>
			<div>
				&nbsp;</div>
		</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">What are &ldquo;genome annotation features&rdquo;?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Genome annotations are discrete regions potentially having functional/regulatory properties.&nbsp;</span><span style="font-family: arial, helvetica, sans-serif;">&nbsp;</span></div>
		<div>
			&nbsp;</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">There are too many genome annotation features! What to choose?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">The genome annotation features are organized by categories mirrored from the UCSC genome browser. Using search box and/or checkboxes, one can select one or more genome annotation categories. Clicking on a genome annotation&rsquo; name will bring up description, if available.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">The ENCODE data are organized by the source/data type, tiers (quality), and by cell types.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;"><strong>Hint:</strong> Several well-known genome annotation features sets are brought forward as &ldquo;default genome annotation features&rdquo;. For <em>Homo Sapiens</em> these include:</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">&nbsp;</span></div>
		<div>
			<table border="1" cellpadding="1" cellspacing="1" style="width: 1000px;">
				<thead>
					<tr>
						<th scope="col">
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">Genome annotation category</span></span></th>
						<th scope="col">
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">Experimental question:&nbsp;Are the SNPs of interest...</span></span></th>
					</tr>
				</thead>
				<tbody>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">coriellVariants (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=coriellDelDup&amp;hgTracksConfigPage=configure" style="font-family: arial, helvetica, sans-serif;" target="_blank">Coriell Cell Line Copy Number Variants</a>, split by cell types)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... enriched in CNVs, and in which cell type?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">dgvVariants (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=dgvPlus&amp;hgTracksConfigPage=configure" style="font-family: arial, helvetica, sans-serif;" target="_blank">Database of Genomic Variants: Structural Variation (CNV, Inversion, In/del)</a>, split by variant type)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... enriched in CNVs, or other types of structural variations?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">repeats (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=rmsk&amp;hgTracksConfigPage=configure" style="font-family: arial, helvetica, sans-serif;" target="_blank">Repeating Elements by RepeatMasker</a>, split by repeat class)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... enriched in regions of low complexity, and in which type?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">gwasCatalog (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=gwasCatalog&amp;hgTracksConfigPage=configure" style="font-family: arial, helvetica, sans-serif;" target="_blank">NHGRI Catalog of Published Genome-Wide Association Studies</a>, split by disease/trait types)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... enriched in known disease-specific SNPs?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">altSplicing (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=knownAlt&amp;hgTracksConfigPage=configure" target="_blank">Alternative Splicing, Alternative Promoter and Similar Events in UCSC Genes</a>, split by splicing type)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... potentially disrupt a specific type of alternative spliced regions?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">ncRNAs (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=wgRna&amp;hgTracksConfigPage=configure" style="font-family: arial, helvetica, sans-serif;" target="_blank">C/D and H/ACA Box snoRNAs, scaRNAs, and microRNAs from snoRNABase and miRBase</a>, split by ncRNA type)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... associated with a class of non-coding elements?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">tfbsEncode (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=wgEncodeRegTfbsClusteredV3&amp;hgTracksConfigPage=configure" target="_blank">Transcription Factor ChIP-seq Clusters V3 (161 targets, 189 antibodies) from ENCODE</a>, split by TFBS name)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... potentially disrupt a specific <u>experimentally defined</u> transcription factor binding site?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">tfbsConserved (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=tfbsConsSites&amp;hgTracksConfigPage=configure">HMR Conserved Transcription Factor Binding Sites</a>, split by TFBS name)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... potentially disrupt a specific&nbsp;<u>computationally predicted</u>&nbsp;transcription factor binding site?</span></span></td>
					</tr>
					<tr>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">chromStates (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=wgEncodeBroadHmm&amp;hgTracksConfigPage=configure" style="font-family: arial, helvetica, sans-serif;" target="_blank">Chromatin State Segmentation by HMM from ENCODE/Broad</a>, Gm12878 cell line, split by chromatin state type)</span></span></td>
						<td>
							<span style="font-size:14px;"><span style="font-family:arial,helvetica,sans-serif;">... preferentially located in certain chromatin states?</span></span></td>
					</tr>
				</tbody>
			</table>
			<p>
				<span style="font-family: arial, helvetica, sans-serif;">Examples of what to choose:</span></p>
		</div>
		<ul>
			<li>
				<span style="font-family:arial,helvetica,sans-serif;">Select &ldquo;tfbsEncode&rdquo; category to get an answer whether the SNPs of interest are enriched in any of the 161 transcription factor binding sites identified by ChIP-seq.</span></li>
			<li>
				<span style="font-family:arial,helvetica,sans-serif;">Select &ldquo;ENCODE/BroadHistone/Tier1/Gm12878&rdquo; category to get an insight whether the SNPs of interest are enriched in histone marks in B lymphocytes.</span></li>
			<li>
				<span style="font-family:arial,helvetica,sans-serif;">Select &ldquo;genes&rdquo; category to answer a question, whether the SNPs of interest are enriched in genes/exons.</span></li>
		</ul>
		<div>
			&nbsp;</div>
		<div>
			&nbsp;</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;"><a name="enrichment"></a>What is Enrichment analysis?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Enrichment analysis answers the question whether a set of SNPs of interest collectively enriched or depleted in regulatory regions, as compared with randomly selected set of SNPs. GenomeRunner Web performs this analysis versus each selected (group of) regulatory tracks, and prioritized over- and underrepresented associations by p-value.</span></div>
		<div>
			&nbsp;</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;"><a name="episimilarity"></a>What is Epigenomic similarity analysis?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Epigenomic similarity analysis visualizes similarity among enrichment profiles for three or more sets of SNPs of interest. It answers the question whether different sets of SNPs are enriched in similar epigenomic elements, hence, may exhibit similar regulatory changes. This analysis may be useful when comparing sets of SNPs from multiple GWA studies, or individual-specific variations. It is particularly useful for the analysis of individual-specific rare variants, to observe whether sets of individual-specific SNPs may affect similar epigenomic elements.</span></div>
		<h3>
			<span style="font-family: arial, helvetica, sans-serif;">I still have questions/suggestions/bug report...</span></h3>
		<p>
			<font face="arial, helvetica, sans-serif">Please, contact&nbsp;</font><br />
			<img height="25%" src="static/images/e-mail.png" width="25%" /></p>
		<div>
			&nbsp;</div>