<!DOCTYPE html>
<html lang="en">
	<head>
		<link href="static/css/bootstrap.min.css" rel="stylesheet">
		<style>
			body {
				padding-top: 60px; /* 60px to make the container go all the way to the bottom of the topbar */
			}
			/* allows the check box tree to be rendered correctly*/
			.ui-icon {float: left;}
			div.tooltip {   
				  position: absolute;           
				  text-align: center;                
				  padding: 5px;             
				  font: 10px sans-serif;        
				  background: #1C1C1C;   
				  border: 0px;      
				  border-radius: 4px;           
				  pointer-events: none;
				  color: #FFFFFF;         
			}
			#content{
			margin:0px auto;
			text-align:left;
			padding:15px;
			}
			#ucsc {
				padding-left: 0px;
			}
			#treeview-inner{
				margin-left: -40px;
			}
			label {
				font-size: 120%;
			}
			.ui-widget-content { background: #ffffff url(images/ui-bg_glass_75_ffffff_1x400.png) 50% 50% repeat-x; color: #404040; }
			#divCheckBox { border: 2px solid #C0C0C0;
						   padding: 3px; }
			/* allows the check box tree to be rendered correctly*/
			.ui-icon {float: left;}
			.popover {background: whiteSmoke;border: 1px solid #5E5E5E;}
			.ui-accordion-content.ui-helper-reset.ui-widget-content.ui-corner-bottom.ui-accordion-content-active {heigh: 100%;}
		</style>

		<!-- For IE6-8 support of HTML5 elements -->
		<!--[if lt IE 9]>
		<script src="http://html5shim.googlecode.com/svn/trunk/html5.js"></script>
		<![endif]-->
		<link type="text/css" href="static/css/demo_table.css" rel="stylesheet" />

		<link type="text/css" href="static/css/ui-lightness/jquery-ui-1.8.16.custom.css" rel="stylesheet" />
		<link type="text/css" href="static/css/main.css"/>
		<link rel="stylesheet" href="static/css/combobox.css" type="text/css" media="screen" title="Test Stylesheet" charset="utf-8" />
		<link type="text/css" href="static/css/TableTools_JUI.css" rel="stylesheet" />

		<script type="text/javascript" src="static/js/jquery.js"></script>
		<script type="text/javascript" src="static/jui/jquery-ui-1.8.12.custom/js/jquery-ui-1.8.12.custom.min.js"></script>
		<script type="text/javascript" src="static/js/jquery.dataTables.min.js"></script>
		<script type="text/javascript" src="static/js/ZeroClipboard.js"></script>

		<script type="text/javascript" src="static/js/TableTools.js"></script>
		<script src="static/js/bootstrap.min.js"></script>
		<script src="static/js/bootstrap-tabs.js"></script>
		<script scr="static/js/main.js"></script>
		<script type="text/javascript" src="static/js/underscore.js"></script>
      	<script type="text/javascript" src="static/js/backbone.js"></script>
      	<script src="static/js/d3.js" charset="utf-8"></script>
      	<script src='static/js/d3.json_heatmap.js' type="text/javascript"></script>
      	<script type="text/javascript" src="static/js/jquery.checkboxtree.js"></script>
      	<script src="static/js/fcbkcomplete.min.js" type="text/javascript" charset="utf-8"></script>
		<script src="static/js/jquery.slimScroll.js" type="text/javascript" charset="utf-8"></script>
		<script type="text/javascript">${script}</script>
</head>
<body bgcolor="#A8D5FF">

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
			</div>
		</div>
	</div>
<div class="well">
	${body}
</div>

</body>
</html>