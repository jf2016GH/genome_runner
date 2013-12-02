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
		</style>

		<!-- For IE6-8 support of HTML5 elements -->
		<!--[if lt IE 9]>
		<script src="http://html5shim.googlecode.com/svn/trunk/html5.js"></script>
		<![endif]-->
		<link type="text/css" href="static/css/ui-lightness/jquery-ui-1.8.16.custom.css" rel="stylesheet" />
		<link type="text/css" href="static/css/main.css"/>
		<script type="text/javascript" src="static/js/jquery.js"></script>
	 	<script src="static/js/bootstrap.min.js"></script>
		<script scr="static/js/main.js"></script>
		
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
					<li><a href="./cite">How to Cite</a></li>
					<li><a href="http://sourceforge.net/projects/genomerunner/">SourceForge</a></li>
					<li><a href="./roadmap">Roadmap</a></li>
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
