<!DOCTYPE html>
<html lang="en">
	<head>
        <script type="text/javascript" src="static/js/jquery.js"></script>
        <script type="text/javascript">${script}</script>
	    <style>
         #shiny_gr {
          background:url("static/images/loading.gif") center center no-repeat;
         }

         .loader {
            position: fixed;
            left: 0px;
            top: 0px;
            width: 100%;
            height: 100%;
            z-index: 9999;
            background: url("static/images/loading.gif") 50% 50% no-repeat rgb(249,249,249);
        }
         </style>
    </head>
<body>
<iframe id="shiny_gr" style="border: none;height: 1000px; width: 100%" src="http://127.0.0.1:5161?job_id=${job_id}"frameborder="0"></iframe>
<!--<iframe id="example1" style="border: none;height: 1000px; width: 100%" src="http://162.216.114.51/shiny-gr?job_id=${job_id}"frameborder="0"></iframe>-->

</body>
<div class="loader"></div>
</html>
