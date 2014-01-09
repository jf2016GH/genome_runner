function getQueryVariable(variable) { 
		var query = window.location.search.substring(1); 
		var vars = query.split("&"); 
		for (var i = 0; i < vars.length; i++) { 
			var pair = vars[i].split("="); 
			if (pair[0] == variable) { 
				return pair[1]; 
			} 
		} 
		return null; 
	} 


function clear_foi_uploads(){
	document.getElementById("inputbedfile").disabled=true;
	var control = $("#inputbedfile");
	control.replaceWith( control = control.clone( true ) );
}

function enable_foi_uploads(){
	document.getElementById("inputbedfile").disabled=false;
}

$(document).ready(function() {




	$("#inputbeddata").text("chr14\t93644378\t93644379\trs1268843\t0\t+\nchr20\t31277093\t31277094\trs210135\t0\t+\nchr17\t41807330\t41807331\trs1513670\t0\t-\nchr7\t12713069\t12713070\trs10488226\t0\t+\nchr22\t27041282\t27041283\trs16982515\t0\t+\nchr19\t51361756\t51361757\trs17632542\t0\t+\nchr16\t53813366\t53813367\trs17817449\t0\t+\nchr10\t96058297\t96058298\trs3765524\t0\t+\nchr7\t16326645\t16326646\trs35681285\t0\t+\nchr13\t38531580\t38531581\trs9548119\t0\t+");
	$("#inputbackgrounddata").text("chr18\t75653522\t75653523\trs11661856\t0\t+\nchr4\t95834033\t95834034\trs1859156\t0\t+\nchr1\t26521139\t26521140\trs11809207\t0\t+\nchr6\t32127476\t32127477\trs3134950\t0\t-\nchr9\t22098573\t22098574\trs4977574\t0\t+\nchr11\t2936951\t2936952\trs16928809\t0\t+\nchr5\t71410355\t71410356\trs2199161\t0\t+\nchr7\t104503812\t104503813\trs10953454\t0\t+\nchr11\t102628051\t102628052\trs7924357\t0\t+\nchr12\t11962572\t11962573\trs7314811\t0\t+\nchr13\t73728138\t73728139\trs9600079\t0\t+\nchr6\t7102083\t7102084\trs675209\t0\t+\nchr2\t54718180\t54718181\trs4557020\t0\t+\nchr5\t38452893\t38452894\trs7715172\t0\t+\nchr18\t11759431\t11759432\trs9947295\t0\t+\nchr5\t88099950\t88099951\trs4521516\t0\t+\nchr1\t34186192\t34186193\trs528059\t0\t+\nchr8\t101841545\t101841546\trs3108919\t0\t+\nchr19\t17590280\t17590281\trs11666579\t0\t+\nchr7\t28189410\t28189411\trs1635852\t0\t+\n");
	
	// Allows the radio buttons for the FOI to be sent as a POST request
	$('div.btn-group[data-toggle-name]').each(function () {
	    var group = $(this);
	    var form = group.parents('form').eq(0);
	    var name = group.attr('data-toggle-name');
	    var hidden = $('input[name="' + name + '"]', form);
	    $('button', group).each(function () {
	        $(this).on('click', function () {
	            hidden.val($(this).data("toggle-value"));
	        });
	    });
	});

	$('#org').change(function() {;
		var org = document.getElementById("org").value;
		var temp = org.split(":");
		var $url="index?organism=" + temp[1];
		window.location = $url;})

	// Set organism value
	if (getQueryVariable('organism') != null){
		$('#org').val("organism:" + getQueryVariable('organism'));}
	else{$('#org').val("organism:${default_organism}");}

	// Crete autocomplete text box 
	$.facebooklist('#gfs', '#preadded', '#grfs-auto',{url:"data/"+$('#org').val().split(":")[1]+"/gfs.php",cache:1}, 10, {userfilter:1,casesensetive:1});

	$(function() {
		$(".accordion").accordion({
			collapsible: true,
			active: false,
			autoHeight: false
		});
	});

	// The tooltips
	/*$('.helptooltip').tipsy({
		gravity: 'e',html: true,
		title: function(){
			var d = $(this).attr("title");
			return '<p>'+d+'</p>';
		}				
	})*/

	// the viewport for the treeview
	$("#treeview-inner").slimScroll({
	    	height: $("#ucsc").height(),
	    	width: '400px',
	    	alwaysVisible: true
	});

	$.post('/gr/get_checkboxtree?organism='+$("select[name='organism']").text(),function(data){
		$('#divCheckBox').html(data)
		$('#ucsc').checkboxTree({
			initializeChecked: 'collapsed',
			initializeUnchecked: 'collapsed',
			collapseEffect: '',
			expandEffect: ''
		});
		// Make the checkbox tree visible
		$('#divCheckBox').css('visibility','visible')
		$('#gfselheader').text('Choose genome annotation features')
	});


	$("#accfoi").bind('accordionchange',
			function () {
				if ($(this).find('.ui-state-active').length == 0){
					enable_foi_uploads();
					$("#btngroup_demo_fois").children().each(function(x){ 
						 $(this).prop('disabled', false);
					});						
					document.getElementById("inputbeddata").disabled=true;
					document.getElementById("inputbeddata").text = ""
					}
				else {
					clear_foi_uploads();
					$('#demo_fois_none').button('toggle')
					$("#btngroup_demo_fois").children().each(function(x){ 
						 $(this).prop('disabled', true);
					});
					document.getElementById("inputbeddata").disabled=false;
					document.getElementById("inputbeddata").text = ""
				}
			});

	$("#accback").bind('accordionchange',
			function () {
				if ($(this).find('.ui-state-active').length == 0){
					document.getElementById("inputbackgroundfile").disabled=false;
					document.getElementById("inputbackgrounddata").disabled=true;
				}
				else {
					document.getElementById("inputbackgroundfile").disabled=true;
					var control = $("#inputbackgroundfile");
					control.replaceWith( control = control.clone( true ) );
					document.getElementById("inputbackgroundfile").disabled=true;
					document.getElementById("inputbackgrounddata").disabled=false;
				}
			});
		
	});

	function viewBoxClick(){
		var minHeight = 300;
		var maxHeight = 500;			
		// resize to the expanded treeview as long as it is less than 500			
		if ($("#ucsc").height() <= maxHeight) {
			var oldHeight = $(".slimScrollDiv").height();
			//$(".slimScrollDiv").css("height", $('#ucsc').height()+'px');
			$(".slimScrollDiv").height($('#ucsc').height());
			$("#treeview-inner").height($("#treeview-inner").height() + $(".slimScrollDiv").height() - oldHeight);

		}
		else if ($("#ucsc").height() > maxHeight){
			//$(".slimScrollDiv").css("height",maxHeight + 'px');
			$(".slimScrollDiv").height(maxHeight);
			$("#treeview-inner").height(maxHeight);
		}

	}
	
	//prevents user from submitting form by accidentally pressing the enter key
	function stopRKey(evt) { 
		var evt = (evt) ? evt : ((event) ? event : null); 
		var node = (evt.target) ? evt.target : ((evt.srcElement) ? evt.srcElement : null); 
		if ((evt.keyCode == 13) && (node.type=="text"))  {return false;} 
	} 

	function changeCheckedStateTreeView(checkAll){
		if (checkAll == true){
			$('#ucsc').checkboxTree('expandAll');
			viewBoxClick();
		}
		else {
			$('#ucsc').checkboxTree('collapseAll');
			viewBoxClick();
		}
	}

	function treeviewCheckAll(){
		$("#ucsc").children().find(':checkbox').prop("checked", true);
	}
	function treeviewUncheckAll(){
		$("#ucsc").children().find(':checkbox').prop("checked", false);
	}

		document.onkeypress = stopRKey; 