 var Heatmap, getConditionNames, getGeneExpressions, isNumber,cur_heatmap,cur_tooltip_matrix, tooltip_matrices = {};
 // ### helper functions
   getGeneExpressions = function(genes, conditionNames) {
      return genes.map(function(gene) {
        return conditionNames.map(function(condition) {

          return +gene[condition];
        });
      });
    };

    getConditionNames = function(genes) {
      return Object.keys(genes[0]).filter(function(columnName) {
        return !columnName.match(/cluster/) && isNumber(genes[1][columnName]);
      });
    };

    isNumber = function(n) {
      return !isNaN(parseFloat(n)) && isFinite(n);
    };
 // ###
 Heatmap = Backbone.View.extend({

    initialize: function() {      
     
      return this.render();
    },
    render: function() {    
      var cell_size, clusterColor, clusters, columns, conditionNames, conditionNamesMargin, extent, geneExpressions, geneNames, geneNamesMargin, getRow, heatmap, heatmapColor, height, margin, rows, textScaleFactor, width, x, y, legend;
      geneExpressions = this.model.get("geneExpressions");
      conditionNames = this.model.get("conditionNames");
      geneNames = this.model.get("geneNames");
      tooltip_matrices[cur_heatmap] = getGeneExpressions(cur_tooltip_matrix, conditionNames); // for tooltips
      extent = this.model.get("extent");
      clusters = this.model.get("clusters");
      clusterColor = this.model.get("clusterColor");
      heatmapColor = d3.scale.linear().domain([-1*color_range, 0, color_range]).range(["#29DA29", "#EAEAEA", "#DF3435"]);
      textScaleFactor = 9;
      conditionNamesMargin = d3.max(conditionNames.map(function(conditionName) {
        return conditionName.length;
      }));
      geneNamesMargin = d3.max(geneNames.map(function(geneName) {
        return geneName.length;
      }));
      margin = {
        top: conditionNamesMargin * textScaleFactor,
        right: 150,
        bottom: conditionNamesMargin * textScaleFactor,
        left: geneNamesMargin * textScaleFactor
      };
      cell_size = 30;
      width = cell_size * geneExpressions[0].length;
      height = cell_size * geneNames.length;
      heatmap = d3.select(this.el).append("svg").attr("width", width + margin.right + margin.left).attr("height", height + margin.top + margin.bottom).attr("id", "heatmap").append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
      x = d3.scale.ordinal().domain(d3.range(geneExpressions[0].length)).rangeBands([0, width]);
      y = d3.scale.ordinal().domain(d3.range(geneNames.length)).rangeBands([0, height]);
      columns = heatmap.selectAll(".column").data(conditionNames).enter().append("g").attr("class", "column").attr("transform", function(d, i) {
        return "translate(" + x(i) + ")rotate(-90)";
      });
      columns.append("text").attr("x", 6).attr("y", x.rangeBand() / 2).attr("dy", "-.5em").attr("dx", ".5em").attr("text-anchor", "start").attr("transform", "rotate(45)").text(function(d, i) {
        return conditionNames[i];
      });
       

      cur_row = -1;
      getRow = function(row) {
        var divtooltip = d3.select("body").append("div")   
          .attr("class", "tooltip")               
          .style("opacity", 0);

        var cell;
        cur_row += 1;
        return cell = d3.select(this).selectAll(".cell").data(row).enter().append("rect").attr("class", "cell").attr("x", function(d, i) {
          return x(i);
        }).attr("width", x.rangeBand()).attr("height", x.rangeBand()).attr("row",matrix[cur_row].gene_name)
        .text(function(d,i) {
          this.setAttribute("tt",tooltip_matrices[cur_heatmap][cur_row][i]);
          this.setAttribute("col",conditionNames[i]);
          return d;
        }).style("fill", function(d) {
          return heatmapColor(d);
      // Tool tips
        }).attr("value",function(d){
          return d
        })
        .on("mouseover", function(d,i) {             
    
            divtooltip.transition()        
                .duration(50)      
                .style("opacity", .9);

            if (log10_tt_value == true) {

               d  = parseFloat(this.getAttribute("value"));
               if (d<0) { d = -1*(1/Math.pow(10,Math.abs(d))); }
               else { d = 1/Math.pow(10,Math.abs(d)); }
               
               if (Math.abs(d) < .05) { 
                  if (d<0) tip1 = "Underrepresented<br>";
                  else tip1 = "Overrepresented<br>";
               }
               else{
                 tip1 = "Not Significant<br>";
               }             
               tip2 = "P-value: "  + d.toExponential(4) + "<br>";               
             }
             else {
                tip1 = parseFloat(this.getAttribute("value")).toPrecision(2);
                tip1 = "Pearsons: " + tip1 + "<br>"
                tip2 = parseFloat(this.getAttribute("tt")).toExponential(4);
                tip2 = "P-value: " + tip2 + "<br>"
             }
            divtooltip .html("<p style=\"color:#C1C1C1; margin-top: 4px; font-size: 16px;\">" + tip1 + tip2 + 
                            "Row: " + this.getAttribute("row")  + "<br>" +
                            "Column: " + this.getAttribute("col") + "</p>")  

                .style("left", (d3.event.pageX+30) + "px")     
                .style("top", (d3.event.pageY - 28) + "px");    
            })                  
          .on("mouseout", function(d) {       
            divtooltip.transition()        
                .duration(50)      
                .style("opacity", 0)});          
      };
      var divtooltip = d3.select("body").append("div")   
                        .attr("class", "tooltip")               
                        .style("opacity", 0);   
      rows = heatmap.selectAll(".row").data(geneExpressions).enter().append("g").attr("class", "row").attr("name", function(d, i) {
        return "gene_" + i;
      }).attr("transform", function(d, i) {
        return "translate(0," + y(i) + ")";
      }).each(getRow);
      return rows.append("text").attr("x", -6).attr("y", x.rangeBand() / 2).attr("dy", ".32em").attr("text-anchor", "end").text(function(d, i) {
        desc = matrix_data_gf_description[i]
        if (desc == "") this.setAttribute("tt", "No Description");
        else this.setAttribute("tt", desc);
        return geneNames[i];
      }).on("mouseover", function(d,i) {
            if (cur_heatmap == "heatmap") {              
              divtooltip.transition()        
                  .duration(50)      
                  .style("opacity", .9);
              
              divtooltip .html("<p style=\"color:#C1C1C1; margin-top: 4px; font-size: 16px;\">" + this.getAttribute("tt") + "</p>")  

                  .style("left", (d3.event.pageX+30) + "px")     
                  .style("top", (d3.event.pageY - 28) + "px");    
            }
          })                  
          .on("mouseout", function(d) {  
           if (cur_heatmap == "heatmap") {          
            divtooltip.transition()        
                .duration(50)      
                .style("opacity", 0)}
              }); 
          }
  });

function generate_heatmaps() {
      // The code that create the heatmaps
      $(document).ready(function() {
        matrix_data_gf_description = matrix_data_gf_description.split("\t");
        var geneExpressionModel, genes, heatmap;
        color_range = 1;
        matrix_cor_pvalues = d3.tsv.parse(matrix_cor_pvalues)
        log10_tt_value = false;
        matrix_cor = d3.tsv.parse(matrix_cor)
        matrix = matrix_cor;
        cur_heatmap = "heatmap_cor";
        cur_tooltip_matrix = matrix_cor_pvalues;    
        create_heatmap("#heatmap_cor", matrix,matrix_cor);
        color_range = 5; 
        matrix_enrich_data = d3.tsv.parse(matrix_data);
        matrix = matrix_enrich_data
        log10_tt_value = true;
        cur_heatmap = "heatmap";
        cur_tooltip_matrix = matrix_enrich_data;
        create_heatmap("#heatmap",matrix);
        // Creates links to download the svg file
        d3.selectAll("#heatmap_download")
            .attr("href", "data:image/svg+xml;charset=utf-8;base64," + 
              btoa(unescape(encodeURIComponent(
                d3.selectAll("#heatmap").selectAll("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg")
               .node().parentNode.innerHTML)
                )
              )
            );
        d3.selectAll("#heatmap_cor_download")
            .attr("href", "data:image/svg+xml;charset=utf-8;base64," + 
              btoa(unescape(encodeURIComponent(
                d3.selectAll("#heatmap_cor").selectAll("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg")
               .node().parentNode.innerHTML)
                )
              )
        );
      });
 

 

}

  create_heatmap = function(target,matrix){
    // Modified to work with GR data
    // matrix_gfs, matrix_fois, matrix are given values in results.html

      if (matrix.length == 0){ 
          heatmap = d3.select("#heatmap").append("svg:heatmap").attr("width", 1195).attr("height", 500);
          heatmap.append("svg:rect")
          .attr("width", 1195)
          .attr("height", 500)
          .attr("fill","white")
          heatmap.append("svg:text")
          .attr("x", 500)
          .attr("y", 200)
          .style("font-size","20px")
          .attr("text-anchor", "middle")
          .text(matrix_data.replace("gene_name",""));
          return; 
      }
      geneExpressionModel = new Backbone.Model;
      geneExpressionModel.set({
        conditionNames: getConditionNames(matrix)
      });
      geneExpressionModel.set({
        geneNames: matrix.map(function(gene) {
          return gene.gene_name;
        })
      });
      geneExpressionModel.set({
        geneExpressions: getGeneExpressions(matrix, geneExpressionModel.get("conditionNames"))
      });
      geneExpressionModel.set({
        extent: d3.extent($.map(geneExpressionModel.get("geneExpressions"), function(item) {
          return item;
        }))
      });

      geneExpressionModel.set({
        clusters: matrix.map(function(gene) {
          return gene.cluster;
        })
      });
      geneExpressionModel.set({
        clusterColor: d3.scale.category20()
      });
      return heatmap = new Heatmap({
        el: target,
        model: geneExpressionModel
      });
    }





   

  // changes the tooltip matrix so that the correct values can be displayed by the tooltips for the 
  // currently active matrix
  change_active_heatmap = function(heatmap_name){
    cur_heatmap = heatmap_name;
    console.log(cur_heatmap);
    if (heatmap_name == "heatmap"){      
      log10_tt_value = true;
    }
    else if (heatmap_name == "heatmap_cor"){
      log10_tt_value = false;
    }
  }

 