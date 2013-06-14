// Generated by CoffeeScript 1.3.1
(function() {
  var Heatmap, getConditionNames, getGeneExpressions, isNumber;

  $(document).ready(function() {
    // Modified to work with GR data
    // matrix_gfs, matrix_fois, genes are given values in results.html
      var geneExpressionModel, genes, heatmap;      
      genes = d3.tsv.parse(matrix_data);
      matrix_gfs = matrix_gfs.split("\t")
      matrix_fois = matrix_fois.split("\t")
      if (genes.length == 0){ 
          heatmap = d3.select("#heatmap").append("svg").attr("width", 1195).attr("height", 500);
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
        conditionNames: getConditionNames(genes)
      });
      geneExpressionModel.set({
        geneNames: genes.map(function(gene) {
          return gene.gene_name;
        })
      });
      geneExpressionModel.set({
        geneExpressions: getGeneExpressions(genes, geneExpressionModel.get("conditionNames"))
      });
      geneExpressionModel.set({
        extent: d3.extent($.map(geneExpressionModel.get("geneExpressions"), function(item) {
          return item;
        }))
      });
      geneExpressionModel.set({
        clusters: genes.map(function(gene) {
          return gene.cluster;
        })
      });
      geneExpressionModel.set({
        clusterColor: d3.scale.category20()
      });
      return heatmap = new Heatmap({
        el: "#heatmap",
        model: geneExpressionModel
      });
  });

  Heatmap = Backbone.View.extend({
    initialize: function() {
      return this.render();
    },
    render: function() {
      var cell_size, clusterColor, clusters, columns, conditionNames, conditionNamesMargin, extent, geneExpressions, geneNames, geneNamesMargin, getRow, heatmap, heatmapColor, height, margin, rows, textScaleFactor, width, x, y, legend;
      geneExpressions = this.model.get("geneExpressions");
      conditionNames = this.model.get("conditionNames");
      geneNames = this.model.get("geneNames");
      extent = this.model.get("extent");
      clusters = this.model.get("clusters");
      clusterColor = this.model.get("clusterColor");
      heatmapColor = d3.scale.linear().domain([-5, 0, 5]).range(["#29DA29", "#EAEAEA", "#DF3435"]);
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
      width = cell_size * geneExpressions[0]v .length;
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
       // Legend stuff including the gradient used
          var legW = 20,legH = 10;
          var gradient = heatmap.append("svg:defs").append("svg:linearGradient").attr("id", "gradient")
                                  .attr("x1", "0%").attr("y1", "0%").attr("x2", "100%").attr("y2", "100%").attr("spreadMethod", "pad");
          gradient.append("stop").attr("offset", "0").attr("stop-color", "#29DA29");
          gradient.append("stop").attr("offset", "0.5").attr("stop-color", "#EAEAEA");
          gradient.append("stop").attr("offset", "1.0").attr("stop-color", "#DF3435");
          legend = d3.select("heatmap").append("rect").attr("class","legend").attr("x",100).attr("y",100)
                          .attr("width",100).attr("height",100).style("fill","url(#gradient)");

      cur_row = -1;
      getRow = function(row) {
        var divtooltip = d3.select("body").append("div")   
          .attr("class", "tooltip")               
          .style("opacity", 0);

        var cell;
        cur_row += 1;
        console.log(cur_row)
        return cell = d3.select(this).selectAll(".cell").data(row).enter().append("rect").attr("class", "cell").attr("x", function(d, i) {
          console.log("i: " + i)
          return x(i);
        }).attr("width", x.rangeBand()).attr("height", x.rangeBand()).text(function(d) {
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

            pvalue = 1/Math.pow(10,Math.abs(d));
            divtooltip .html("<p style=\"color:#C1C1C1; margin-top: 4px; font-size: 16px;\">p-value: " + pvalue.toExponential(4) + "<br>"+ "log: " + d +
              "<br>GF: " + matrix_gfs[i] + "<br>FOI: " + matrix_fois[cur_row] + "</p>")  
                .style("left", (d3.event.pageX) + "px")     
                .style("top", (d3.event.pageY - 28) + "px");    
            })                  
          .on("mouseout", function(d) {       
            divtooltip.transition()        
                .duration(50)      
                .style("opacity", 0)});          
      };

      rows = heatmap.selectAll(".row").data(geneExpressions).enter().append("g").attr("class", "row").attr("name", function(d, i) {
        return "gene_" + i;
      }).attr("transform", function(d, i) {
        return "translate(0," + y(i) + ")";
      }).each(getRow);
      return rows.append("text").attr("x", -6).attr("y", x.rangeBand() / 2).attr("dy", ".32em").attr("text-anchor", "end").text(function(d, i) {
        return geneNames[i];
      });
     
    }
  });

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

}).call(this);