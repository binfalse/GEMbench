// based on https://bl.ocks.org/kerryrodden/766f8f6d31f645c39f488a0befa1e3c8

// Dimensions of sunburst.
var width = 750;
var height = 600;
var radius = Math.min(width, height) / 2;

// Breadcrumb dimensions: width, height, spacing, width of tip/tail.
var b = {
  w: 10, h: 30, s: 3, t: 10
};

// Mapping of step names to colors.
var colors = {
  "Patient": "#004697",
  
  "Cell-line": "#d7191c",
  
  "Microarray (EMTAB-37)": "#d7301f",
  "RNA-Seq (HPA)": "#fc8d59",
  "Mass Spec (ProteomeNCI60)": "#ffd288",
  
  "Microarray (GSE2109)": "#0570b0",
  "RNA-Seq (TCGA)": "#74a9cf",
  "Mass Spec (ProteomePatients)": "#bdc9e1"
};
var fgColors = {
  "Patient": "#fff",
  
  "Cell-line": "#fff",
  
  "Microarray (EMTAB-37)": "#fff",
  "RNA-Seq (HPA)": "#000",
  "Mass Spec (ProteomeNCI60)": "#000",
  
  "Microarray (GSE2109)": "#fff",
  "RNA-Seq (TCGA)": "#000",
  "Mass Spec (ProteomePatients)": "#bdc9e1"
};



function colorMap (node) {
  if (node.data.name == "root")
    return "#000";
  if (node.data.name == "Proteome") {
    if (node.parent.data.name == "Patient")
      return "#1e91ca"
    else
      return "#ff9a78"
  }
  if (node.data.name == "Transcriptome") {
    if (node.parent.data.name == "Patient")
      return "#1e63ca"
    else
      return "#e61500"
  }
  //console.log (colors[node.data.name])
  if (colors[node.data.name])
    return colors[node.data.name];
  //console.log ("not found: ", node.data.name)
  return colors[node.parent.data.name];
}
function fgColorMap (node) {
  if (node.data.name == "root")
    return "#000";
  if (node.data.name == "Proteome") {
    if (node.parent.data.name == "Patient")
      return "#fff"
    else
      return "#000"
  }
  if (node.data.name == "Transcriptome") {
    if (node.parent.data.name == "Patient")
      return "#fff"
    else
      return "#fff"
  }
  //console.log (colors[node.data.name])
  if (fgColors[node.data.name])
    return fgColors[node.data.name];
  //console.log ("not found: ", node.data.name)
  return fgColors[node.parent.data.name];
}


// Total size of all segments; we set this later, after loading the data.
var totalSize = 0; 

var vis = d3.select("#sunburst").append("svg:svg")
    .attr("width", width)
    .attr("height", height)
    .append("svg:g")
    .attr("id", "container")
    .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

var partition = d3.partition()
    .size([2 * Math.PI, radius * radius]);

var arc = d3.arc()
    .startAngle(function(d) { return d.x0; })
    .endAngle(function(d) { return d.x1 - 0.002; })
    .innerRadius(function(d) { return Math.sqrt(d.y0) + 1; })
    .outerRadius(function(d) { return Math.sqrt(d.y1) - 1; });

  console.log ("test");
d3.dsv(";","/data/data_viz_sunburst.csv", function(d) {
    //return //Object.keys(d).map(function(key){return d[key];});
        return [
          d.Source,
          d.Location,
          d.Type,
          d.Cancer,
          d.Samples
        ];
   }).then (function(text) {
  console.log ("test2");
  console.log (text);
  //var dsv = d3.dsvFormat(';')
  //csv = dsv.parseRows(text)
  ////var csv = d3.csvParseRows(text);
  //console.log (csv);
  var json = buildHierarchy(text);
  console.log (json);
  createVisualization(json);
})
  .catch(error => console.log(error));;

// Main function to draw and set up the visualization, once we have the data.
function createVisualization(json) {

  // Basic setup of page elements.
  //initializeBreadcrumbTrail();
  //drawLegend();
  //d3.select("#togglelegend").on("click", toggleLegend);

  // Bounding circle underneath the sunburst, to make it easier to detect
  // when the mouse leaves the parent g.
  vis.append("svg:circle")
      .attr("r", radius)
      .style("opacity", 0);

  // Turn the data into a d3 hierarchy and calculate the sums.
  var root = d3.hierarchy(json)
      .sum(function(d) { return d.size; })
      .sort(function(a, b) { return b.value - a.value; });
  
  // For efficiency, filter nodes to keep only those large enough to see.
  var nodes = partition(root).descendants()
      .filter(function(d) {
          return (d.x1 - d.x0 > 0.005); // 0.005 radians = 0.29 degrees
      });
      
      console.log (nodes);

  var path = vis.data([json]).selectAll("path")
      .data(nodes)
      .enter().append("svg:path")
      .attr("display", function(d) { return d.depth ? null : "none"; })
      .attr("d", arc)
      .attr("fill-rule", "evenodd")
      //.style("fill", function(d) { return colors[d.data.name]; })
      .style("fill", function(d) { return colorMap (d); })
      .style("opacity", 1)
      .on("mouseover", mouseover);

  // Add the mouseleave handler to the bounding circle.
  //d3.select("#container").on("mouseleave", mouseleave);

  // Get total size of the tree = value of root node from partition.
  totalSize = path.datum().value;
 };


// Fade all but the current sequence, and show it in the breadcrumb trail.
function mouseover(d) {

  var percentage = (100 * d.value / totalSize).toPrecision(3);
  var percentageString = percentage + "%";
  if (percentage < 0.1) {
    percentageString = "< 0.1%";
  }

  d3.select("#percentage")
      .text(percentageString);

  d3.select("#sunburst-explanation")
      .style("visibility", "");

  var sequenceArray = d.ancestors().reverse();
  sequenceArray.shift(); // remove root node from the array
  updateBreadcrumbs(sequenceArray, percentageString);

  // Fade all the segments.
  d3.selectAll("path")
      .style("opacity", 0.3);

  // Then highlight only those that are an ancestor of the current segment.
  vis.selectAll("path")
      .filter(function(node) {
                return (sequenceArray.indexOf(node) >= 0);
              })
      .style("opacity", 1);
}

// Restore everything to full opacity when moving off the visualization.
/*function mouseleave(d) {
$("#sunburst-bc").empty ();
  d3.select("#sunburst-explanation")
      .style("visibility", "hidden");
}*/
/*
function initializeBreadcrumbTrail() {
  // Add the svg area.
  var trail = d3.select("#sequence").append("svg:svg")
      .attr("width", width)
      .attr("height", 50)
      .attr("id", "trail");
  // Add the label at the end, for the percentage.
  trail.append("svg:text")
    .attr("id", "endlabel")
    .style("fill", "#000");
}*/


// Update the breadcrumb trail to show the current sequence and percentage.
function updateBreadcrumbs(nodeArray, percentageString) {
console.log (nodeArray);

$("#sunburst-bc").empty ();
for (var i = 0; i < nodeArray.length; i++) {
  var arrow = "<span class='sbbcnoarrow' ></span>";
  var ml = "";
  if (i < nodeArray.length - 1)
    arrow = "<span class='sbbcarrow' style='border-left-color:"+colorMap(nodeArray[i])+";' ></span>";
  
  if (i > 0)
    ml = ";margin-left:-1em;"
  $("#sunburst-bc").append ("<div class='my-2 sbbc' style='background-color:"+colorMap(nodeArray[i])+";color:"+fgColorMap(nodeArray[i])+ml+"'>" + nodeArray[i].data.name + arrow + "</div>");
}



}


function buildHierarchy(csv) {
  var root = {"name": "root", "children": []};
  for (var i = 0; i < csv.length; i++) {
    var l = csv[i].length;
    var parts = csv[i].slice (0, l - 1);
    
    //console.log ("row ", i)
    //console.log (csv[i])
    //console.log (csv[i].slice (0, l - 1))
    //console.log (csv[i][l - 1])
    //break;
    //var sequence = csv[i][0];
    var size = +csv[i][l - 1];
    if (isNaN(size)) { // e.g. if this is a header row
      continue;
    }
    //var parts = sequence.split("-");
    var currentNode = root;
    for (var j = 0; j < parts.length; j++) {
      var children = currentNode["children"];
      var nodeName = parts[j];
      var childNode;
      if (j + 1 < parts.length) {
   // Not yet at the end of the sequence; move down the tree.
 	var foundChild = false;
 	for (var k = 0; k < children.length; k++) {
 	  if (children[k]["name"] == nodeName) {
 	    childNode = children[k];
 	    foundChild = true;
 	    break;
 	  }
 	}
  // If we don't already have a child node for this branch, create it.
 	if (!foundChild) {
 	  childNode = {"name": nodeName, "children": []};
 	  children.push(childNode);
 	}
 	currentNode = childNode;
      } else {
 	// Reached the end of the sequence; create a leaf node.
 	childNode = {"name": nodeName, "size": size};
 	children.push(childNode);
      }
    }
  }
  return root;
};
