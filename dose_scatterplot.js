//function to move svg elements  to the front
//from https://github.com/wbkd/d3-extended
d3.selection.prototype.moveToFront = function() {  
  return this.each(function(){
	this.parentNode.appendChild(this);
  });
};

function DoseScatterPlot(data){
	this.data = data;
	this.selectedColor = '#e41a1c';
}

DoseScatterPlot.prototype.draw = function(target, selectedPerson = null){
	div = document.getElementById(target);
	this.width = div.clientWidth;
	this.height = .5*this.width;
	this.xMargin = .04*this.width;
	this.yMargin = .03*this.width;
	this.clusterMargin = 20;
	d3.select("#"+target).selectAll('.scatterSvg').remove();
	this.svg = d3.select("#"+target).insert('svg',':first-child')
					.attr('class', 'scatterSvg')
					.attr('width', this.width)
					.attr('height', this.height)
					.append('g');
	this.tooltip = d3.select("div.tooltip")
			.attr('class','tooltip')
			.style('visibility','hidden');
	this.getColor = this.getColorMap();
	this.getXAxis = function(d){ return d.dose_pca[0]; };
	this.getYAxis = function(d){ return d.dose_pca[1]; };
	this.drawCircles();
	if (selectedPatient != null){
		this.highlightSelectedPatients(selectedPerson);
	}
	this.drawClusterCircles(this.clusterMargin);
	this.setupTooltip();
	this.setupSwitchButtons();
}

DoseScatterPlot.prototype.setupSwitchButtons = function(){
	var setSelected = function(element){
		element.style.opacity = 1;
		element.style.background = 'hsl(240, 90%, 50%)';
	}
	var setUnselected = function(element){
		element.style.opacity = .4;
		element.style.background = 'hsl(240, 50%, 60%)';
	}
	var doseButton = document.getElementById('doseScatterButton');
	var distanceButton = document.getElementById('distanceScatterButton');
	var stagingButton = document.getElementById('stagingScatterButton');
	setSelected(doseButton);
	setUnselected(distanceButton);
	setUnselected(stagingButton);
	var self = this;
	doseButton.addEventListener('click', function(){
		self.switchAxisVariable('dose');
		setSelected(doseButton);
		setUnselected(distanceButton);
		setUnselected(stagingButton);
	});
	distanceButton.addEventListener('click', function(){
		self.switchAxisVariable('distance');
		setSelected(distanceButton);
		setUnselected(doseButton);
		setUnselected(stagingButton);
	});
	stagingButton.addEventListener('click', function(){
		self.switchAxisVariable('staging');
		setSelected(stagingButton);
		setUnselected(doseButton);
		setUnselected(distanceButton);
	});
}

DoseScatterPlot.prototype.getAxisScales = function(){
	self = this;
	var xDomain = d3.extent(this.data, function(d){return self.getXAxis(d);});
	var yDomain = d3.extent(this.data, function(d){return self.getYAxis(d);})
	
	var xScale = d3.scaleLinear()
		.domain(xDomain)
		.range([this.xMargin, this.width - this.xMargin]);
	var yScale = d3.scaleLinear()
		.domain(yDomain)
		.range([this.height - this.yMargin, this.yMargin]);
	return([xScale, yScale])
}

DoseScatterPlot.prototype.getFeatureScales = function(){
	var sizeScale = d3.scalePow().exponent(2)
		.domain( d3.extent(this.data, function(d){ return d.mean_error; }) )
		.range([3, 8]);
	return sizeScale;
}

DoseScatterPlot.prototype.getColorMap = function(){
	//get a cluster value based color map.  assumes clusters index 1-n_clusters
	//fir
	var data = this.data;
	//make this better somehow
	var classColors = ['orange', 'cyan', 'magenta', 'blue', 'deeppink', 'purple', 'green', 'goldenrod', 'steelblue', 'brown', 'silver', 'burlywood', 'greenyellow', 'darkslategray']
	var errorPrinted = false;
	//scale luminosity based on best score?
	var luminosityScale = d3.scaleLinear()
		.domain(d3.extent(this.data, function(d){ return d.scores_ssim[1]; }))
		.range([.4, .6]);
	var colorMap = function(dataPoint){
		var cluster = dataPoint.cluster;
		if(cluster >= classColors.length){
			var color = 'black';
		}
		else{
			var color = classColors[cluster-1];
		}
		return color;
	}
	return colorMap;
}

DoseScatterPlot.prototype.drawCircles = function(){
	var self = this;
	var [xScale, yScale] = this.getAxisScales();
	var sizeScale = this.getFeatureScales();
	var getColor = this.getColor;
	this.circles = this.svg.selectAll('.dot')
		.data(this.data).enter()
		.append('circle')
		.attr('class', 'dot')
		.attr('id', function(d){ return 'scatterDot' + d.ID;})
		.attr('r', function(d){
			return sizeScale(d.mean_error)})
		.attr('cx', function(d) {
			return xScale(self.getXAxis(d));})
		.attr('cy', function(d) {
			return yScale(self.getYAxis(d));});
	this.circles
		.attr('fill', function(d){ 
			return getColor(d);})
		.attr('stroke', 'black')
		.attr('opacity', .8)
		.attr('stroke-width', 1)
		.on('click', function(d){
			switchPatient(d.ID_internal);//from camprt.js
		});
}

DoseScatterPlot.prototype.drawClusterCircles = function(margin){
	var self = this;
	var clusters = new Map();
	var [xScale, yScale] = this.getAxisScales();
	var toPoint = function(d){
		var x = xScale(self.getXAxis(d));
		var y = yScale(self.getYAxis(d));
		return [x,y];
	}
	this.data.forEach(function(d){
		var cluster = d.cluster;
		if(!clusters.has(cluster)){
			clusters.set(cluster, []);
		}
		var current = clusters.get(cluster);
		current.push(toPoint(d))
		clusters.set(cluster, current);
	}, clusters);
	var interpolateLine = function(x0, x1){
		var magnitude = Math.sqrt((x1[1] - x0[1])**2 + (x1[0] - x0[0])**2);
		var vect = [ (x1[0] - x0[0])/magnitude, (x1[1] - x0[1])/magnitude ];
		var point = [x1[0] + vect[0]*margin, x1[1] + vect[1]*margin];
		return point;
	}
	var offsetHulls = [];
	for (var [key, value] of clusters.entries()) {
		var hull = []
		var convexHull = d3.polygonHull(value);
		var centroid = d3.polygonCentroid(convexHull);
		convexHull.forEach(function(point){
			var offsetPoint = interpolateLine(centroid, point);
			hull.push(offsetPoint);
		});
		hull.cluster = +key;
		offsetHulls.push(hull);
	}
	var arcPath = d3.line()
		.x(function(d){ return d[0];})
		.y(function(d){ return d[1];})
		.curve(d3.curveCardinalClosed);
	var arc = this.svg.selectAll('.clusterCurves')
		.data(offsetHulls);
	arc.exit().remove();
	arc.enter()
		.append('path')
		.attr('class','clusterCurves')
		.attr('fill', 'none')
		.attr('stroke', function(d) {return self.getColor(d);})
		.attr('stroke-width', margin/10)
		.attr('opacity', .1)
		.merge(arc).transition().duration(800)
		.attr('d', function(d){return arcPath(d);})
		.attr('opacity', .9);
	console.log(offsetHulls);

}

DoseScatterPlot.prototype.setAxisVariable = function(axisFunction, axis){
	axis = +axis
	if(axis != 1 || axis != 0){
		console.log('invalid axis to scatterplot setAxisVariable.  Value mut be 1 or 0');
	}
	if(axis == 1){
		this.getXAxis = axisFunction;
	}else{
		this.getYAxis = axisFunction;
	}
}

DoseScatterPlot.prototype.animateAxisChange = function(){
	var [xScale, yScale] = this.getAxisScales();
	this.drawClusterCircles(this.clusterMargin);
	this.circles.transition().duration(800)
		.attr('cx', function(d) {
			return xScale(self.getXAxis(d));})
		.attr('cy', function(d) {
			return yScale(self.getYAxis(d));}); 
}

DoseScatterPlot.prototype.switchAxisVariable = function(type){
	if(type == 'distance'){
		this.getXAxis = function(d){ return d.distance_pca[0]; };
		this.getYAxis = function(d){ return d.distance_pca[1]; };
	} else if(type == 'staging'){
		this.getXAxis = function(d){ return (d.gtvp_volume > Math.E)? Math.log(d.gtvp_volume): d.gtvp_volume; };
		this.getYAxis = function(d){ return (d.gtvn_volume > Math.E)? Math.log(d.gtvn_volume): d.gtvn_volume; };
	}
	else{
		this.getXAxis = function(d){ return d.dose_pca[0]; };
		this.getYAxis = function(d){ return d.dose_pca[1]; };
	}
	this.animateAxisChange();
}

DoseScatterPlot.prototype.highlightSelectedPatients = function(selectedPerson){
	//recolor everything first
	var getColor = this.getColor;
	this.circles
		.attr('fill', function(d){ 
			return getColor(d);})
		.attr('stroke', 'black')
		.attr('stroke-width', 1)
		.attr('opacity', .5);
	//get entries linked to the person
	var selectedData = this.data[selectedPerson - 1];
	this.subsetData = new Array();
	selectedData.similarity_ssim.forEach(function(x){
		this.subsetData.push(this.data[x - 1]);
	}, this);
	//recolor people matched with selected patient, make them darker and colorfuler
	this.subsetData.forEach(function(x){
		var dataPoint = d3.select('#scatterDot' + x.ID)
			.attr('opacity', 1)
			.attr('stroke-width', 2)
			.moveToFront();
	});
	//make main patient red and stuff
	d3.select('#scatterDot' + selectedData.ID)
		.attr('opacity', 1)
		.attr('stroke-width', 2)
		.attr('fill', 'FireBrick')
		.moveToFront();
}

DoseScatterPlot.prototype.setupTooltip = function(d){
	var tooltip = this.tooltip;
	this.circles.on('mouseover', function(d){
		tooltip.html(d.name + '</br>' 
			+ 'Error: ' + Math.round(d.mean_error*100)/100 + '</br>'
			+ 'Cluster: ' + d.cluster)
			.style('left', d3.event.pageX + 10 + 'px')
			.style('top', d3.event.pageY - 30 + 'px');
		tooltip.transition().duration(50).style('visibility','visible');
	}).on('mouseout', function(d){
		tooltip.transition().duration(50).style('visibility', 'hidden');
	});
}