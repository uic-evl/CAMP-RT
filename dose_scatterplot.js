//function to move svg elements  to the front
//from https://github.com/wbkd/d3-extended
d3.selection.prototype.moveToFront = function() {  
  return this.each(function(){
	this.parentNode.appendChild(this);
  });
};

function DoseScatterPlot(data){
	this.data = data;
	this.ids = data.getInternalIdList();
	this.selectedColor = '#e41a1c';
}

DoseScatterPlot.prototype.draw = function(target, selectedPerson = null){
	div = document.getElementById(target);
	this.width = div.clientWidth;
	this.height = div.clientHeight;
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
	this.getXAxis = function(d){ return self.data.getDosePCA(d, 1); };
	this.getYAxis = function(d){ return self.data.getDosePCA(d, 2); };
	this.drawCircles();
	if (selectedPatient != null){
		this.highlightSelectedPatients(selectedPerson);
	}
	this.drawClusterCircles(this.clusterMargin);
	this.setupTooltip(selectedPerson);
	this.setupSwitchButtons();
}

DoseScatterPlot.prototype.setupSwitchButtons = function(){
	var setSelected = function(element){
		element.style.opacity = 1;
		element.style.background = 'hsl(240, 90%, 30%)';
	}
	var setUnselected = function(element){
		element.style.opacity = .4;
		element.style.background = 'hsl(240, 50%, 30%)';
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
	var xDomain = d3.extent(this.ids, function(d){return self.getXAxis(d);});
	var yDomain = d3.extent(this.ids, function(d){return self.getYAxis(d);})
	
	var xScale = d3.scaleLinear()
		.domain(xDomain)
		.range([this.xMargin, this.width - this.xMargin]);
	var yScale = d3.scaleLinear()
		.domain(yDomain)
		.range([this.height - this.yMargin, this.yMargin]);
	return([xScale, yScale])
}

DoseScatterPlot.prototype.getFeatureScales = function(){
	var self = this;
	var sizeScale = d3.scalePow().exponent(2)
		.domain( d3.extent(this.ids, function(d){ return data.getPatientMeanError(d); }) )
		.range([3, 8]);
	return sizeScale;
}

DoseScatterPlot.prototype.getColorMap = function(){
	//get a cluster value based color map.  assumes clusters index 1-n_clusters
	//fir
	var data = this.data;
	var self = this;
	//make this better somehow
	var classColors = ['orange', 'cyan', 'magenta', 'blue', 'deeppink', 'purple', 'green', 'goldenrod', 'steelblue', 'brown', 'silver', 'burlywood', 'greenyellow', 'darkslategray']
	//scale luminosity based on best score?
	var luminosityScale = d3.scaleLinear()
		.domain(d3.extent(this.ids, function(d){ return self.data.getPatientSimilarityScores(d)[1]; }))
		.range([.4, .6]);
	var colorMap = function(dataPoint, id_input = true){
		var cluster = self.data.getCluster(dataPoint);
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
		.data(this.ids).enter()
		.append('circle')
		.attr('class', 'dot')
		.attr('id', function(d){ return 'scatterDot' + d;})
		.attr('r', function(d){
			return sizeScale(self.data.getPatientMeanError(d))})
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
			switchPatient(d);//from camprt.js
		});
}

DoseScatterPlot.prototype.drawClusterCircles = function(margin){
	var self = this;
	var clusters = new Map();
	var clusterColors = new Map();
	var [xScale, yScale] = this.getAxisScales();
	var toPoint = function(d){
		var x = xScale(self.getXAxis(d));
		var y = yScale(self.getYAxis(d));
		return [x,y];
	}
	this.ids.forEach(function(d){
		var cluster = self.data.getCluster(d);
		if(!clusters.has(cluster)){
			clusters.set(cluster, []);
			var color = self.getColor(d);
			clusterColors.set(cluster, color);
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
		var hull = [];
		var convexHull = d3.polygonHull(value);
		var centroid = d3.polygonCentroid(convexHull);
		convexHull.forEach(function(point){
			var offsetPoint = interpolateLine(centroid, point);
			hull.push(offsetPoint);
		});
		hull.color = clusterColors.get(key);
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
		.attr('stroke', function(d) {return d.color;})
		.attr('stroke-width', margin/3)
		.attr('opacity',.3)
		.merge(arc).transition().duration(800)
		.attr('d', function(d){return arcPath(d);})
	this.setupCurveTooltip();
}

DoseScatterPlot.prototype.setupCurveTooltip = function(){
	var clusterStats = new Map();
	var self = this;
	this.ids.forEach(function(d){
		var cluster = self.data.getCluster(d);
		if(!clusterStats.has(cluster)){
			var base = new Object();
			base.numPoints = 0;
			base.meanDose = 0;
			base.meanError = 0;
			clusterStats.set(cluster, base);
		}
		var current = clusterStats.get(cluster);
		current.numPoints += 1;
		current.meanError += self.data.getPatientMeanError(d);
		current.meanDose += self.data.getPatientMeanDose(d);
		clusterStats.set(cluster, current);
	});
	for(var stats of clusterStats.values()){
		stats.meanError = stats.meanError/stats.numPoints;
		stats.meanDose = stats.meanDose/stats.numPoints;
	}
	d3.selectAll('path').filter('.clusterCurves')
		.on('mouseover', function(d){
			d3.select(this).attr('opacity', 1);
			var stats = clusterStats.get(d.cluster);
			self.tooltip.html('Cluster ' + d.cluster + '</br>'
			+ 'Size: ' + stats.numPoints + '</br>'
			+ 'Mean Dose: ' + stats.meanDose.toFixed(1) + 'Gy </br>'
			+ 'Mean Prediction Error: ' + stats.meanError.toFixed(1) + 'Gy' )
				.style('left', d3.event.pageX + 8 + 'px')
				.style('top', d3.event.pageY - 20 + 'px');
			self.tooltip.transition().duration(50).style('visibility','visible');
		}).on('mouseout', function(d){
			d3.select(this).attr('opacity', .3);
			self.tooltip.transition().duration(50).style('visibility', 'hidden');
		});
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
		this.getXAxis = function(d){ return self.data.getDistancePCA(d, 1); };
		this.getYAxis = function(d){ return self.data.getDistancePCA(d, 2); };
	} else if(type == 'staging'){
		this.getXAxis = function(d){ 
			var volume = self.data.gtvpVol(d);
			return (volume > Math.E)? Math.log(volume): volume; };
		this.getYAxis = function(d){ 
			var volume = self.data.gtvnVol(d);
			return (volume > Math.E)? Math.log(volume): volume; };
	}
	else{
		this.getXAxis = function(d){ return self.data.getDosePCA(d, 1); };
		this.getYAxis = function(d){ return self.data.getDosePCA(d, 2); };
	}
	this.animateAxisChange();
}

DoseScatterPlot.prototype.highlightSelectedPatients = function(selectedPerson){
	//recolor everything first
	var getColor = this.getColor;
	var self = this;
	this.circles
		.attr('fill', function(d){ 
			return getColor(d);})
		.attr('stroke', 'black')
		.attr('stroke-width', 1)
		.attr('opacity', .5);
	//get entries linked to the person
	var selectedMatches = new Array();
	self.data.getPatientMatches(selectedPerson).forEach(function(id){
		selectedMatches.push(id);
	}, this);
	//recolor people matched with selected patient, make them darker and colorfuler
	selectedMatches.forEach(function(x){
		d3.select('#scatterDot' + x)
			.attr('opacity', 1)
			.attr('stroke-width', 2)
			.moveToFront();
	});
	//make main patient red and stuff
	d3.select('#scatterDot' + selectedPerson)
		.attr('opacity', 1)
		.attr('stroke-width', 2)
		.attr('fill', 'FireBrick')
		.moveToFront();
}

DoseScatterPlot.prototype.setupTooltip = function(selectedPatient){
	var tooltip = this.tooltip;
	var self = this;

	this.circles.on('mouseover', function(id){
		tooltip.html(self.data.getPatientName(id) + '</br>' 
			+ 'Dose: ' + self.data.getPatientMeanDose(id).toFixed(3) + ' Gy</br>'
			+ 'Error: ' + self.data.getPatientMeanError(id).toFixed(3) + ' Gy</br>'
			+ 'Cluster: ' + self.data.getCluster(id) + '</br>'
			+ 'x Value: ' + self.getXAxis(id).toFixed(3) + '</br>'
			+ 'y Value: ' + self.getYAxis(id).toFixed(3))
			.style('left', d3.event.pageX + 10 + 'px')
			.style('top', d3.event.pageY - 30 + 'px');
		tooltip.transition().duration(50).style('visibility','visible');
		Controller.brush(id);
	}).on('mouseout', function(d){
		tooltip.transition().duration(50).style('visibility', 'hidden');
		Controller.unbrush(d);
	});
}