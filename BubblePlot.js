var OrganBubblePlot = (function(){
	"use strict";
	var instance;
	var num_organs;
	var self = this;
	function Graph(){
		if (instance){
			return instance;
		}
		instance = this
	}
	Graph.getInstance = function(){
		return instance || new Graph();
	}
	Graph.init = function(target, patientInternalId, patientData){
		console.log('drawing');
		div = document.getElementById(target);
		this.setPatient(patientInternalId, patientData);
		
		this.width = div.clientWidth;
		this.height = div.clientHeight
		this.xMargin = .025*this.width;
		this.yMargin = .05*this.height;
		this.binWidth = (this.width - 2*this.xMargin)/num_organs;
		
		d3.select("#"+target).selectAll('.bubbleSvg').remove();
		this.svg = d3.select("#"+target).append('svg')
						.attr('class', 'bubbleSvg')
						.attr('width', this.width)
						.attr('height', this.height)
						.append('g');
		this.tooltip = d3.select("div.bubbletooltip")
				.attr('class','tooltip')
				.style('visibility','hidden');
		console.log(this.organList);
		
		this.drawOrgans();
	}
	
	Graph.setPatient = function(patientInternalId, patientData){
		this.data = patientData;
		this.patient = patientData[patientInternalId - 1];
		
		var organList = Object.keys(this.patient.organData)
		num_organs = organList.length;
		if('GTVp' in this.patient.organData){ num_organs = num_organs - 1; }
		if('GTVn' in this.patient.organData){ num_organs = num_organs - 1; }
		this.organList = organList.slice(0, num_organs);
		
		this.patientMatches = [];
		for(var i = 1; i < this.patient.similarity_ssim.length; i++){
			let match = this.patient.similarity_ssim[i];
			this.patientMatches.push( this.data[match - 1] );
		}
	}
	
	Graph.drawOrgans = function(){
		var [xAxis, yAxis] = this.getAxes();
		var rAxis = this.getRScale(this.patient, this.binWidth);
		this.organList.forEach(function(d){
			this.drawOrganBubbles(d, xAxis, yAxis, rAxis);
		},this);
		this.drawDoseRects(xAxis, yAxis);
	}
	
	Graph.drawDoseRects = function(xScale, yScale){
		var patient = this.patient;
		var width = .8*this.binWidth;
		var height = .2*this.binWidth;
		var self = this;
		var drawDose = function(variable, color){
			self.svg.selectAll('.'+ variable + 'Rect')
				.data(self.organList).enter()
				.append('rect')
				.attr('class', 'meanDoseRect')
				.attr('x', function(o,index){ return xScale(index) - width/2;})
				.attr('y', function(o) { return yScale(patient.organData[o][variable]) - height/2; })
				.attr('fill', color)
				.attr('width', width)
				.attr('height', height)
				.attr('opacity', .8)
				.attr('stroke-width', .5)
				.attr('stroke', 'black');
		}
		drawDose('meanDose', 'red');
		drawDose('estimatedDose', 'blue');
	}
	
	Graph.drawOrganBubbles = function(organName, xScale, yScale, rScale){
		var position = this.organList.indexOf(organName);
		var x = xScale(position);
		var patient = this.patient;
		var getRadius = function(x){
			var index = patient.similarity_ssim.indexOf(+x.ID_internal);
			if(index == 1){ console.log("bad stuff");}
			if(index != -1){
				return rScale(patient.scores_ssim[index]);
			} else{
				return 0;
				console.log('error in getRadius?');
			}
		}
		this.svg.selectAll('.organDot').filter('.' + organName).remove();
		var bubbles = this.svg.selectAll('.organDot').filter('.' + organName)
			.data(this.patientMatches).enter()
			.append('circle')
			.attr('class', 'organDot ' + organName)
			.attr('cx', x)
			.attr('cy', function(d){ return yScale(d.organData[organName].meanDose);})
			.attr('r', getRadius)
			.attr('opacity', .4)
			.attr('fill', 'gold')
			.attr('stroke', 'black')
			.attr('stroke-width', .1);
	}
	
	Graph.getRScale = function(patient, width){
		var ratio = .4;
		var scoreExtent = d3.extent( patient.scores_ssim.slice(1), function(d) {return d;});
		console.log([ratio*width, width*(1-ratio)]);
		var rScale = d3.scalePow().exponent(2)
			.domain(scoreExtent)
			.range([ratio*width/2, width*(1-ratio)/2]);
		return(rScale)
	}
	
	Graph.getAxes = function(){
		var xScale = d3.scaleLinear()
			.domain([0, this.organList.length])
			.range([this.binWidth/2 + this.xMargin, this.width - this.xMargin - this.binWidth/2]);
		//convoluted way of getting ths max and min dose for everything
		var allPatients = [...this.patientMatches];
		allPatients.push(this.patient);
		let maxDose = 0;
		let minDose = 0;
		allPatients.forEach(function(d){
			this.organList.forEach(function(x){
				var newDose = d.organData[x].meanDose;
				maxDose = (maxDose > newDose)? maxDose: newDose;
				minDose = (minDose < newDose)? minDose: newDose;
			});
		}, this);
		var yScale = d3.scaleLinear()
			.domain( [minDose, maxDose])
			.range([this.height - this.yMargin, this.yMargin]);
		return [xScale, yScale];
	}
	return Graph;
}());