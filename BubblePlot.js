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
		this.height = .6*this.width;
		this.xMargin = .04*this.width;
		this.yMargin = .03*this.width;
		
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
		
		this.getAxes();
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
		for(var i = 0; i < this.patient.similarity_ssim.length; i++){
			let match = this.patient.similarity_ssim[i];
			this.patientMatches.push( this.data[match - 1] );
		}
	}
	
	Graph.getAxes = function(){
		var binWidth = (this.width - 2*this.xMargin)/num_organs;
		var xScale = d3.scaleLinear()
			.domain([0, this.organList.length])
			.range([binWidth/2, this.width - this.xMargin - binWidth/2]);
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
	}
	return Graph;
}());