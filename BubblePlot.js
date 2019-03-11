var OrganBubblePlot = (function(){
	"use strict";
	var instance;
	function Graph(){
		if (instance){
			return instance;
		}
		instance = this
	}
	Graph.getInstance = function(){
		return instance || new Graph();
	}
	return Graph;
}());