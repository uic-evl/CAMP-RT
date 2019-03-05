'use strict';

if (!Detector.webgl) {
    Detector.addGetWebGLMessage();
}

var parent = document.getElementById("content");
var nodeDetails = document.getElementById("details");

var detailsOffsetX = 15;
var detailsOffsetY = 15;

var organName = document.getElementById("details_organName"),
    dosePerVolume = document.getElementById("details_dosePerVolume"),
    lineSeparator = document.getElementById("details_line"),
    volumeVal = document.getElementById("details_volume_val"),
    meanDoseVal = document.getElementById("details_meanDose_val"),
    minDoseVal = document.getElementById("details_minDose_val"),
    maxDoseVal = document.getElementById("details_maxDose_val"),
    pDisplayed = document.getElementById("pDisplayed"),
    pTarget_chart = document.getElementById("pTarget_chart"),
    pPrediction_chart = document.getElementById("pPrediction_chart"),
    pDifference_chart = document.getElementById("pDifference_chart");


var scenes = [],
    scenesRP = [],
    renderer, renderer2;

	
var selectedPatient = 1;
//patients shown on load screen?
var patientsToShow = 5;
//patients to use in the macching algorithm
var numMatches = 5;

var totalModelCount;

var syncCameras = true,
    syncCamerasInterval,
    detailsOnRotate = true;

var raycaster;

var mouse = new THREE.Vector2(-500, -500);

var mouseNorm = new THREE.Vector2(-500, -500),
    INTERSECTED = null,
    nodeHover = null;

var cameraDistZ = 500;

// 36 steps

var domainColorScale = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105];
var rangeColorScale = ['#ffffe0', '#fff8d2', '#fff0c4', '#ffe9b8', '#ffe2ae', '#ffdaa3', '#ffd39a', '#ffcb91', '#ffc389', '#ffbb82', '#ffb27c', '#ffab77', '#ffa272', '#ff986e', '#fe906a', '#fb8768', '#f98065', '#f67762', '#f26f60', '#ee675d', '#eb5f5b', '#e75758', '#e25055', '#dd4852', '#d8404e', '#d3394a', '#cc3146', '#c62a41', '#c0223b', '#b91c35', '#b3152f', '#ab0e28', '#a40820', '#9b0317', '#93010e', '#8b0000'];

var color = d3.scaleLinear()
    .domain(domainColorScale)
    .range(rangeColorScale);

var domainColorScale2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
var rangeColorScale2 = ['#999999', '#98949f', '#968fa5', '#958aaa', '#9384b0', '#9180b5', '#8e7aba', '#8c75bf', '#8971c5', '#856bca', '#8166d0', '#7d61d5', '#795bda', '#7356e0', '#6e50e4', '#674bea', '#5f44ef', '#553ff5', '#4b38fa', '#3d32ff'];

var color2 = d3.scaleLinear()
    .domain(domainColorScale2)
    .range(rangeColorScale2);

// data
var organs, oAtlas, links, patients;
var organModels = new THREE.Group();

var partitions = ["Oral Cavity & Jaw", "Throat", "Salivary Glands", "Eyes", "Brainstem & Spinal Cord", "Other"];

var organRef = [];

var pRankingOrder, pScores;

var listItems, arrayOfDivs = [],
    currScene;

var files = ["data/organAtlas.json", "PYTHON/data/patient_dataset.json"];
var promises = [];

var master = document.getElementById("masterList");

var materialArray;
var materialArray2;

var canvas = document.getElementById("c");
var canvas2 = document.getElementById("c2");
var template = document.getElementById("template").text;

var data = {
    y: "Dose (GY)",
    series: [],
    dates: []
};

var manager = new THREE.LoadingManager();

manager.onLoad = function () {
	//this may break this because I moved it to the front
    console.log('Loading complete!');
    //updateOrder(selectedPatient);
    document.getElementById("loadScreen").style.display = "none";
};

var loadProgress = document.getElementById("loadProgress");

manager.onProgress = function (url, itemsLoaded, itemsTotal) {
    loadProgress.innerHTML = parseInt(itemsLoaded / itemsTotal * 100) + " %"
};

files.forEach(function (url) {
    promises.push(d3.json(url));
});

Promise.all(promises).then(function (values) {
    start(values[0], values[1]);
});

function start(organAtlas, patientsData) {

    console.log(patientsData);

    oAtlas = organAtlas[0];
    patients = patientsData;

    selectedPatient = populateDropDownMenu();

    pRankingOrder = patients[selectedPatient - 1].similarity_ssim;

    pScores = patients[selectedPatient - 1].scores_ssim;

    populateColorScale();

    flipGraph(); // fixes orientation of organs
    computeCenterOfGraphAndShift(); // compute center of graph and shift to origin

    init(); // initialize

    populateOrganMasterList();

    listItems = document.getElementsByClassName("list-item");

    for (var i = 0, ref = arrayOfDivs.length = listItems.length; i < ref; i++) {
        arrayOfDivs[i] = listItems[i];
    }

    currScene = scenes[0];

    document.addEventListener("mousedown", onMouseDown, false);
    document.addEventListener("mouseup", onMouseUp, false);

    document.addEventListener("touchstart", onTouchStart, false);
    document.addEventListener("touchend", onTouchEnd, false);

    document.addEventListener("mousemove", onDocumentMouseMove, false);

    animate(); // render
}

// ----------------------------------------------------------------

function populateColorScale() {

    var parentDiv = document.getElementById("colorScale");
	//appears to make a color scale made from just like, a bunch of 10px wide divs for the hard-coded color scales
    rangeColorScale.forEach(function (color, index) {

        var tempDiv = document.createElement("div");

        tempDiv.style.height = "100%";
        tempDiv.style.width = "10px";
        tempDiv.style["backgroundColor"] = color;
        tempDiv.style.display = "inline-block";
        tempDiv.style.transition = "0.3s";

        tempDiv.value = domainColorScale[index];
        tempDiv.className = "colorScaleEntry";

        parentDiv.appendChild(tempDiv);
    });

    var colorScaleEntries = document.getElementsByClassName("colorScaleEntry");
	//adds mouseover tooltip
    for (var i = 0, ref = colorScaleEntries.length; i < ref; i++) {
        colorScaleEntries[i].addEventListener('mouseover', showColorScaleLabel, false);
        colorScaleEntries[i].addEventListener('mouseout', hideColorScaleLabel, false);
    }
	//does it again for the blue one
    parentDiv = document.getElementById("colorScale2");

    rangeColorScale2.forEach(function (color, index) {

        var tempDiv = document.createElement("div");

        tempDiv.style.height = "100%";
        tempDiv.style.width = "17px";
        tempDiv.style["backgroundColor"] = color;
        tempDiv.style.display = "inline-block";
        tempDiv.style.transition = "0.3s";

        tempDiv.value = domainColorScale2[index];
        tempDiv.className = "colorScaleEntry2";

        parentDiv.appendChild(tempDiv);
    });

    var colorScaleEntries2 = document.getElementsByClassName("colorScaleEntry2");

    for (var i = 0, ref = colorScaleEntries2.length; i < ref; i++) {
        colorScaleEntries2[i].addEventListener('mouseover', showColorScaleLabel2, false);
        colorScaleEntries2[i].addEventListener('mouseout', hideColorScaleLabel2, false);
    }
}

function showColorScaleLabel(event) {
	
    var details = document.getElementById("colorScaleDetails");

    details.style.left = event.target.offsetLeft + "px";
    details.innerHTML = "" + event.target.value;

    details.style.display = "block";
}

function hideColorScaleLabel(event) {
	
    var details = document.getElementById("colorScaleDetails");

    details.style.display = "none";
}

function showColorScaleLabel2(event) {
	
    var details = document.getElementById("colorScaleDetails2");

    details.style.left = event.target.offsetLeft + "px";
    details.innerHTML = "" + event.target.value;

    details.style.display = "block";
}

function hideColorScaleLabel2(event) {

    var details = document.getElementById("colorScaleDetails2");

    details.style.display = "none";
}

function compareID(a, b) {

    var a_ID = a.ID_int;
    var b_ID = b.ID_int;

    var comparison = 0;

    if (a_ID > b_ID) {
        comparison = 1;
    } else if (a_ID < b_ID) {
        comparison = -1;
    }
    return comparison;
}

function populateDropDownMenu() {
	//holds an array of patient internal ids
    var menu = document.getElementById("patientMenu");

    // copy of patients sorted
    var patients_sorted = patients.concat().sort(compareID);

    patients_sorted.forEach(function (patient, index) {

        var tempOption = document.createElement("option");

        tempOption.value = patient.ID_internal;

        tempOption.innerHTML = patient.name;

        menu.appendChild(tempOption);
    });

    // first patient 
    var firstPatient = patients_sorted[0].ID_internal;

    //THis appears to look at the url to see if a patient is selected there
	//if so, sets "first patient" to this guy? otherwise, uses the lowest id (above)
    var patientURL = getQueryVariable("id");

    if (patientURL != false) {

        var convertedPatient = getInternalID(patientURL);

        if (convertedPatient != false) {
            firstPatient = convertedPatient;
            menu.value = convertedPatient;
        }
    }

    patients_sorted.length = 0;

    return firstPatient;
}

function getInternalID(searchString) {

    for (var i = 0; i < patients.length; i++) {

        if (patients[i].ID == searchString)
            return patients[i].ID_internal;

    }

    return false;
}

function getQueryVariable(variable) {
    var query = window.location.search.substring(1);
    var vars = query.split("&");
    for (var i = 0; i < vars.length; i++) {
        var pair = vars[i].split("=");
        if (pair[0] == variable) {
            return pair[1];
        }
    }
    return (false);
}


function flipGraph() {
	//coordinate rotation and scaling for organ positions
    for (var i = 0; i < patients.length; i++) {

        var patientOrganList = patients[i].organData;

        for (var pOrgan in patientOrganList) {

            var tOrganX = (patientOrganList[pOrgan].x * -1);
            var tOrganY = (patientOrganList[pOrgan].y * -1);
            var tOrganZ = (patientOrganList[pOrgan].z * -1);

            patientOrganList[pOrgan].x = tOrganY * 1.3;
            patientOrganList[pOrgan].y = tOrganZ * 2.5;
            patientOrganList[pOrgan].z = tOrganX * 1.1;
        }
    }

    for (var pOrgan in oAtlas) {

        var tOrganX = (oAtlas[pOrgan].x * -1);
        var tOrganY = (oAtlas[pOrgan].y * -1);
        var tOrganZ = (oAtlas[pOrgan].z * -1);

        oAtlas[pOrgan].x = tOrganY * 1.3;
        oAtlas[pOrgan].y = tOrganZ * 2.5;
        oAtlas[pOrgan].z = tOrganX * 1.1;
    }
}

function computeCenterOfGraphAndShift() {
	//coordiante translation so that the organs are centered, I think
    for (var i = 0; i < patients.length; i++) {

        var sceneCenter = [0.0, 0.0, 0.0];

        var xyzMin = new Array(3);
        var xyzMax = new Array(3);

        var positions = [];

        var patientOrganList = patients[i].organData;

        for (var pOrgan in patientOrganList) {

            var xyz = {
                x: patientOrganList[pOrgan].x,
                y: patientOrganList[pOrgan].y,
                z: patientOrganList[pOrgan].z
            };

            positions.push(xyz);
        }

        xyzMin = getMin(positions);
        xyzMax = getMax(positions);

        sceneCenter = [
            ((xyzMin[0] + xyzMax[0]) / 2),
            ((xyzMin[1] + xyzMax[1]) / 2),
            ((xyzMin[2] + xyzMax[2]) / 2)
        ];

        for (var pOrgan in patientOrganList) {

            patientOrganList[pOrgan].x = (patientOrganList[pOrgan].x - sceneCenter[0]);
            patientOrganList[pOrgan].y = (patientOrganList[pOrgan].y - sceneCenter[1]);
            patientOrganList[pOrgan].z = (patientOrganList[pOrgan].z - sceneCenter[2]);
        }

    }

    var sceneCenter = [0.0, 0.0, 0.0];

    var xyzMin = new Array(3);
    var xyzMax = new Array(3);

    var positions = [];
	//shifts organs in organ atlas also?  couldn't this just be hard coded?
    for (var pOrgan in oAtlas) {

        var xyz = {
            x: oAtlas[pOrgan].x,
            y: oAtlas[pOrgan].y,
            z: oAtlas[pOrgan].z
        };

        positions.push(xyz);
    }

    xyzMin = getMin(positions);
    xyzMax = getMax(positions);

    sceneCenter = [
        ((xyzMin[0] + xyzMax[0]) / 2),
        ((xyzMin[1] + xyzMax[1]) / 2),
        ((xyzMin[2] + xyzMax[2]) / 2)
        ];

    for (var pOrgan in oAtlas) {

        oAtlas[pOrgan].x = (oAtlas[pOrgan].x - sceneCenter[0]);
        oAtlas[pOrgan].y = (oAtlas[pOrgan].y - sceneCenter[1]);
        oAtlas[pOrgan].z = (oAtlas[pOrgan].z - sceneCenter[2]);
    }
}

function getMin(pos) {

    var x = pos.reduce(function (min, obj) {
        return obj.x < min ? obj.x : min;
    }, Infinity);

    var y = pos.reduce(function (min, obj) {
        return obj.y < min ? obj.y : min;
    }, Infinity);

    var z = pos.reduce(function (min, obj) {
        return obj.z < min ? obj.z : min;
    }, Infinity);

    if (x == Infinity || y == Infinity || z == Infinity)
        return [0.0, 0.0, 0.0];

    return [x, y, z];
}

function getMax(pos) {

    var x = pos.reduce(function (max, obj) {
        return obj.x > max ? obj.x : max;
    }, -Infinity);

    var y = pos.reduce(function (max, obj) {
        return obj.y > max ? obj.y : max;
    }, -Infinity);

    var z = pos.reduce(function (max, obj) {
        return obj.z > max ? obj.z : max;
    }, -Infinity);

    if (x == -Infinity || y == -Infinity || z == -Infinity)
        return [0.0, 0.0, 0.0];

    return [x, y, z];
}

function handleCheckBoxSingle(event) {

    if (event.checked) {

        scenes.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);
            var model = scene.getObjectByName(String(event.value) + "_model");

            //console.log(event.value);

            if (node && model) {
                node.visible = true;
                model.visible = true;
            }
        });

        scenesRP.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);
            var model = scene.getObjectByName(String(event.value) + "_model");

            //console.log(event.value);

            if (node && model) {
                node.visible = true;
                model.visible = true;
            }
        });

        //d3.select("#line_" + event.value).style("opacity", 1.0);
        d3.select("#line_" + event.value).attr("display", null);

    } else {

        scenes.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);
            var model = scene.getObjectByName(String(event.value) + "_model");

            if (node && model) {
                node.visible = false;
                model.visible = false;
            }
        });

        scenesRP.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);
            var model = scene.getObjectByName(String(event.value) + "_model");

            if (node && model) {
                node.visible = false;
                model.visible = false;
            }
        });

        //d3.select("#line_" + event.value).style("opacity", 0.0);
        d3.select("#line_" + event.value).attr("display", "none");
    }


}

function handleCheckBoxGroup(event) {

    var children = master.getElementsByClassName(event.id[0] + "_GroupChildren");

    if (event.checked) {

        for (var i = 0; i < children.length; i++) {

            children[i].checked = true;

            d3.select("#line_" + children[i].value).attr("display", null);
        }

        scenes.forEach(function (scene, index) {

            for (var i = 0; i < children.length; i++) {

                var node = scene.getObjectByName(children[i].value);
                var model = scene.getObjectByName(String(children[i].value) + "_model");

                if (node && model) {
                    node.visible = true;
                    model.visible = true;
                }

            }
        });

        scenesRP.forEach(function (scene, index) {

            for (var i = 0; i < children.length; i++) {

                var node = scene.getObjectByName(children[i].value);
                var model = scene.getObjectByName(String(children[i].value) + "_model");

                if (node && model) {
                    node.visible = true;
                    model.visible = true;
                }
            }
        });

    } else {

        for (var i = 0; i < children.length; i++) {

            children[i].checked = false;
            
            d3.select("#line_" + children[i].value).attr("display", "none");
        }

        scenes.forEach(function (scene, index) {

            for (var i = 0; i < children.length; i++) {

                var node = scene.getObjectByName(children[i].value);
                var model = scene.getObjectByName(String(children[i].value) + "_model");
				
                if (node && model) {
                    node.visible = false;
                    model.visible = false;
                }

            }
        });

        scenesRP.forEach(function (scene, index) {

            for (var i = 0; i < children.length; i++) {

                var node = scene.getObjectByName(children[i].value);
                var model = scene.getObjectByName(String(children[i].value) + "_model");

                if (node && model) {
                    node.visible = false;
                    model.visible = false;
                }
            }
        });
    }


}

function populateOrganMasterList() {

    // make partition input first
    partitions.forEach(function (group, i) {

        var tempDiv = document.createElement("div");

        tempDiv.setAttribute("class", "checkbox_group");
        tempDiv.setAttribute("id", String(i + 1) + "_group_container");

        var tempInput = document.createElement("INPUT");

        tempInput.setAttribute("type", "checkbox");
        tempInput.setAttribute("id", String(i + 1) + "_title");
        tempInput.setAttribute("value", group);
 
        tempInput.setAttribute("onchange", "handleCheckBoxGroup(this)");

        tempInput.setAttribute("checked", true);


        var tempLabel = document.createElement("label");
        tempLabel.setAttribute("for", String(i + 1) + "_title");

        tempLabel.style.fontSize = "16px";
        tempLabel.style.fontWeight = "bold";
        tempLabel.style.paddingLeft = "15px";
        tempLabel.innerHTML = group;

        // ----------

        tempDiv.appendChild(tempInput);
        tempDiv.appendChild(tempLabel);


        // ----------

        var tempDiv2 = document.createElement("div");

        tempDiv2.setAttribute("id", String(i + 1) + "_single_container");



        master.appendChild(tempDiv);
        master.appendChild(tempDiv2);

        var tempDiv3 = document.createElement("div");

        tempDiv3.setAttribute("class", "dummy");

        master.appendChild(tempDiv3);

    });

	////unhelpful documentation:
    // for loop bad, iterates in an unspecified order!!!!!!!!!!!!!!!!!!!!!!

    // individual organs
    for (var organ in oAtlas) {

        if (organ != "GTVn" && organ != "GTVp") {

            var parent = document.getElementById(String(oAtlas[organ].partition) + "_single_container");

            var tempDiv = document.createElement("div");
            tempDiv.style.marginRight = "20px";
            tempDiv.style.marginLeft = "50px";
            tempDiv.style.paddingLeft = "10px";
            tempDiv.style.paddingTop = "1px";
            tempDiv.style.paddingBottom = "1px";
            tempDiv.setAttribute("class", "GroupChildren");

            var tempInput = document.createElement("INPUT");

            tempInput.setAttribute("type", "checkbox");
            tempInput.setAttribute("id", organ + "_checkList");
            tempInput.setAttribute("class", String(oAtlas[organ].partition) + "_GroupChildren");
            tempInput.setAttribute("value", organ);
            tempInput.setAttribute("onchange", "handleCheckBoxSingle(this)");

            tempInput.setAttribute("checked", true);


            var tempLabel = document.createElement("label");
            tempLabel.setAttribute("for", organ + "_checkList");
            tempLabel.style.paddingLeft = "15px";

            tempLabel.innerHTML = organ;

            // ----------

            tempDiv.appendChild(tempInput);
            tempDiv.appendChild(tempLabel);

            parent.appendChild(tempDiv);

        }

    }
}

function checkOrganMasterList() {

    organs.forEach(function (organ, index) {

        var tempItem = document.getElementById(organ.name);

        if (scenes[selectedPatient - 1].getObjectByName(organ.name).userData.dosePerVolume != null) {

            if (tempItem.checked != true)
                tempItem.setAttribute("checked", true);

        } else {

            if (tempItem.checked != false)
                tempItem.setAttribute("checked", false);
        }
    });

}

function formatOrganMasterList() {
    var organList = master.children;

    for (var i = 0; i < organList.length; i++) {

        if (i % 2 == 0)
            organList[i].style["backgroundColor"] = "#3a3a3a";
        else
            organList[i].style["backgroundColor"] = "#444444";

    }
}

function init() {
	
	//renderer for main views?
	var getRenderer = function(canvas){
		var r = new THREE.WebGLRenderer({
			canvas: canvas,
			antialias: true
		});
		r.setClearColor(0xffffff, 1);
		r.setPixelRatio(window.devicePixelRatio);
		r.sortObjects = true;
		return r
	}
    
	renderer = getRenderer(canvas);

	//renderer for dose estimation views
    renderer2 = getRenderer(canvas2);

    raycaster = new THREE.Raycaster();

    var maxAnisotropy = renderer.getMaxAnisotropy();

    var textureLoader = new THREE.TextureLoader();

    var texture0 = textureLoader.load('resources/anterior.png'), // xpos, Right
        texture1 = textureLoader.load('resources/posterior.png'), // xneg, Left
        texture2 = textureLoader.load('resources/superior.png'), // ypos, Top
        texture3 = textureLoader.load('resources/inferior.png'), // yneg, Bottom
        texture4 = textureLoader.load('resources/right.png'), // zpos, Back
        texture5 = textureLoader.load('resources/left.png'); // zneg, Front

    texture0.anisotropy = maxAnisotropy;
    texture1.anisotropy = maxAnisotropy;
    texture2.anisotropy = maxAnisotropy;
    texture3.anisotropy = maxAnisotropy;
    texture4.anisotropy = maxAnisotropy;
    texture5.anisotropy = maxAnisotropy;
    //getcontext2d. draw image
    materialArray = [
            new THREE.MeshBasicMaterial({
            map: texture0
        }),
            new THREE.MeshBasicMaterial({
            map: texture1
        }),
            new THREE.MeshBasicMaterial({
            map: texture2
        }),
            new THREE.MeshBasicMaterial({
            map: texture3
        }),
            new THREE.MeshBasicMaterial({
            map: texture4
        }),
            new THREE.MeshBasicMaterial({
            map: texture5
        })
    ];
    //

    //
    var maxAnisotropy2 = renderer2.getMaxAnisotropy();

	//Do we even use this?
    var textureLoader2 = new THREE.TextureLoader();

    var texture0_2 = textureLoader2.load('resources/anterior.png'), // xpos, Right
        texture1_2 = textureLoader2.load('resources/posterior.png'), // xneg, Left
        texture2_2 = textureLoader2.load('resources/superior.png'), // ypos, Top
        texture3_2 = textureLoader2.load('resources/inferior.png'), // yneg, Bottom
        texture4_2 = textureLoader2.load('resources/right.png'), // zpos, Back
        texture5_2 = textureLoader2.load('resources/left.png'); // zneg, Front

    texture0_2.anisotropy = maxAnisotropy2;
    texture1_2.anisotropy = maxAnisotropy2;
    texture2_2.anisotropy = maxAnisotropy2;
    texture3_2.anisotropy = maxAnisotropy2;
    texture4_2.anisotropy = maxAnisotropy2;
    texture5_2.anisotropy = maxAnisotropy2;
    
    //getcontext2d. draw image
    materialArray2 = [
            new THREE.MeshBasicMaterial({
            map: texture0_2
        }),
            new THREE.MeshBasicMaterial({
            map: texture1_2
        }),
            new THREE.MeshBasicMaterial({
            map: texture2_2
        }),
            new THREE.MeshBasicMaterial({
            map: texture3_2
        }),
            new THREE.MeshBasicMaterial({
            map: texture4_2
        }),
            new THREE.MeshBasicMaterial({
            map: texture5_2
        })
    ];
    
    updateOrder(selectedPatient);
}

function updateScenes(selectedPatient, material){
	var scenes = [] //scenes is a wonderful global for now
	for (var i = 0; i < patients[selectedPatient - 1].similarity_ssim.length; i++) {
		var id = patients[selectedPatient-1].similarity_ssim[i]
		var newScene = showPatient(patients, material, id);
		scenes.push(newScene);
	}
	return scenes
}

function placeOrganModels(pOrgan, organProperties, scene, nodeColor) {

    let loader = new THREE.VTKLoader(manager);

    if (!(pOrgan == "GTVn" || pOrgan == "GTVp")) {

        loader.load('resources/models/' + pOrgan + '.vtk', function (geometry) {

            geometry.computeVertexNormals();
            geometry.center();

            let material = new THREE.MeshBasicMaterial({
                color: nodeColor,
                opacity: 0.2,
                transparent: true,
                //side: THREE.DoubleSide,
                depthTest: true,
                depthWrite: true,
                depthFunc: THREE.LessEqualDepth
            });

            let mesh = new THREE.Mesh(geometry, material);
            mesh.name = (String(pOrgan) + "_model");
            mesh.userData.type = "node_model";

            mesh.position.x = organProperties.x;
            mesh.position.y = organProperties.y;
            mesh.position.z = organProperties.z;

            mesh.rotation.x = -Math.PI / 2.0;
            mesh.rotation.z = -Math.PI / 2;


            // oral cavity
            if (pOrgan == "Tongue")
                mesh.renderOrder = -10;
            else if (pOrgan == "Genioglossus_M")
                mesh.renderOrder = -10;
            else if (pOrgan == "Lt_Ant_Digastric_M")
                mesh.renderOrder = -10;
            else if (pOrgan == "Mylogeniohyoid_M")
                mesh.renderOrder = -10;
            else if (pOrgan == "Rt_Ant_Digastric_M")
                mesh.renderOrder = -10;

            else if (pOrgan == "Extended_Oral_Cavity")
                mesh.renderOrder = -9;

            // throat
            else if (pOrgan == "Larynx")
                mesh.renderOrder = -10;
            else if (pOrgan == "Supraglottic_Larynx")
                mesh.renderOrder = -9;


            scene.add(mesh);
        });
    }
}

function showPatient(patients, materialArray, id){
	var scene = new THREE.Scene();
	var patient = patients[id-1];
	var patientOrganList = patient.organData;

	// make a list item
	var element = document.createElement("div");
	element.className = "list-item";
	element.id = id;
	element.innerHTML = template.replace('$', patient.name);
	
	var totDoseElement = element.querySelector(".totDose");
	totDoseElement.innerHTML = "Total Dose: " + "<b>" + patient.total_Dose + "</b>" + " GY";

	var tVolumeElement = element.querySelector(".tVolume");
	tVolumeElement.innerHTML = "GTV: " + "<b>" + patient.tumorVolume + "</b>" + " cc";

	var lateralityElement = element.querySelector(".laterality");
	lateralityElement.innerHTML = "<b>(" + patient.laterality + ")</b> " + " " + patient.tumorSubsite;

	// Look up the element that represents the area
	// we want to render the scene
	scene.userData.element = element.querySelector(".scene");
	//var parent = document.getElementById("content");
	
	if(!document.getElementById( element.id )){
		parent.appendChild(element);
	}

	var scalarVal = 2.4; //4.1

	var camera = new THREE.OrthographicCamera(scene.userData.element.offsetWidth / -scalarVal, 
		scene.userData.element.offsetWidth / scalarVal, 
		scene.userData.element.offsetHeight / scalarVal, 
		scene.userData.element.offsetHeight / -scalarVal, 
		1, 100000);
		
	camera.position.z = cameraDistZ;

	camera.updateProjectionMatrix();
	scene.userData.camera = camera;

	// orientation marker, patient coordinate system
	var MovingCubeMat = new THREE.MultiMaterial(materialArray);
	var MovingCubeGeom = new THREE.CubeGeometry(25, 25, 25, 1, 1, 1, materialArray);
	var MovingCube = new THREE.Mesh(MovingCubeGeom, MovingCubeMat);

	camera.add(MovingCube);
	MovingCube.position.set(121, -121, -250);

	//
	var controls = new THREE.OrbitControls(scene.userData.camera, scene.userData.element);
	controls.minDistance = 2;
	controls.maxDistance = 5000;
	controls.enablePan = false;
	controls.enableZoom = false;

	scene.userData.controls = controls;

	var geometry = new THREE.SphereGeometry(4, 16, 16);

	var material = new THREE.MeshStandardMaterial({
		color: new THREE.Color().setHex(0xa0a0a0),
		roughness: 0.5,
		metalness: 0,
		shading: THREE.FlatShading
	});

	var outlineMaterial = new THREE.MeshBasicMaterial({
		color: 0x3d3d3d,
		side: THREE.BackSide
	});

	var linkMaterial = new THREE.LineBasicMaterial({
		color: 0x3d3d3d,
		opacity: 1,
		linewidth: 3
	});
	
	for (var pOrgan in patientOrganList) {
		//this looks like it draws the organs in each patient?
		
		var organSphere = new THREE.Mesh(geometry, material.clone());

		organSphere.position.x = (patientOrganList[pOrgan].x);
		organSphere.position.y = (patientOrganList[pOrgan].y);
		organSphere.position.z = (patientOrganList[pOrgan].z);

		organSphere.name = pOrgan;
		organSphere.userData.type = "node";

		// outline
		var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

		outlineMesh.name = pOrgan + "_outline";

		if (organSphere.name == "GTVp")
			outlineMesh.scale.multiplyScalar(1.6);
		else if (organSphere.name == "GTVn")
			outlineMesh.scale.multiplyScalar(1.5);
		else
			outlineMesh.scale.multiplyScalar(1.3);

		// color
		var nodeColor;

		organSphere.userData.volume = patientOrganList[pOrgan].volume;
		organSphere.userData.minDose = patientOrganList[pOrgan].minDose;
		organSphere.userData.meanDose = patientOrganList[pOrgan].meanDose;
		organSphere.userData.maxDose = patientOrganList[pOrgan].maxDose;
		
		organSphere.userData.estimatedDose = patientOrganList[pOrgan].estimatedDose;

		// do this in python script maybe
		//grays are already in joules per kilogram?!?!? I might want to delete this because it's misleading to users
		organSphere.userData.dosePerVolume = (patientOrganList[pOrgan].meanDose / patientOrganList[pOrgan].volume).toFixed(3);

		if (organSphere.userData.meanDose >= 0.0) //null == -1 in json, pearson problems
			nodeColor = color(organSphere.userData.meanDose);
		else {
			nodeColor = '#a0a0a0'; //new THREE.Color().setHex(0xa0a0a0)
			organSphere.userData.meanDose = undefined;
			organSphere.userData.dosePerVolume = undefined;
		}

		organSphere.material.color.setStyle(nodeColor);

		scene.add(organSphere);
		organSphere.add(outlineMesh);

		placeOrganModels(pOrgan, patientOrganList[pOrgan], scene, nodeColor);

	}

	var tmp_geo = new THREE.Geometry();

	var source = scene.getObjectByName("GTVp");
	var target = scene.getObjectByName("GTVn");

	if (source != null && target != null) {
		//draws a line between the gtvp and gtvn is they are there
		tmp_geo.vertices.push(source.position);
		tmp_geo.vertices.push(target.position);

		var line = new THREE.LineSegments(tmp_geo, linkMaterial);
		line.scale.x = line.scale.y = line.scale.z = 1;
		line.originalScale = 1;

		scene.add(line);
	}

	// check for missing data
	for (var organ in oAtlas) {

		if (!patientOrganList.hasOwnProperty(organ)) {

			// node
			var organSphere = new THREE.Mesh(geometry, material.clone());

			organSphere.position.x = (oAtlas[organ].x);
			organSphere.position.y = (oAtlas[organ].y);
			organSphere.position.z = (oAtlas[organ].z);

			organSphere.name = organ;
			organSphere.userData.type = "node";

			// outline
			var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

			outlineMesh.name = organ + "_outline";

			if (organSphere.name == "GTVp")
				outlineMesh.scale.multiplyScalar(1.6);
			else if (organSphere.name == "GTVn")
				outlineMesh.scale.multiplyScalar(1.5);
			else
				outlineMesh.scale.multiplyScalar(1.3);

			// color
			var nodeColor = '#a0a0a0';

			organSphere.userData.volume = undefined;
			organSphere.userData.minDose = undefined;
			organSphere.userData.meanDose = undefined;
			organSphere.userData.maxDose = undefined;
			
			organSphere.userData.estimatedDose = undefined;
			
			organSphere.userData.dosePerVolume = undefined;

			organSphere.material.color.setStyle(nodeColor);

			scene.add(organSphere);
			organSphere.add(outlineMesh);

			placeOrganModels(organ, oAtlas[organ], scene, nodeColor);

		}

	}

	scene.add(camera);

	// light
	var light = new THREE.AmbientLight(0xffffff, 1.0); // white light   

	scene.add(light);
	return scene;

}

function removeOldViews(selectedPatientObject){
	//remove list-items not matched to the patient
	var matches = selectedPatientObject.similarity_ssim;
	var parentNode = document.getElementById('content');
	var patientViews = parentNode.getElementsByClassName('list-item');
	var element;
	console.log(patientViews);
	for(var i = patientViews.length - 1; i >= 0; i--){
		element = patientViews[i];
		if( !matches.includes( +element.id ) ){
			parentNode.removeChild(element);
		}
	}
}

function updateOrder(updatedPatient) {

    var pNames = [];

	selectedPatient = updatedPatient;
	var patientObject = patients[selectedPatient - 1];
    pRankingOrder = patientObject.similarity_ssim;
    pScores = patientObject.scores_ssim;
	scenes = updateScenes(selectedPatient, materialArray);//populates required views
	removeOldViews(patientObject); //removes old views

    var lastPatient = document.getElementById(pRankingOrder[pRankingOrder.length - 1]);
    var firstPatient = document.getElementById(pRankingOrder[0]);
    firstPatient.style.display = "none";

    //insert last element from pRankingOrder in last place (before null)
    parent.insertBefore(lastPatient, null);

    // first patient always has score of 1, clear it
    var pScoreElement = firstPatient.querySelector(".pScore");
    pScoreElement.innerHTML = "";


    for (var i = (pRankingOrder.length - 2); i >= 0; i--) {

        var first = document.getElementById(pRankingOrder[i]);
        var second = document.getElementById(pRankingOrder[i + 1]);
        // order div elements
        parent.insertBefore(first, second);

        pScoreElement = second.querySelector(".pScore");

        // update patient score
        pScoreElement.innerHTML = pScores[i + 1].toFixed(5);

        // hide patients
     
        second.style.display = "none";
    }

    var trueFirst = document.getElementById(selectedPatient);

    // if there is a tie, fix first patient
    parent.insertBefore(trueFirst, firstPatient);

    var pScoreElement1 = document.getElementById(selectedPatient).querySelector(".pScore");
    var pScoreElement2 = firstPatient.querySelector(".pScore");

    pScoreElement2.innerHTML = pScoreElement1.innerHTML;
    pScoreElement1.innerHTML = "";

    for (var j = 0; j < patientsToShow; j++) {
        var p = document.getElementById(pRankingOrder[j]);
        p.style.display = "inline-block";
        if (j <= 5) {

            var pName = p.querySelector(".description");
            pNames.push(String(j) + ": " + pName.innerHTML);
        }
    }

    pNames[0] = "0: Estimation";


    //initializeRiskPrediction(trueFirst, selectedPatient, pNames);

}

function initializeRiskPrediction(firstPatient, rank, pNames) {

    data.dates = pNames;
    data.series = [];
    
    var simScores = patients[rank - 1].scores_ssim;

    // remove scenes
    scenesRP.length = 0;

    // -----------------------------------------
    // -----------------------------------------

    var cln = firstPatient.cloneNode(true);

    cln.className = "list-item";

    cln.removeAttribute("class");
    cln.removeAttribute("id");

    cln.setAttribute("class", "list-item-RP");
    cln.value = 1;

    if (pTarget_chart.childNodes.length > 0) {
        pTarget_chart.removeChild(pTarget_chart.childNodes[0]);
    }

    pTarget_chart.appendChild(cln);

    // -----------------------------------------

    var source_scene = scenes[rank - 1];
    var target_scene = new THREE.Scene();

    target_scene.userData.element = pTarget_chart.childNodes[0].querySelector(".scene");

    for (var i = 0; i < source_scene.children.length; i++) {

        if (source_scene.children[i].userData.type == "node" ||
            source_scene.children[i].userData.type == "node_model") {

            var organ = source_scene.children[i].clone();
            organ.material.color.setStyle(source_scene.children[i].material.color);
            target_scene.add(organ);
        }
    }

    var scalarVal = 2.4; //4.1

    //var camera = new THREE.PerspectiveCamera(35, scene.userData.element.offsetWidth / scene.userData.element.offsetHeight, 1, 100000);
    var target_camera = new THREE.OrthographicCamera(target_scene.userData.element.offsetWidth / -scalarVal, target_scene.userData.element.offsetWidth / scalarVal, target_scene.userData.element.offsetHeight / scalarVal, target_scene.userData.element.offsetHeight / -scalarVal, 1, 100000);

    target_camera.position.x = scenes[0].userData.camera.position.x;
    target_camera.position.y = scenes[0].userData.camera.position.y;
    target_camera.position.z = scenes[0].userData.camera.position.z;

    target_camera.updateProjectionMatrix();
    target_scene.userData.camera = target_camera;


    var MovingCubeMat2 = new THREE.MultiMaterial(materialArray2);
    var MovingCubeGeom2 = new THREE.CubeGeometry(25, 25, 25, 1, 1, 1, materialArray2);
    var MovingCube2 = new THREE.Mesh(MovingCubeGeom2, MovingCubeMat2);

    //target_camera.add(MovingCube2);
    //MovingCube2.position.set(121, -121, -250);


    var target_controls = new THREE.OrbitControls(target_scene.userData.camera, target_scene.userData.element);
    target_controls.minDistance = 2;
    target_controls.maxDistance = 5000;
    target_controls.enablePan = false;
    target_controls.enableZoom = false;

    target_scene.userData.controls = target_controls;


    var light = new THREE.AmbientLight(0xffffff, 1.0); // white light
    target_scene.add(light);

    scenesRP.push(target_scene);


    //

    // -----------------------------------------
    // -----------------------------------------
    // make a list item
    var element = document.createElement("div");
    element.className = "list-item-RP";
    element.innerHTML = template.replace('$', "Estimation").replace('!', "");
    
    var totDoseElement = element.querySelector(".totDose");
    totDoseElement.innerHTML = "";

    var tVolumeElement = element.querySelector(".tVolume");
    tVolumeElement.innerHTML = "";

    var lateralityElement = element.querySelector(".laterality");
    lateralityElement.innerHTML = "";

    element.value = 2;

    if (pPrediction_chart.childNodes.length > 0) {
        pPrediction_chart.removeChild(pPrediction_chart.childNodes[0]);
    }

    pPrediction_chart.appendChild(element);

    // -----------------------------------------

   
    //source_scene = scenes[rank - 1];
    target_scene = new THREE.Scene();

    target_scene.userData.element = pPrediction_chart.childNodes[0].querySelector(".scene");

    for (var i = 0; i < source_scene.children.length; i++) {

        var organ = source_scene.children[i].clone();

        if (organ.userData.type == "node") {

            var organ_lc = {
                name: organ.name,
                values: []
            };

            var organTotal= 0;
			var scoreTotal = 0;
	
            for (var x = 1; x <= numMatches; x++) {

                var scene = scenes[pRankingOrder[x] - 1];

                var node = scene.getObjectByName(organ.name);

                if (node) {

                    var meanDose = node.userData.meanDose;
                    if (isNaN(meanDose)) {
                        console.log(meanDose);
                        meanDose = 0.0
                    }
                    organTotal += (meanDose * simScores[x]);
					scoreTotal += simScores[x]
                    //organAverage += (meanDose);
                    
                    //why is this storing both the actaul mean dose and the doses?
					if(x >= 5){
						organ_lc.values.push(meanDose);
					}
                }
            }

            var averageOrganDose = (organTotal / scoreTotal).toFixed(3);
            organ_lc.values.unshift(averageOrganDose);

            data.series.push(organ_lc);


            organ.userData.volume = undefined;
            organ.userData.minDose = undefined;
            organ.userData.meanDose = averageOrganDose;
            organ.userData.maxDose = undefined;

            organ.userData.dosePerVolume = undefined;

            var nodeColor = color(organ.userData.meanDose);
            organ.material = source_scene.children[i].material.clone();
            organ.material.color.setStyle(nodeColor);

            target_scene.add(organ);

            var source_model = source_scene.getObjectByName(String(organ.name) + "_model");

            if (source_model != null) {
                var target_model = source_model.clone();
                target_model.material = source_model.material.clone();
                target_model.material.color.setStyle(nodeColor);
                target_scene.add(target_model);
            }

        } 
        else {
            organ = undefined
        }
    }

    var scalarVal = 2.4; //4.1

    var target_camera = new THREE.OrthographicCamera(target_scene.userData.element.offsetWidth / -scalarVal, target_scene.userData.element.offsetWidth / scalarVal, target_scene.userData.element.offsetHeight / scalarVal, target_scene.userData.element.offsetHeight / -scalarVal, 1, 100000);

    target_camera.position.x = scenes[0].userData.camera.position.x;
    target_camera.position.y = scenes[0].userData.camera.position.y;
    target_camera.position.z = scenes[0].userData.camera.position.z;

    target_camera.updateProjectionMatrix();
    target_scene.userData.camera = target_camera;

    var MovingCubeMat2 = new THREE.MultiMaterial(materialArray2);
    var MovingCubeGeom2 = new THREE.CubeGeometry(25, 25, 25, 1, 1, 1, materialArray2);
    var MovingCube2 = new THREE.Mesh(MovingCubeGeom2, MovingCubeMat2);

    //target_camera.add(MovingCube2);
    //MovingCube2.position.set(121, -121, -250);

    var target_controls = new THREE.OrbitControls(target_scene.userData.camera, target_scene.userData.element);
    target_controls.minDistance = 2;
    target_controls.maxDistance = 5000;
    target_controls.enablePan = false;
    target_controls.enableZoom = false;

    target_scene.userData.controls = target_controls;


    var light = new THREE.AmbientLight(0xffffff, 1.0); // white light
    target_scene.add(light);

    scenesRP.push(target_scene);


    // -----------------------------------------
    // -----------------------------------------
    // make a list item
    var element = document.createElement("div");
    element.className = "list-item-RP";
    element.innerHTML = template.replace('$', "Difference").replace('!', "");
    
    var totDoseElement = element.querySelector(".totDose");
    totDoseElement.innerHTML = "";

    var tVolumeElement = element.querySelector(".tVolume");
    tVolumeElement.innerHTML = "";

    var lateralityElement = element.querySelector(".laterality");
    lateralityElement.innerHTML = "";

    element.value = 3;

    if (pDifference_chart.childNodes.length > 0) {
        pDifference_chart.removeChild(pDifference_chart.childNodes[0]);
    }

    pDifference_chart.appendChild(element);

    // -----------------------------------------

    //source_scene = scenes[rank - 1];
    target_scene = new THREE.Scene();

    target_scene.userData.element = pDifference_chart.childNodes[0].querySelector(".scene");

    for (var i = 0; i < source_scene.children.length; i++) {

        var organ = source_scene.children[i].clone();

        if (organ.userData.type == "node") {

            var organSum = 0;


            var scene1 = scenesRP[0];
            var scene2 = scenesRP[1];

            var node1 = scene1.getObjectByName(organ.name);
            var node2 = scene2.getObjectByName(organ.name);


            organSum = Math.abs(node1.userData.meanDose - node2.userData.meanDose);


            organ.userData.volume = undefined;
            organ.userData.minDose = undefined;
            organ.userData.meanDose = organSum.toFixed(3);
            organ.userData.maxDose = undefined;

            organ.userData.dosePerVolume = undefined;

            var nodeColor = color2(organ.userData.meanDose);
            organ.material = source_scene.children[i].material.clone();
            organ.material.color.setStyle(nodeColor);

            target_scene.add(organ);

            var source_model = source_scene.getObjectByName(String(organ.name) + "_model");

            if (source_model != null) {
                var target_model = source_model.clone();
                target_model.material = source_model.material.clone();
                target_model.material.color.setStyle(nodeColor);
                target_scene.add(target_model);
            }

        } else {
            organ = undefined;
        }

    }

    var scalarVal = 2.4; //4.1

    target_camera = new THREE.OrthographicCamera(target_scene.userData.element.offsetWidth / -scalarVal, target_scene.userData.element.offsetWidth / scalarVal, target_scene.userData.element.offsetHeight / scalarVal, target_scene.userData.element.offsetHeight / -scalarVal, 1, 100000);

    target_camera.position.x = scenes[0].userData.camera.position.x;
    target_camera.position.y = scenes[0].userData.camera.position.y;
    target_camera.position.z = scenes[0].userData.camera.position.z;

    target_camera.updateProjectionMatrix();
    target_scene.userData.camera = target_camera;

    var MovingCubeMat2 = new THREE.MultiMaterial(materialArray2);
    var MovingCubeGeom2 = new THREE.CubeGeometry(25, 25, 25, 1, 1, 1, materialArray2);
    var MovingCube2 = new THREE.Mesh(MovingCubeGeom2, MovingCubeMat2);

    var target_controls = new THREE.OrbitControls(target_scene.userData.camera, target_scene.userData.element);
    target_controls.minDistance = 2;
    target_controls.maxDistance = 5000;
    target_controls.enablePan = false;
    target_controls.enableZoom = false;

    target_scene.userData.controls = target_controls;


    var light = new THREE.AmbientLight(0xffffff, 1.0); // white light
    target_scene.add(light);

    scenesRP.push(target_scene);


    multiLineChart();


}


function multiLineChart() {

    d3.selectAll("#mLineChart > *").remove();

    var svg = d3.select("#mLineChart"),
        margin = {
            top: 20,
            right: 60,
            bottom: 50,
            left: 60
        },
        width = svg.attr("width") - margin.left - margin.right,
        height = svg.attr("height") - margin.top - margin.bottom,
        g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var x = d3.scalePoint().range([0, width]).padding([.2]);

    var y = d3.scaleLinear()
        .range([height, 0]);

    var line = d3.line()
        .x(function (d, i) {
            //return x(d.date);
            return x(data.dates[i]);
        })
        .y(function (d) {
            //return y(d.value);
            return y(d);
        });

    // --------------------------------

    x.domain(data.dates);
    x.invert = (function () {
        var domain = x.domain();
        var range = x.range();
        var scale = d3.scaleQuantize().domain(range).range(domain);

        return function (x) {
            return scale(x);
        }
    })();
  
    y.domain([0, d3.max(data.series, d => d3.max(d.values))]).nice();

    g.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x).tickSizeOuter([0]))
        .append("text")
        .attr("x", width + 5)
        .attr("y", 0.5)
        .attr("dy", "0.32em")
        .style("text-anchor", "start")
        .style("fill", "#000")
        .style("font-weight", "bold")
        .text("Patients");

    svg.append("text")
        .attr("transform", "translate(" + (width / 2 + 50) + " ," + (20) + ")")
        .style("text-anchor", "middle")
        .text("Organ Dose Across Most Similar Patients");

    g.append("g")
        .attr("class", "axis axis--y")
        //.call(d3.axisLeft(y).ticks(10, "%"))
        .call(d3.axisLeft(y).ticks(10))
        .append("text")
        .attr("x", 4)
        .attr("y", 0.5)
        .attr("dy", "0.32em")
        .style("text-anchor", "start")
        .style("fill", "#000")
        .style("font-weight", "bold")
        .text(data.y);

    const path = g.append("g")
        .attr("class", "cities")
        .selectAll("path")
        //.data(data_lc)
        .data(data.series)
        .enter().append("path")
        .attr("d", function (d) {
            d.line = this;
            return line(d.values);
        })
        .attr("id", function (d) {
            return ("line_" + d.name);
        })
        .attr("display", function (d) {

            var el = document.getElementById(d.name + "_checkList")

            if (d.name == "GTVp" || d.name == "GTVn")
                return (null);
            else if (el && !el.checked)
                return ("none");
            else
                return (null);
        });

    svg.call(hover, path, x, y, data);

}

function hover(svg, path, x, y) {

    svg
        .style("position", "relative");

    if ("ontouchstart" in document) svg
        .style("-webkit-tap-highlight-color", "transparent")
        .on("touchmove", moved)
        .on("touchstart", entered)
        .on("touchend", left)
    else svg
        .on("mousemove", moved)
        .on("mouseenter", entered)
        .on("mouseleave", left);

    const dot = svg.append("g")
        .attr("display", "none")
        .attr("class", "dot");

    dot.append("circle")
        .attr("r", 3.5)
        .style("fill", "#00e4ff");

    dot.append("text")
        .style("font", "15px sans-serif")
        .attr("text-anchor", "middle")
        .attr("y", -8);

    function moved() {
        d3.event.preventDefault();

        var xAxisOrder = [...Array(5).keys()];

        const ym = y.invert(d3.event.offsetY - 20);
        const xm = x.invert(d3.event.offsetX - 60);
        const i1 = d3.bisectRight(data.dates, xm, 1);
        const i0 = i1 - 1;
        const i = xm - data.dates[i0] > data.dates[i1] - xm ? i1 : i0;
        const s = data.series.reduce((a, b) => Math.abs(a.values[i] - ym) < Math.abs(b.values[i] - ym) ? a : b);
        path.attr("stroke", d => d === s ? "#00e4ff" : "#d1d1d1").attr("stroke-width", d => d === s ? "4px" : "2px").filter(d => d === s).raise();
        
        dot.attr("transform", `translate(${x(data.dates[i])+60},${y(s.values[i])+20})`);
        dot.select("text").text(s.name + " (" + s.values[i] + " GY)");

        scenesRP.forEach(function (scene, index) {

            for (var i = 0; i < scene.children.length; i++) {


                for (var j = 0; j < scene.children[i].children.length; j++) {

                    if (scene.children[i].children[j].name.includes("_outline")) {

                        var organ = scene.children[i].children[j];
                        organ.material.color.setHex(0x3d3d3d);
                    }
                }
            }

        });

        scenesRP.forEach(function (scene, index) {

            var tempObject = scene.getObjectByName(s.name + "_outline");

            tempObject.material.color.setHex(0x00e4ff);

        });

    }

    function entered() {

        path.attr("stroke", null).attr("stroke-width", "4px");
        dot.attr("display", null);
    }

    function left() {

        path.attr("stroke", "#d1d1d1").attr("stroke-width", "2px");
        dot.attr("display", "none");

        scenesRP.forEach(function (scene, index) {

            for (var i = 0; i < scene.children.length; i++) {

                for (var j = 0; j < scene.children[i].children.length; j++) {

                    if (scene.children[i].children[j].name.includes("_outline")) {

                        var organ = scene.children[i].children[j];
                        organ.material.color.setHex(0x3d3d3d);
                    }
                }
            }
        });

    }

}

function animate() {
    render();
    requestAnimationFrame(animate);
}

function render() {

    updateSize();

    renderer.setClearColor(0xffffff);
    renderer.setScissorTest(false);
    renderer.clear();

    renderer.setClearColor(0xa5a5a5);
    renderer.setScissorTest(true);

    updateMainView();

    renderer2.setClearColor(0x444444, 0);
    renderer2.setScissorTest(false);
    renderer2.clear();

    renderer2.setClearColor(0xa5a5a5);
    renderer2.setScissorTest(true);

    updateRiskPView();

}

function updateMainView(rotMatrix) {

    var rotMatrix = new THREE.Matrix4();
	//scenes = updateScenes(selectedPatient, materialArray);
    pRankingOrder = patients[selectedPatient - 1].similarity_ssim;

	for (var index = 0; index < scenes.length; index++) {

		var scene = scenes[index];
		var controls = scene.userData.controls;
		var camera = scene.userData.camera;

		var orientMarkerCube = camera.children[0];

		// get the element that is a place holder for where we want to
		// draw the scene
		var element = scene.userData.element;

		// get its position relative to the page's viewport
		var rect = element.getBoundingClientRect();

		// check if it's offscreen. If so skip it
		if (rect.bottom < 0 || rect.top > renderer.domElement.clientHeight ||
			rect.right < 0 || rect.left > renderer.domElement.clientWidth) {

			return; // it's off screen
		}

		// update orientation marker
		rotMatrix.extractRotation(controls.object.matrix);
		orientMarkerCube.rotation.setFromRotationMatrix(rotMatrix.transpose());

		// set the viewport
		var width = rect.right - rect.left;
		var height = rect.bottom - rect.top;
		var left = rect.left;
		var bottom = renderer.domElement.clientHeight - rect.bottom;

		renderer.setViewport(left, bottom, width, height);
		renderer.setScissor(left, bottom, width, height);

		// raycaster
		raycaster.setFromCamera(mouseNorm, currScene.userData.camera);

		var intersects = raycaster.intersectObjects(currScene.children);

		if (intersects.length >= 1 && detailsOnRotate) {

			for (var i = intersects.length - 1; i >= 0; i--) {

				if (intersects[i].object.userData.type == "node") {

					nodeHover = intersects[i].object;
					var tempObject = scene.getObjectByName(nodeHover.name + "_outline");

					if (INTERSECTED != tempObject) {

						if (INTERSECTED) {
							INTERSECTED.material.color.setHex(INTERSECTED.currentHex);
						}

						INTERSECTED = tempObject;

						if (INTERSECTED) {
							INTERSECTED.currentHex = INTERSECTED.material.color.getHex();
							INTERSECTED.material.color.setHex(0x00e4ff);
						}

						// details

						populateAndPlaceDetails("SHOW");

					}

					break;

				} else {
					populateAndPlaceDetails("HIDE");
				}


			}
		} else {

			if (INTERSECTED) {
				INTERSECTED.material.color.setHex(INTERSECTED.currentHex);
				// details
				populateAndPlaceDetails("HIDE");

			}

			INTERSECTED = null;
		}

		renderer.render(scene, camera);
	}


}

function updateRiskPView(rotMatrix) {

    var rotMatrix = new THREE.Matrix4();

    scenesRP.forEach(function (scene, index) {

        var scene = scene;
        var controls = scene.userData.controls;
        var camera = scene.userData.camera;

        var orientMarkerCube = camera.children[0];

        // get the element that is a place holder for where we want to
        // draw the scene
        var element = scene.userData.element;

        // get its position relative to the page's viewport
        var rect = element.getBoundingClientRect();

        // check if it's offscreen. If so skip it
        if (rect.bottom < 0 || rect.top > renderer2.domElement.clientHeight ||
            rect.right < 0 || rect.left > renderer2.domElement.clientWidth) {

            return; // it's off screen
        }

        // update orientation marker
        rotMatrix.extractRotation(controls.object.matrix);

        // set the viewport
        var width = rect.right - rect.left;
        var height = rect.bottom - rect.top;
        var left = rect.left;
        var bottom = renderer2.domElement.clientHeight - rect.bottom;

        renderer2.setViewport(left, bottom, width, height);
        renderer2.setScissor(left, bottom, width, height);

        // raycaster
        raycaster.setFromCamera(mouseNorm, currScene.userData.camera);

        var intersects = raycaster.intersectObjects(currScene.children);

        if (intersects.length >= 1 && detailsOnRotate) {

            for (var i = intersects.length - 1; i >= 0; i--) {

                if (intersects[i].object.userData.type == "node") {
                    nodeHover = intersects[i].object;

                    var tempObject = scene.getObjectByName(nodeHover.name + "_outline");

                    if (INTERSECTED != tempObject) {

                        if (INTERSECTED) {
                            INTERSECTED.material.color.setHex(INTERSECTED.currentHex);

                            for (var organ in oAtlas) {

                                d3.select("#line_" + organ).attr("stroke", "#d1d1d1");
                                d3.select("#line_" + organ).style("stroke", null).style("stroke-width", null);
                            }
                        }

                        INTERSECTED = tempObject;

                        if (INTERSECTED) {
                            INTERSECTED.currentHex = INTERSECTED.material.color.getHex();
                            INTERSECTED.material.color.setHex(0x00e4ff);
                            d3.select("#line_" + nodeHover.name).style("stroke", "#00e4ff").style("stroke-width", "4px").raise();
                     
                        }
                    }

                    break;


                } 

            }
        } else {

            if (INTERSECTED) {
                INTERSECTED.material.color.setHex(INTERSECTED.currentHex);

            }

            INTERSECTED = null;
            for (var organ in oAtlas) {

                d3.select("#line_" + organ).style("stroke", null).style("stroke-width", null);
            }
            
        }

        renderer2.render(scene, camera);

    });


}

function updateSize() {

    var width = canvas.clientWidth;
    var height = canvas.clientHeight;

    if (canvas.width !== width || canvas.height != height)
        renderer.setSize(width, height, false);

    width = canvas2.clientWidth;
    height = canvas2.clientHeight;

    if (canvas2.width !== width || canvas2.height != height)
        renderer2.setSize(width, height, false);
}

function populateAndPlaceDetails(state) {

    if (state == "SHOW") {

        nodeDetails.style.display = "block";
        // PLACEMENT
        // check if details are offscreen, then shift appropriately
        // X, add 10 pixels for buffer, since width is dynamic
        if (mouse.x + detailsOffsetX + nodeDetails.offsetWidth + 5 >= canvas.clientWidth) {
            nodeDetails.style.left = (mouse.x - detailsOffsetX - nodeDetails.offsetWidth) + "px";
        } else {
            nodeDetails.style.left = (mouse.x + detailsOffsetX) + "px";
        }

        // Y
        if (mouse.y + detailsOffsetY + nodeDetails.offsetHeight + 5 >= canvas.clientHeight) {
            nodeDetails.style.top = (mouse.y - detailsOffsetY - nodeDetails.offsetHeight) + "px";
        } else {
            nodeDetails.style.top = (mouse.y + detailsOffsetY) + "px";
        }
       
        // POPULATE

        // Organ name
        organName.innerHTML = nodeHover.name;

        // Dose Per Volume
        dosePerVolume.innerHTML = nodeHover.userData.dosePerVolume;

        // line separator
        lineSeparator.style["borderColor"] = "#" + nodeHover.material.color.getHexString();

        // Volume
        volumeVal.innerHTML = nodeHover.userData.volume + "";

        // Mean Dose
        meanDoseVal.innerHTML = nodeHover.userData.meanDose + "  GY";

        // Min Dose
        minDoseVal.innerHTML = nodeHover.userData.minDose + "";

        // Max Dose
        maxDoseVal.innerHTML = nodeHover.userData.maxDose + "";

    } else if (state == "HIDE") {

        nodeDetails.style.display = "none";
        nodeDetails.style.top = -500 + "px";
        nodeDetails.style.left = -500 + "px";
    }
}

function onMouseDown(event) {

    if (event.target) {

        detailsOnRotate = false;
        handleInputRotate(event);
    }
}

function onTouchStart(event) {

    if (event.target) {

        detailsOnRotate = false;
        handleInputRotate(event);
    }
}

function handleInputRotate(event) {

    var targ = event.target,
        cameraToCopy;

    if (targ.className == "scene" && syncCameras == true) {
		var index;
        if (targ.parentNode.hasAttribute("id")) {
			index = patients[selectedPatient - 1].similarity_ssim.indexOf( +targ.parentNode.id );
            cameraToCopy = scenes[index].userData.camera;

        } else {
            cameraToCopy = scenesRP[targ.parentNode.value - 1].userData.camera;
        }

        // 20 milliseconds interval => 50 FPS
        syncCamerasInterval = setInterval(syncAllCameras, 20, cameraToCopy);
    }
}

function syncAllCameras(cameraToCopy) {

    for( var i = 0; i < scenes.length; i++) {

        var scene = scenes[i];
        var camera = scene.userData.camera;
        var controls = scene.userData.controls;

        camera.position.subVectors(cameraToCopy.position, controls.target);
        camera.position.setLength(cameraDistZ);
        camera.lookAt(scene.position);
    };

    scenesRP.forEach(function (scene, index) {
        var camera = scene.userData.camera;
        var controls = scene.userData.controls;

        camera.position.subVectors(cameraToCopy.position, controls.target);
        camera.position.setLength(cameraDistZ);
        camera.lookAt(scene.position);
    });
}

function onMouseUp(event) {

    detailsOnRotate = true;
    clearInterval(syncCamerasInterval);
}

function onTouchEnd(event) {

    detailsOnRotate = true;
    clearInterval(syncCamerasInterval);
}

function onDocumentMouseMove(event) {

    //event.preventDefault();

    if (event.target) {

        var targ = event.target;

        if (targ.className == "scene") {

            if (targ.parentNode.hasAttribute("id")) {
				let index = patients[selectedPatient - 1].similarity_ssim.indexOf( +targ.parentNode.id );
                currScene = scenes[index];
            } else {
                currScene = scenesRP[targ.parentNode.value - 1];
            }

            mouse.x = event.clientX;
            mouse.y = event.clientY;

            mouseNorm.x = (event.offsetX / targ.offsetWidth) * 2 - 1;
            mouseNorm.y = -(event.offsetY / targ.offsetHeight) * 2 + 1;
        }
    }
}

var opacSlider = document.getElementById("opacSlider");

opacSlider.oninput = function () {

    var opac = (this.value / 100.0);

    scenes.forEach(function (scene, index) {

        for (var pOrgan in oAtlas) {

            var tempObject = scene.getObjectByName(pOrgan + "_model");

            if (tempObject)
                tempObject.material.opacity = opac;

        }


    });

    scenesRP.forEach(function (scene, index) {

        for (var pOrgan in oAtlas) {

            var tempObject = scene.getObjectByName(pOrgan + "_model");

            if (tempObject)
                tempObject.material.opacity = opac;

        }


    });
}

var onAfterRender = function (renderer, scene, camera, geometry, material, group) {};
