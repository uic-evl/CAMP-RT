'use strict';

// SOME QUICK, TEMPORARY NOTES:
//      update color scale/range
//      top navbar slides down to show extended view of selected patient?
//          shows volume rendered detailed models of organs
//          exploded view
//          on left show list of all organs, highlights organ when hovering over node

if (!Detector.webgl) {
    Detector.addGetWebGLMessage();
}

var canvas;

//console.log(JSON.stringify());

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
    maxDoseVal = document.getElementById("details_maxDose_val");

var scenes = [],
    renderer;

var selectedPatient = 1;

var syncCameras = true,
    syncCamerasInterval,
    detailsOnRotate = true;

var raycaster;

var mouse = new THREE.Vector2(-500, -500);

var mouseNorm = new THREE.Vector2(-500, -500),
    INTERSECTED = null,
    nodeHover;

var width, height;

var cameraDistZ = 300;

// 36 steps

//var domainColorScale = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0];
var domainColorScale = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105];
var rangeColorScale = ['#ffffe0', '#fff8d2', '#fff0c4', '#ffe9b8', '#ffe2ae', '#ffdaa3', '#ffd39a', '#ffcb91', '#ffc389', '#ffbb82', '#ffb27c', '#ffab77', '#ffa272', '#ff986e', '#fe906a', '#fb8768', '#f98065', '#f67762', '#f26f60', '#ee675d', '#eb5f5b', '#e75758', '#e25055', '#dd4852', '#d8404e', '#d3394a', '#cc3146', '#c62a41', '#c0223b', '#b91c35', '#b3152f', '#ab0e28', '#a40820', '#9b0317', '#93010e', '#8b0000'];

var color = d3.scaleLinear()
    .domain(domainColorScale)
    .range(rangeColorScale);


// data
var organs, links, patients;

var organRef = [];

var pRankingOrder, pScores;

var listItems, arrayOfDivs = [],
    currScene;

d3.queue()
    .defer(d3.json, "data/organs.json")
    .defer(d3.json, "data/links.json")
    .defer(d3.json, "data/patients_V5.json")
    .await(start);

function start(error, organsData, linksData, patientsData) {
    if (error) return alert("Data invalid: " + error);

    organs = organsData;
    links = linksData;
    patients = patientsData;

    pRankingOrder = patients[selectedPatient - 1].similarity;
    pScores = patients[selectedPatient - 1].scores;

    populateDropDownMenu();
    populateColorScale();

    flipGraph(); // fixes orientation of organs
    computeCenterOfGraphAndShift(); // compute center of graph and shift to origin
    //shiftGraphToOrigin(); // center graph to origin
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

    updateOrder(selectedPatient); // update order in GUI

    document.getElementById("loadScreen").style.display = "none";

    animate(); // render
}

// ----------------------------------------------------------------

function populateColorScale() {

    var parentDiv = document.getElementById("colorScale");

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

    for (var i = 0, ref = colorScaleEntries.length; i < ref; i++) {
        colorScaleEntries[i].addEventListener('mouseover', showColorScaleLabel, false);
        colorScaleEntries[i].addEventListener('mouseout', hideColorScaleLabel, false);
    }
}

function showColorScaleLabel(event) {

    //console.log(event.target.value);

    var details = document.getElementById("colorScaleDetails");


    details.style.left = event.target.offsetLeft + "px";
    details.innerHTML = "" + event.target.value;


    details.style.display = "block";

}

function hideColorScaleLabel(event) {

    //console.log(event.target.value);

    var details = document.getElementById("colorScaleDetails");

    details.style.display = "none";

}

function populateDropDownMenu() {

    var menu = document.getElementById("patientMenu");

    patients.forEach(function (patient, index) {

        var tempOption = document.createElement("option");

        //tempOption.value = index + 1;
        tempOption.value = patient.ID_internal;

        tempOption.innerHTML = patient.name;

        menu.appendChild(tempOption);
    });

}

function populateOrganMasterList() {

    var master = document.getElementById("masterList");

    for (var i = 0; i < patients.length; i++) {

        var patientOrganList = patients[i].organData;

        for (var pOrgan in patientOrganList) {

            // pOrgan == string name of organ
            // patientOrganList[pOrgan] == the properties of current object

            if (!organRef.includes(pOrgan))
                organRef.push(pOrgan);

        }

        //organRef.sort();
        console.log(organRef);
    }

    organs.forEach(function (organ, index) {

        var tempDiv = document.createElement("div");

        tempDiv.setAttribute("class", "checkboxContainer");
        tempDiv.setAttribute("id", organ.name + "_Master");


        var tempInput = document.createElement("INPUT");

        tempInput.setAttribute("type", "checkbox");
        tempInput.setAttribute("id", organ.name);
        tempInput.setAttribute("value", organ.name);
        tempInput.setAttribute("name", "organMasterList");
        tempInput.setAttribute("onchange", "handleCheckBox(this)");


        var tempLabel = document.createElement("label");

        tempLabel.setAttribute("for", organ.name);
        tempLabel.innerHTML = organ.name;

        // ----------

        tempDiv.appendChild(tempInput);
        tempDiv.appendChild(tempLabel);

        master.appendChild(tempDiv);

    });

    //checkOrganMasterList(); // to see which organs in the list should be checked initially
    formatOrganMasterList(); // alternate background color for better reading
}

function checkOrganMasterList() {

    //var master = document.getElementById("masterList");

    //var organList = master.getElementsByClassName("list-item");
    //var organList = master.children;

    //console.log(organList.length);

    organs.forEach(function (organ, index) {

        var tempItem = document.getElementById(organ.name);

        //if (scenes[selectedPatient - 1].getObjectByName(organ.name).userData.dosePerVolume != null) {
        if (scenes[selectedPatient - 1].getObjectByName(organ.name).userData.meanDose != null) {

            if (tempItem.checked != true)
                tempItem.setAttribute("checked", true);

        } else {

            if (tempItem.checked != false)
                tempItem.setAttribute("checked", false);
        }
    });

}

function formatOrganMasterList() {
    var master = document.getElementById("masterList");
    var organList = master.children;

    for (var i = 0; i < organList.length; i++) {

        if (i % 2 == 0)
            organList[i].style["backgroundColor"] = "#3a3a3a";
        else
            organList[i].style["backgroundColor"] = "#444444";

    }

}

function handleCheckBox(event) {

    if (event.checked) {

        scenes.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);

            node.visible = true;
        });

    } else {

        scenes.forEach(function (scene, index) {

            var node = scene.getObjectByName(event.value);

            node.visible = false;
        });
    }
}

function flipGraph() {

    for (var i = 0; i < patients.length; i++) {

        var patientOrganList = patients[i].organData;

        for (var pOrgan in patientOrganList) {

            // pOrgan == string name of organ
            // patientOrganList[pOrgan] == the properties of current object

            var tOrganX = (patientOrganList[pOrgan].x * -1);
            var tOrganY = (patientOrganList[pOrgan].y * -1);
            var tOrganZ = (patientOrganList[pOrgan].z * -1);

            patientOrganList[pOrgan].x = tOrganY;
            patientOrganList[pOrgan].y = tOrganZ;
            patientOrganList[pOrgan].z = tOrganX;
        }
    }
}

function computeCenterOfGraphAndShift() {

    for (var i = 0; i < patients.length; i++) {

        var sceneCenter = [0.0, 0.0, 0.0];

        var xyzMin = new Array(3);
        var xyzMax = new Array(3);

        var positions = [];

        var patientOrganList = patients[i].organData;

        for (var pOrgan in patientOrganList) {

            // pOrgan == string name of organ
            // patientOrganList[pOrgan] == the properties of current object

            var xyz = {
                x: patientOrganList[pOrgan].x,
                y: patientOrganList[pOrgan].y,
                z: patientOrganList[pOrgan].z
            };

            positions.push(xyz);
        }

        xyzMin = getMin(positions);
        xyzMax = getMax(positions);

        console.log(patients[i].ID);
        console.log(xyzMin);
        console.log(xyzMax);


        sceneCenter = [
            ((xyzMin[0] + xyzMax[0]) / 2),
            ((xyzMin[1] + xyzMax[1]) / 2),
            ((xyzMin[2] + xyzMax[2]) / 2)
        ];

        for (var pOrgan in patientOrganList) {

            // pOrgan == string name of organ
            // patientOrganList[pOrgan] == the properties of current object

            patientOrganList[pOrgan].x = (patientOrganList[pOrgan].x - sceneCenter[0]);
            patientOrganList[pOrgan].y = (patientOrganList[pOrgan].y - sceneCenter[1]);
            patientOrganList[pOrgan].z = (patientOrganList[pOrgan].z - sceneCenter[2]);
        }

        //console.log(positions);
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

function shiftGraphToOrigin() {

    //organs.forEach(function (organ, index) {

    //    organ.x = (organ.x - sceneCenter[0]);
    //    organ.y = (organ.y - sceneCenter[1]);
    //    organ.z = (organ.z - sceneCenter[2]);
    //});
}

function init() {

    canvas = document.getElementById("c");

    var template = document.getElementById("template").text;

    renderer = new THREE.WebGLRenderer({
        canvas: canvas,
        antialias: true
    });

    renderer.setClearColor(0xffffff, 1);
    renderer.setPixelRatio(window.devicePixelRatio);
    //renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.sortObjects = false;

    raycaster = new THREE.Raycaster();

    //
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
    var materialArray = [
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

    for (var i = 0; i < patients.length; i++) {

        var scene = new THREE.Scene();

        var patientOrganList = patients[i].organData;


        // make a list item
        var element = document.createElement("div");
        element.className = "list-item";
        //element.id = patients[i].id;
        element.id = patients[i].ID_internal;
        element.innerHTML = template.replace('$', patients[i].name);

        // Look up the element that represents the area
        // we want to render the scene
        scene.userData.element = element.querySelector(".scene");

        parent.appendChild(element);

        var scalarVal = 4.1;

        var camera = new THREE.PerspectiveCamera(35, scene.userData.element.offsetWidth / scene.userData.element.offsetHeight, 1, 100000);
        //var camera = new THREE.OrthographicCamera(scene.userData.element.offsetWidth / -scalarVal, scene.userData.element.offsetWidth / scalarVal, scene.userData.element.offsetHeight / scalarVal, scene.userData.element.offsetHeight / -scalarVal, 1, 100000);
        camera.position.z = cameraDistZ;
        //camera.position.z = cameraDistZ; // 2000

        camera.updateProjectionMatrix();
        scene.userData.camera = camera;

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

        //console.log(typeof patientOrganList);

        for (var pOrgan in patientOrganList) {

            // pOrgan == string name of organ
            // patientOrganList[pOrgan] == the properties of current object

            //console.log(patientOrganList[pOrgan]);

            // node
            var organSphere = new THREE.Mesh(geometry, material.clone());

            organSphere.position.x = (patientOrganList[pOrgan].x);
            organSphere.position.y = (patientOrganList[pOrgan].y);
            organSphere.position.z = (patientOrganList[pOrgan].z);

            organSphere.name = pOrgan;
            organSphere.userData.type = "node";

            // outline
            var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

            outlineMesh.name = pOrgan + "_outline";

            if (organSphere.name == "GTVn" || organSphere.name == "GTVp")
                outlineMesh.scale.multiplyScalar(1.4);
            else
                outlineMesh.scale.multiplyScalar(1.15);

            //outlineMesh.scale.multiplyScalar(1.15);

            // color
            var nodeColor;

            //if (patientOrganList[organ.name] != null) {
            //console.log(patientOrganList[organ.name]);

            organSphere.userData.volume = undefined;
            organSphere.userData.minDose = undefined;
            organSphere.userData.meanDose = patientOrganList[pOrgan].meanDose;
            organSphere.userData.maxDose = undefined;
            organSphere.userData.dosePerVolume = undefined;

            if (organSphere.userData.meanDose >= 0.0) //null == -1 in json, pearson problems
                nodeColor = color(organSphere.userData.meanDose);
            else
                nodeColor = '0xa0a0a0'; //new THREE.Color().setHex(0xa0a0a0)

            //if (organSphere.name == "GTVn" || organSphere.name == "GTVp")
            //    nodeColor = '#000000';

            organSphere.material.color.setStyle(nodeColor);

            //} else {
            //console.log("no");

            //organSphere.userData.dosePerVolume = null;
            //nodeColor = "rgb(131, 131, 131)";
            //organSphere.material.color.setStyle(nodeColor);

            //organSphere.visible = false;
            //}

            scene.add(organSphere);
            organSphere.add(outlineMesh);


        }

        var tmp_geo = new THREE.Geometry();

        var source = scene.getObjectByName("GTVp");
        var target = scene.getObjectByName("GTVn");

        if (source != null && target != null) {

            tmp_geo.vertices.push(source.position);
            tmp_geo.vertices.push(target.position);

            var line = new THREE.LineSegments(tmp_geo, linkMaterial);
            line.scale.x = line.scale.y = line.scale.z = 1;
            line.originalScale = 1;

            //line.frustumCulled = false;

            scene.add(line);
        }


        //for (var z = 0; z < patientOrganList.length; z++) {
        //var organ = patientOrganList[i];
        //}


        /*
                // nodes, node outline, and color
                organs.forEach(function (organ, index) {

                    // node
                    var organSphere = new THREE.Mesh(geometry, material.clone());

                    organSphere.position.x = (organ.x);
                    organSphere.position.y = (organ.y);
                    organSphere.position.z = (organ.z);

                    organSphere.name = organ.name;
                    organSphere.userData.type = "node";

                    // outline
                    var outlineMesh = new THREE.Mesh(geometry, outlineMaterial.clone());

                    outlineMesh.name = organ.name + "_outline";

                    outlineMesh.scale.multiplyScalar(1.1);

                    // color
                    var nodeColor;

                    if (patientOrganList[organ.name] != null) {
                        //console.log(patientOrganList[organ.name]);

                        organSphere.userData.volume = patientOrganList[organ.name].volume;
                        organSphere.userData.minDose = patientOrganList[organ.name].minDose;
                        organSphere.userData.meanDose = patientOrganList[organ.name].meanDose;
                        organSphere.userData.maxDose = patientOrganList[organ.name].maxDose;
                        organSphere.userData.dosePerVolume = patientOrganList[organ.name].dosePerVolume;

                        nodeColor = color(organSphere.userData.dosePerVolume);
                        organSphere.material.color.setStyle(nodeColor);

                    } else {
                        //console.log("no");

                        organSphere.userData.dosePerVolume = null;
                        nodeColor = "rgb(131, 131, 131)";
                        organSphere.material.color.setStyle(nodeColor);

                        organSphere.visible = false;
                    }

                    scene.add(organSphere);
                    organSphere.add(outlineMesh);

                });
                
                */


        // links
        /*
            links.forEach(function (link, index) {

                var tmp_geo = new THREE.Geometry();

                var source = scene.getObjectByName(link.source);
                var target = scene.getObjectByName(link.target);

                tmp_geo.vertices.push(source.position);
                tmp_geo.vertices.push(target.position);

                var line = new THREE.LineSegments(tmp_geo, linkMaterial);
                line.scale.x = line.scale.y = line.scale.z = 1;
                line.originalScale = 1;

                //line.frustumCulled = false;

                scene.add(line);
            });
        
            */


        scene.add(camera);

        // orientation marker, patient coordinate system
        var MovingCubeMat = new THREE.MultiMaterial(materialArray);
        var MovingCubeGeom = new THREE.CubeGeometry(15, 15, 15, 1, 1, 1, materialArray);
        var MovingCube = new THREE.Mesh(MovingCubeGeom, MovingCubeMat);

        camera.add(MovingCube);
        //MovingCube.position.set(65, -65, -100);
        MovingCube.position.set(100, -60, -250);

        // light
        var light = new THREE.AmbientLight(0xffffff, 1.0); // white light
        scene.add(light);


        scenes.push(scene);
    }




    ///







    ///

    ///models
    // Load STL model

    //var tempObject = scenes[0].getObjectByName("Upper_Lip");

    //console.log("loading model");
    //var loaderSTL = new THREE.STLLoader();
    //loaderSTL.load('resources/models/Upper_Lip.stl',
    //    function (geometry) {
    //        var material = new THREE.MeshPhongMaterial({
    //            color: 0xF44336,
    //            specular: 0x111111,
    //            shininess: 200
    //       });
    //        var mesh = new THREE.Mesh(geometry, material);
    // to LPS space
    /*
                var RASToLPS = new THREE.Matrix4();
                RASToLPS.set(-1, 0, 0, 0,
                    0, -1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
                mesh.applyMatrix(RASToLPS);
    */

    //mesh.position.x = organs[23].x;
    //mesh.position.y = organs[23].y;
    //mesh.position.z = organs[23].z;

    //mesh.position.set(tempObject.position.x, tempObject.position.y, tempObject.position.z);

    //mesh.position.set(0, 0, 0);



    //scenes[0].add(mesh);
    //  });

    ///


    //renderer.autoClear = false;
}

function updateOrder(updatedPatient) {

    selectedPatient = updatedPatient;
    pRankingOrder = patients[selectedPatient - 1].similarity;
    pScores = patients[selectedPatient - 1].scores;

    var lastPatient = document.getElementById(pRankingOrder[pRankingOrder.length - 1]);
    var firstPatient = document.getElementById(pRankingOrder[0]);

    //insert last element from pRankingOrder in last place (before null)
    parent.insertBefore(lastPatient, null);

    var pScoreElement = firstPatient.querySelector(".pScore");

    // first patient always has score of 1, clear it
    pScoreElement.innerHTML = "";

    for (var i = (pRankingOrder.length - 2); i >= 0; i--) {

        var first = document.getElementById(pRankingOrder[i]);
        var second = document.getElementById(pRankingOrder[i + 1]);

        // order div elements
        parent.insertBefore(first, second);

        pScoreElement = second.querySelector(".pScore");

        // update patient score
        pScoreElement.innerHTML = pScores[i + 1].toFixed(5);
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

    pRankingOrder = patients[selectedPatient - 1].similarity;

    var rotMatrix = new THREE.Matrix4();

    pRankingOrder.forEach(function (rank, index) {

        var scene = scenes[rank - 1];
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

        //controls.update();

        // raycaster
        //raycaster.setFromCamera(mouseNorm, camera);
        raycaster.setFromCamera(mouseNorm, currScene.userData.camera);

        //var intersects = raycaster.intersectObjects(scene.children);
        var intersects = raycaster.intersectObjects(currScene.children);

        if (intersects.length >= 1 && intersects[0].object.userData.type == "node" && detailsOnRotate) {

            nodeHover = intersects[0].object;
            var tempObject = scene.getObjectByName(nodeHover.name + "_outline");
            //var tempObject = nodeHover.children[0]; // this breaks something with details?

            if (INTERSECTED != tempObject) {

                if (INTERSECTED)
                    INTERSECTED.material.color.setHex(INTERSECTED.currentHex);

                INTERSECTED = tempObject;

                INTERSECTED.currentHex = INTERSECTED.material.color.getHex();
                INTERSECTED.material.color.setHex(0x00e4ff);

                // details

                populateAndPlaceDetails("SHOW");

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
    });

    //renderer.render( currScene, currScene.userData.camera );
}

function updateSize() {

    width = canvas.clientWidth;
    height = canvas.clientHeight;

    if (canvas.width !== width || canvas.height != height)
        renderer.setSize(width, height, false);
}

function populateAndPlaceDetails(state) {

    if (state == "SHOW") {
        // PLACEMENT
        // check if details are offscreen, then shift appropriately
        // X, add 10 pixels for buffer, since width is dynamic
        if (mouse.x + detailsOffsetX + nodeDetails.offsetWidth + 10 > canvas.clientWidth) {
            nodeDetails.style.left = (mouse.x - detailsOffsetX - nodeDetails.offsetWidth) + "px";
        } else {
            nodeDetails.style.left = (mouse.x + detailsOffsetX) + "px";
        }

        // Y
        if (mouse.y + detailsOffsetY + nodeDetails.offsetHeight > canvas.clientHeight) {
            nodeDetails.style.top = (mouse.y - detailsOffsetY - nodeDetails.offsetHeight) + "px";
        } else {
            nodeDetails.style.top = (mouse.y + detailsOffsetY) + "px";
        }

        nodeDetails.style.display = "block";
        //nodeDetails.style.opacity = .95;
        //nodeDetails.style.opacity = "1.0";

        //nodeDetails.style["borderColor"] = "#" + nodeHover.material.color.getHexString();

        // POPULATE

        // Organ name
        organName.innerHTML = nodeHover.name;

        // Dose Per Volume
        dosePerVolume.innerHTML = nodeHover.userData.dosePerVolume + " (GY/cc)";

        // line separator
        lineSeparator.style["borderColor"] = "#" + nodeHover.material.color.getHexString();

        // Volume
        volumeVal.innerHTML = nodeHover.userData.volume;

        // Mean Dose
        meanDoseVal.innerHTML = nodeHover.userData.meanDose;

        // Min Dose
        minDoseVal.innerHTML = nodeHover.userData.minDose;

        // Max Dose
        maxDoseVal.innerHTML = nodeHover.userData.maxDose;

    } else if (state == "HIDE") {

        nodeDetails.style.display = "none";
        //nodeDetails.style.opacity = .3;
        //nodeDetails.style.opacity = "0.0";

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

        cameraToCopy = scenes[targ.parentNode.id - 1].userData.camera;

        // 20 milliseconds interval => 50 FPS
        syncCamerasInterval = setInterval(syncAllCameras, 20, cameraToCopy);
    }
}

function syncAllCameras(cameraToCopy) {

    pRankingOrder.forEach(function (rank, index) {

        var scene = scenes[rank - 1];
        var camera = scene.userData.camera;
        var controls = scene.userData.controls;

        camera.position.subVectors(cameraToCopy.position, controls.target);
        camera.position.setLength(cameraDistZ);
        //camera.position.setLength(2000);
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

            currScene = scenes[targ.parentNode.id - 1];

            mouse.x = event.clientX;
            mouse.y = event.clientY;

            mouseNorm.x = (event.offsetX / targ.offsetWidth) * 2 - 1;
            mouseNorm.y = -(event.offsetY / targ.offsetHeight) * 2 + 1;
        }
    }
}
