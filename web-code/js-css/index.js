var listNode=[], listLink=[], /*rList, */geoList, ctryList, listVirus=[], listSite, listChange={}, outGrp='EPI_ISL_402131';
var wMap=504, stdR=4.5, cR1=4, mapR=3.5, cR2=8, distance=12//, lnkdist=1;
var annoDataW, geoDataW, countryData, map_ctry, refData,
    ids=[], annoDataC, geoDataC, pubData, source, symData, orfData, compData, seqData, dndData;
var wTree=100, tabW, unit=30, marginL=1, wYear=140,
    hUnit=19.5, yEdge=22, hTree, hSvg,
    orgColor = {1:"DodgerBlue", 2:"magenta", 3:"yellowGreen", 4:"orange"},
    geoColor, dateScale, shipColor='navy',
    ratio = 36.4, hORF=6, transBar, orfLine, sub_site=[],
    unitW=6.6, showAA;
var ratioC, cid_comp=[1,2,3,4,5,6,8,9,10], compLine, seqCid, tabNum, seqScroll;
var dbLink ={"pubmed":"http://www.ncbi.nlm.nih.gov/pubmed/","ncbi":"http://www.ncbi.nlm.nih.gov/nuccore/"};
var transTime=2000, timeT=500;

$(document).ready(function(){
    $("#mainTab, #tabs").tabs();
    $('#mainTab').tabs({ active:1 });
    $("#treeContainer").css("width", marginL*2 + wTree);
    
    $.ajax({
        async:false, dataType:"json", url:"js-css/anno.json",
        success: function(data) {
            annoDataW = data.annoW;
            geoDataW = data.geoW;
            countryData = data.country;
            map_ctry = data.map_geo;
            refData = data.ref;
            annoDataC = data.annoC;
            geoDataC = data.geoC;
            pubData = data.pub;
            source = data.source;
            symData = data.sym;

            $.ajax({
                async:false, dataType:"json", url:"js-css/net.json",
                success: function(data) {
                    readHaps(data);
                    asignColor();
                    drawChart();
                    drawGeoW();
                    drawAdmin();
                    drawDateChart();
                    drawRef()
                    drawLg();
                }
            })
        }
    });

    $(':radio').change(function(){changeColor(Number($(':radio:checked').val()))});

    $.ajax({
        async:false, dataType:"text", url:"js-css/tree.dnd",
        success: function(data){
            var obj = parseNewick(data);
            drawTree(obj);
            addIds();
            drawYear();
            drawGeoC()

            $.ajax({
                async:false, dataType:"json", url:"js-css/seq.json",
                success: function(data) {seqData = data} });
        }
    });
    $('#showWhat').click(function(){
        $(this).html(showAA? "Show Nt" : "Show AA");
        if (showAA){
            d3.select("#seqNt").transition().duration(transTime).style("opacity", 0);
            if (!$('#seqAA').length) { aaSeqTbl() }
            d3.select("#seqAA").transition().duration(transTime*1.25).delay(transTime/4).style("opacity", 1);
            showAA = 0
        } else {
            d3.select("#seqAA").transition().duration(transTime).style("opacity", 0)
            d3.select("#seqNt").transition().duration(transTime*1.25).delay(transTime/4).style("opacity", 1);
            showAA = 1
        }
    });

    var gvW = unitW*3-1;
    $("#guideV").css("top", 0+"px").css("width", gvW + "px");
    $("#guideH").css("height", (unitW*3-1) + "px");
	$("#seqContainer").mousemove(function(e) {
        var xMouse = e.pageX - $(this).offset().left + $("#seqContainer").scrollLeft(),
            yMouse = e.pageY - 23;
        if (xMouse>0 && xMouse < svgSeqW){$("#guideV").css("left",(xMouse-gvW)+"px")}
	   	if (yMouse > 127 && yMouse<580) { $("#guideH").css("top", (yMouse-15) + "px") }

        if (!seqCid){return}
        var px = xMouse/unitW/(ratioC+0.02),
            wSeqCont = $("#seqContainer").width();
        var x = px+seqScroll,
            scroll = px<wSeqCont/2 ? seqScroll : x-wSeqCont/2;
        $("#dnaComp").scrollLeft(scroll);
        compLine.attr("x1", x).attr("x2", x)
	});

    $("#dnaComp").mousemove(function(e){
        var x0 = e.pageX - $(this).offset().left,
            xMouse = x0 + $("#dnaComp").scrollLeft();
        compLine.attr("x1", xMouse).attr("x2", xMouse);
        var px;
        if (tabNum==1){
            px = xMouse*ratioC/ratio;
            orfLine.attr("x1",px).attr("x2",px)
        } else {
            px = (xMouse-seqScroll)*ratioC*unitW;
            $("#seqContainer").scrollLeft(px-x0+unitW*2.2);
            $("#guideV").css("left", px + "px")
        }
	});
    
	$("#genome").mousemove(function(e){
		var xMouse = e.pageX - $(this).offset().left;
	    orfLine.attr("x1",xMouse).attr("x2",xMouse);
        var px = xMouse*ratio/ratioC;
        $("#dnaComp").scrollLeft(px-xMouse);
        compLine.attr("x1", px).attr("x2", px)
    });

    $('#closeAck').on("click", function(){ $('#ack').hide(); return false});

var client = algoliasearch('0F0AD3F6TP', '5139b826dbcf021c70db5bb3be8abec1');
var index = client.initIndex('isolate');

    $('#searchAcc').autocomplete({hint: false}, [
        {source: $.fn.autocomplete.sources.hits(index, { hitsPerPage:8 }),
         displayKey: 'name',
         templates: {
            suggestion: function(suggestion){
                var res = suggestion._highlightResult;
    //            $('#suggestion').html(res.acc.value + ' ' +  res.name.value)
                return res.acc.value + ' | ' +  res.name.value
            }
          }
        }
      ]).on('autocomplete:selected', function(event, sug){
        currAcc = [sug.hap, sug.geo, [sug.acc]]
        overM()
      });
});

function asignColor(){
    var geo={};
    listVirus.forEach(function(d){
        var an = annoDataW[d];
        if (an.geoId<500) geo[an.geoId]=1
    });
    geoList = Object.keys(geo).map(d=>d*1);
    
    var ctr = {};
    geoList.forEach(function(d){ ctr[geoDataW[d].ctry]=1 });
    ctryList = Object.keys(ctr).map(d=>d*1);
    
    var l = ctryList.length+5;
    geoColor = d3.scaleOrdinal().domain(ctryList).range(d3.range(1,l).map(function(d){return d3.hsl(360/l * d, 1, 0.45)}));
    geoList.push(501,502,503)
    ctryList.push(701)
}

var chartW=655, node_gid=[];
var currAcc=[]
function drawChart(){
    var width = chartW, height=592,
        linewidth = 1;

    var pie = d3.pie().sort(null).value(function(d){return 1});

    var svg = d3.select("#netGraph").attr("width", width).attr("height", height);
    svg.append("rect").attr("width", width).attr("height", height)
                    .attr("class", 'whiteFill')
                    .on('click',function(d){ ulNodes(); hideNodeInf(); ulLink(); ulGeo() })


    for (var i=0; i<listLink.length; i++){
        var lk = listLink[i];
        svg .append('defs')
            .append('marker')
            .attr('id', 'arrow'+ i)
            .attr('viewBox', [0, -5, 10, 10])
            .attr('refX', lk.target.radius*0.6+5)
            .attr('refY', -0.5)
            .attr('markerWidth', 8)
            .attr('markerHeight', 7)
            .attr('orient', 'auto')
            .append('path')
            .attr("d", "M0,-3L7,0L0,3")
    
        if (!lk.change) continue;
        var cc = lk.change;
        for (var j=0; j<cc.length; j++){
            var ch = cc[j];
            if (!ch.src_aa || ch.src_aa==ch.tgt_aa) continue;
            var ll=[];
            if (listChange[ch.label]) ll = listChange[ch.label];
            ll.push(i)
            listChange[ch.label] = ll
        }
    }

    var chrg=-100;
    var simulation = d3.forceSimulation()
        .force("link", d3.forceLink()
               .id(function(d){ return d.index })
               .distance(function(d){ return d.source.haps? d.ldist : d.ldist*2.4})
               .iterations(35))
        .force("charge", d3.forceManyBody().strength(chrg))
        .force("center", d3.forceCenter(width/2, height/2))
        .force("y", d3.forceY(0))
        .force("x", d3.forceX(0));

    var link = svg.selectAll(".link").data(listLink).enter()
        .append("line").attr('class', function(d){
                return 'link ' + (d.change? 'finger': 'dashStroke')
            })
        .attr("id", function(d,i){return 'll_'+i})
        .attr('marker-end', function(d,i){return 'url(#arrow' + i + ')'})
        .on('click',function(d,i){
            if (!d.change){return}
            ulLink();

            $('#chgTbl td').remove();
            var trs = [];
            d.change.forEach(function(c){
                var noSyn = c.src_aa && c.src_aa != c.tgt_aa;
                d3.select('#site_'+c.site)
                    .style("stroke", noSyn? "red" : "DodgerBlue")
                var tds = [];
                tds.push(c.site, c.label, c.pos);
                if (c.src_aa){
                    tds.push(Math.ceil(c.pos/3), hlCodonP(c.src_codon,c.cd_pos, noSyn) + '(' + c.src_aa + ')', hlCodonP(c.tgt_codon,c.cd_pos, noSyn) + '(' + c.tgt_aa + ')')
                } else {
                    tds.push('', c.src_base, c.tgt_base)
                }
                trs.push('<tr><td>' + tds.join('</td><td>') + '</td></tr>')
            });
            $('#chgTbl').append(trs.join('')).show()

            d3.select('#ll_'+i).classed("orangeStroke",true);
            d3.select('#arrow'+i+' path').classed("orangeFill",true);

            var sel = d.change.map(x=>'#site_'+x.site).join(',');
            d3.selectAll(sel).classed("hidden",false)
        })

    var node = svg.selectAll(".node").data(listNode).enter()
        .append("g")
        .attr("class", "node")
        .attr('id', function(d){ return d.id})
        .call(d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended))

    simulation.nodes(listNode).on("tick", function(){
        link.attr("x1", function(d) { return d.source.x })
            .attr("y1", function(d) { return d.source.y })
            .attr("x2", function(d) { return d.target.x })
            .attr("y2", function(d) { return d.target.y });
        node.attr("x", function(d){return d.x= Math.max(d.radius, Math.min(width - d.radius, d.x))})
            .attr("y", function(d){return d.y= Math.max(d.radius, Math.min(height - d.radius, d.y))})
            .attr("transform", function(d){return "translate(" + d.x + "," + d.y + ")"})
        });
    simulation.force("link").links(listLink);

    var nLab=10;
    listNode.forEach(function(d){
        var n = d3.select('#'+d.id);
        n.append("circle").attr("r", d.radius).style('fill',"white")

        if (d.haps && d.haps.length>1){
            var sector = d3.arc().innerRadius(0).outerRadius(d.radius);
            
            n.selectAll('path')
                .data(pie(d.haps)).enter()
                .append('path')
                .attr("class", function(n){return 'nn_' + n.data[0] + " iso crosshair"})
                .attr('d', sector)
                .style('fill', function(n){
                    var gid = n.data[0];
                    return gid<500? geoColor(geoDataW[gid].ctry) : shipColor
                })
                .on('click',function(n){
                    $('#searchAcc').val('')
                    currAcc = [d.id, n.data[0], n.data[1]]
                    overM()
                })
                .on("mouseover", function(n){
                    var num =  n.data[1].length
                    tipNet(num + ' pateint' + (num==1? '' : 's'))
                })
                .on('mouseout',function(){$("#tipNet").hide()})
            
            var center = n.append("circle").attr("r", (d.haps.length>nLab? 7 : 2)).style('fill',"white")
            if (d.haps.length>nLab){
                n.append("text").attr("class","nodeId").attr("dy",3).text(d.id)
                    .on("mouseover", function(){tipNet(d.id, 20)})
                    .on('mouseout',function(){$("#tipNet").hide()})
            } else {
                center.style("fill-opacity",0.1).style("stroke","none")
                    .on("mouseover", function(){tipNet(d.id, 20)})
                    .on('mouseout',function(){$("#tipNet").hide()})
            }
        } else {
            if (!d.haps){
                n.append("circle").attr("r", d.radius)
                    .style('fill',"white").style('stroke','#999').style('stroke-width','1px')
                    .on('mouseover',function(d){tipNet('Hypothetical ' + (d.id=='H0'? 'root' : 'haplotype'))})
                    .on('mouseout',function(){$("#tipNet").hide()})
            } else {
                var sl = d.haps[0],
                    gid = sl[0],
                    ac = sl[1][0],
                    an = annoDataW[ac]
                var circle;
                if (ac!=outGrp){
                    circle = n.append("circle")
                        .attr("r", d.radius)
                        .style('fill', gid<500? geoColor(geoDataW[gid].ctry) : shipColor);
//                n.append("circle").attr("r", 8).style('fill',"white")
//                n.append("text").attr("class","nodeId").attr("dy",3).text(d.id.replace(/ST/,'H'))
                } else {
                    $("#bat").appendTo(('#'+d.id));
                    circle = n.append("svg")
                            .style("fill","gray")
                            .attr("x",-40).attr("y",-7)
                            .attr("width",48).attr("height",23.996)
                            .attr("viewBox","562 485.829 50 23.996")
                            .append("path").attr("d",batPath)
/*                    circle = n.append("text").attr("id", 'nn_'+ac).attr("class",'bat finger')
                        .attr("dy","0.4em").attr("dx","10px").html('&#129415;')*/
                }

                circle.attr("class",'nn_' + gid + " iso crosshair")
                    .on('click',function(){
                    $('#searchAcc').val('')
                    currAcc = [d.id, gid, [ac]]
                    overM()
                })
                .on("mouseover", function(n){if(n.haps[0][1][0]==outGrp) return; tipNet('1 pateint (' + d.id + ')')})
                .on('mouseout',function(){$("#tipNet").hide()})
            }
        }
    })

    function hlCodonP(ch,p,n){
        var arr = ch.split('');
        arr[p-1] = '<b class="' + (n? "redFill":"blueColor") + '">' + arr[p-1] + '</u></b>'
        return arr.join('')
    }
    
    function dragstarted(d){
        if (!d3.event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }
    function dragged(d) {
        d.fx = d3.event.x;
        d.fy = d3.event.y;
    }
    function dragended(d) {
        if (!d3.event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }      
}

function hlNodes(sel){
    if (!Number($(':radio:checked').val())){
        d3.selectAll('.iso').transition().duration(timeT).style('fill-opacity',0.1)
        sel.transition().duration(timeT).style('fill-opacity',1)
    } else {
        d3.selectAll('.iso').transition().duration(timeT).style("fill","#999")
        sel.transition().duration(timeT).style('fill','purple')
    }
}

function overM(){
    var accs = currAcc[2],
        gid = currAcc[1]

    var sel = d3.select("#"+currAcc[0] + ' .nn_' + currAcc[1])
    hlNodes(sel)

    var rg=[]
    if (accs.length>1){
        rg = d3.extent(accs, d=>tParser(annoDataW[d].col));
    } else {
        var ac = accs[0]
        if (ac!=outGrp) rg = [tParser(annoDataW[ac].col)]
        
    }
    if (rg[0]) hlDate(rg)
    
    hlGeo(gid)

    $("#nPatient").html(rg[0]? (accs.length + ' patient' + (accs.length==1? '' : 's')) : 'Outgroup').show();
    var tbls=[]
    for (var i=0; i<accs.length; i++){
        var ac = accs[i],
            an = annoDataW[ac],
            lab = an.author + '<em> et al</em><br>Submit: ' + an.subdate + '<br>' + an.lab

        var host, patient = '';
        if (ac==outGrp){
            host = 'Host'
            patient = an.patient
            d3.select("#" + currAcc[0] + " .nn_"+gid).style("fill", geoColor(geoDataW[gid].ctry))
        } else if (/Environment/i.test(an.patient)){
            host = 'Environment';
            patient = an.patient.split(":")[1]
        } else {
            host = 'Patient';
            if (an.patient){patient = an.patient}
        }

        tbls.push('<table>' + '<tr><th>Name</th><td>' + an.iso + '</td><th>Accession</th><td>' + ac + '</td></tr><tr><th>Collect</th><td>' + an.col + '</td><th>' + (an.city? 'City' : '') + '</th><td>' + (an.city? an.city : '') + '</td></tr><tr><th>' + host + '</th><td colspan="3">' + patient + '</td></tr><tr><th>Citation</th><td colspan="3">' + lab + '</td></tr></table>')
    }
    $('#nodeInf').html(tbls.join('<hr>')).scrollTop(0)
}

function hlGeo(geoId){
    d3.selectAll('#mapW circle:not(#mapW_'+geoId+')').transition().duration(timeT)
        .attr("r", mapR/transK).style("fill-opacity",0.1).style("stroke","none");
    d3.select('#mapW_'+geoId).transition().duration(timeT)
        .attr("r", stdR*2/transK).style("fill-opacity",0.8);
    $("#siteName").html(geoDataW[geoId].name).show()
}
function hlGeos(gids){
    var arr = gids.map(function(d){return '#mapW_'+d}),
        sel = d3.selectAll(arr.join(','))
    d3.selectAll('#mapW circle').transition().duration(timeT)
        .attr("r", mapR/transK).style("fill-opacity",0).style("stroke","none");
    sel.transition().duration(timeT).attr("r", mapR/transK).style("fill-opacity",0.8).style("stroke","white");
}

function ulGeo(){
    d3.selectAll('#mapW circle').transition().duration(timeT)
        .attr("r", mapR/transK).style("fill-opacity",0.8).style("stroke","white");
    $("#siteName").hide()
}

function hlDate(arr){
    d3.select('#dateRg').classed("hidden",false)
        .transition().duration(transTime)
        .attr("x", ddScale(arr[0])).attr("width",arr[1]? (ddScale(arr[0])==ddScale(arr[1])? 1 : ddScale(arr[1])-ddScale(arr[0])) : 1);
}

function ulLink(){
    d3.selectAll('#netGraph .link').classed("orangeStroke",false)//.style("stroke","#999");
    d3.selectAll('#netGraph marker path').classed("orangeFill",false)//.style("fill","#999")
    $("#chgTbl").hide()
    d3.selectAll('.siteLine').classed("hidden",true)
}

function hideNodeInf(){
    if (!currAcc[0]) return
    $("#nodeInf table, #nodeInf hr, #nPatient").hide()
    d3.select('#dateRg').classed("hidden",true)
    $('#searchAcc').val('')
    currAcc=[]
}
function geo2node(gids){
    var arr = gids.map(function(d){return '.nn_'+d}),
        sel = d3.selectAll(arr.join(','));
    hlNodes(sel)
    hideNodeInf()
}

function ulNodes(){
    if (!Number($(':radio:checked').val())){
        d3.selectAll('.iso').transition().duration(timeT).style("fill-opacity",1)
        d3.select('#'+outId + ' path').style("fill","gray")
    } else {
        d3.selectAll('.iso').transition().duration(timeT).style("fill","purple")
    }
}

function drawLg(){
    var width=wMap, height=24,
        margin_left=7, margin_top=2,
        radius=9, wpic=radius*2+3;
    var svg = d3.select("#graphLg").append("svg").attr("width", width).attr("height", height),
        chart = svg.append("g").attr("transform", "translate("+[margin_left, margin_top]+")"),
        pic = chart.append("g").attr("class", "blackLine whiteFill")
                .attr("transform", "translate("+[radius, radius]+")"),
        txt1 = chart.append("g").attr("class","lgTitle").attr("transform", "translate("+[0, 9]+")"),
        txt2 = chart.append("g").attr("class","lgExp").attr("transform", "translate("+[0, 22]+")"),
        txt3 = chart.append("g").attr("class","lgExp grayFill").attr("transform", "translate("+[0,9]+")")
 
    pic.append("circle").attr("r", radius)
        .on("mouseover", function(){tipNet('Node size proportional to geographic range')})
        .on('mouseout',function(){$("#tipNet").hide()})
    txt1.append("text").text("Haplotype").attr("x",wpic)
    txt2.append("text").html("a group of similar genomes").attr("x",wpic);
    txt3.append("text").html("N=212").attr("x",wpic+56)

    ww = 152
    pic.append("circle").attr("r", radius).attr("cx", ww).style("stroke","#ccc")
    
    var pie = d3.pie().startAngle(0.5*Math.PI).endAngle(0.7*Math.PI),
        arc = d3.arc().innerRadius(0).outerRadius(radius);
    pic.append("g").attr("transform", "translate("+[ww, 0]+")")
        .selectAll("g").data(pie([10])).enter()
        .append("path").attr("d", arc)
        .on("mouseover", function(){tipNet('One or more isolates from one location')})
        .on('mouseout',function(){$("#tipNet").hide()})
    txt1.append("text").text("Isolate(s)").attr("x",ww+wpic);
    txt2.append("text").html("from the same location").attr("x",ww+wpic);
    txt3.append("text").html('N=2334 from <a href="https://www.gisaid.org" target="_blank">GISAID</a>').attr("x",ww+wpic+52)
    
    ww += 173//249
    svg .append('defs').append('marker').attr('id', 'arrow')
        .attr('viewBox', [0, -5, 10, 10])
        .attr('refX', 7).attr('refY', -0.5).attr('markerWidth', 7).attr('markerHeight', 7)
        .attr('orient', 'auto').append('path').attr("d", "M0,-3L7,0L0,3");
    pic.append("line").attr("class","link")
        .attr("transform", "translate("+[ww-radius+3, 0]+")")
        .attr("y1",radius-3).attr("x2",2*radius-6).attr("y2", 3-radius)
        .attr('marker-end', 'url(#arrow)')
        .on("mouseover", function(){tipNet('Changes between two haplotypes')})
        .on('mouseout',function(){$("#tipNet").hide()})

    txt1.append("text").attr("x",ww+wpic-3).text("Mutation(s)");
    txt2.append("text").attr("x",ww+wpic-3).html("genetic changes");
    txt3.append("text").html("at 146 genome sites").attr("x",ww+wpic+61)
}

function tipNet(txt,x){
    $("#tipNet").css("left", (d3.event.pageX-(x? x : 50))+"px").css("top", (d3.event.pageY-77)+"px").html(txt).show()
}

function drawRef(){
    var width = wMap-6,
        height=52,
        margin=2, margin_top=12,
        chartW = width-2*margin,
        base = 16,
        halfH = 5.5;
    var seqL = refData.len;
    var svg = d3.select("#refOrf svg").attr("width", width).attr("height", height),
        chart = svg.append("g")
            .attr("width", chartW).attr("height", height)
            .attr("transform", "translate("+[margin, margin_top]+")");

    svg.append("text").attr("class","infoTitle").attr("y",13).html('&nbsp; Genome map & mutations');

    var scale = d3.scaleLinear().domain([1, seqL]).range([0, chartW]);
    
    var xAxis = d3.axisBottom().scale(scale).tickSize(-3, 0).ticks(seqL/1000)
                .tickFormat(function(f){ return f<2000? '1k' : f/1000 });
    chart.append("g").attr("class", "xaxis").attr("transform", "translate(0," + (base+9.5) + ")").call(xAxis);

    var orfChart = chart.append("g").attr("transform", "translate("+[0, base]+")");
    orfChart.append("line").attr("id","refLine").attr("x2", scale(seqL));

    var orfs = refData.orf;
    orfChart.append("g").attr("class","finger")
            .selectAll("rect").data(orfs).enter()
            .append("rect")
            .attr("x", function(d){return scale(d[0])})
            .attr("y", -halfH)
            .attr("width", function(d){return scale(d[1]) - scale(d[0])})
            .attr("height", halfH*2)
            .on('mouseover',function(d){
                d3.select(this).style("fill", "orange");
                ulLink()

                var gene = d[2];
                var lks = listChange[gene],
                    tipTxt = '<em>' + gene + (lks? '<br><span>red arrows - presence of<br>non-synonymous mutations</span>' : '') + '</em>'
                $("#tipNet").css("left", (d3.event.pageX-(lks?120:16))+"px").css("top", (d3.event.pageY-(lks?103:73))+"px").html(tipTxt).show()

                if (!lks){return}
                var ll = lks.map(function(l){return '#ll_'+l}).join(','),
                    ar = lks.map(function(l){return '#arrow'+l+' path'}).join(',');
                if (!Number($(':radio:checked').val())){
                    d3.selectAll('.iso').transition().duration(timeT*2).style('opacity',0.15)
                } else {
                    d3.selectAll('.iso').transition().duration(timeT*2).style('fill',"#ccc")
                }
                d3.selectAll('#netGraph .link').style("stroke","#eee")
                d3.selectAll(ll).style('stroke','red')
                d3.selectAll('#netGraph marker path').style("fill","#eee")
                d3.selectAll(ar).style("fill",'red')
            })
            .on('mouseout', function(){
                 d3.select(this).style("fill", "#ffe4b2")
                $("#tipNet").hide()
                if (!Number($(':radio:checked').val())){
                    d3.selectAll('.iso').transition().duration(timeT*2).style('opacity',1)
                } else {
                    d3.selectAll('.iso').transition().duration(timeT*2).style('fill',"purple")
                }
                d3.selectAll('#netGraph .link').style("stroke","#999");
                d3.selectAll('#netGraph marker path').style("fill","#999")
            })

    var toShow = {"orf1ab":1, "S":1}
    orfChart.append("g").attr("class","refSym")
            .selectAll(".refSym").data(orfs).enter()
            .append("text")
            .attr("x", function(d){return (scale(d[0])+scale(d[1]))/2})
            .attr("y", halfH-2)
            .text(function(d){return toShow[d[2]]? d[2] : ''})
    
    
//    base -= 9+halfH;
    chart.append("g")//.attr("transform", "translate("+[0, base]+")")
            .selectAll(".siteLine").data(listSite).enter()
            .append("line")
            .attr("class","siteLine hidden")
            .attr("id", function(d){return "site_"+d })
            .attr("x1", function(d){ return scale(d)})
            .attr("x2", function(d){return scale(d)})
            .attr("y1", 6)
            .attr("y2", 30)
}


var transK=1, wAdmin=95;
function drawGeoW(){
    var width = wMap - wAdmin,
        height = $('#mapW').height() - 45;
    var projection = d3.geoMercator().translate([width/2, height/2]).scale(70).center([18,22]),
        mapPath = d3.geoPath().projection(projection);

    var active = d3.select(null);
    var zoom = d3.zoom().scaleExtent([1,15]).on("zoom", zoomed);

    var svg = d3.select("#mapW svg").attr("width",width).attr("height",height);

    var svgMap = svg.append("g");
    svg.call(zoom);

    svgMap.append("rect").attr("class","mapBg").attr("width",width).attr("height",height).on("click",reset);

    d3.json("js-css/110m.json", function(error, world){
        var countries = topojson.feature(world, world.objects.countries).features;
        svgMap.selectAll(".country")
            .data(countries).enter()
            .append("path")
            .attr("class", "country")
            .attr("id",function(d){return 'ctry_'+map_ctry[d.id*1]})
            .attr("d", mapPath)
            .on('mouseover', function(d){
                var id = map_ctry[d.id*1]
                if (id){
                    $("#admin").scrollTop(ctry.indexOf(id) * $('#admin td').height() - $('#admin').height()/2);
                    hlCtry(id)
                } else {
                     d3.select(this).classed("country_sel",true)
                }
            })
            .on('mouseout', function(d){
                if (map_ctry[d.id*1]){ulCtry()}
                else {d3.select(this).classed("country_sel",false)}
            })
            .on("click", clicked);
 
        svgMap.selectAll("circle")
            .data(geoList).enter()
            .append("circle")
//            .attr("class","crosshair")
            .attr("id", function(d){ return 'mapW_'+d })
            .attr("fill", function(d){return d<500? geoColor(geoDataW[d].ctry) : shipColor})
//            .style("fill-opacity",0.8)
            .attr("r", mapR)
            .attr("cx", function(d){
                var site = geoDataW[d].locate,
                    coords = projection([site[1], site[0]]);
                return coords[0]
            })
            .attr("cy", function(d){
                var site = geoDataW[d].locate,
                    coords = projection([site[1], site[0]]);
                return coords[1]
            })
            .on('mouseover', function(d){
                hlGeo(d)
                geo2node([d])
            })
            .on('mouseout', function(){
                ulGeo()
                ulNodes()
        })
    })

    function clicked(d) {
//        if (active.node() === this) return reset();
        active = d3.select(this).classed("active", true);

        var bounds = mapPath.bounds(d),
            dx = bounds[1][0] - bounds[0][0],
            dy = bounds[1][1] - bounds[0][1],
            x = (bounds[0][0] + bounds[1][0]) / 2,
            y = (bounds[0][1] + bounds[1][1]) / 2,
            scale = Math.max(1, Math.min(8, 0.9 / Math.max(dx/width, dy/height))),
            translate = [width/2 - scale*x, height/2 - scale*y];

        svgMap.transition()
            .duration(750)
            .call(zoom.transform, d3.zoomIdentity.translate(translate[0],translate[1]).scale(scale));
    }

    function reset(){
        active.classed("active", false);
        active = d3.select(null);
        transK=1;
        svgMap.transition().duration(1000).call(zoom.transform, d3.zoomIdentity)
    }

    function zoomed() {
        var trans = d3.event.transform;
        transK=trans.k;
        d3.selectAll('.country, circle').style("stroke-width", 1/transK + "px");
        svgMap.attr("transform", trans);
        d3.selectAll('#mapW circle').attr("r",mapR/transK)
    }
}

var ctry;
function drawAdmin(){
    ctry = ctryList.sort((a,b) => countryData[a] > countryData[b] ? 1 : -1);
    
    var trs = ctry.map(function(d){
        return '<tr><td id="adm_'+ d + '">' + countryData[d] + '</td></tr>'
    })
    $('#admin table').append(trs.join(''))
 
    d3.selectAll("#admin td")
        .on('mouseover', function(){
            var id = this.id.split('_')[1];
            hlCtry(id)
        })
        .on('mouseout', ulCtry)
}

function hlCtry(id){
    ulGeo()
    $('#adm_'+id).css("background", id<700? geoColor(id) : shipColor).css("color", "white")
    d3.select('#ctry_'+id).classed('country_sel',true)

    var gids = geoList.filter(function(d){return geoDataW[d].ctry==id})
    hlGeos(gids)
    geo2node(gids)
}

function ulCtry(){
    $('#admin td').css("background", "none").css("color", "inherit")
    d3.selectAll('#mapW .country').classed("country_sel",false)
    d3.selectAll('#mapW circle').transition().duration(timeT)
        .style("fill-opacity",0.8).style("stroke","white");
    ulNodes()
}

var tParser = d3.timeParse("%Y-%m-%d"), ddScale;
function drawDateChart(){
    var width = wMap, height=56;
    var margin = {top:43, left:0, bottom:0, right:0},
        chartWidth = width - (margin.left+margin.right),
        chartHeight = height - (margin.top+margin.bottom);

    var dateObj={};
    listVirus.forEach(function(ac){
        var an = annoDataW[ac];
        if (ac!=outGrp) dateObj[an.col]=1
    });

    var dateArr = Object.keys(dateObj),
        rg = d3.extent(dateArr, d=>tParser(d));
    ddScale = d3.scaleLinear().domain(rg).range([stdR*2,chartWidth-stdR*2]);

    var svg = d3.select("#dateChart").append("svg").attr("width",width).attr("height",height),
        chart = svg.append("g")
            .attr("width", chartWidth).attr("height", chartHeight)
            .attr("transform", "translate("+[margin.left, margin.top]+")");

    var opRange = [1,0.05];
    dateScale = d3.scaleLinear().domain(rg).range(opRange);

    var linearGrad = chart.append("defs").append("linearGradient").attr("id", "dateGrad");
    linearGrad.append("stop").attr("offset",  "0%").attr("stop-color","purple").attr("stop-opacity", opRange[0]);
    linearGrad.append("stop").attr("offset","100%").attr("stop-color","white").attr("stop-opacity", opRange[1]);

    chart.append("rect").attr("id","dateGradient")
        .attr("width", ddScale(rg[1]) - ddScale(rg[0])).attr("height", 8)
        .attr("x",ddScale(rg[0])).attr("y",-5)
        .attr("class","grayLine")
        .style("fill", "url(#dateGrad)")
        .style("opacity", "url(#dateGrad)")

    chart.append("rect").attr("id","dateRg").attr("class","hidden grayFill")
        .attr("y",3).attr("height",8).attr("width", 200)
        
    var xAxis = d3.axisTop().scale(ddScale).ticks(6).tickSize(-12, 0).tickFormat(d3.timeFormat("%m/%d/%y"));
    chart.append("g").attr("class", "xaxis")
        .attr("transform", "translate(" + 0 + "," + (-8) + ")")
        .call(xAxis);
    svg.append("text").attr("class","infoTitle").attr("y",14).html("&nbsp; Collection Date");
}

function changeColor(val){
    for (var i=0; i<node_gid.length; i++){
        var oo = node_gid[i],
            gid = oo[1],
            sel = d3.select('#'+oo[0] + ' .nn_'+gid);
        if (val){
            sel.style("fill-opacity", dateScale(oo[2]))
        } else {
            var ctr = geoDataW[gid].ctry;
            sel.style("fill", oo[0]==outId? 'gray' : (ctr<700? geoColor(ctr) : shipColor))
        }
    }
    if (val){
        d3.selectAll('.iso').transition().duration(timeT).style("fill","purple").style("stroke","none");
    } else {
        d3.selectAll('.iso').transition().duration(timeT).style("stroke","white").style("fill-opacity",1)
    }
    if (currAcc[0]){hlNodes(d3.select("#"+currAcc[0] + ' .nn_' + currAcc[1]))}
}

function drawTree(obj){
    var root = d3.hierarchy(obj, function(d){return d.children});
    root.leaves().forEach(function(d){ids.push(d.data.name)});
    hTree = ids.length*hUnit;
    hSvg = hTree + yEdge;

    var cluster = d3.cluster().size([hTree, 0]).separation(function(a, b) { return 1 });
    cluster(root);

    resetY(root, root.data.length = 0, wTree/maxLength(root));

    var svg=d3.select("#treeContainer svg")
            .attr("width", marginL*2 + wTree)
            .attr("height", hTree)
            .append("g").attr("id","tree").attr("transform", "translate(" + marginL + ",0)");

    svg.append("g")
        .attr("class", "blackLine")
        .selectAll("path").data(root.links()).enter()
        .append("path")
        .attr("d", link)

    svg.append("g")
        .attr("class", "link-extensions")
        .selectAll("path").data(root.links().filter(function(d){ return !d.target.children})).enter()
        .append("path")
        .attr("d", linkExtension);

//draw scale:
    var svgTreeScale = d3.select("#treeScale").append("svg").attr("width",130).attr("height",59),
        treeScale = svgTreeScale.append("g").attr("transform", "translate(5,45)"),
        scaleLine = treeScale.append("g").attr("class",'grayLine');
    scaleLine.append("line").attr("x2", unit);
    scaleLine.append("line").attr("y1", -2).attr("y2", 2);
    scaleLine.append("line").attr("x1", unit).attr("x2", unit).attr("y1", -2).attr("y2", 2);

    sub_site[0] = d3.format(".2g")(unit*maxLength(root)/wTree);
    treeScale.append("text").attr("x",unit+5).attr("y",4)
        .html('<tspan></tspan> sub/site');
    $('#treeScale tspan').html(sub_site[0])
}

function drawTreeC(cid){
    var obj = parseNewick(dndData[cid]),
        root = d3.hierarchy(obj, function(d){return d.children});
    id=[];
    root.leaves().forEach(function(d){id.push(d.data.name)});
    var hTr = id.length*hUnit,
        hSv = hTr + yEdge;

    var cluster = d3.cluster().size([hTr, 0]).separation(function(a, b) { return 1 });
    cluster(root);

    resetY(root, root.data.length = 0, wTree/maxLength(root));

    var svg=d3.select("#treeContainer svg").append("g")
        .attr("id","tree_cid").attr("transform", "translate(" + marginL + ",0)");

    svg.append("g")
        .attr("class", "blackLine")
        .selectAll("path").data(root.links()).enter()
        .append("path")
        .attr("d", link)

    svg.append("g")
        .attr("class", "link-extensions")
        .selectAll("path").data(root.links().filter(function(d){ return !d.target.children})).enter()
        .append("path")
        .attr("d", linkExtension);

    sub_site[1] = d3.format(".2g")(unit*maxLength(root)/wTree);
    $('#treeScale tspan').html(sub_site[1])
}

function mOver(d){
    var an = annoDataC[d];
    d3.selectAll('#o_'+d + ', #lgDot'+an.orgId).attr("r", cR2+2);
    d3.select('#map'+an.geoId+'_'+an.orgId).attr("r", cR2*2);
    $('#id_'+d + ', #co_'+d).css("background", orgColor[an.orgId]).css("color", "white")
}
function mOut(d){
    d3.selectAll('#col_year circle, #lg_org circle').attr("r", cR1);
    d3.selectAll('#mapC circle').attr("r", cR2);
    $('#ids td, #collect td').css("background", "inherit");
    $("#collect td").css("color", "black");
    if (d.length==3){
        var acc = ids.filter(function(ac){ var an=annoDataC[ac]; return an.geoId==d[0] && an.orgId==d[1]});
        acc.forEach(function(ac){ $('#id_'+ac).css("color", orgColor[annoDataC[ac].orgId]) })
    }
    else {
        $('#id_'+d).css("color", orgColor[annoDataC[d].orgId])
    }
    $('#affiliation').html('')
}

var tdId={}, tdCol={};
function addIds(){
    var tdPh=[];
    ids.forEach(function(d){
        var an = annoDataC[d],
            pId = an.pubId,
            pu = pubData[pId];
        var ii = '<td id="id_' + d + '" style="color:' + orgColor[an.orgId] + '"><a href="' + dbLink.ncbi + d + '" target="_blank">' + an.iso + '</a></td>',
            cc = '<td id="co_' + d + '">' + an.col + ' | ' + geoDataC[an.geoId].name + '</td>',
            pp = '<td id="ph_' + d + '">' + (pId>1000? '<a href="' + dbLink.pubmed + pId + '" target="_blank">' : '') + pu[0] + ' <em>et al</em>. '  + pu[1] + (pId>1000? '</a>' : '') + '</td>';
        var tr1 = '<tr class="p' + an.geoId + 'o' + an.orgId + '">',
            tr2 = '</tr>';

        tdId[d] = tr1 + ii + tr2;
        tdCol[d] = tr1 + cc + tr2;
        tdPh.push('<tr>' + pp + tr2);
    });
    
    appendPhe(ids);
    $('#citation').append(tdPh.join(''));
    $('.pheno th').css("height", (yEdge-2)+'px');
    addEv()

    d3.csv("js-css/orf.csv", function(data){ orfData = data; drawORF()});
}

function drawYear(){
    var year = ids.map(function(d){return annoDataC[d].col}),
        minYear = d3.min(year) - 1,
        scaleYear = d3.scaleLinear().domain([minYear, d3.max(year)]).range([marginL+4, wYear+marginL]),
        start = scaleYear(minYear);
    
    var svg = d3.select("#col_year").append("svg")
                .attr("width", 14 + wYear)
                .attr("height", hSvg+17);

    var xAxis = d3.axisBottom().scale(scaleYear).tickSize(-hTree-5, 0)
                .tickFormat(function(f){ var y=f-2000; return y<10? '0'+y : y });
    svg.append("g").attr("class", "xaxis").attr("transform", "translate(" + 0 + "," + (hSvg+6) + ")").call(xAxis);

    svg.append("text").attr("x",27).attr("y",13).attr("class","lg").text("Collect Year")

    var lineFun = d3.line()
        .x(function(d) { return scaleYear(annoDataC[d].col); })
        .y(function(d) { return hUnit*(ids.indexOf(d)+0.5); })
        .curve(d3.curveLinear);

    svg.append("path")
        .attr("transform", "translate(0," + yEdge + ")")
        .attr("d", lineFun(ids))
        .attr("class", "blackLine")

    svg.append("g")
        .attr("transform", "translate(0," + yEdge + ")")
        .selectAll("circle").data(ids).enter()
        .append("circle")
        .attr("id", function(d){return 'o_' + d})
        .attr("class", function(d){var obj = annoDataC[d]; return 'g' + obj.geoId + 'o' + obj.orgId})
        .attr("transform", function(d,i){ return "translate(" + scaleYear(annoDataC[d].col) + "," + hUnit*(i+0.5) + ")" })
        .attr("r", cR1)
        .attr("fill", function(d){return orgColor[annoDataC[d].orgId]})
        .on("mouseover", mOver)
        .on("mouseout", mOut)
}

function drawGeoC(){
    var hMap=333,
        projection = d3.geoMercator().translate([wMap/2, hMap/2]).scale(420).center([104.2,38]),
        mapPath = d3.geoPath().projection(projection);

    var svg = d3.select("#mapC").append("svg").attr("width",wMap).attr("height",hMap);

    d3.json("js-css/china.geo.json", function(error, china){
        svg.selectAll(".country")
            .data(china.features).enter()
            .append("path")
            .attr("class", "country")
            .attr("d", mapPath)
            .on('mouseover', function(d){d3.select(this).classed("country_sel",true)})
            .on('mouseout', function(d){ d3.select(this).classed("country_sel",false)});

        var geo_org = {};
        ids.forEach(function(d){
            var an = annoDataC[d];
            var obj = {};
            if (geo_org[an.geoId]){ obj = geo_org[an.geoId] }
            obj[an.orgId] = 1;
            geo_org[an.geoId]=obj
        });

        var g_o_n = [];
        Object.keys(geo_org).forEach(function(d){
            var obj = geo_org[d],
                nG = Object.keys(obj).length,
                offset = 0.45;
            Object.keys(obj).forEach(function(g,i){ g_o_n.push([d*1,g*1, nG==1? 0 : (i? offset : -offset)]) })
        });
        
        svg.selectAll("circle")
            .data(g_o_n).enter()
            .append("circle")
            .attr("id", function(d){ return 'map'+d[0]+'_'+d[1]})
            .attr("fill", function(d){return orgColor[d[1]]})
            .attr("r", function(d){return cR2})
            .attr("cx", function(d){
                var site = geoDataC[d[0]].locate,
                    coords = projection([site[1] + d[2], site[0]]);
                return coords[0]
            })
            .attr("cy", function(d){
                var site = geoDataC[d[0]].locate,
                    coords = projection([site[1], site[0]]);
                return coords[1]
            })
            .on('mouseover', geoTip)
            .on('mouseout', mOut);

        var geoIds = ids.map(function(d){return annoDataC[d].geoId}),
            geoId = $.unique(geoIds.sort(function(a,b){return b-a}));
        svg.selectAll("text")
            .data(geoId).enter()
            .append("text")
            .attr("x", function(d){
                var site = geoDataC[d].locate,
                    coords = projection([site[1], site[0]]);
                return coords[0]
            })
            .attr("y", function(d){
                var site = geoDataC[d].locate,
                    coords = projection([site[1], site[0] + (d==71 || d==76? -2.4 : 1.4)]);
                return coords[1]
            })
            .text(function(d){return geoDataC[d].name})
    });

// draw legend
    var svgLg=d3.select("#lg_org").append("svg")
            .attr("width", wMap).attr("height", hUnit*6-5)
            .append("g").attr("transform", "translate(130,10)");
    svgLg.append("text").attr("x",15).attr("class","lg").text("Organism")
    var lgOrg = svgLg.append("g").attr("class","lgOrg")
        .attr("transform", "translate(18,13)")
        .selectAll("g").data([2,3,4,1]).enter()
        .append("g").attr("transform", function(d,i) { return "translate(0," + (i*(hUnit+1)) + ")" });
    lgOrg.append("circle").attr("id", function(d){ return 'lgDot' + d})
        .attr("cy",3.6).attr("r", cR1).attr("fill", function(d){return orgColor[d]});
    lgOrg.append("text").attr("x", 12).attr("dy", ".7em")
        .text(function(d) { return source[d] });
    var orgChinese = {2:'新冠病毒 (人)',3:'新冠病毒 (蝙蝠)',4:'SARS (人)',1:'SARS (蝙蝠)'};
    lgOrg.append("text").attr("x", 117).attr("dy", ".7em")
        .text(function(d) { return orgChinese[d] })
}

function geoTip(d){
    d3.select('#map'+d[0]+'_'+d[1]).attr("r", cR2*2);
    d3.selectAll('.g' + d[0] + 'o' + d[1] + ', #lgDot' + d[1]).attr("r", cR2);
    $('#ids .p'+d[0]+'o'+d[1]+' td' + ', #collect .p'+d[0]+'o'+d[1]+' td').css("background", orgColor[d[1]]).css("color", "white");
}

var base0 = yEdge + hUnit/2;
function drawORF(){
    var maxSeqL = d3.max(ids.map(function(d){return annoDataC[d].len}));
    var scale = d3.scaleLinear().domain([1, maxSeqL]).range([0, maxSeqL/ratio]);
    
    var svg = d3.select("#genome").append("svg")
                .attr("width", scale(maxSeqL)).attr("height", hSvg+1);

    var xAxis = d3.axisTop().scale(scale).tickSize(-3, 0).ticks(maxSeqL/1000)
                .tickFormat(function(f){ return f/1000 + 'k' });
    svg.append("g").attr("class", "xaxis").attr("transform", "translate(" + 0 + ",11)").call(xAxis);

    svg.append("text").attr("x",-50).attr("y",yEdge+2).attr("id","symbol");

// set orf colors by cdhit:
    var cidList = orfData.map(function(d){return d.cid? d.cid*1 : 0}).sort(function(a,b){return a-b});
    cidList = $.unique(cidList);
    cidList = jQuery.grep(cidList, function(d){ return d && d != -1 });

    color_class = d3.scaleOrdinal()
        .domain(cidList)
        .range(d3.range(1,cidList.length).map(function(d){
            return d3.hsl(360/(cidList.length+1) * (d-1), 1, 0.45)}));

    var base = base0,
        lTri=8,
        p1 = scale(1);
    ids.forEach(function(ac){
        var myStrain = svg.append("g").attr("class", 'con'+ac),
            line = myStrain.append("g").attr("class", "grayLine"),
            ll = annoDataC[ac].len;
        myStrain.attr('len', ll);

        //draw base line
        var p2 = scale(ll);
        line.append("line").attr("x1", p1).attr("x2", p2).attr("y1", base).attr("y2", base);
        line.append("line").attr("x1", p1).attr("x2", p1).attr("y1", base-hORF).attr("y2", base+hORF);
        line.append("line").attr("x1", p2).attr("x2", p2).attr("y1", base-hORF).attr("y2", base+hORF);
   
        var oo = orfData.filter(function(d){ return d.acc==ac});
            
        oo.forEach(function(obj){
            var start = Math.round(scale(obj.start)*10)/10,
                end = Math.round(scale(obj.stop)*10)/10;

            var xv = [start, end-lTri, end, end-lTri, start],
                yv = [base-hORF/2, base-hORF/2, base, base+hORF/2, base+hORF/2],
                points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');

            var mycolor = obj.cid>0? color_class(obj.cid) : '#bbb';
            var thisOrf = myStrain.append("polygon")
                    .attr("points", points)
                    .attr("fill", mycolor).attr("stroke", mycolor)
               
            if (!obj.orth) {
                var myId = obj.cid>0? obj.cid : obj.orf_id,
                    myClass = ac + '_' + myId; 
                thisOrf.attr("class", myClass);
                thisOrf.classed('whiteFill',true);
                if (!obj.cid) { thisOrf.classed('dotStroke',true)}
            }

            var tipText;
            if (obj.cid>0){
                tipText = '<em>' + symData[obj.cid] + '</em>'
            } else if (obj.sym){ tipText = '<em>' + obj.sym + '</em>' }
            else { tipText = '[' + obj.orf_id + ']' }
            thisOrf.on("mouseover", function(){
                $("#tipORF").css("left", (d3.event.pageX-22)+"px").css("top", (d3.event.pageY-127)+"px").show().html(tipText);
                mOver(ac)
                })
                .on("mouseout", function(){ $("#tipORF").hide(); mOut(ac)})

            if (!obj.cid || obj.cid<0){ return }

            thisOrf.classed('finger',true)
                   .on("dblclick", function(){ drawSeq(obj.cid) });
            
            if (obj.orth){
                thisOrf.attr("id", ac+'-'+obj.cid + '-' + obj.orth)
                        .attr('end', obj.stop)
                        .attr('len', obj.stop - obj.start + 1)
                        .on("click", function(){ alignORF(obj.cid, obj.orth) })
            } else {
                thisOrf.on("click", function(){ locateORF(obj.cid>0? obj.cid : obj.orf_id)})
            }
        });
        base += hUnit
    });

    orfLine = svg.append("line").attr("class","redStroke").attr("x1",-50).attr("x2",-50).attr("y2",hSvg);
    
    tabW = $('#tabs').width();
    $("#dnaComp").css("width", (tabW-5)+"px");
    $("#seqContainer").css("width", (tabW-10)+"px");

    $.ajax({
        async:false, dataType:"json", url:"js-css/dnacomp.json",
        success: function(data) {compData = data; drawComp()}
    }) 
}

function alignORF(cid, oid){
    var accs = orfData.filter(function(d){ return d.cid==cid && d.orth==oid})
            .map(function(d){return d.acc})
            .sort(function(a,b){return ids.indexOf(a) - ids.indexOf(b)});

    d3.selectAll('#genome svg .dotStroke').remove();
    var end_len=[], endMax=0, remainMax=0;
    
    $.each(accs, function(i,acc){
        var myid = '#'+acc+'-'+cid + '-' + oid,
            myorf = $(myid)[0],
            end =  myorf.attributes['end'].value * 1,
            len =  myorf.attributes['len'].value * 1;
        end_len.push([acc, end, len]);

        if (end > endMax) { endMax = end }
        var remain = $('.con' + acc)[0].attributes['len'].value - end;
        if (remain > remainMax) { remainMax = remain }
    });
    lMax = endMax + remainMax;

    d3.select("#genome svg").attr("width", lMax/ratio+1);

    transBar = {};
    var shadeX=[], shadeY=[];
    end_len.forEach(function(d,i){
        var acc = d[0],
            end = d[1],
            len = d[2],
            n = ids.indexOf(acc);

        var toTrans = Math.round((endMax-end)/ratio*10)/10;
        d3.select('.con' + acc).attr("transform", "translate(" + toTrans + ")");

        if (!i) { d3.select("genome .xaxis").attr("transform", "translate(" + toTrans + ',' + 20 + ")") }
        transBar[acc] = toTrans;

        //shade points
		var base = base0 + hUnit * n;
	    if (!shadeX.length){
            shadeX.push(endMax, (endMax-len));
            shadeY.push(base-hUnit/2, base-hUnit/2)
        }
	    shadeX.push(endMax-len);
        shadeY.push(base)
    });

    var lastX = shadeX[shadeX.length-1], lastY = shadeY[shadeY.length-1];
    shadeX.push(lastX, shadeX[0]); shadeY.push(lastY+hUnit/2, lastY+hUnit/2);

    //  draw shade
    var points = $.map(shadeX, function(vx,i){ return vx/ratio + ',' + shadeY[i]}).join(' ');
    d3.selectAll('#genome svg .shade').remove();
    d3.select('#genome svg').append("polygon").attr("id", "shadeC").attr("class","shade").attr("points", points);

    // gene symbol
    d3.select('#symbol').attr("x", (shadeX[0]+lastX)/2/ratio).text(symData[cid])
}

function locateORF(cid){
    d3.selectAll('#genome svg .dotStroke').remove();
    
    for (var i=0; i<ids.length; i++){
        var acc = ids[i],
            base = base0 + hUnit*i,
            myorfs = $('#genome .' + acc + '_'  + cid);
        
        if (!myorfs[0]){ continue }
        
        $.each(myorfs, function(n,myorf){
            var points = myorf.attributes.points.value.split(' ');
            var p1 = points[0].split(',')[0] * 1,
                p2 = points[2].split(',')[0],
                px = transBar[acc]? transBar[acc] : 0;
            
            var shd = d3.select('#genome svg').append("rect")
                .attr("class", 'grayLine dotStroke')
                .attr("x", p1+px).attr("y", base-7)
                .attr("width", p2-p1).attr("height", 14);
        })
    }
}

function appendPhe(myIds){
    $('#ids tr:not(.firstRow), #collect tr:not(.firstRow)').remove();
    $('#ids').append(myIds.map(function(x){return tdId[x]}).join(''));
    $('#collect').append(myIds.map(function(x){return tdCol[x]}).join(''));
    addEv()
}
function addEv(){
    $(".pheno td").css("height", hUnit+"px");
    d3.selectAll(".pheno td")
        .on('mouseover', function(){
            var arr = this.id.split('_');
            mOver(arr[1]);
            if (arr[0]=='ph'){$('#affiliation').html(pubData[annoDataC[arr[1]].pubId][2])}})
        .on('mouseout', function(){var arr = this.id.split('_'); mOut(arr[1])})
}

var cdhit, seqNt, firstNt, firstId, firstAA, svgSeqW, id;
function drawSeq(cid){
    tabNum = 2;
    $('#tabs').tabs({ active:tabNum });
    $('#tabs li:nth-of-type(3)').removeClass("disabled");
    $("#geneName em").html(symData[cid]);
    $('#geneName').show();
    
    if (cid && cdhit && cdhit==cid){showIns(2); return }
    cdhit = cid;

    if (!dndData){
        $.ajax({
            async:false, dataType:"json", url:"js-css/dndData.json",
            success: function(data) {dndData = data}
        })
    }
    d3.select('#tree_cid').remove();
    drawTreeC(cid);
    
    showIns(2);

    var orfs = orfData.filter(function(d){ return d.cid == cid});
    seqNt = {};
    orfs.forEach(function(d){ seqNt[d.acc] = seqData[d.orf_id] });

    var seqLen = seqNt[orfs[0].acc].length;
    svgSeqW = seqLen*unitW;

    d3.select('#seq svg').remove();
    var svg = d3.select('#seq').append("svg").attr("width", svgSeqW+unitW/2).attr("height", id.length*hUnit+yEdge);

    var codonBg = svg.append("g").attr("class", 'codonBg');
    for (var i=0; i<seqLen/6; i++){
        codonBg.append("rect").attr("x",i*6*unitW).attr("y", yEdge).attr("width",unitW*3).attr("height",id.length*hUnit)
    }

//    firstNt = '';
    svg.append("g").attr("id", "seqNt")
        .attr("transform", "translate(0," + yEdge + ")")
        .selectAll("text").data(id).enter().append("text")
        .attr("transform", function(d,i) { return "translate(0," + (i*hUnit) + ")" })
        .attr("dy", "1em")
        .html(function(d,i) {
            var ss = seqNt[d];
            if (!i){ firstNt = ss.split(''); firstId=d }
            else{ ss = matchSeq(ss)}
            return ss.replace(/-/g, '&nbsp;')
        });

    var axis=svg.append("g").attr("transform", "translate(-3," + (yEdge-hUnit+3) + ")"),
        tick = axis.append("g").attr("class", "xaxis"),
    	mark = axis.append("g").attr("class", "mark");
    for (var i=1; i<seqLen/30; i++){
        var px = i*30*unitW;
        tick.append("line").attr("x1",px).attr("x2",px).attr("y1",10).attr("y2",16);
        mark.append("text").attr("x",px+3).attr("y",6).text(i*30)
    }
    $("#seqContainer").scrollLeft(0);
    $('#showWhat').html("Show AA").fadeIn(3000);
    showAA=1;
    $("#guideH").css("left",(wTree+12)+'px').css("width", ($('#mainTab').width()-wTree-18)+"px").css("top", "-50px");
    $("#guideV").css("height",(hSvg)+"px").css("left", "-50px");

    if ($.inArray(cid*1, cid_comp)==-1){seqCid=0; $('#compContainer').hide()}
    else { seqCid=cid; showComp(2) }
}
function aaSeqTbl(){
    var seqAA = {};
    Object.keys(seqNt).forEach(function(d){seqAA[d] = transAA(seqNt[d])})

    d3.select("#seq svg").append("g").attr("id", "seqAA")
        .attr("transform", "translate(" + unitW + "," + yEdge + ")")
        .selectAll("text").data(id).enter().append("text")
        .attr("transform", function(d,i) { return "translate(0," + (i*hUnit) + ")" })
        .attr("dy", "1em")
        .html(function(d,i) {
            var ss = seqAA[d];
            if (d==firstId){ firstAA = ss.split('')}
            var arr = matchSeqAA(ss,d);
            return arr.join('&nbsp;&nbsp;').replace(/-/g, '&nbsp;')
        });

    d3.select('#seqAA').classed("opc0", true)
}
function matchSeq(d){
    return arr = d.split('').map(function(a,i){ return a==firstNt[i] && a!='-'? '.' : a }).join('')
}
function matchSeqAA(d,acc){
    var arr = d.split('').map(function(a,i){
        if (acc!=firstId && a==firstAA[i] && a!='-'){ a = '.' }
        return a
    })
    return arr
}
var codon = {AAA:"K",AAG:"K",AAC:"N",AAT:"N",GAT:"D",GAC:"D",GAG:"E",GAA:"E",CAA:"Q",CAG:"Q",CAC:"H",CAT:"H",TAT:"Y",TAC:"Y",TAG:"*",TAA:"*",TGA:"*",TGG:"W",TGC:"C",TGT:"C",CGT:"R",CGC:"R",CGG:"R",CGA:"R",GGA:"G",GGG:"G",GGC:"G",GGT:"G",AGT:"S",AGC:"S",AGG:"R",AGA:"R",ACA:"T",ACG:"T",ACC:"T",ACT:"T",GCT:"A",GCC:"A",GCG:"A",GCA:"A",CCA:"P",CCG:"P",CCC:"P",CCT:"P",TCT:"S",TCC:"S",TCG:"S",TCA:"S",TTA:"L",TTG:"L",TTC:"F",TTT:"F",CTT:"L",CTC:"L",CTG:"L",CTA:"L",GTA:"V",GTG:"V",GTC:"V",GTT:"V",ATT:"I",ATC:"I",ATG:"M",ATA:"I"};
function transAA(s){
    var arr = s.split('');
    var aa = '';
    for (var i=0; i<arr.length; i+=3){
        var c = arr[i] + arr[i+1] + arr[i+2];
        aa += c=='---'? '-' : codon[c.toUpperCase()]
    }
    return aa
}

function drawComp(){
    var wLg=300;
//maxComp
    var maxComp = d3.max(Object.keys(compData).map(function(d){
            var arr = compData[d];
            return max = d3.max(arr.map(function(j){ return d3.max(j)}))
        }));
    
    var stdAcc='MN996527',
        stdOrf = orfData.filter(function(d){ return d.acc==stdAcc }),
        stdCid1 = stdOrf.filter(function(d){return d.cid==1})[0];
    ratioC = d3.format(".3n")((stdCid1.stop-stdCid1.start+1)/(tabW-10)) * 1;

    var stdSeqL = annoDataC[stdAcc].len;
    var compH=100,
        marginB = 17,
        chartH = compH - marginB;
    var scaleX = d3.scaleLinear().domain([1, stdSeqL]).range([0, stdSeqL/ratioC]),
        scaleY = d3.scaleLinear().domain([0,maxComp]).range([chartH,0]),
        axisX = d3.axisBottom().tickSize(-3, 0).tickFormat(function(f){ return f/1000 + 'k' }),
        axisY = d3.axisLeft().scale(scaleY).tickSize(-3, 0).ticks(maxComp/2);

    //draw legend
    var svgLegend = d3.select("#compY").append("svg").attr("width",wLg).attr("height",chartH+24);
    svgLegend.append("line").attr("class","grayLine").attr("x2",wLg-7);
    svgLegend.append("g").attr("class","xaxis").attr("transform", "translate("+(wLg-3)+",0)").call(axisY);
    svgLegend.append("text").attr("transform", "translate(270,82)rotate(-90)").text("NT Changes");

    var svgLg = svgLegend.append("g").attr("transform", "translate(134,22)");
    svgLg.append("text").text("Codon Position");
    var lgComp = svgLg.append("g").attr("class","lgOrg")
        .attr("transform", "translate(10,17)")
        .selectAll("g").data([1,2,3]).enter()
        .append("g").attr("transform", function(d,i) { return "translate(0," + (i*(hUnit+1)) + ")" }),
        lenL = 40;
    lgComp.append("line").attr("x2",lenL).attr("class", function(d){return 'path'+(d-1)});
    lgComp.append("text").attr("x", lenL+4).attr("dy", ".35em")
        .text(function(d) { return d==1? "1st" : (d==2? "2nd" : "3rd") });

    var svg = d3.select("#dnaComp").append("svg").attr("width", scaleX(stdSeqL)).attr("height",compH);

    svg.append("line").attr("class","grayLine").attr("transform", "translate(0," + chartH/2 + ")").attr("x1", scaleX(1)).attr("x2", scaleX(stdSeqL));

    var orfRect = svg.append('g').attr("class","grayLine whiteFill")
    var arrLine = [];
    for (var i=1; i<=parseInt(maxComp); i++){arrLine.push(i)}

    stdOrf.forEach(function(obj){
        var start = obj.start,
            end = obj.stop,
            len = end - start +1;
        orfRect.append("rect").attr("x", scaleX(start)).attr("width", len/ratioC).attr("height", chartH);
        svg.append("rect").attr("class","lowOpc").attr("x", scaleX(start)).attr("width", len/ratioC).attr("height", chartH).attr("fill", color_class(obj.cid));
        
        if ($.inArray(obj.cid*1, cid_comp)!=-1){
            drawC(obj.cid, start, len);
        }
    });
    compLine = svg.append("line").attr("class","redStroke")
    			.attr("x1","-50px").attr("x2","-50px")
    			.attr("y2",chartH+3);

    function drawC(cid, start, len){
        var scale = d3.scaleLinear().domain([1, len]).range([0, len/ratioC]);
        axisX.scale(scale).ticks(len/200);
        var chart = svg.append("g").attr("id", "comp"+cid).attr("class", "charts")
                .attr("transform", "translate(" + scaleX(start) + ",0)");
        chart.append("g").attr("class", "xaxis")
            .attr("transform", "translate(0,"+(chartH+6)+")").call(axisX);
        chart.append("g").attr("class", "compBgLine")
            .selectAll("line").data(arrLine).enter()
            .append("line").attr("x2", scaleX(len))
            .attr("y1", function(d){return scaleY(d)}).attr("y2", function(d){return scaleY(d)});

        var funLine=[];
        for (var i=0; i<3; i++){
            funLine[i] = d3.line()
                    .x(function(d,j){return scale(j*3+5+i)})
                    .y(function(d){return scaleY(d)})
                    .curve(d3.curveLinear)
        }
        
        var data = compData[cid];
        for (var i=0; i<3; i++){
            chart.append("path").attr("class", 'path'+i).attr("d", funLine[i](data[i]))
        }
    }
}

//tree functions:
function parseNewick(a){
    if (!/^[\s\n\r]*\(.+\);[\s\n\r]*$/.test(a)){return 0}
    for(var e=[],r={},s=a.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){var n=s[t];switch(n){case"(":var c={};r.children=[c],e.push(r),r=c;break;case",":var c={};e[e.length-1].children.push(c),r=c;break;case")":r=e.pop();break;case":":break;default:var h=s[t-1];")"==h||"("==h||","==h?r.name=n:":"==h&&(r.length=parseFloat(n))}}return r
}
function maxLength(d) {return d.data.length + (d.children ? d3.max(d.children, maxLength) : 0)}
function resetY(d, x0, k){
    d.y = (x0 += d.data.length) * k;
    if (d.children) d.children.forEach(function(d) { resetY(d, x0, k); });
}
function link(d) {return linkStep(d.source.y, d.source.x, d.target.y, d.target.x)}
function linkExtension(d) {return linkStep(d.target.y, d.target.x, wTree, d.target.x)}
var lineFunction = d3.line().x(function(d){return d[0]}).y(function(d){return d[1]}).curve(d3.curveLinear);
function linkStep(x1, y1, x2, y2) {
    var arr = [[x1,y1]];
    if (y2!=y1){ arr.push([x1,y2]) }
    arr.push([x2,y2]);
    return lineFunction(arr)
}

function showIns(d, cid){
    tabNum = d;
    if (d==1){$('#instruction').show()}
    else {
        if (d==2 && $('#tabs li:nth-of-type(3)').hasClass("disabled")){ return }
        $('#instruction').hide()
    }
    if (d==2){
        $("#tree").fadeOut(1500);
        $('#tree_cid').fadeIn(1500);
        appendPhe(id);
        $('#treeScale tspan').html(sub_site[1]);
        $('#geneName, #showWhat, #guideH').show()
    } else {
        $('#tree_cid').fadeOut(1500);
        $('#tree').fadeIn(1500);
        $('#treeScale tspan').html(sub_site[0]);
        appendPhe(ids);
        $('#geneName, #showWhat, #guideH').hide()
    }
    if (d==1 || (d==2 && seqCid)){showComp(d)} else {$('#compContainer').hide()}
}

function showComp(d){
    $('#compContainer, .charts').show();
    if (d==2){
        $('.charts:not(#comp' + seqCid + ')').hide();
        seqScroll = $('#comp'+seqCid)[0].attributes.transform.value.split('(')[1].split(',')[0]*1;
        $("#dnaComp").scrollLeft(seqScroll)        
    }
}

function showAck(){$('#ack').show()}

var outId;
function readHaps(data){
    var allNodes = data[0];
    var area = Math.PI*Math.pow(stdR,2)

    for (var i=0; i<allNodes.length; i++){
        var node = allNodes[i],
            radius = stdR;

        if (!node.haps){
            listNode.push({id:node.id, radius:radius})
            continue
        }

        var nn = node.haps
        
        listVirus = listVirus.concat(nn)

        var grp={};
        nn.forEach(function(ac){
            var gid = annoDataW[ac].geoId;
            var gg = [];
            if (grp[gid]){gg = grp[gid]}
            gg.push(ac);
            grp[gid] = gg
        })
        var grps = Object.keys(grp).map(function(id){return [id*1, grp[id]]})
                .sort(function(a,b){return geoDataW[a[0]].ctry - geoDataW[b[0]].ctry}),
            ng = grps.map(function(oo){
                    var min = d3.min(oo[1].map(function(ac){ return tParser(annoDataW[ac].col)}));
                    return [node.id, oo[0], min]})
        node_gid = node_gid.concat(ng)

        if (grps.length>1) radius = Math.sqrt(grps.length*area/Math.PI);
        listNode.push({id:node.id, haps:grps, radius:radius})
            
        if (nn[0]==outGrp){outId=node.id}
    }

    var ro = {id:"H0", radius:stdR};
    var listEdge = data[1],
        listChg = [];
//    var rObj={};
    for (var i=0; i<listEdge.length; i++){
        var edge = listEdge[i];
        var so = listNode.filter(function(n){ return n.id==edge.source})[0],
            ta = listNode.filter(function(n){ return n.id==edge.target})[0],
            ch = edge.change;

        if (so.id==outId){
            listLink.push({
                source: ro,
                target: so,
                ldist: ro.radius + so.radius + distance,
            });
            listLink.push({
                source: ro,
                target: ta,
                ldist: ro.radius + ta.radius + distance,
                change:ch
            })
        } else {
            listLink.push({
                source: so,
                target: ta,
                ldist: so.radius + ta.radius + distance,// + lnkdist,
                change:ch
            })
        }
        listChg = listChg.concat(ch)
    }
    listNode.push(ro);
    
    var site = {};
    for (var i=0; i<listChg.length; i++){
        site[listChg[i].site]=1
    }
    listSite = Object.keys(site).map(d=>d*1)
}

var batPath="M598.935,485.829c-2.286,5.292-8.146,7.917-8.146,7.917c0.025-0.684-0.068-1.365-0.279-2.016c1.077-1.284,1.363-3.057,0.745-4.614c-1.251,0.396-2.317,1.229-3.006,2.345c-0.419-0.118-0.853-0.175-1.287-0.169c-0.426,0.016-0.848,0.09-1.253,0.22c-0.711-1.098-1.788-1.908-3.04-2.286c-0.587,1.574-0.262,3.344,0.847,4.606c-0.229,0.639-0.346,1.312-0.347,1.99c0,0-5.928-2.54-8.299-7.799c-3.796,4.166-8.123,7.815-12.87,10.855c0,0,5.69-0.991,7.553,3.073c0,0,6.029-3.167,7.553,2.092c0,0,5.826-5.927,10,7.781c3.963-13.76,9.873-7.942,9.873-7.942c1.432-5.283,7.52-2.21,7.52-2.21c1.795-4.09,7.502-3.192,7.502-3.192C607.195,493.521,602.802,489.939,598.935,485.829zM574.016,490.63c-1.797,2.231-3.075,4.834-3.743,7.621c-0.024,0.147-0.154,0.256-0.305,0.254H569.9c-0.173-0.008-0.307-0.154-0.3-0.327c0.002-0.033,0.009-0.065,0.021-0.097c0.683-2.876,2.009-5.56,3.878-7.85c0.102-0.14,0.298-0.169,0.438-0.067c0.14,0.103,0.17,0.299,0.067,0.438c-0.007,0.01-0.015,0.019-0.022,0.027H574.016zM578.748,494.212c-0.708,1.744-1.022,3.623-0.923,5.504c0.021,0.172-0.1,0.329-0.271,0.355h-0.042c-0.159,0.005-0.296-0.113-0.313-0.271c-0.123-1.983,0.204-3.969,0.957-5.809c0.052-0.16,0.225-0.247,0.385-0.194c0.007,0.003,0.014,0.005,0.021,0.008c0.162,0.06,0.247,0.238,0.188,0.401C578.75,494.208,578.749,494.21,578.748,494.212L578.748,494.212z M586.217,492.561c-0.102,0.22-0.627,0.187-1.177-0.068c-0.551-0.254-0.923-0.635-0.847-0.847c0.076-0.211,0.618-0.186,1.177,0.068S586.31,492.349,586.217,492.561z M588.901,492.459c-0.542,0.263-1.076,0.305-1.178,0.093s0.263-0.601,0.847-0.847c0.585-0.245,1.076-0.305,1.178-0.085C589.85,491.841,589.451,492.196,588.901,492.459z M596.793,499.673c-0.024,0.158-0.154,0.278-0.313,0.288l0,0c-0.172-0.022-0.295-0.175-0.279-0.348c0.068-1.883-0.276-3.759-1.008-5.495c-0.063-0.161,0.016-0.343,0.177-0.406h0.001c0.163-0.063,0.349,0.016,0.415,0.178C596.545,495.719,596.889,497.694,596.793,499.673z M604.058,498.217h-0.06c-0.148,0-0.277-0.102-0.313-0.246c-0.718-2.761-2.039-5.329-3.869-7.519c-0.117-0.127-0.117-0.322,0-0.449c0.136-0.11,0.336-0.092,0.448,0.043c1.926,2.248,3.318,4.902,4.073,7.765c0.029,0.175-0.089,0.341-0.265,0.37c-0.005,0.001-0.01,0.002-0.015,0.002V498.217z"

//function numberWithCommas(x) { return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") }


function seen(s, arr){
    var see={}, yes=0;
    for (var i=0; i<arr.length; i++){
        if (s==arr[i]){yes=1; break}
    }
    return yes
}