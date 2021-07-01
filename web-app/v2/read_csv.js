function get_materials_object(callback) {
    var materials = new Object()
    var base_phrase = []
    var myQuestions = []
    $.ajax({
	type: 'get',
	url: 'read_local.php',
	data: '',
	success: function(csv) {
	    var objects = Papa.parse(csv)
	    // remove the header row
	    // and also the \n ending row
	    data = objects.data
	    data.shift()
	    data.length--
	    // data[n] = the nth sentence object
	    // data[n][0] = the base sentence
	    // data[n][1] = the set (which base sentence object)
	    // data[n][2] = the transformation number: 1 = base, 2-5 = transformations
	    // 1-4 because we now want to present each base phrase with just one transformation
	    i=0	
	    for (var key in data) {
		if (data[key][2] == 1) {
		    base_phrase.push(data[key][0])
		    materials[i] = []
		    i=i+1
		}
	    }
    
	    // create a convention and stick to it
	    // condA = meaning preserve
	    // condB = noun change
	    // condC = preposition change
	    // condD = adjective change
	    // list 1: j = 0
	    // list 2: j = 1
	    // list 3: j = 2
	    // list 4: j = 3
	    var j = 3
	    for (var key in materials) {
		i = 0
		while (i<data.length) {
		    if (base_phrase[key] == data[i][0]) {
			sub_materials = new Object()
			j = j+1
			sub_materials[data[i+j][0]]=data[i+j][2]
			materials[key] = sub_materials
			if (j==4) {
			    j = 0
			}
			break;
		    }
		    i = i + 1
		}
	    }
	    console.log(materials)
	    for (var key in materials) {
		temp = new Object()                                                                 
		temp["question"] = base_phrase[key]
		temp["answers"] = materials[key]
		myQuestions.push(temp)
	    }
	    return callback(null, myQuestions)
	}
    });
}
