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
	    i=0
	    for (var key in data) {
		if (data[key][2] == 1) {
		    base_phrase.push(data[key][0])
		    materials[i] = []
		    i=i+1
		}
	    }
	    // each base sentence is now a key in materials
	    for (var key in materials) {
		i = 0
		while (i<data.length) {
		    if (base_phrase[key] == data[i][0]) {
			sub_materials = new Object()
			sub_materials[data[i+1][0]]=data[i+1][2]
			sub_materials[data[i+2][0]]=data[i+2][2]
			sub_materials[data[i+3][0]]=data[i+3][2]
			sub_materials[data[i+4][0]]=data[i+4][2]
			materials[key] = sub_materials
			break;
		    }
		    i = i + 1
		}
	    }
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
