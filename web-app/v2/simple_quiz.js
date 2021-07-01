$(document).ready(function() {

    //Below is the question/answer structure
    var transform_set = [] 

    function makeRandom20() {
	var text = "";
	var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

	for (var i = 0; i < 20; i++)
	{
	    text += possible.charAt(Math.floor(Math.random() * possible.length));
	}
	return text;
    }
    
    var shuffle = function (array) {
        var currentIndex = array.length;
        var temporaryValue, randomIndex;
        // While there remain elements to shuffle...                                                                                                                       
        while (0 !== currentIndex) {
                // Pick a remaining element...                                                                                                                             
                randomIndex = Math.floor(Math.random() * currentIndex);
                currentIndex -= 1;
                // And swap it with the current element.                                                                                                                   
                temporaryValue = array[currentIndex];
                array[currentIndex] = array[randomIndex];
                array[randomIndex] = temporaryValue;
        }
        return array;
    };

    function buildQuiz(myQuestions, quizContainer) {
	const output = [];		
	// for each question...
	myQuestions.forEach((currentQuestion, questionNumber) => {
	    // we'll want to store the list of answer choices
	    const answers = [];
	    var counter = 1;
	    var ansNumber = 0;
	    var currentAnswers = []
	    for (key in currentQuestion.answers)
	    {
		currentAnswers.push(key);
	    }	    
	    var transformations = []
	    for (key in currentAnswers)
	    {
		transformations.push(currentQuestion.answers[currentAnswers[key]])
	    }
	    transform_set.push(transformations)
	    for (val in currentAnswers) {
		currentAnswer = currentAnswers[val]
		ansNumber = ansNumber + 1;
		// and for each available answer...
		// currentQuestion.answers.forEach((currentAnswer, ansNumber) => {
		// for (num in currentQuestion.answers) {
		var q_choices = [];
		q_choices[questionNumber] = [];
		//console.log(num)
		q_choices[questionNumber].push (
		    `${currentAnswer} \t` )
		// ...add 1-4 option choices
		for (i = 1; i <= 4; i++) { 
		    if (i==1) {
			q_choices[questionNumber].push (
			    `<div class="qa${questionNumber}" name="choice${ansNumber}"> <font size="4px">VERY SIMILAR</font> <input type="radio" name="answersT${ansNumber}" value="${i}" >`
			);
			counter+=1
			continue
		    }
		    if (i==4) {
			q_choices[questionNumber].push (
			    `<input type="radio" name="answersT${ansNumber}"  value="${i}"> <font size="4px">VERY DISSIMILAR</font> </div>`
			);
			counter+=1
			break
		    }
		    q_choices[questionNumber].push (
			`<input type="radio" name="answersT${ansNumber}"  value="${i}">`
		    );
		    counter+=1
		    continue
		}
		//console.log(q_choices[0].join(""))
		// ...add the HTML rating per question
		answers.push(
		    `<label>${q_choices[questionNumber].join("")}</label>`
		);
	    }

	    // add this question and its answers to the output
	    output.push(
		`<div class="slide">
           <div class="question"> ${currentQuestion.question} </div>
           <div class="answers"> ${answers.join("")} </div>
         </div>`
	    );
	});		    
	// finally combine our output list into one string of HTML and put it on the page
	quizContainer.innerHTML = output.join("");
    }
    
    function showResults() {
	quizContainer = this.quizContainer
	myQuestions = this.myQuestions
	validated = validateAllButtons(this.nextButton)
	if (validated == true){
	    document.getElementById(this.id).disabled = true;
            document.getElementById(this.id).style.backgroundColor = "#B7B7B7"
            document.getElementById(this.id).textContent="Submitted";
	    var active = $("div[class='slide active-slide']")
            const confirmation_code = makeRandom20()
            active[0].innerText = 'Confirmation code: ' + confirmation_code
	    answers = this.answers 
	    const header = "set,transformation,ranking"
	    var line = header + "\n"
	    iter = 0
	    for (var i = 0; i < answers.length; i++)
	    {
		iter = iter + 1
		line += iter+","+transform_set[i][0]+","+answers[i][0]+"\n"
		if (iter == 32)
		{
		    iter = 0
		}
	    }
	    // form data and then php to save to server
	    var data = new FormData();
	    data.append("data" , line);
	    data.append('id', confirmation_code)
	    var xhr = new XMLHttpRequest();
	    xhr.open( 'post', 'save_id.php', true);
	    xhr.send(data);
	}
    }
  /*  
  function disableOtherButtons(button) {
    //Grab slide active-slide
    var active = $("div[class='slide active-slide']")
    //console.log(active)
    //Grab radios within the active slide
    var radios = $("input[name!='"+ button.name +"']",active)
    //console.log(radios)
    for(var i = 0, max = radios.length; i < max; i++) {
        r_ques = radios[i].name
        r_ans = radios[i].value
        if (r_ans == button.value){
          console.log(r_ans)
          radios[i].disabled = true;
        }
        //console.log(radios[i]);
    }
  }

  function preDisableOtherButtons(n){
    //Grab slide active-slide
    var active = $("div[class='slide active-slide']")
    //Grab radios within the active slide
    var radios = $("input[name*='answersT']",active)
    for(var i = 0, max = radios.length; i < max; i++) {
    radios[i].onclick = function() {
        console.log(this);
        disableOtherButtons(this)
      }
    }
  }
*/
  function validateAllButtons(nextButton) {
      //Grab slide active-slide
      var active = $("div[class='slide active-slide']")
      //Grab radios within the active slide
      var radios = $("input[name*='answersT']:checked",active)
      // debugging
      // return true
      if (radios.length < 1) {
	alert("Please rate all options")
	  return false
      } 
      var temp = []
      for (i=0; i<1; i++) {
	  temp.push(radios[i].value)
      }
      nextButton.submitButton.answers.push(temp)      
      return true
  }

    function showSlide(slides, currentSlide, nextButton, submitButton, n) {
	slides[currentSlide].classList.remove("active-slide");
	slides[n].classList.add("active-slide");
	nextButton.currentSlide = n;
	if (nextButton.currentSlide === slides.length - 1) {
	    nextButton.style.display = "none";
	    submitButton.style.display = "inline-block";
	} else {
	    nextButton.style.display = "inline-block";
	    submitButton.style.display = "none";
	}	
  }

    function showNextSlide(slides, currentSlide, nextButton, submitButton) {
	validated = validateAllButtons(this)
	if (validated == true){
	    showSlide(this.slides, this.currentSlide, this, this.submitButton, this.currentSlide + 1);
      //Check that all buttons are pressed
    }
  }
   // start the "main" function 
   get_materials_object(function(err,results){
       if (err){console.log(err)}
       else { 
	   const myQuestions = results
	   const quizContainer = document.getElementById("quiz");
	   const resultsContainer = document.getElementById("results");
	   const submitButton = document.getElementById("submit");

	   // display quiz right away
	   buildQuiz(myQuestions, quizContainer);

	   const nextButton = document.getElementById("next");
	   const slides = document.querySelectorAll(".slide");
	   let currentSlide = 0;

	   // textbook object oriented programming 
	   nextButton.slides = slides;
	   nextButton.currentSlide = currentSlide
	   nextButton.submitButton = submitButton
	   
	   submitButton.quizContainer = quizContainer
	   submitButton.myQuestions = myQuestions
	   submitButton.answers = []
	   submitButton.nextButton = nextButton
	   // start the sequence -- show first slide
	   showSlide(slides, currentSlide, nextButton, submitButton, 0);

	   // on submit, show results
	   submitButton.addEventListener("click", showResults);
	   nextButton.addEventListener("click", showNextSlide);
       }
   });
})
