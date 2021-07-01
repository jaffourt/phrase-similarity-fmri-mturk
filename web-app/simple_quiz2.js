$(document).ready(function() {
  //Below is the question/answer structure

    function buildQuiz() {

	get_materials_object(function(err,results){
	    if (err){console.log(err)}
	    else {
		const output = [];
		const myQuestions = results
		
		// for each question...
		myQuestions.forEach((currentQuestion, questionNumber) => {
		    // we'll want to store the list of answer choices
		    const answers = [];
		    var counter = 1;
		    //console.log(questionNumber);
		    // console.log(currentQuestion);
		    // console.log(currentQuestion.answers);

		    // and for each available answer...
		    currentQuestion.answers.forEach((currentAnswer, ansNumber) => {
			// for (num in currentQuestion.answers) {
			var q_choices = [];
			q_choices[questionNumber] = [];
			//console.log(num)
			q_choices[questionNumber].push (
			    `${currentQuestion.answers[ansNumber]}: \t`
			);
			
			// ...add 1-4 option choices
			for (i = 1; i <= 4; i++) { 
			    if (i==1) {
				q_choices[questionNumber].push (
				    `<div class="qa${questionNumber}" name="choice${ansNumber}" > ${i}<input type="radio" name="answersT${ansNumber}" value="${i}">`
				);
				counter+=1
				continue
			    }
			    if (i==4) {
				q_choices[questionNumber].push (
				    `<input type="radio" name="answersT${ansNumber}"  value="${i}">${i} </div>`
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
			//console.log(q_choices[num].join(""))
			// ...add the HTML rating per question
			answers.push(
			    `<label>${q_choices[questionNumber].join("")}</label>`
			);

		    });

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
	});
  }

  //Make a button listener to block answers  
  
/*
  function showResults() {
    // gather answer containers from our quiz
    const answerContainers = quizContainer.querySelectorAll(".answers");

    // keep track of user's answers
    let numCorrect = 0;

    // for each question...
    myQuestions.forEach((currentQuestion, questionNumber) => {
      // find selected answer
      const answerContainer = answerContainers[questionNumber];
      const selector = `input[name=question${questionNumber}]:checked`;
      const userAnswer = (answerContainer.querySelector(selector) || {}).value;

      // if answer is correct
      if (userAnswer === currentQuestion.correctAnswer) {
        // add to the number of correct answers
        numCorrect++;

        // color the answers green
        answerContainers[questionNumber].style.color = "lightgreen";
      } else {
        // if answer is wrong or blank
        // color the answers red
        answerContainers[questionNumber].style.color = "red";
      }
    });

    // show number of correct answers out of total
    resultsContainer.innerHTML = `${numCorrect} out of ${myQuestions.length}`;
  }
*/
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

  function preDisbaleOtherButtons(n){
    //Grab slide active-slide
    var active = $("div[class='slide active-slide']")
    console.log(active)
    //Grab radios within the active slide
    var radios = $("input[name*='answersT']",active)
    console.log(radios)
    for(var i = 0, max = radios.length; i < max; i++) {
    radios[i].onclick = function() {
        console.log(this);
        disableOtherButtons(this)
      }
    }
  }

  function validateAllButtons() {
    //Grab slide active-slide
    var active = $("div[class='slide active-slide']")
    console.log(active)
    //Grab radios within the active slide
    var radios = $("input[name*='answersT']:checked",active)
    console.log(radios.length)
    if (radios.length < 4) {
      alert("Please rate all options")
      return false
    } 
    //Check all values are unique
    var uniq = [];
    for(var i = 0, max = radios.length; i < max; i++) {
      uniq.push(radios[i].value)
    }

    if (( new Set(uniq)).size !== uniq.length ) {
      alert("Please rate each option differently")
      return false
    }
    
    return true
  }

  function showSlide(n) {
    slides[currentSlide].classList.remove("active-slide");
    slides[n].classList.add("active-slide");
    currentSlide = n;
    
    if (currentSlide === 0) {
      previousButton.style.display = "none";
    } else {
      previousButton.style.display = "inline-block";
    }
    
    if (currentSlide === slides.length - 1) {
      nextButton.style.display = "none";
      submitButton.style.display = "inline-block";
    } else {
      nextButton.style.display = "inline-block";
      submitButton.style.display = "none";
    }

  }

  function showNextSlide() {
    validated = validateAllButtons()

    if (validated == true){
      showSlide(currentSlide + 1);
      //Check that all buttons are pressed
    }
  }

  function showPreviousSlide() {
    showSlide(currentSlide - 1);
  }
    
  const quizContainer = document.getElementById("quiz");
  const resultsContainer = document.getElementById("results");
  const submitButton = document.getElementById("submit");

  // display quiz right away
  buildQuiz();

  const previousButton = document.getElementById("previous");
  const nextButton = document.getElementById("next");
  const slides = document.querySelectorAll(".slide");
  let currentSlide = 0;

  showSlide(0);

  // on submit, show results
  submitButton.addEventListener("click", showResults);
  previousButton.addEventListener("click", showPreviousSlide);
  nextButton.addEventListener("click", showNextSlide);
})
