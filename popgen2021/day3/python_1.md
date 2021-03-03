### Paleo-Programing Facts 

| 1843:  [Ada Lovelace](https://en.wikipedia.org/wiki/Ada_Lovelace) wrote [the first computer program](https://twobithistory.org/2018/08/18/ada-lovelace-note-g.html) calculating the Bernoulli numbers. |
|:--:| 
| <img src="../img/ada.jpg" alt="drawing" width="400"/> | 
|*> In 1842, the Italian mathematician Luigi Federico Menabrea published a description of the Analytical Engine, a proposed mechanical general-purpose computer, based on a lecture by Babbage in French. In 1843, the description was translated into English and extensively annotated by Ada Lovelace, who had become interested in the engine eight years earlier. In recognition of her additions to Menabrea's paper, which included a way to calculate Bernoulli numbers using the machine (widely considered to be the first complete computer program), she has been described as the first computer programmer.  (from https://en.wikipedia.org/wiki/Analytical_Engine* )|



| 1950s: IBM develops the [FORTRAN](https://en.wikipedia.org/wiki/Dorothy_Vaughan) for scientific and engineering applications, FORTRAN came to dominate this area of programming early on and has been in continuous use for over half a century in computationally intensive areas. |
|:--:| 
| <img src="../img/FortranCardPROJ039.agr.jpg" alt="drawing" width="300"/> | 


| 1949: [Dorothy Vaughan](https://en.wikipedia.org/wiki/Dorothy_Vaughan) became acting supervisor of the West Area Computers, the first African-American woman to supervise a group of staff at the center. |
|:--:| 
| <img src="../img/doroth_vaughan_nasa.jpg" alt="drawing" width="300"/> | 
|*> During her 28-year career, Vaughan prepared for the introduction of machine computers in the early 1960s by teaching herself and her staff the programming language of FORTRAN; she later headed the programming section of the Analysis and Computation Division (ACD) at Langley. (source https://en.wikipedia.org/wiki/Dorothy_Vaughan)* |

| 1950s: [Grace Hopper](https://en.wikipedia.org/wiki/Grace_Hopper) popularized the idea of machine-independent programming languages, which led to the development of COBOL, an early high-level programming language still in use today. |
|:--:| 
| <img src="../img/grace-hopper-3-600x403-1.jpg" alt="drawing" width="300"/> |
|*>While she was working on a Mark II Computer at Harvard University in 1947, her associates discovered a moth that was stuck in a relay; the moth impeded the operation of the relay. While neither Hopper nor her crew mentioned the phrase "debugging" in their logs, the case was held as an instance of literal "debugging." For many years, the term bug had been in use in engineering. **The remains of the moth can be found in the group's log book at the Smithsonian Institution's National Museum of American History in Washington, D.C.**  (from https://en.wikipedia.org/wiki/Grace_Hopper)*| 
| <img src="../img/bug.jpg" alt="drawing" width="300"/> |


***
### Programing languages evolution 

| <img src="../img/programingtree.jpg" alt="drawing" width="400"/> | The "phylogenetic" tree of languages | 
|:--:|:--:|
| *(source http://www.cs.sjsu.edu/faculty/pearce/modules/lectures/languages3/history/evolution.htm)* |


***
### A quick history of Python 

|<img src="../img/pythonlogo.jpg" alt="drawing" width="300"/> | <img src="../img/DO6GvRhi.gif" alt="drawing" width="300"/>|
|:--:|:--:|
| |https://en.wikipedia.org/wiki/Guido_van_Rossum | 


Python2.x will be dismissed in 2020. Current relase is Python3.x. Read more about [current Python release](https://wiki.python.org/moin/Python2orPython3)

The name Python was inspired by the [Monty Python](http://www.montypython.com/python_The_Pythons/14) comedy series. 


| | |
|----|----|
| Python tutorial | https://docs.python.org/3/tutorial/  |
|Python for mathematics, science, and engineering | https://scipy.org/ |
| Python official website | https://www.python.org/  |


*** 
### Hands-on Python 

#### Setup 
0. Make sure you have a plain text editor ready like notepad, geany, notepad++... do not use word for editing files 
1. Connect to the server as described [here](../WiFi-SSHinstruction.md) **Make sure you connect using the -Y flag**
2. In your home make a folder that you will use for all the python. Use a meaningful name *be consistent and descriptive*. See suggestions for naming [here](https://library.stanford.edu/research/data-management-services/data-best-practices/best-practices-file-naming)

#### Software Carpentry lessons 

We will use the [Software Carpentry lessons](https://software-carpentry.org/lessons/) to learn Python. 

3. Download to your home directory on the server the toy-dataset and the code from the SC lesson repository [here](http://swcarpentry.github.io/python-novice-inflammation/setup/). Use `wget`. 

```
wget http://swcarpentry.github.io/python-novice-inflammation/data/python-novice-inflammation-data.zip

wget http://swcarpentry.github.io/python-novice-inflammation/code/python-novice-inflammation-code.zip
```
After that you will see two new folders named  `data` and  `code` in your home directory. 

We will now go through the SC Python lessons that can be found [here](http://swcarpentry.github.io/python-novice-inflammation/). 

#### Links to modules used during lessons 
| | |
|----|----|
|numpy |http://www.numpy.org/ |
|time | https://docs.python.org/3.7/library/time.html | 
|matplotlib | https://matplotlib.org/ | 

#### Objectives of today's lesson 

- Learn about libraries 
- `import` libraries and use the functions in the library 
- Learn about variables
- Read tabular data from a file into a program and access slices of data 
- Plot simple graphs from data
