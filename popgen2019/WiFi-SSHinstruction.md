

##  WiFi

network name: embopopgen2019

password: embopopgen2019

##  SSH
For this course we will work on the High Performance Cluster of the [Istituto Nazionale di Fisica Nucleare](https://www.ba.infn.it/index.php/en/)

The filesystem has been generated according to this structure

<img src="./img/directoryscheme.png" alt="yay">

Each has writing permission only in the folder(s) with their name and reading permission to any folder.


To connect to the INFN machine using SSH from a Linux terminal use:

```
ssh -Y yourusername@elixirschool.recas.ba.infn.it
```
the -Y tag enable graphical options

Username and password were given to you at the reception


## Transfer files from the server to local computers

#### I have a linux computer
open a terminal and scp to your local machine

```
scp username@elixirschool.recas.ba.infn.it:/home/username/myfiletotransfer.extension .

```

#### I have a mac
open a terminal and scp to your local machine

```
scp username@elixirschool.recas.ba.infn.it:/home/username/myfiletotransfer.extension .

```

in alternative download and use [cyberduck](https://cyberduck.io/)

#### I have a PC

For window user we recommend using MobaXterm [Download the program here](http://mobaxterm.mobatek.net/download-home-edition.html)

- choose Portable edition
- unzip the file
- open and run the 'exe' file.
- a press start a new session. Enter server, login and password


## Open pdfs on the INFN machine

To open pdfs use [evince](https://en.wikipedia.org/wiki/Evince)

```
evince  mypdf.pdf
```


 #### thanks to Stefano Nicotri, Giacinto Donvito
