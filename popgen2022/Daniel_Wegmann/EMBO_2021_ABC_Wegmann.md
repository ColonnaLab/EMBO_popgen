Material
--------

While this lecture will take place on a virtual black board, here is
some accompaning material:

-   The script [The Art of Statistical Modeling and
    Inference](https://bitbucket.org/wegmannlab/statistical_modelling_and_inference/raw/master/Statistical_Modelling_and_Inference.pdf),
    with a detailed chapter on ABC.
-   The book-chapter [A Guide to General-Purpose Approximate Bayesian
    ComputationSoftware](https://arxiv.org/pdf/1806.08320.pdf).
-   A video on [Bayesian
    statistics](https://www.youtube.com/watch?v=S873W1RwSn4&t=902s).
-   A video on [Approximate Bayesian
    computation](https://www.youtube.com/watch?v=xtPHTSVD4T8&t=2s).
-   A video on
    [ABC-GLM](https://www.youtube.com/watch?v=r3iAqCjODrk&t=16s).
-   The R script
    [ABC\_demo.r](https://bitbucket.org/wegmannlab/statistical_modelling_and_inference/raw/master/Scripts/ABC_demo.r)
    used to illustrate ABC.

Exercises A (Morning): Getting started with ABC
-----------------------------------------------

### Step A1: Compile the latest version of ABCtoolbox and check executables

Let us begin by creating a directory “bin” for all executables
(binaries) of programs you will use. Then enter that directory.

    mkdir bin
    cd bin

Let’s now compile the latest version of ABCtoolbox. For this you first
have to download the latesest code from the git repository. You can do
this as follows:

    git clone --depth 1 https://bitbucket.org/wegmannlab/abctoolbox.git

Here, the argument “–depth 1” implies tat you only get the latest
version, as opposed to the whole history of the program. The code is now
cloned into a diectory called abctoolbox.

Now lets compile ABCtoolbox! This is done easily by changing into the
directory of the code and using the provided make file by simply typing

    cd abctoolbox
    make

Now let’s move the `ABCtoolbox` binary into the bin folder, and then
change into that directory.

    mv ABCtoolbox ../
    cd ..

Next make sure ABCtoolbox is executable.

    chmod +x ABCtoolbox

Next lets download also two other binaries: the latest versions of
`fastsimcoal2`, a coalescent simulator, and `arlsumstat`, a program to
calculate summary statistics on fastsimcoal2 output. Starting with
former:

    wget http://cmpg.unibe.ch/software/fastsimcoal2/downloads/fsc26_linux64.zip
    unzip fsc26_linux64.zip

Then copy the fastsimcoal esecutable into the bin folder and make sure
it is executable. When doing so, we will also rename the binary so the
exercises below work regardless of the version you have.

    cp fsc2[0-9]_linux64/fsc2[0-9] fastsimcoal2
    chmod +x fastsimcoal2

Now we do the same for `arlsumstat`.

    wget http://cmpg.unibe.ch/software/arlequin35/linux/arlsumstat_linux.zip
    unzip arlsumstat_linux.zip
    cp arlsumstat_linux/arlsumstat*64bit arlsumstat
    chmod +x arlsumstat

Finally, and to have easy access to these executables, add the bin
folder to your PATH.

    PATH="`pwd`:${PATH}"

You can easily test if that worked. Simply try to launch one of the
programs from another directory:

    cd ..
    ABCtoolbox
    fastsimcoal2

### Step A2: Generate simulations with known parameters as truth set

To get a feeling for ABC, we will play around with a simple model of
constant population size. You can can get all necessary files as a zip
file from the git repo:

    wget https://github.com/ColonnaLab/EMBO_popgen/raw/main/popgen2021/Daniel_Wegmann/ABCDemo.tar.gz
    tar -xf ABCDemo.tar.gz

Then enter the folder `ABDdemo`.

    cd ABCDemo

Let us begin by generating some pseudo-observed data. To generate
simulations, we use `fastsimcoal2`, which requires an input file
specifying the model. Here we will use the file `constsize_obs.par`.
Have a look at that file and locate the two important parameters:

1.  the population size and
2.  the mutation rate.

You may also consult the manual of fastsimcoal, which is available
[here](http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal26.pdf).

Use the following command to generate the simulations (run fastsimcoal2
without arguments to get some help):

    fastsimcoal2 -i constsize_obs.par -n 1

Here, the option `-n 1` implies that a single simulation is conducted.

`fastsimcoal2` will generate a directory with the same name as the input
will. In that directory you will find a file with the ending .arp, which
contains the simulated data.

### Step A3: Calculate summary statistics from the observed data

We will next summarize the pseudo-observed data with summary statistics.
For that we will use the command line version of Arlequin (called
`arlsumstat`), which requires two setting files to run correctly:

1.  the file arl\_run.ars that tells arlsumstats what statistics are to
    be calculated
2.  ssdefs.txt, which specified what summary statistics to be printed.

The files provided here are set to calculate essentially everything
possible. Have a look at those file, but note that both files can be
generated and modified with the graphical version of Arlequin.

To run `arlsumstat` on the generated data, use the following command
line (run `arlsumstat` without arguments first to understand the
arguments):

    arlsumstat constsize_obs/constsize_obs_1_1.arp constsize.obs 0 1

The results are written to the file constsize.obs. Have a look at it!

We will now use ABC to estimate the parameters we used to genrate this
data.

### Step A4: Generate simulations with ABCtoolbox

`ABCtoolbox` allows you to use external programs to conduct simulations
with parameters taken from prior distributions and to calculate summary
statistics from the resulting data. It is configured by means of an
input file and a file defining the model parameters and their prior
distributions. Here we will use the files `constsize.input` and
`constsize.est`, respectively. Have a look at these files and try to
figure out how the use of external programs is specified. Note that
`fastimcoal2` is to be run with an input file on its own. Here we will
use the file `constsize.par`, which contains tags for the parameters
that are to be sampled from the prior distributions.

Also, have a look at how the prior are specified. Why are we not using
the logunif prior distribution provided by `ABCtoolbox` directly on
`N_NOW`?

To run ABCtoolbox with these input files, simply run

    ABCtoolbox constsize.input

The resulting simulations are found in the file
`sims_constsize_sampling1.txt`. Also have a look at the file
`childOutput.txt`, which contains the output of all the `fastsimcoal2`
and `arlsumstat` calls.

Since generating these simulations may take some time, I already
prepared a file containing about 50,000 of them. You can now add those
to the 100 simulations you generated with the following command:

    tail -n+2 sims_constsize_50K.txt >> sims_constsize_sampling1.txt

### Step A5: Estimate posterior distributions with ABCtoolbox

Let’s next use ABC-GLM to infer posterior distributions. Again, the
settings telling `ABCtoolbox` to perform parameter estimation are
provided in an input file. Here we will use the file `estimate.input` -
have a look at it and try to understand what ABCtoolbox will do when
using it.

Then launch ABCtoolbox with this input file as follows:

    ABCtoolbox estimate.input 

ABCtoolbox will now complain that there are highly correlated statistics
and not perform an estimation. In order to automatically prune one of
any pair of statistics that are highly correlated, add the option
`pruneCorrelatedStats` to the input file and rerun the estimation. The
estimation should now have finished successfully. If that was the case,
the output was written to a series of files with tag
`ABC_estimation_constsize_`, as was specified in the input file.

### Step A6: Plot posterior distributions in R

`ABCtoolbox` ships with a bunch of R scripts to plot posterior
distributions and other metrics. However, I recommend to plot them
yourself so that you get a feel for the structure of the output. For
this, simply open an R terminal in your working directory and follow the
steps below.

Should your computer not allow X forwarding from the server, you can
also output all plots to a pdf file. To open pdf output in R, simply
type `pdf("filename", width=5, height=5)`, where filename is the name of
the outputfile (e.g. `posteriors.pdf`). After plotting, you need to
close the pdf output using `dev.off()`.

Begin by plotting the marginal posterior distributions of the population
size (`N_NOW`) and the mutation rate (`MUTATION_RATE`) as follows.

    post <- read.table("ABC_estimation_constsize_model0_MarginalPosteriorDensities_Obs0.txt", header=T)
    par(mfrow=c(1,2))
    plot(post$LOG10_N_NOW, post$LOG10_N_NOW.density, type='l')
    plot(post$LOG10_MUTATION_RATE, post$LOG10_MUTATION_RATE.density, type='l')

Now let’s compare these marginal estimates to the 2D joint posterior we
also estimated

    twoD <- read.table("ABC_estimation_constsize_model0_jointPosterior_1_2_Obs0.txt", header=T)
    x <- unique(twoD$LOG10_N_NOW)
    y <- unique(twoD$LOG10_MUTATION_RATE)
    z <- matrix(twoD$density, nrow=length(x), byrow=T)
    contour(x,y,z, xlab="Log10(N_Now)", ylab="Log10(MutRate)")

Do you notice something? Very surprisingly, we can not estimate N and mu
together! So let’s try to estimate *θ* = 2*N**μ* instead. (wonder why
not *θ* = 4*N**μ*? Well, `fastsimcoal2` takes the haploid population
size).

### Step A7: Generate simulations with ABCtoolbox for theta

The model for theta is defined in the est file `constsize_theta.est`.
Note that the same par file can be used!

Can you generate an input file constsize\_theta.input for ABCtoolbox to
generate simulations under this model? Hint: copy the file from the N
and Mu model and modify it. Make sure you change the output tag to
`sims_constsize_theta_`. Then, you can generate simulations as follows:

    ABCtoolbox constsize_theta.input

Since generating these simulations may take some time, I already
prepared a file containing about 50,000 of them. You can now add those
to the 100 simulations you generated with the following command:

    tail -n+2 sims_constsize_theta_50K.txt >> sims_constsize_theta_sampling1.txt

Step A8: Estimate the posterior distributionon theta with ABCtoolbox
--------------------------------------------------------------------

Now estimate *θ* using `ABCtoolbox`. For this, create an input file
`estimate_theta.input` by copying and then modifying the input file we
used for the previous model. When doing so, make sure to use the option
`pruneCorrelatedStats`. Also, change the output prefix to
`ABC_estimation_constsize_theta_` to avoid overwriting your results.
Then launch `ABCtoolbox` with this input file as follows:

    ABCtoolbox estimate_theta.input 

You can now plot the marginal posterior for *θ* in R. How well does the
estimate fit the parameters we used to generate the observed data? Check
the file `dna_singlepop_constsize_obs.par` to see what parameters we
used and remember that *θ* = 2*N**μ* (for haploid species).

Change the parameters in `constsize_obs.par` and regenerate the observed
data, re-calculate summary statistics from them and rerun the estimation
(no need to generate new simulations!). Can you accurately estimate *θ*?

Exercises B (Afternoon): Choice of Summary Statistics & Validation
------------------------------------------------------------------

### Step B1: Finding appropriate summary statistics using PLS

Some of the summary statistics are highly correlated, suggesting that we
may reduce the summary statistics space while retaining all the
information. One method to do so is via a Partial Least Squares (PLS)
transformation. For such a transformation you may use the R script
find\_pls.r, which is provided in the exercise directory. Simply call it
as follows:

    Rscript find_pls.r sims_constsize_theta_sampling1.txt

This script will create two files:

1.  one called `PLSdef_sims_dna_singlepop_constsize_theta_sampling1.txt`
    containing the PLS definitions, and
2.  one called
    `RMSE_sims_dna_singlepop_constsize_theta_sampling1.txt.pdf`, which
    contains the RMSE plot to determine the lowest number of PLS
    components that are required to retain most of the information.

Have a look at it. You will probably agree that a single PLS component
is capturing all the info (in fact, it is known that the number of
segregating sites S is a sufficient statistics for *θ*).

### Step B2: Transforming statistics via PLS definitions

`ABCtoolbox` offers the possibility to transform the statistics into PLS
components. While you may again write an input file for that, it is time
to learn how to use `ABCtoolbox` via the command line. To transform the
statistics using the PLS definitions you generate before, use ABCtoolbox
as follows:

    ABCtoolbox task=transform linearComb=PLSdef_sims_constsize_theta_sampling1.txt numLinearComb=1 input=sims_constsize_theta_sampling1.txt output=sims_constsize_theta_sampling1.txt.pls boxcox verbose

This will generate a new simulation file containing all parameters and,
instead of the summary statistics, PLS components. To now rerun the
estimation, you also need to transform the observed data the same way.
I’m sure you’ll manage to do that!

Once the observed data has also been transformed, rerun the estimation.
Note that you need to modify the input file to use the transformed files
now! How does the posterior look like? Any visible changes?

### Step B3: Validation of parameter inference: are the stats reproduced?

Validation is an important part of any ABC application. `ABCtoolbox`
offers plenty of validation tools - let’s look at some of them!

First, let’s check if the model is capable of reproducing the observed
summary statistics. For this purpose, `ABCtoolbox` offers two statistics
to compare the observed statistics with those simulated: the marginal
density and the Tukey depth. To perform tests using these two
statistics, you need to add the following two arguments to the input
file (note: the 1000 refer to the number of retained simulations to be
used when calculating p-values):

    marDensPValue 1000
    tukeyPValue 1000

Then, simply rerun the estimation and check either the output written to
screen or the file with tag “modelFit.txt”. Are these tests passed? Did
using PLS make a difference?

### Step B4: Validation of parameter inference: power and bias

Next, let’s ask `ABCtoolbox` to perform some estimation on
pseudo-observed data sets. This is easily done: simply add the following
argument to the input file:

    randomValidation 1000

When rerunning the estimation with this argument, `ABCtoolbox` will
select 1000 random simulations as pseudo-observed data sets and then
perform estimation on these. The results are written to the file with
tag `RandomValidation.txt`. Does the content of this file make sense to
you?

You can now use the content of this file to contrast the true value of
the simulation with the estimated one and compute correlations between
them:

    d <- read.table("ABC_estimation_constsize_theta_model0_RandomValidation.txt", header=T)
    plot(d$LOG10_THETA, d$LOG10_THETA_mode)
    cor(d$LOG10_THETA, d$LOG10_THETA_mode)

Additionally, you can check for biases in the posterior distributions by
using either the quantile to HDI. As discussed in class, these measures
should be distributed uniformly.

    hist(d$LOG10_THETA_quantile)
    hist(d$LOG10_THETA_HDI)

Are they? What may affect these distributions and in what way?

Exercises C (Afternoon): Moderl Choice
--------------------------------------

### Step C1: Setting up a competing model with a bottleneck

Let’s next generate some simulations for a model with a bottleneck. The
goal is to have three parameters: the current size, the magnitude of the
population size change and when that change happened. As for the model
of constant size, it does make sense to define the population size in
units of *θ* and to use priors on the log<sub>10</sub> scale.
Specifically, let us have a uniform priors on

1.  the current theta log10(THETA\_CUR) ~ U\[-4, -1\],
2.  the old population size before the bottleneck relative to the
    current size log10(OLD\_SIZE\_RELATIVE) ~ U\[-2,2\] and
3.  the time when the bottleneck occurred log10(T\_BOT) ~ U\[1,3\].

Can you generate an `est` and `par` file (named `bottleneck.est` and
`bottleneck.par`, respectively) for this model? Hint: copy the files of
the *θ* example and modify them.

Then, prepare an input file `bottleneck.input` to generate 100
simulations under this model and run it using `ABCtoolbox`. Choose
`sims_bottleneck` as the output file name.

Since generating enough simulations may take some time, I already
prepared a file containing about 50,000 of them. You can now add those
to the 100 simulations you generated with the following command:

    tail -n+2 sims_bottleneck_50K.txt >> sims_bottleneck_sampling1.txt

### Step C2: Finding statistics for model choice

Let us now check if ABC is capable of distinguishing between these two
models. To do so, we first need to find summary statistics that are
informative about these models. Unfortunately, PLS does not work for
this. However, `ABCtoolbox` implements a greedy search to find such
summary statistics for model choice. To invoke this greedy search, you
can use the prepared input file `findStatsModelChoice.input` by
launching `ABCtoolbox` with it:

    ABCtoolbox findStatsModelChoice.input 

This run will take quite some time. But once it finished, you should
find a file called
`ABC_findStats_greedySearchForBestStatisticsForModelChoice.txt`. This
file lists the power of all tested statistic combinations. Among the
combinations with highest power, choose a combinations with few
statistics and low correlation among them. Which is your combination?

Prepare an observed file containing only these statistics for the
observed data.

### Step C3: Performing model choice

Running model choice with `ABCtoolbox` is straight forward: simply add
additional models to the arguments `simName` and `params`. In our case,
first copy the file `estimate_theta.input` and add the bottleneck model
to simName and to params and change the output prefix. It shoudl like
this;

    simName sims_constsize_sampling1.txt;sims_bottleneck_sampling1.txt.old
    params 2;2-4 
    outputPrefix ABC_modelchoice_

Finally, run `ABCtoolbox` on this input file and check both the output
written to screen as well as the model fit file. Which model is
preferred?

### Step C4: Model choice validation

Model choice validation is invoked with the argument
`modelChoiceValidation`, followed by the number of pseudo-observed data
sets to be used. Add the following to your input file:

     modelChoiceValidation 1000

And then run `ABCtoolbox` on it. You should now find the two files
`ABC_modelchoice_confusionMatrix.txt` and
`ABC_modelchoice_modelChoiceValidation.txt` containing the model choice
results.

First have a look at the confusion matrix: is there power to distinguish
between these models?

The file `ABC_modelchoice_modelChoiceValidation.txt` can be used to
infer biases in the Bayes factor calculation. This is a rather complex
analysis for which an R script is provided. Run it as follows:

    Rscript Make_model_choice_power_plot.r

This scripts produces a plot called `model_choice_power_fig.pdf`
comparing the posterior probability as estimated via ABC to those
obtained empirically from the validation results. Can you trust the ABC
posterior probabilities?
