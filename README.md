# Authors: Zhang & Hermans

# Dowser 

This program aims to search for a water molecule's energy positions in cavities of proteins having minimum energy. Read Zhang and Hermans' paper for more information:
  [Zhang & Hermans](https://www.ncbi.nlm.nih.gov/pubmed/9162944).

Dowser's oficial link: http://danger.med.unc.edu/hermans/dowser/dowser.htm 

However, it is hardly ever available, which is why I decided to upload the Dowser code to ensure program availability.

## Installation

Prior to installing Dowser, type:

    $ sudo apt update
    $ sudo apt install build-essential

You need gfotran compiler in linux to install Dowser. However, gfortran 77 is no longer available as a g77. Of course you could find it, but it is advisable to use the latest gfortran compiler, like Gfortran e.g.: 

    $ sudo apt-get install gfortran

You have to name your gfortran compiler absolute path as g77 in other Dowser versions, e.g.:

    $ sudo cp /usr/bin/gfortran /usr/bin/g77
    
The gfortran compiler has to be named gfortran for Dowser in this GitHub repository version (https://github.com/caraortizmah/dowser.git).

You have to install csh prior to install Dowser, e.g.:

     $ sudo apt update
     $ sudo apt install csh

Download Dowser by typing:

     $ git clone https://github.com/caraortizmah/dowser.git

Once you have downloaded the code, type:

     $ cd $PATH_dowser/
     $ chmod u+x Install
     $ chmod u+x bin/dowser
     $ ./Install
     
When installation ends, you need to export the Dowser program location to your PATH. For BASH, edit your .bashrc or .bash_profile file by typing:

     export DOWSER=   # write the Dowser absolute installation path here
     export DOW_MACH=linux
     export PATH=$PATH:$DOWSER/bin:$DOWSER/bin/$DOW_MACH

For CSH type:

     setenv DOWSER    # write the Dowser absolute installation path here
     setenv DOW_MACH linux
     set path = ( $path $DOWSER/bin $DOWSER/bin/$DOW_MACH )

And then, source the edited file (whichever in your case):

    $ source .bashrc

or

    $ source .bash_profile
    
That's all there is to it.

For further information about the installation go to: http://somagliablog.blogspot.com/2013/11/install-dowser-program.html
