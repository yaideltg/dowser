# Authors: Zhang & Hermans

# Dowser 

This program aims to search energy positions of a water molecule in cavities of proteins having minimum energy. Read Zhang and Hermans' paper for more information:
  [Zhang & Hermans](https://www.ncbi.nlm.nih.gov/pubmed/9162944).

Dowser is oficially in the following link: http://danger.med.unc.edu/hermans/dowser/dowser.htm however, it is hardly ever 
   available. That is why I decided to uplodad the Dowser code to mantain the availability of the program.

## Installation

Prior to install Dowser type:

    $ sudo apt update
    $ sudo apt install build-essential

In linux, you need gfotran compiler to intall Dowser. However gfortran 77 is no longer available as a g77. Of course you could find it but it is advisable to use the latest gfortran compiler like Gfortran e.g.: 

    $ sudo apt-get install gfortran

In other Dowser versions, you have to name your gfortran compiler absolute path as g77 e.g.:

    $ sudo cp /usr/bin/gfortran /usr/bin/g77
    
For Dowser in this (https://github.com/caraortizmah/dowser.git) GitHub repository version, gfortran compiler has to be named as gfortran.

Download Dowser typing:

     $ git clone https://github.com/caraortizmah/dowser.git

Once you download the code, type:

     $ cd $PATH_dowser/
     $ chmod u+x Install
     $ ./Install
     
When installation ends type:
     
     $ source dowserinit
     
However, If you have some issues launching Dowser please check your Dowser program location:
     
You need to export the Dowser program location to your PATH. For BASH, edit your .bashrc or .bash_profile file typing:

     export DOWSER=   # put here the Dowser absolute installation path
     export DOW_MACH=linux
     export PATH=$PATH:$DOWSER/bin:$DOWSER/bin/$DOW_MACH

For CSH type:

     setenv DOWSER    # put here the Dowser absolute installation path
     setenv DOW_MACH linux
     set path = ( $path $DOWSER/bin $DOWSER/bin/$DOW_MACH )

And then, source the edited file (in any case):

    $ source .bashrc

or

    $ source .bashrc
    
That's all.

Further information about the installation: http://somagliablog.blogspot.com/2013/11/install-dowser-program.html
