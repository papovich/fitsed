# fitsed

add to your .bashrc something like:   

IDL_PATH=$IDL_PATH:+$IDL_DIR/lib:$IDL_DIR/examples:+$HOME/idl/astron/pro:+$HOME/Source/fitsed

it does still need the idl/astron library:

https://idlastro.gsfc.nasa.gov/

but that should be it other than the basic IDL stuff (I tried hard to get all our hand-made routines into the fitsed source files!) 

3. get the bc files (they have to use *my* file name convention
   because I haven’t fixed that yet).  Contact me if you want them.

4. cd to the example_phot directory, and edit the fitsed.param file to point to all the right files and directories (I’m hoping I don’t have to explain it). 

5. type at the command line: 

fitsed fitsed.param

you can also run it from within IDL>  fitsed,'fitsed.param' (should work).

6. from IDL you can run the 'play.pro' script to see best fit
   results. 
