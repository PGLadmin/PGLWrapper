After a "near-death" experience, I decided to include Main&Menu, FlashSubs, and MinPack in the repository. 
A seven-year old laptop is not a reliable place to store the code that tests how all the other codes work together.
These codes are even more of an embarrassment than SourceCodeUA. They are sometimes redundant and full of spaghetti. 
I always dreamed of upgrading FlashSubs to implement the methodologies of Michelsen & Mollerup, but never did.
My best suggestion is that you ignore these codes entirely. Nevertheless, you can construct a working project by inserting: 
Main&Menu, FlashSubs, MinPack, SourceCodeUa, SourceCode2a. 
Running requires the input and output directories and having them accessible. The default location is in the same
folder as the SourceCode__ dir's. The PGLInputDir is an element of GlobConst so it will be available everywhere 
once you set it however you want. It may be helpful to start with a (crudely) working project, then systematically 
upgrade the components that are most egregious. 

If your only interest is the PGLDLL.DLL, you only need: SourceCode2a, SourceCodeUa, Input, Output.

Good luck!
JRE 20250425

I added a dir for Python adaptations of the PGL code. The hope is that all methods evaluated in PGL 
will be accessible through Python. This may take time so check progress occasionally. 
JRE 20250520