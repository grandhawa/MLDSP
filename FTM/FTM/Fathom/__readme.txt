FATHOM TOOLBOX FOR MATLAB
by David L. Jones, PhD
djones14@mail.usf.edu

This a collection of MATLAB functions I've written during my daily work as a fish ecologist and fisheries oceanographer. I am releasing it to the public under the GNU General Public License (version 2) to encourage the sharing of code and prevent duplication of effort.

If you use this Toolbox to produce a report or publication, PLEASE include the citation listed below to help justify maintenance and further development of the code. I'd also appreciate bug reports and suggestions for improvements. While I've made every attempt to write code that provides accurate and precise results, the functions in this Toolbox are provided as is, with no guarantees.

CITATION:
Jones, D. L. 2014. The Fathom Toolbox for Matlab: multivariate ecological and oceanographic data analysis. College of Marine Science, University of South Florida, St. Petersburg, Florida, USA. Available from: http://www.marine.usf.edu/user/djones/

INSTALLATION:
1) unzip the compressed FTM.zip file
2) place the 'Fathom' folder on your computer where your Matlab toolboxes are stored
3) Add the 'Fathom' folder to your Matlab path: Set Path... > Add Folder...
Use the Matlab 'help' command to learn what arguments are required by each function. For example:

==========================================================================
>> help f_stnd
  - standardize values of a matrix, column-wise (= z-scores)
 
  Usage: Xstd = f_stnd(x,{y});
 
  x    = matrix to standardize
  y    = optional, when present X will be standardized according to Y
 
  Xstd = standardized matrix
 
  SEE ALSO: f_center, f_ranging, f_transform
==========================================================================

Most functions have detailed documentation and references added as text comments at the beginning of each file. To view this information, load a *.m file in a text editor. Check out the 'examples' folder located within the main 'Fathom' folder to see demonstrations of how to use many of the functions. Additional information is also available in the 'doc' folder.
