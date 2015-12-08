###The MIT License

Copyright (c) 2015 University of California, San Francisco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[http://www.opensource.org/licenses/mit-license.php](http://www.opensource.org/licenses/mit-license.php)

##Overview

**Extract from publication:** 

Using [caDNAno](https://github.com/sdouglas/cadnano2), a computeraided design tool for DNA origami, we first routed the scaffold to approximate a 3D shape and exported the design from caDNAno as a text file in JSON format. Next, we used a custom Python script (provided here) and input the JSON file and desired number of unique staples (e.g. 10). We included the option to specify specific staple sequences as input, or to input a range of desired staple lengths and automatically generate sequences. The script determines each custom scaffold sequence by generating a random staple layout, repeatedly assigning a set of staple sequences to that layout, and then assigning the appropriate complementary bases to the scaffold. Because highly repetitive scaffold sequences can be difficult to synthesize, we used 10,000 iterations to create a library of custom scaffold sequences. We ranked the sequences according to total fraction of nucleotides that appear in a “repetitive” motif, defined as a 12 base window that appears morethan once in the scaffold. Then we selected the leastrepetitive scaffold sequence for synthesis.

####Userguide:





