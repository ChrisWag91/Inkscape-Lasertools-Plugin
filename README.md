# Inkscape-Lasertools-Plugin
A simple plugin to convert Inkscape vector graphics to Gcode for DIY laser engravers

## Version 0.9
Tested with Inkscape 0.92
This project is based on the Gcodetools extension.

## How to use Lasertools
Watch this Video for a quick introduction
https://youtu.be/NhUvRJsa4D0

step by step workflow:

1. import picture to Inkscape
2. set the size of the page to the size you want your engraving to be [under File/Document Properties/Page]
3. IMPORTANT: set *Display Units* and *Units* to *mm* and set *Scale x* to 1 [under File/Document Properties/Page]
4. convert your picture into a vector graphics [Path/Trace Bitmap]
5. remove the bitmap so that only the path is left
6. open Lasertools and set your parameters [under Extensions/Lasertools]

Make sure the specified directory exists.
Only use the Live Preview for smaller engravings.
Calculating the Gcode can take quite some time on larger engravings. I would advice to test with a small engraving.

## Installation
Copy the .py and .inx file to the Inkscape Extentions folder.

under linux:    /usr/share/inkscape/extensions/ 
under windows   C:\Program Files\Inkscape\share\extensions 

## Warning
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY Without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

License text: 
https://github.com/nevir/readable-licenses/blob/master/markdown/GPLv2-LICENSE.md

## Open Bugs
- Scaling factor needs to be Set to 1
- displayed units need to be set to mm
- performance on larger engravings very poor [in some cases the calculation of the toolpaths can take hours]

## Features to be implemented / tested
- make offsets work
- check if *Passes* feature works
- make a help section in the ui with the step by step instructions

## Help appreciated
If you find bugs or want to implement a feature or function please feel free to commit a change to the project.
